import pysam
import itertools
import os
import io
import sys
from collections import Counter
from firecloud import fiss

# Get the Google billing project name and workspace name
billing_project = os.environ['WORKSPACE_NAMESPACE']
workspace=os.path.basename(os.path.dirname(os.getcwd()))
output_bucket = os.environ['WORKSPACE_BUCKET'] + "/"
bucket = 'gs://' # insert Google Cloud address here

# Verify that we've captured the environment variables
sys.stderr.write("Billing project: " + billing_project)
sys.stderr.write("Workspace: " + workspace)
sys.stderr.write("Workspace storage bucket: " + bucket)

workspace_id = bucket[5:]
bucket_in_console = "https://console.cloud.google.com/storage/browser/{}".format(workspace_id)

# print(bucket_in_console)
# print(workspace)
# print(bucket)

with open('bams.tsv', 'r') as f:
	bams = [line.strip() for line in f.readlines()]

variants_path= '{}vcf/vcf_439_output.tsv'.format(output_bucket)
!gsutil cp $variants_path .
!ls

# hard-coded coordinates - this script was run on the GRCh37 BAMS before re-alignment to hg38
orf_text = 'chr20:4,679,867-4,680,628'
orf_start_coord = 4679867
orf_end_coord = 4680628

class Variant:
	def __init__(self, chrom, pos, ref, alt, refname='', altname='', name=''):
		self.chrom = chrom
		self.pos = pos
		self.ref = ref
		self.alt = alt
		self.refname = refname if refname != '' else str(pos) + '-' + ref
		self.altname = altname if altname != '' else str(pos) + '-' + alt
		if name == '':
			self.name = self.chrom + '-' + str(self.pos) + '-' + self.ref + '-' + self.alt
		else:
			self.name = name
	def __str__(self):
		return self.chrom + '-' + str(self.pos) + '-' + self.ref + '-' + self.alt + ' (' + self.refname + '/' + self.altname + ')'

class ReadSupport:
	def __init__(self, variant):
		self.variant = variant
		self.ref_qnames = []
		self.alt_qnames = []

def getReadSupport(samfile, variant):
	readSup = ReadSupport(variant)
	iter = samfile.pileup(variant.chrom, variant.pos, variant.pos+1)
	for pileupcolumn in iter:
		if pileupcolumn.pos == variant.pos:
			for pileupread in pileupcolumn.pileups:
				if pileupread.query_position is None:
					continue
				base = pileupread.alignment.query_sequence[pileupread.query_position-1]
				qname = pileupread.alignment.query_name
				if base == variant.ref:
					readSup.ref_qnames.append(qname)
				elif base == variant.alt:
					readSup.alt_qnames.append(qname)
	return readSup

class PhaseSupport:
	def __init__(self, readSupA, readSupB):
		self.variant_a = readSupA.variant
		self.variant_b = readSupB.variant
		self.refa_refb = len(set(readSupA.ref_qnames) & set(readSupB.ref_qnames))
		self.refa_altb = len(set(readSupA.ref_qnames) & set(readSupB.alt_qnames))
		self.alta_refb = len(set(readSupA.alt_qnames) & set(readSupB.ref_qnames))
		self.alta_altb = len(set(readSupA.alt_qnames) & set(readSupB.alt_qnames))
	#
	def summarize(self):
		return  self.variant_a.refname + '/' + self.variant_b.refname + ': ' + str(self.refa_refb) + '\n' + \
				self.variant_a.refname + '/' + self.variant_b.altname + ': ' + str(self.refa_altb) + '\n' + \
				self.variant_a.altname + '/' + self.variant_b.refname + ': ' + str(self.alta_refb) + '\n' + \
				self.variant_a.altname + '/' + self.variant_b.altname + ': ' + str(self.alta_altb)
	def output(self, samplename):
		return  samplename + '\t' + self.variant_a.refname + '\t' + self.variant_b.refname + '\t' + str(self.refa_refb) + '\n' + \
				samplename + '\t' + self.variant_a.refname + '\t' + self.variant_b.altname + '\t' + str(self.refa_altb) + '\n' + \
				samplename + '\t' + self.variant_a.altname + '\t' + self.variant_b.refname + '\t' + str(self.alta_refb) + '\n' + \
				samplename + '\t' + self.variant_a.altname + '\t' + self.variant_b.altname + '\t' + str(self.alta_altb) + '\n'
	def __str__(self):
		return self.summarize()

def phaseVariants(samfile, variant_a, variant_b, samplename='', outfile=sys.stdout):
	aSupport = getReadSupport(samfile, variant_a)
	bSupport = getReadSupport(samfile, variant_b)
	phase = PhaseSupport(aSupport,bSupport)
	outfile.write(phase.output(samplename))

def parseVariants(varpath):
	variants = []
	with open(varpath) as f:
		colnames = f.readline().strip().split('\t')
		for line in f.readlines():
			cells = line.strip().split('\t')
			chrom = cells[colnames.index('chrom')]
			pos = int(cells[colnames.index('pos')])
			ref = cells[colnames.index('ref')]
			alt = cells[colnames.index('alt')]
			variants.append(Variant(chrom, pos, ref, alt))
	return(variants)

def indivVariants(varpath, sampleid):
	variants = []
	with open(varpath) as f:
		colnames = f.readline().strip().split('\t')
		for line in f.readlines():
			cells = line.strip().split('\t')
			sample = cells[colnames.index('sample')]
			if sample == sampleid:            
				chrom = cells[colnames.index('chrom')]
				pos = int(cells[colnames.index('pos')])
				ref = cells[colnames.index('ref')]
				alt = cells[colnames.index('alt')]
				variants.append(Variant(chrom, pos, ref, alt))
	return(variants)

def phaseAll(bampath, variants, outfile):
	samfile = pysam.AlignmentFile(bampath, "rb")
	vrnts = variants
	i = 1
	for pair in itertools.combinations(vrnts,2):
		sys.stderr.write('\rOn '+bampath+'... phasing pair '+str(i)+'...')
		phaseVariants(samfile, *pair, samplename=str.replace(bampath,'.bam',''), outfile=outfile)
		i += 1
	sys.stderr.write('...done!\n')

def subsetVariants(variants, minpos, maxpos):
	return ([x for x in variants if x.pos >= minpos and x.pos <= maxpos])

def phaseOneBucketBam(sampleid, variantfile, outfile):
	with open(outfile, 'a') as f:
		#sys.stderr.write('Copying bam...\n')
		#bucket_bam_path = ("{}/"+sampleid+".bam").format(bucket)
		#bucket_bai_path = ("{}/"+sampleid+".bai").format(bucket)
		#!gsutil cp $bucket_bam_path .
		#!gsutil cp $bucket_bai_path .
		codon129 = Variant('20',4680251,'A','G')
		sys.stderr.write('Parsing variants...\n')
		if '/' in sampleid: # strip subdirectories
			sampleid = sampleid[(sampleid.index('/')+1):]
		indiv_variants = indivVariants(variantfile,sampleid)
		indiv_variants = subsetVariants(indiv_variants, orf_start_coord-2500, orf_end_coord+2500)
		coords = [v.pos for v in indiv_variants]
		if not codon129.pos in coords:
			indiv_variants.append(codon129)
		sys.stderr.write('Phasing...\n')
		with open(outfile, 'a') as f:
			phaseAll(sampleid+'.bam', indiv_variants, f)

with open('phase_output.tsv', 'w') as f:
	f.write('sample\tallele1\tallele2\tsupporting_pairs\n')
for sampleid in bams:
	phaseOneBucketBam(sampleid,'vcf_439_output.tsv','phase_output.tsv')
output_bucket_dest = '{}output/phase_output.tsv'.format(output_bucket)
!gsutil cp phase_output.tsv $output_bucket_dest

