library(tidyverse)

vcf = read_tsv('vcf_439_output.tsv') # tableized VCF file of genotypes
phs = read_tsv('phase_output.tsv') # output of Python script
vnames = read_tsv('variant-names.tsv') # descriptive names for each variant, e.g. M129V

vnames$altid = paste0(vnames$pos, '-', vnames$alt)

allele_names = read_tsv('allele_names.tsv') # descriptive names for each allele, e.g. 129V
allele_names %>%
  mutate(refid = paste0(pos,'-',ref)) %>%
  select(id=refid, name=refname) -> ref_names
allele_names %>%
  mutate(altid = paste0(pos,'-',alt)) %>%
  select(id=altid, name=altname) -> alt_names
anames = rbind(ref_names, alt_names)

phs$a1name = anames$name[match(phs$allele1, anames$id)]
phs$a2name = anames$name[match(phs$allele2, anames$id)]

phs$pos1 = as.integer(gsub('-.*','',phs$allele1))
phs$pos2 = as.integer(gsub('-.*','',phs$allele2))

phs %>%
  filter(supporting_pairs > 0) %>%
  arrange(sample, pos1, pos2, desc(supporting_pairs)) %>%
  group_by(sample, pos1, pos2) %>%
  slice(1:2) %>%
  ungroup() %>% 
  select(sample, allele1, allele2, a1name, a2name, supporting_pairs) -> phs_browse

phs_browse %>%
  inner_join(vnames, by=c('allele1'='altid')) %>%
  filter(grepl('129',a2name)) %>%
  mutate(cis129=a2name) %>%
  select(sample, aa_change, cis129) -> phs1

phs_browse %>%
  inner_join(vnames, by=c('allele2'='altid')) %>%
  filter(grepl('129',a1name))  %>%
  mutate(cis129=a1name) %>% 
  select(sample, aa_change, cis129) -> phs2

rbind(phs1, phs2) %>%
  filter(aa_change %in% pathogenic_variants_in_mgh_study) %>%
  mutate(indiv = substr(iv,1,4)) -> phs_output

write_tsv(phs_output, 'mgh_c129_phases.tsv')
