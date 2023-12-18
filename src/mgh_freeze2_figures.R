##############
# STARTUP
##############

overall_start_time = Sys.time()
tell_user = function(x) { cat(file=stderr(), x); flush.console() }

tell_user('Loading required packages...')

options(stringsAsFactors=F)
if(interactive()) setwd('~/d/sci/src/mgh_prnp_freeze2/')
suppressMessages(library(tidyverse))
suppressMessages(library(janitor))
suppressMessages(library(weights))
suppressMessages(library(openxlsx))
suppressMessages(library(Hmisc))
suppressMessages(library(MASS)); select = dplyr::select; summarize = dplyr::summarize
suppressMessages(library(pdftools))
suppressMessages(library(magick))


##############
# FUNCTIONS & CONSTANTS
##############

cat(file=stderr(), 'done.\nSetting constants and functions...')

percent = function(x, digits=0, signed=F) gsub(' ','',paste0(ifelse(x <= 0, '', ifelse(signed, '+', '')),formatC(100*x,format='f',digits=digits),'%'))

upper = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) + sds*sd(x)/sqrt(sum(!is.na(x)))
}
lower = function(x, ci=0.95) { 
  alpha = 1 - ci
  sds = qnorm(1-alpha/2)
  mean(x) - sds*sd(x)/sqrt(sum(!is.na(x)))
}

alpha = function(rgb_hexcolor, proportion) {
  hex_proportion = sprintf("%02x",round(proportion*255))
  rgba = paste(rgb_hexcolor,hex_proportion,sep='')
  return (rgba)
}
ci_alpha = 0.35 # degree of transparency for shading confidence intervals in plot

meansd = function(x, digits=1) {
  paste0(formatC(mean(x, na.rm=T),format='f',digits=digits),'±',formatC(sd(x, na.rm=T),format='f',digits=digits))
}

rangestring = function(x, digits=1) {
  paste0(formatC(min(x, na.rm=T),format='f',digits=digits),'-',formatC(max(x, na.rm=T),format='f',digits=digits))
}

summary_stats = function(x, digits=1) {
  paste0(meansd(x, digits=digits),' (',rangestring(x, digits=digits),')')
}

fraction_positive = function(x) {
  paste0(sum(x=='+',na.rm=T),'/',sum(x %in% c('+','-'), na.rm=T))
}

clipdist = function(x, minx, maxx) {
  return (pmin(maxx,pmax(minx,x)))
}

bin_age = function(age) {
  case_when(age < 20 ~ '<20',
            age >= 20 & age < 90 ~ paste0(floor(age/5)*5,'-',floor(age/5)*5+4),
            age >= 90 ~ '≥90')
}


rbind_files = function(path, grepstring) {
  all_files = list.files(path, full.names=T)
  these_files = all_files[grepl(grepstring,all_files)]
  if (exists('rbound_table')) rm('rbound_table')
  for (this_file in these_files) {
    this_tbl = read_delim(this_file, col_types=cols()) %>% clean_names()
    this_tbl$file = gsub('.*\\/','',gsub('\\.[tc]sv','',this_file))
    if (exists('rbound_table')) {
      rbound_table = rbind(rbound_table, this_tbl)
    } else {
      rbound_table = this_tbl
    }
  }
  return (rbound_table)
}


clipcopy = function(tbl) {
  clip = pipe("pbcopy", "w")  
  write.table(tbl, file=clip, sep = '\t', quote=F, row.names = F, na='')
  close(clip)
}


dviz = function(tbl,
                xlims,
                ylims,
                xvar,
                yvar,
                colorvar=NULL,
                pchvar=NULL,
                pchbg='#FFFFFF',
                pchcex=1,
                xcols=character(0), # columns that group with x TBD
                xats=xbigs,
                xbigs,
                xlwds=1,
                xbiglabs=xbigs,
                xaxcex=1,
                xlabcex=1,
                xlabline=1.5,
                yats=ybigs,
                ybigs,
                ylwds=1,
                ybiglabs=ybigs,
                yaxcex=1,
                ylabcex=1,
                ylabline=2.25,
                log,
                mar=c(3,3,1,1),
                jitamt=0.1,
                randseed=1,
                boxwidth=NA,
                barwidth=NA,
                polygon=NA,
                arrowlength=0.05,
                test=NA,
                control_group=NA,
                xlab='',
                ylab='',
                legtext=NULL,
                legcol='#000000',
                legtextcol=legcol,
                leglty=1,
                legpch=20,
                leglwd=1,
                legcex=1,
                crosshairs=F
) {
  
  if (is.null(pchvar)) {
    pchvar='pch'
    tbl$pch = 19
  }
  
  if (is.null(colorvar)) {
    colorvar='color'
    tbl$color = '#000000'
  }
  
  tbl %>%
    mutate(x=!!as.name(xvar), y=!!as.name(yvar), color=!!as.name(colorvar), pch=!!as.name(pchvar)) %>%
    select(x, y, color, pch, all_of(xcols)) -> tbl
  
  if (!crosshairs) {
    xcols = c('x', xcols)
  }
  
  tbl %>%
    group_by(color, across(all_of(xcols))) %>%
    summarize(.groups='keep',
              n = n(),
              mean = mean(y),
              l95 = mean(y) - 1.96 * sd(y) / sqrt(n()),
              u95 = mean(y) + 1.96 * sd(y) / sqrt(n()),
              median = median(y),
              q25 = quantile(y, .25),
              q75 = quantile(y, .75),
              cv = sd(y) / mean(y),
              x_mean = mean(x),
              x_l95  = mean(x) - 1.96 * sd(x) / sqrt(n()),
              x_u95  = mean(x) + 1.96 * sd(x) / sqrt(n())) %>%
    ungroup() %>%
    mutate(l95 = ifelse(l95 < 0 & log %in% c('xy','y'), min(ylims), l95)) %>%
    arrange(x) -> tbl_smry
  
  if (crosshairs) {
    tbl_smry$x = tbl_smry$x_mean
  }
  
  par(mar=mar)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log=log)
  axis(side=1, at=xlims, labels=NA, lwd.ticks=0)
  if (!is.null(xats)) {
    axis(side=1, at=xats, tck=-0.025, lwd.ticks=xlwds, labels=NA)
  }
  if (!is.null(xbigs)) {
    axis(side=1, at=xbigs, tck=-0.05, lwd.ticks=xlwds, labels=NA)
    axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5, labels=xbiglabs, cex.axis=xaxcex)
  }
  mtext(side=1, line=xlabline, text=xlab, cex=xlabcex)
  if (!is.null(yats)) {
    axis(side=2, at=yats, tck=-0.025, labels=NA)
  }
  if (!is.null(ybigs)) {
    axis(side=2, at=ybigs, tck=-0.05, labels=NA)
    axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=-0.3, labels=ybiglabs, cex.axis=yaxcex)
  }
  mtext(side=2, line=ylabline, text=ylab, cex=ylabcex)
  if (crosshairs) {
    jitamt = 0
  }
  set.seed(randseed)
  points(x=jitter(tbl$x,amount=jitamt), y=tbl$y, col=alpha(tbl$color,ci_alpha), pch=tbl$pch, bg=pchbg)
  
  if (crosshairs) {
    segments(x0=tbl_smry$x_l95, x1=tbl_smry$x_u95, y0=tbl_smry$mean, col=tbl_smry$color, lwd=1.5)
    segments(x0=tbl_smry$x_mean, y0=tbl_smry$l95, y1=tbl_smry$u95, col=tbl_smry$color, lwd=1.5)
  }
  if (!is.na(barwidth)) {
    segments(x0=tbl_smry$x-barwidth, x1=tbl_smry$x+barwidth, y0=tbl_smry$mean, col=tbl_smry$color, lwd=1.5)
    arrows(x0=tbl_smry$x, y0=tbl_smry$l95, y1=tbl_smry$u95, code=3, angle=90, length=arrowlength, col=tbl_smry$color, lwd=1.5)
  }
  if (!is.na(boxwidth)) {
    rect(xleft=tbl_smry$x-boxwidth, xright=tbl_smry$x+boxwidth, ybottom=tbl_smry$q25, ytop=tbl_smry$q75, border=tbl_smry$color, lwd=1.5, col=NA)
    segments(x0=tbl_smry$x-boxwidth, x1=tbl_smry$x+boxwidth, y0=tbl_smry$median, col=tbl_smry$color, lwd=1.5)
  }
  
  if (!is.na(polygon)) {
    for (clr in unique(tbl_smry$color)) {
      if (polygon=='iqr') {
        subs = subset(tbl_smry, color==clr & !is.na(q25) & !is.na(q75))
        points( x=  subs$x, y=  subs$median, type='l', lwd=2, col=subs$color)
        polygon(x=c(subs$x, rev(subs$x)), y=c(subs$q25, rev(subs$q75)), col=alpha(subs$color, ci_alpha), border=NA)
      } else if (polygon=='ci') {
        subs = subset(tbl_smry, color==clr & !is.na(l95) & !is.na(u95))
        points( x=  subs$x, y=  subs$mean, type='l', lwd=2, col=subs$color)
        polygon(x=c(subs$x, rev(subs$x)), y=c(subs$l95, rev(subs$u95)), col=alpha(subs$color, ci_alpha), border=NA)
      }
    }
  }
  
  if (!is.na(test)) {
    testfun = get(test) # e.g. ks.test
    control_color = control_group
    tbl_smry$p = as.numeric(NA)
    for (i in 1:nrow(tbl_smry)) {
      this_x = tbl_smry$x[i]
      this_rows = tbl$x == this_x & tbl$color == tbl_smry$color[i]
      ctrl_rows = tbl$x == this_x & tbl$color == control_group
      test_obj = suppressWarnings(testfun(tbl$y[this_rows],  tbl$y[ctrl_rows]))
      tbl_smry$p[i] = test_obj$p.value
    }
    tbl_smry$p_symb = ''
    tbl_smry$p_symb[!is.na(tbl_smry$p) & tbl_smry$p < 0.05] = '*'
    tbl_smry$p_symb[!is.na(tbl_smry$p) & tbl_smry$p < 0.01] = '**'
    tbl_smry$p_symb[!is.na(tbl_smry$p) & tbl_smry$p < 0.001] = '***'
    text(x=tbl_smry$x[tbl_smry$color != control_color], y=max(ylims)*.95, labels=tbl_smry$p_symb[tbl_smry$color != control_color])
  }
  
  if (!is.null(legtext)) {
    par(xpd=T)
    legend(x=max(xlims),y=max(ylims),legtext,col=legcol,text.col=legtextcol,pch=legpch,lwd=leglwd, cex=legcex,lty=leglty, bty='n')
    par(xpd=F)
  }
  
  return(tbl_smry)
  
}


smooth_p = function(p) {
  if (length(p)<2) {
    return (p)
  }
  while (any(cummin(p) != p)) { # while p is non-monotonic
    for (i in 2:length(p)) { # starting from the second dilution and going up,
      if (p[i] > p[i-1]) { # if the current one has higher p than the last one...
        # select all dilutions below the current one that have lower p...
        indices_to_average = c(which(1:length(p) < i & p < p[i]), i)
        # and average them all together with the current value.
        p[indices_to_average] = mean(p[indices_to_average])
      }
    }
  }
  return (p)
}

cmode = function(x) {
  as.integer(table(x)[which.max(table(x))])
}

spearman_karber = function(x, n, d) {
  # to formally define LLQ and ULQ, assume all samples would be positive at one dilution higher
  # and all negative at one dilution lower
  if (length(x)==1) {
    dd = 1
  } else {
    dd = mean(abs(diff(d))) # find dd, the difference between dilutions
    if (any(abs(diff(d)) - dd > .05)) { # .05 provides some numeric tolerance, esp in CSF dilution where rounding error is substantial
      stop("x are not uniformly distributed\n")
    }
  }
  n_mode = cmode(n) # the most common number of replicates
  x = c(n_mode, x, 0) # to define llq & ulq, assume it'd be 100% pos at 1 higher dilution and 0% pos at 1 lower dilution
  n = c(n_mode, n, n_mode) # assume the higher and lower dilution would've been the modal N
  d = c(d[1]+dd, d, d[length(d)]-dd) # assume the higher and lower dilution would have followed the same dilution factor
  llq = mean(d[1:2]) # define the LLQ as what you'd get if you had 100% pos at 1 higher dilution but 0% at highest dilution tested
  ulq = mean(d[(length(d)-1):length(d)]) # define the ULQ as what you'd get if you had 0% pos at 1 lower dilution but 100% at lowest dilution tested
  # if a non-zero positive is observed at a dilution higher than where zeroes are observed, call this noise and remove it
  for (j in 2:length(x)) {
    if (x[j-1]==0 & x[j] > 0) {
      x[j] = 0
    }
  }
  p = smooth_p(x/n) # smooth the proportion positive to be non-increasing
  nd = length(x) # number of dilutions
  # find i, the dilutest dilution at which 100% positive replicates
  if (all(p < 1)) {
    i = 1
  }  
  else  {
    i = max(which(p==1))
  }
  sd50 = (d[i]+dd) - dd*(1/2 + sum(p[i:nd])) # calculate log10(sd50)
  se = sqrt(dd^2 * sum(p*(1-p)/(n-1))) # calculate standard error
  # invert estimate and limits of quantification, in order to be in titer space rather than dilution space:
  return_value = cbind(estimate=-sd50, se=se, llq=-llq, ulq=-ulq)
  return (return_value)
}




##############
# OUTPUT STREAMS
##############



tell_user('done.\nCreating output streams...')

text_stats_path = 'display_items/stats_for_text.txt'
write(paste('Last updated: ',Sys.Date(),'\n',sep=''),text_stats_path,append=F) # start anew - but all subsequent writings will be append=T
write_stats = function(...) {
  write(paste(list(...),collapse='',sep=''),text_stats_path,append=T)
  write('\n',text_stats_path,append=T)
}

supplement_path = 'display_items/supplement.xlsx'
supplement = createWorkbook()
# options("openxlsx.numFmt" = "0.00") # this looks better for residuals but terrible for p values and weeks post-dose
supplement_directory = tibble(name=character(0), title=character(0))
write_supp_table = function(tbl, title='') {
  # write Excel sheet for supplement
  table_number = length(names(supplement)) + 1
  table_name = paste0('s',formatC(table_number,'d',digits=0,width=2,flag='0'))
  addWorksheet(supplement,table_name)
  bold_style = createStyle(textDecoration = "Bold")
  writeData(supplement,table_name,format(as.data.frame(tbl),digits=3,na.encode=F),na.string='',keepNA=T,headerStyle=bold_style,withFilter=T)
  freezePane(supplement,table_name,firstRow=T)
  saveWorkbook(supplement,supplement_path,overwrite = TRUE)
  
  # also write tab-sep version for GitHub repo
  write_tsv(as_tibble(format(as.data.frame(tbl),digits=3)),paste0('display_items/table-',table_name,'.tsv'), na='')
  
  # and save the title in the directory tibble for later
  assign('supplement_directory',
         supplement_directory %>% add_row(name=table_name, title=title),
         envir = .GlobalEnv)
}




##############
# ANALYTIC DATASET
##############



tell_user('done.\nReading analytic dataset...')

params = read_tsv('data/mgh_study_params.tsv', col_types=cols())
leg = read_tsv('data/legend.tsv', col_types=cols())
conv_table = read_tsv('../mgh_freeze2/nosync/analytic/conv_table.tsv', col_types=cols()) %>%
  mutate(symptomatic = months_from_onset > 0)
gt_dem = read_tsv('../mgh_freeze2/nosync/analytic/gt_dem.tsv', col_types=cols())
all_biomarkers = read_tsv('../mgh_freeze2/nosync/analytic/all_biomarkers.tsv', col_types=cols())
endpoint_smry = read_tsv('../mgh_freeze2/nosync/analytic/endpoint_smry.tsv', col_types=cols())
bsyn_chaps = read_tsv('../mgh_freeze2/nosync/analytic/bsyn_chaps.tsv', col_types=cols())
bsqc = read_tsv('../mgh_freeze2/nosync/analytic/bsqc.tsv', col_types=cols())
ttau_chaps = read_tsv('../mgh_freeze2/nosync/analytic/ttau_chaps.tsv', col_types=cols())
two_taus = read_tsv('../mgh_freeze2/nosync/analytic/two_taus.tsv', col_types=cols())
tau_cmpr = read_tsv('../mgh_freeze2/nosync/analytic/ttau_study_v_poscon.tsv', col_types=cols())
qcs = read_tsv('../mgh_freeze2/nosync/input/qcs.tsv', col_types=cols())
qcleg = read_tsv('../mgh_freeze2/nosync/input/qc_legend.tsv', col_types=cols())
rtq_avkin_plot = read_tsv('../mgh_freeze2/nosync/analytic/rtq_avkin_plot.tsv', col_types=cols())

##############
# Table 1
##############

tell_user('done.\nCreating Table 1...')

gt_dem %>%
  group_by(mut) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() %>%
  inner_join(params, by=c('mut'='gt')) %>%
  arrange(x) %>%
  mutate(count = paste0(n,' ',mut)) %>%
  group_by(cat) %>%
  summarize(.groups='keep', muts = toString(count)) %>%
  ungroup() %>%
  mutate(status = case_when(cat=='-' ~ 'control',
                            cat=='+' ~ 'carrier')) -> gt_smry
gt_dem %>%
  group_by(status) %>%
  summarize(.groups='keep',
            n = n(),
            m = sum(sex=='M'),
            f = sum(sex=='F'),
            mean_age = mean(age_last_seen),
            sd_age = sd(age_last_seen),
            mean_years_elapsed = mean(years_elapsed),
            sd_years_elapsed = sd(years_elapsed),
            total_visits = sum(n_visits),
            total_csf = sum(n_csf),
            total_plasma = sum(n_plasma)) %>%
  ungroup() %>%
  mutate(age = paste0(formatC(mean_age,format='f',digits=1),'±',formatC(sd_age,format='f',digits=1)),
         sex = paste0(m,'M / ',f,'F'),
         elapsed = paste0(formatC(mean_years_elapsed,format='f',digits=1),'±',formatC(sd_years_elapsed,format='f',digits=1))) %>%
  inner_join(gt_smry, by=c('status')) -> table1_raw

table1_raw %>%
  select(status, n, sex, age, elapsed, total_visits, total_csf, total_plasma, muts) -> table1

write_tsv(table1, 'display_items/table_1.tsv')


####
# Supplementary tables
####



tell_user('done.\nCreating supplementary tables...')

conv_table %>%
  distinct(indiv) %>%
  mutate(individual_code = LETTERS[row_number()]) -> conv_codes

conv_table %>%
  select(-path_variant) %>%
  inner_join(all_biomarkers, by=c('iv','indiv','visit')) %>%
  inner_join(conv_codes, by='indiv') %>%
  select(-indiv, -iv) %>%
  relocate(individual_code) %>%
  mutate(age_of_onset_bin = bin_age(age_of_onset)) %>%
  select(individual_code, visit, path_variant, codon129, age_of_onset_bin, months_from_onset, 
         csf_rtquic, csf_rtquic_replicates, csf_tau, csf_nfl, csf_bsyn, csf_prp, plasma_nfl, plasma_gfap,
         moca, mrc_scale) -> all_converter_markers

write_supp_table(all_converter_markers, 'All biomarker values from all visits by individuals who developed active disease.')

all_biomarkers %>%
  filter(!(indiv %in% conv_table$indiv)) %>%
  mutate(group = case_when(cat=='+' ~ 'non-converting mutation carrier',
                           cat=='-' ~ 'mutation-negative control')) %>%
  group_by(group) %>%
  summarize(.groups='keep',
            n=n(),
            rtq_positive=fraction_positive(csf_rtquic),
            csf_tau = summary_stats(csf_tau, digits=0),
            csf_nfl = summary_stats(csf_nfl, digits=0),
            csf_bsyn = summary_stats(csf_bsyn, digits=0),
            csf_prp = summary_stats(csf_prp, digits=0),
            plasma_nfl = summary_stats(plasma_nfl, digits=1),
            plasma_gfap = summary_stats(plasma_gfap, digits=1),
            moca = summary_stats(moca, digits=0),
            mrc = summary_stats(mrc_scale, digits=0)) %>%
  ungroup() -> carrier_control_stats

write_supp_table(carrier_control_stats, 'Means, standard deviations, and ranges of biomarker values from all visits by participants who did not develop active disease, by mutation status.')

all_biomarkers %>%
  filter(!is.na(csf_prp)) %>% 
  group_by(indiv, mut) %>%
  summarize(.groups='keep', mean_prp = mean(csf_prp, na.rm=T)) %>%
  ungroup() %>%
  group_by(mut) %>%
  summarize(.groups='keep', n=n(), mean=mean(mean_prp, na.rm=T), sd=sd(mean_prp, na.rm=T)) %>%
  ungroup() %>%
  mutate(rel = mean/max(mean)) %>%
  inner_join(select(params, gt, x), by=c('mut'='gt')) %>%
  arrange(x) %>%
  select(-x) %>%
  rename(mutation=mut) -> csf_prp_by_mutation

all_biomarkers %>%
  filter(!is.na(csf_prp)) %>% 
  group_by(indiv, mut) %>%
  summarize(.groups='keep', mean_prp = mean(csf_prp, na.rm=T)) %>%
  ungroup() %>%
  mutate(mut = factor(mut, levels=params$gt)) -> csf_prp_lmdata
m = lm(mean_prp ~ mut, data=csf_prp_lmdata)
summary(m)$coefficients %>%
  as_tibble(rownames='group') %>%
  mutate(mut = case_when(group=='(Intercept)' ~ 'none',
                         TRUE ~ gsub('mut','',group))) %>%
  mutate(lm_pval = case_when(mut=='none' ~ 1,
                             TRUE ~ `Pr(>|t|)`)) -> lm_coefs
csf_prp_by_mutation$lm_pvalue = lm_coefs$lm_pval[match(csf_prp_by_mutation$mutation, lm_coefs$mut)]

write_supp_table(csf_prp_by_mutation, 'Mean CSF PrP concentration (ng/mL) by mutation.')


all_biomarkers %>%
  filter(!is.na(csf_prp)) %>%
  select(indiv, visit, csf_prp) %>%
  arrange(indiv, visit) %>%
  group_by(indiv) %>%
  slice(1) %>%
  rename(baseline=csf_prp) -> baselines

all_biomarkers %>%
  filter(!is.na(csf_prp)) %>%
  filter(indiv %in% conv_table$indiv | years_on >= 3) %>%
  distinct(indiv) -> fu3 # ppl with at least 3y followup

all_biomarkers %>%
  filter(!is.na(csf_prp)) %>%
  select(indiv, iv, visit, disp_years_on, mut, csf_prp, csf_rtquic) %>%
  inner_join(fu3, by='indiv') %>%
  inner_join(baselines, by='indiv') %>%
  mutate(rel = csf_prp/baseline) %>%
  left_join(conv_table, by=c('indiv','iv')) %>%
  mutate(converter = !is.na(age_of_onset)) %>%
  inner_join(select(params, gt, color, cat), by=c('mut'='gt')) %>%
  mutate(plot_color = case_when(converter ~ color,
                                mut=='none' ~ '#000000',
                                TRUE ~ '#D782D7')) %>%
  mutate(plot_lwd = case_when(converter ~ 4,
                              mut=='none' ~ 0.5,
                              TRUE ~ 0.5)) %>%
  mutate(plot_pch = case_when(converter ~ 19,
                              mut=='none' ~ 20,
                              TRUE ~ 20)) %>%
  mutate(plot_cex = case_when(converter ~ 1.25,
                              mut=='none' ~ 0.5,
                              TRUE ~ 0.5)) %>%
  select(indiv, iv, converter, disp_years_on, age_of_onset, csf_prp, csf_rtquic, rel, plot_color, plot_lwd, plot_pch, plot_cex, mut, cat) -> cselisa_plotdata


cselisa_plotdata %>%
  filter(!converter) %>%
  group_by(indiv, cat) %>%
  summarize(.groups='keep', n=n(), mean=mean(csf_prp), sd=sd(csf_prp), cv=sd(csf_prp)/mean(csf_prp)) %>%
  ungroup() %>%
  filter(n > 1) %>%
  group_by(cat) %>%
  summarize(.groups='keep', n_indiv = length(unique(indiv)), n_samples_total=sum(n), mean_cv = mean(cv)) %>%
  ungroup() %>%
  mutate(disp_cv = percent(mean_cv,digits=1)) -> csf_prp_cvs

csf_prp_cvs %>%
  mutate(group = case_when(cat=='-' ~ 'control',
                           cat=='+' ~ 'non-converting carrier')) %>%
  select(group, n_indiv, n_samples_total, mean_cv) -> csf_prp_cvs_out

write_supp_table(csf_prp_cvs_out, 'Long-term test-retest reliability of CSF PrP.')



####
# Figure S1
####

tell_user('done.\nCreating Figure S1...')
consort = image_read_pdf("data/consort_diagram.pdf")
resx = 300
png('display_items/figure-s1.png',width=6.5*resx,height=8*resx)
plot(consort)
unnecessary_message = dev.off()


####
# Figure S2
####

tell_user('done.\nCreating Figure S2...')

plot_csf_prp_bymut = function() {

  csf_prp_lmdata %>%
    inner_join(params, by=c('mut'='gt')) %>%
    mutate(rel = mean_prp / mean(mean_prp[mut=='none'])) %>%
    mutate(y = rel) -> csf_prp_indivs
  csf_prp_mut_summary = dviz(csf_prp_indivs, xlims=range(csf_prp_indivs$x) + c(-0.5, 0.5), ylims=c(0,1.75), log='',
                             xvar='x', yvar='rel', colorvar='color', xcols=c('mut'), mar=c(3,4,1,1),
                             xbigs=NULL, xlwds=0, 
                             yats=0:8/4, ybigs=0:4/2, ybiglabs=percent(0:4/2),
                             ylab='CSF PrP (% control mean)', jitamt=0.1, xlab='', ylabline = 3,
                             barwidth=.33)
  abline(h=1, lty=3, lwd=0.5)
  use_params = params %>% filter(gt != 'poscon')
  mtext(side=1, line=0.25, at=use_params$x, text=use_params$gt, col=use_params$color, cex=0.8, font=2)
  
  csf_prp_by_mutation %>%
    inner_join(params, by=c('mutation'='gt')) %>%
    filter(mutation!='none') -> csf_prp_modeldata
  
  mtext(side=3, line=0, at=csf_prp_modeldata$x, col=csf_prp_modeldata$color, cex=0.5, text=paste0('P=',formatC(csf_prp_modeldata$lm_pvalue, format='g', digits=2)))
}

resx=600
png('display_items/figure-s2.png',width=3.25*resx,height=3.25*resx,res=resx)
plot_csf_prp_bymut()
unnecessary_message = dev.off()


#####
# Functions to produce Figure 1 panels
#####

plot_csf_prp_tr = function() {
  par(mar=c(3,4,3,1))
  xlims = c(0,6)
  ylims = c(0.25, 6)
  xats = 0:20/2
  xbigs = 0:10
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
  
  par(xpd=T)
  mtext(side=3, line=0.5, text='CSF PrP')
  par(xpd=F)
  
  yats = 0:60/10
  ybigs = c(0.25, 0.5, 1, 2, 3, 4, 5, 6)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5)
  mtext(side=1, line=1.5, text='years from baseline visit', cex=0.8)
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, labels=percent(ybigs-1,signed=T), tck=-0.05, las=2, lwd=0, line=-0.25, cex.axis=.8)
  mtext(side=2, line=2.0, text='∆', cex=0.8, las=2)
  
  abline(h=1, lty=3)
  cselisa_plotdata %>%
    arrange(converter, age_of_onset, indiv) %>%
    distinct(indiv) -> cselisa_indivs
  for (this_indiv in cselisa_indivs$indiv) {
    subs = subset(cselisa_plotdata, indiv==this_indiv) %>%
      left_join(conv_table, by='iv') %>%
      mutate(symptomatic = replace_na(symptomatic, F))
    points(subs$disp_years_on, subs$rel, col=subs$plot_color, lwd=subs$plot_lwd/2, lty=3, type='l')
    points(subs$disp_years_on[!subs$symptomatic], subs$rel[!subs$symptomatic], type='l', lty=1, lwd=subs$plot_lwd, col=subs$plot_color)
    points(subs$disp_years_on[subs$symptomatic],  subs$rel[subs$symptomatic],  type='l', lty=1, lwd=subs$plot_lwd, col=subs$plot_color)
    
    points(subs$disp_years_on, subs$rel, col=subs$plot_color, pch=subs$plot_pch, cex=subs$plot_cex)
    if (this_indiv %in% conv_table$indiv) {
      points(subs$disp_years_on[subs$csf_rtquic=='+'], subs$rel[subs$csf_rtquic=='+'], col=leg$color[leg$disp=='RT-QuIC positive timepoint'], pch=leg$pch[leg$disp=='RT-QuIC positive timepoint'], lwd=leg$ptlwd[leg$disp=='RT-QuIC positive timepoint'])
      points(subs$disp_years_on[subs$months_from_onset > 0], subs$rel[subs$months_from_onset > 0], col=leg$color[leg$disp=='symptomatic timepoint'], pch=leg$pch[leg$disp=='symptomatic timepoint'], lwd=leg$ptlwd[leg$disp=='symptomatic timepoint'], cex=1.25)
    }
    
  }

}


resx=300
png(paste0('for_slides/figure-csf-prp-tr.png'),width=3.25*resx,height=3.25*resx,res=resx)
plot_csf_prp_tr()
unnecessary_message = dev.off()



bioms_meta = tibble(disp = c('CSF PrP','CSF NfL','CSF T-tau','CSF β-syn','plasma NfL','plasma GFAP', 'plasma beta-syn'),
                    varname = c('csf_prp','csf_nfl','csf_tau','csf_bsyn','plasma_nfl','plasma_gfap', 'plasma_bsyn'),
                    llq = c(2, 5.4, 1.68, 15.9, 5.4, 2.744, 31.8),
                    ymax = c(150, 4000, 525, 1500, 60, 800, 100),
                    unit = c('ng/mL', 'pg/mL', 'pg/mL', 'pg/mL', 'pg/mL', 'pg/mL', 'pg/mL')) %>%
  filter(!(varname %in% c('csf_prp','plasma_bsyn')))
bioms_meta %>%
  add_column(n_individuals = numeric(nrow(bioms_meta)),
             n_samples = numeric(nrow(bioms_meta)),
             n_samples_carrier = numeric(nrow(bioms_meta)),
             n_samples_control = numeric(nrow(bioms_meta)),
             test_retest_mean_cv_without_converters = numeric(nrow(bioms_meta)),
             overall_y_intercept = numeric(nrow(bioms_meta)),
             overall_annual_change = numeric(nrow(bioms_meta)),
             overall_annual_change_pval = numeric(nrow(bioms_meta)),
             overall_carrier_diff = numeric(nrow(bioms_meta)),
             overall_carrier_pval = numeric(nrow(bioms_meta)),
             nccarriers_y_intercept = numeric(nrow(bioms_meta)),
             nccarriers_annual_change = numeric(nrow(bioms_meta)),
             nccarriers_annual_change_pval = numeric(nrow(bioms_meta)),
             controls_y_intercept = numeric(nrow(bioms_meta)),
             controls_annual_change = numeric(nrow(bioms_meta)),
             controls_annual_change_pval = numeric(nrow(bioms_meta))) -> bioms_meta


plot_biom_age = function(biom) {
  
  i = which(bioms_meta$varname==biom)
  
  par(mar=c(3,4,3,1))
  xlims = c(20, 80)
  ylims = c(0, bioms_meta$ymax[i])
  xats = 0:20*5
  xbigs = 0:10*10
  if (max(ylims) > 550) {
    yats = 0:1000 * 100
    ybigs = 0:200 * 500
  } else {
    yats = 0:1000 * 10
    ybigs = 0:200 * 50
  }
  
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
  par(xpd=T)
  mtext(side=3, line = 0.5, at=max(xlims), text=bioms_meta$disp[i])
  par(xpd=F)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5)
  mtext(side=1, line=1.5, text='age', cex=0.8)
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=-0.5)
  mtext(side=2, line=2.5, text=bioms_meta$unit[i], cex=0.8)
  all_biomarkers %>%
    select(-mut, -path_variant) %>%
    inner_join(gt_dem, by='indiv') %>%
    mutate(value = !!as.name(biom)) %>%
    filter(!is.na(value)) %>%
    rename(rtq = csf_rtquic) %>%
    mutate(age = disp_age + disp_years_on) %>%
    mutate(true_age = age_at_first_visit + years_on) %>%
    inner_join(select(params, gt, color), by=c('mut'='gt')) %>%
    mutate(carrier = mut != 'none') %>%
    mutate(converter = indiv %in% conv_table$indiv) %>%
    left_join(select(conv_table, iv, symptomatic), by=c('iv')) %>%
    mutate(symptomatic = replace_na(symptomatic, F)) %>%
    mutate(lwd = case_when(converter ~ 2, TRUE ~ 0.5)) %>%
    mutate(kcolor = case_when(mut == 'none' ~ '#000000',
                              converter ~ color,
                              TRUE ~ '#D782D7')) %>%
    select(indiv, age, true_age, value, kcolor, lwd, mut, carrier, converter, rtq, symptomatic) %>% 
    arrange(converter, indiv, age) -> this_biom
  par(xpd=T)
  for (ind in unique(this_biom$indiv)) {
    subs = subset(this_biom, indiv==ind)
    points(x=subs$age, y=subs$value, type='l', lty=3, lwd=subs$lwd/2, col=subs$kcolor) # dashed line default (will show up when connecting pre-symp and symp timepoints)
    points(x=subs$age[!subs$symptomatic], y=subs$value[!subs$symptomatic], type='l', lty=1, lwd=subs$lwd[!subs$symptomatic], col=subs$kcolor[!subs$symptomatic]) # solid line between pre-symp timepoints
    points(x=subs$age[subs$symptomatic], y=subs$value[subs$symptomatic], type='l', lty=1, lwd=subs$lwd[!subs$symptomatic], col=subs$kcolor[!subs$symptomatic]) # solid line between symp timepoints
    points(x=subs$age, y=subs$value, pch=20, cex=0.5, col=subs$kcolor)
    if (subs$converter[1]) {
      points(x=subs$age, y=subs$value, pch=19, cex=1.5, col=subs$kcolor)
      points(x=subs$age[subs$rtq=='+'], y=subs$value[subs$rtq=='+'], col=leg$color[leg$disp=='RT-QuIC positive timepoint'], pch=leg$pch[leg$disp=='RT-QuIC positive timepoint'], lwd=leg$ptlwd[leg$disp=='RT-QuIC positive timepoint'], cex=1.25)
      points(x=subs$age[subs$symptomatic], y=subs$value[subs$symptomatic], col=leg$color[leg$disp=='symptomatic timepoint'], pch=leg$pch[leg$disp=='symptomatic timepoint'], lwd=leg$ptlwd[leg$disp=='symptomatic timepoint'], cex=0.9)
    }
  }
  par(xpd=F)
  abline(h=bioms_meta$llq[i], lty=3, lwd=.25, col='blue')
  mtext(side=4, cex=0.6, line=0.25, at=bioms_meta$llq[i], text='LLQ', col='blue', las=2)
  
  # plot the single best fit curve for non-converting carriers
  m = lm(log(value) ~ true_age, data=subset(this_biom, mut != 'none' & !converter))
  x = seq(20,80,.1)
  y = exp(predict(m, list(true_age=x)))
  points(x, y, type='l', lwd=1.5, col='#D782D7')

  # plot the single best fit curve for controls
  m = lm(log(value) ~ true_age, data=subset(this_biom, mut == 'none' & !converter))
  x = seq(20,80,.1)
  y = exp(predict(m, list(true_age=x)))
  points(x, y, type='l', lwd=1.5, col='#000000')

}


plot_biom_delta = function(biom) {
  
  i = which(bioms_meta$varname==biom)
  
  all_biomarkers %>%
    select(-mut, -path_variant) %>%
    inner_join(gt_dem, by='indiv') %>%
    mutate(value = !!as.name(biom)) %>%
    filter(!is.na(value)) %>%
    rename(rtq = csf_rtquic) %>%
    inner_join(select(params, gt, color), by=c('mut'='gt')) %>%
    mutate(carrier = mut != 'none') %>%
    mutate(converter = indiv %in% conv_table$indiv) %>%
    mutate(lwd = case_when(converter ~ 2, TRUE ~ 0.5)) %>%
    mutate(pch = case_when(converter ~ 19, TRUE ~ 20)) %>%
    mutate(kcolor = case_when(mut == 'none' ~ '#000000',
                              converter ~ color,
                              TRUE ~ '#D782D7')) %>%
    left_join(select(conv_table, iv, symptomatic), by='iv') %>%
    mutate(symptomatic = replace_na(symptomatic, F)) %>%
    select(iv, indiv, visit, value, kcolor, pch, lwd, mut, carrier, converter, rtq, disp_years_on, symptomatic) %>% 
    arrange(converter, indiv, visit) -> this_biom
  
  this_biom %>%
    group_by(indiv) %>%
    slice(1) %>%
    select(indiv, baseline=value) -> baselines
  
  this_biom %>%
    inner_join(baselines, by='indiv') %>%
    mutate(change = value/baseline) %>%
    left_join(select(conv_table, iv, months_from_onset), by='iv') %>%
    mutate(months_from_onset = case_when(converter ~ months_from_onset,
                                         TRUE ~ disp_years_on*12/6 - 12*5.5)) -> this_biom
  
  
  par(mar=c(3,3,3,2))
  xlims = c(-6, 1)
  xats = -48:12/12
  xbigs = -4:1
  ylims = c(.25, 6)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
  yats = 0:60/10
  ybigs = c(0.25, 0.5, 1, 2, 3, 4, 5, 6)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5)
  mtext(side=1, line=1.5, text='years from onset', cex=0.8)
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, labels=percent(ybigs-1,signed=T), tck=-0.05, las=2, lwd=0, line=-0.5, cex.axis=.8)
  mtext(side=2, line=2.0, text='∆', cex=0.8, las=2)
  abline(h=1, lty=1, lwd=.25)
  abline(v=0, col='red',  lwd=0.5)
  mtext(side=3, line=0.1, cex=0.5, col='red', at=0, text='onset', font=3)
  par(xpd=T)
  for (this_indiv in unique(this_biom$indiv)) {
    subs = subset(this_biom, indiv==this_indiv)
    points(subs$months_from_onset/12, subs$change, type='l', lty=3, lwd=subs$lwd/2, col=subs$kcolor)
    points(subs$months_from_onset[!subs$symptomatic]/12, subs$change[!subs$symptomatic], type='l', lty=1, lwd=subs$lwd, col=subs$kcolor)
    points(subs$months_from_onset[subs$symptomatic]/12, subs$change[subs$symptomatic], type='l', lty=1, lwd=subs$lwd, col=subs$kcolor)
    points(subs$months_from_onset/12, subs$change, pch=subs$pch, col=subs$kcolor, bg='#FFFFFF')
    
    if (subs$converter[1]) {
      points(x=subs$months_from_onset[subs$rtq=='+']/12,    y=subs$change[subs$rtq=='+'], col=leg$color[leg$disp=='RT-QuIC positive timepoint'], pch=leg$pch[leg$disp=='RT-QuIC positive timepoint'], lwd=leg$ptlwd[leg$disp=='RT-QuIC positive timepoint'], cex=1.25)
      points(x=subs$months_from_onset[subs$symptomatic]/12, y=subs$change[subs$symptomatic], col=leg$color[leg$disp=='symptomatic timepoint'], pch=leg$pch[leg$disp=='symptomatic timepoint'], lwd=leg$ptlwd[leg$disp=='symptomatic timepoint'], cex=0.9)
      
    }
    
  }
  par(xpd=F)
  
}

for (i in 1:nrow(bioms_meta)) {
  biom = bioms_meta$varname[i]
  
  all_biomarkers %>%
    select(indiv, iv, cat, all_of(biom)) %>%
    rename(value = all_of(biom)) %>%
    filter(!is.na(value)) -> biom_values
  
  bioms_meta$n_individuals[i] = length(unique(biom_values$indiv))
  bioms_meta$n_samples[i] = nrow(biom_values)
  bioms_meta$n_samples_carrier[i] = length(unique(biom_values$iv[biom_values$cat=='+']))
  bioms_meta$n_samples_control[i] = length(unique(biom_values$iv[biom_values$cat=='-']))
  
  biom_values %>%
    filter(!indiv %in% conv_table$indiv) -> biom_values_for_cv
  biom_values_for_cv %>%
    group_by(indiv) %>%
    summarize(.groups='keep', n=n(), cv=sd(value)/mean(value)) %>%
    ungroup() %>%
    filter(!is.na(cv)) %>%
    summarize(mean_cv = mean(cv), n=n()) -> biom_cv
  bioms_meta$test_retest_mean_cv_without_converters[i] = biom_cv$mean_cv

  
  resx=300
  png(paste0('for_slides/biom_age_',bioms_meta$varname[i],'.png'),width=3.25*resx,height=3.25*resx,res=resx)
  plot_biom_age(bioms_meta$varname[i])
  unnecessary_message = dev.off()
  
  all_biomarkers %>%
    select(-mut, -path_variant) %>%
    inner_join(gt_dem, by='indiv') %>%
    mutate(value = !!as.name(biom)) %>%
    filter(!is.na(value)) %>%
    rename(rtq = csf_rtquic) %>%
    mutate(age = disp_age + disp_years_on) %>%
    mutate(true_age = age_at_first_visit + years_on) %>%
    inner_join(select(params, gt, color), by=c('mut'='gt')) %>%
    mutate(carrier = mut != 'none') %>%
    mutate(converter = indiv %in% conv_table$indiv) %>%
    mutate(lwd = case_when(converter ~ 2, TRUE ~ 0.5)) %>%
    mutate(kcolor = case_when(mut == 'none' ~ '#000000',
                              converter ~ color,
                              TRUE ~ '#D782D7')) %>%
    select(indiv, age, true_age, value, kcolor, lwd, mut, carrier, converter, rtq) %>% 
    arrange(converter, indiv, age) -> this_biom
 
  m = lm(log(value) ~ true_age + carrier, data=subset(this_biom, !converter))
  summary(m)
  bioms_meta$overall_y_intercept[i] = exp(summary(m)$coefficients['(Intercept)','Estimate'])
  bioms_meta$overall_annual_change[i] = exp(summary(m)$coefficients['true_age','Estimate'])-1
  bioms_meta$overall_annual_change_pval[i] = summary(m)$coefficients['true_age','Pr(>|t|)']
  bioms_meta$overall_carrier_diff[i] = exp(summary(m)$coefficients['carrierTRUE','Estimate'])-1
  bioms_meta$overall_carrier_pval[i] = summary(m)$coefficients['carrierTRUE','Pr(>|t|)']
  
  # plot the single best fit curve for non-converting carriers
  m = lm(log(value) ~ true_age, data=subset(this_biom, mut != 'none' & !converter))
  bioms_meta$nccarriers_y_intercept[i] = exp(summary(m)$coefficients['(Intercept)','Estimate'])
  bioms_meta$nccarriers_annual_change[i] = exp(summary(m)$coefficients['true_age','Estimate'])-1
  bioms_meta$nccarriers_annual_change_pval[i] = summary(m)$coefficients['true_age','Pr(>|t|)']
  
  # plot the single best fit curve for controls
  m = lm(log(value) ~ true_age, data=subset(this_biom, mut == 'none' & !converter))
  bioms_meta$controls_y_intercept[i] = exp(summary(m)$coefficients['(Intercept)','Estimate'])
  bioms_meta$controls_annual_change[i] = exp(summary(m)$coefficients['true_age','Estimate'])-1
  bioms_meta$controls_annual_change_pval[i] = summary(m)$coefficients['true_age','Pr(>|t|)']
  
  
  all_biomarkers %>%
    select(-mut, -path_variant) %>%
    inner_join(gt_dem, by='indiv') %>%
    mutate(value = !!as.name(biom)) %>%
    filter(!is.na(value)) %>%
    rename(rtq = csf_rtquic) %>%
    inner_join(select(params, gt, color), by=c('mut'='gt')) %>%
    mutate(carrier = mut != 'none') %>%
    mutate(converter = indiv %in% conv_table$indiv) %>%
    mutate(lwd = case_when(converter ~ 2, TRUE ~ 0.5)) %>%
    mutate(pch = case_when(converter ~ 19, TRUE ~ 20)) %>%
    mutate(kcolor = case_when(mut == 'none' ~ '#000000',
                              converter ~ color,
                              TRUE ~ '#D782D7')) %>%
    select(iv, indiv, visit, value, kcolor, pch, lwd, mut, carrier, converter, rtq, disp_years_on) %>% 
    arrange(converter, indiv, visit) -> this_biom
  
  this_biom %>%
    group_by(indiv) %>%
    slice(1) %>%
    select(indiv, baseline=value) -> baselines
  
  this_biom %>%
    inner_join(baselines, by='indiv') %>%
    mutate(change = value/baseline) %>%
    left_join(select(conv_table, iv, months_from_onset), by='iv') %>%
    mutate(months_from_onset = case_when(converter ~ months_from_onset,
                                         TRUE ~ disp_years_on*12/6 - 12*5.5)) -> this_biom
  
  resx=300
  png(paste0('for_slides/biom_delta_',bioms_meta$varname[i],'.png'),width=3.25*resx,height=3.25*resx,res=resx)
  plot_biom_delta(bioms_meta$varname[i])
  unnecessary_message = dev.off()
 
}

write_supp_table(bioms_meta, 'Descriptive statistics and log-linear model fits on CSF and plasma biomarkers.')



plot_rtq_endpoint = function() {
  par(mar=c(3,4,3,1))
  xlims = c(-4, 1.5)
  xats = -48:12/12
  xbigs = -4:1
  ylims = c(0.01, 3)
  plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='y')
  mtext(side=3, line=0.5, text='RT-QuIC endpoint')
  yats = rep(1:9,4)  * 10^rep(-2:1,each=9)
  ybigs = 10^(-2:1)
  axis(side=1, at=xats, tck=-0.025, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, labels=NA)
  axis(side=1, at=xbigs, tck=-0.05, lwd=0, line=-0.5)
  mtext(side=1, line=1.5, text='years from onset', cex=0.8)
  axis(side=2, at=yats, tck=-0.025, labels=NA)
  axis(side=2, at=ybigs, tck=-0.05, labels=NA)
  axis(side=2, at=ybigs, labels=ybigs, tck=-0.05, las=2, lwd=0, line=-0.25)
  mtext(side=2, line=2.0, text='seeds / µL', cex=0.8)
  llq = max(endpoint_smry$llq_per_ul)
  ulq = max(endpoint_smry$ulq_per_ul)
  #abline(h=c(llq, ulq), lty=3, col='blue')
  segments(x0=min(xats), x1=max(xats), y0=llq, lty=3, col='blue')
  text(x=min(xlims) + diff(xlims)*0.08, y=llq, pos=1, labels='LLQ', cex=.8, font=3, col='blue')
  abline(v=0, col='red')
  endpoint_smry$seeds_per_ul = pmax(endpoint_smry$seeds_per_ul, llq)
  endpoint_smry$pch = ifelse(endpoint_smry$seeds_per_ul <= llq, 21, 19)
  par(xpd=T)
  for (this_indiv in unique(conv_table$indiv)) {
    subs = subset(endpoint_smry, indiv==this_indiv)
    arrows(x0=subs$months_from_onset/12, y0=subs$l95, y1=subs$u95, code=3, angle=90, length=0.05, col=subs$color, lwd=0.5)
    points(subs$months_from_onset/12, subs$seeds_per_ul, type='l', lwd=3, col=subs$color)
    points(subs$months_from_onset/12, subs$seeds_per_ul, pch=subs$pch, col=subs$color, bg='#FFFFFF')
  }
  par(xpd=F)
  endpoint_smry$codon129 = conv_table$codon129[match(endpoint_smry$indiv, conv_table$indiv)]
  endpoint_smry %>%
    group_by(indiv) %>%
    slice_max(months_from_onset) -> last_sample
  par(xpd=T)
  text(x=last_sample$months_from_onset/12, y=last_sample$seeds_per_ul, pos=4, labels=last_sample$codon129, col=last_sample$color, cex=0.8)
  par(xpd=F)
}

resx=300
png('for_slides/csf_rtq_endpoint.png',width=3.25*resx,height=3.25*resx,res=resx)
plot_rtq_endpoint()
unnecessary_message = dev.off()





####
# Figure S3
####

tell_user('done.\nCreating Figure S3...')

resx=600
png('display_items/figure-s3.png', width=5.75*resx, height=3.25*resx, res=resx)

layout_matrix = matrix(c(1,2,3,
                         1,2,4), nrow=2, byrow=T)
layout(layout_matrix, widths=c(0.8,.25,.5))
panel = 1

par(mar=c(4,4,3,1))
lims = c(0,800)
ats = 0:8*100
plot(NA, NA, xlim=lims, ylim=lims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=ats, labels=NA)
axis(side=1, at=ats, line=-0.4, lwd=0, cex=0.6)
mtext(side=1, text='T-tau by Fujirebio ELISA (pg/mL)', line=2.25)
axis(side=2, at=ats, labels=NA)
axis(side=2, at=ats, line=-0.4, las=2, lwd=0, cex=0.6)
mtext(side=2, text='T-tau by Ella (pg/mL)', line=2.2)
points(two_taus$pgml_av, two_taus$ttau, pch=20)
abline(a=0, b=1, lty=3)
m = lm(ttau ~ pgml_av, data=two_taus)
abline(m, col='red', lwd=0.5)
intercept = coefficients(m)['(Intercept)']
slope = coefficients(m)['pgml_av']
pearson_obj = cor.test(two_taus$pgml_av, two_taus$ttau)
write_stats('Pearson correlation Ella vs. Fujirebio T-tau: r = ',
            formatC(pearson_obj$estimate, format='f',digits=2),', P = ',
            formatC(pearson_obj$p.value, format='e', digits=2),' * elisa')
write_stats('Best linear regression fit for Ella vs. Fujirebio T-tau: ella = ',
            formatC(intercept, format='f',digits=1),' + ',
            formatC(slope, format='f', digits=2),' * elisa')
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5); panel = panel + 1


tau_cmpr %>%
  mutate(diagnosis = case_when(id %in% qcs$sample[qcs$diagnosis=='suspected prion disease'] ~ 'suspected prion disease',
                         TRUE ~ 'study samples')) %>%
  inner_join(qcleg, by='diagnosis') %>%
  mutate(indiv = gsub('-.*','',id)) %>%
  group_by(x, color, indiv) %>%
  summarize(.groups='keep', ttau = mean(ttau)) %>%
  ungroup() -> tau_cmpr_plot

par(mar=c(10,4,3,1))
xlims = c(0.5, 2.5)
ylims = c(0, 1500)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
mtext(side=1, line=0.25, at=qcleg$x[qcleg$x %in% tau_cmpr_plot$x], text=qcleg$diagnosis[qcleg$x %in% tau_cmpr_plot$x], las=2, cex=0.7, col=qcleg$color[qcleg$x %in% tau_cmpr_plot$x])
axis(side=2, at=0:3*500, tck=-0.08, labels=NA)
axis(side=2, at=0:3*500, lwd=0, line=-0.5, labels=0:3*500, las=2)
axis(side=2, at=0:15*100, tck=-0.03, labels=NA)
mtext(side=2, text='CSF T-tau (pg/mL)', line=2.5)
set.seed(1)
points(jitter(tau_cmpr_plot$x, amount=0.25), tau_cmpr_plot$ttau, col=alpha(tau_cmpr_plot$color, ci_alpha), pch=19)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5); panel = panel + 1


two_taus %>%
  mutate(indiv = substr(iv,1,4)) %>%
  group_by(indiv) %>%
  summarize(.groups='keep', 
            ella_mean=mean(ttau), 
            ella_cv=sd(ttau)/mean(ttau),
            fuj_mean = mean(pgml_av),
            fuj_cv = sd(pgml_av)/mean(pgml_av)) %>%
  ungroup() -> ella_ttau_stats
ella_ttau_stats %>%
  summarize(ella_mean_cv = mean(ella_cv, na.rm=T),
            fuj_mean_cv = mean(fuj_cv, na.rm=T)) %>%
  ungroup() %>%
  pivot_longer(ella_mean_cv:fuj_mean_cv) %>%
  mutate(x = row_number()) %>%
  rename(mean_cv = value) %>%
  mutate(assay = case_when(grepl('ella',name) ~ 'Ella',
                           grepl('fuj',name) ~ 'ELISA')) -> ttau_cvs

write_stats('Test-retest CVs for CSF T-tau: Ella = ',
            percent(ttau_cvs$mean_cv[ttau_cvs$assay=='Ella'],digits=1),
            ' ELISA = ',
            percent(ttau_cvs$mean_cv[ttau_cvs$assay=='ELISA'],digits=1))

par(mar=c(3,4,3,1))
xlims = c(0.5, 2.5) 
ylims = c(0, .15)
yats = seq(0, .15, .01)
ybigs = seq(0, .15, .05)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='')
barwidth = 0.33
rect(xleft=ttau_cvs$x-barwidth, xright=ttau_cvs$x+barwidth, ybottom=rep(0, nrow(ttau_cvs)), ytop=ttau_cvs$mean_cv, col='#A9A990', border=NA)
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
mtext(side=1, line=0.25, text=ttau_cvs$assay, at=ttau_cvs$x, cex=0.8)
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=-0.5, labels=percent(ybigs), cex.axis=0.8)
mtext(side=2, line=2.0, text='mean test-retest %CV', cex=0.8)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5); panel = panel + 1

par(mar=c(3,4,3,1))
xlims = c(0.5, 2.5)
ylims = c(0,200)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, yaxs='i', xaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
mtext(side=1, at=c(1,2), line=0.25, text=c('neat','CHAPS'), cex=0.8)
for (i in 1:nrow(ttau_chaps)) {
  points(1:2, c(ttau_chaps$neat[i], ttau_chaps$chaps[i]), pch=20)
  points(1:2, c(ttau_chaps$neat[i], ttau_chaps$chaps[i]), type='l')
}
axis(side=2, at=0:2*100, las=2, lwd=0, line=-0.5)
axis(side=2, at=0:2*100, tck=-0.05, labels=NA)
axis(side=2, at=0:20*10, labels=NA, tck=-0.02)
mtext(side=2, line=2.0, text='CSF T-tau (pg/mL)', cex=0.8)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5); panel = panel + 1

tobj = t.test(ttau_chaps$neat, ttau_chaps$chaps, paired=T)
chaps_ratio = mean(ttau_chaps$chaps / ttau_chaps$neat)

write_stats('CSF T-tau with vs. without CHAPS: mean ',percent(chaps_ratio-1,signed=T,digits=1),', T test P value = ',tobj$p.value)


unnecessary_message = dev.off()





####
# Figure S4
####

tell_user('done.\nCreating Figure S4...')

resx=600
png('display_items/figure-s4.png', width=6.5*resx, height=3.25*resx, res=resx)

layout_matrix = matrix(c(1,2,
                  1,3), nrow=2, byrow=T)
layout(layout_matrix, widths=c(1,.5))
panel = 1

bsyn_llq = 7.95

bsqc %>%
  filter(!grepl('qc', sample_name)) %>%
  filter(!is.na(bsyn)) %>%
  group_by(biofluid, sample_name, dilution_factor) %>%
  summarize(.groups='keep',
            n = n(),
            cv = sd(bsyn)/mean(bsyn),
            mn = mean(bsyn)) %>%
  ungroup() %>%
  arrange(sample_name, dilution_factor) -> cv1

cv1 %>%
  filter(!is.na(mn) & !is.na(cv)) %>%
  group_by(biofluid, dilution_factor) %>%
  summarize(.groups='keep', mean_cv = mean(cv)) -> mean_cvs

all_biomarkers %>%
  select(iv, plasma_bsyn) %>%
  mutate(dilution = plasma_bsyn / bsyn_llq) -> instudy_plasma_bsyn

par(mar=c(3,4,3,9))
xlims = c(1.5, 10)
ylims = c(10, 10000)
xats = 1:10
xbigs = c(2, 4, 8)
yats = rep(1:9, 4) * 10^rep(1:4, each = 9)
ybigs = 10^(1:4)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='xy')
axis(side=1, at=xats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=xbigs, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-0.25)
mtext(side=1, line=2, text='dilution factor', cex=0.8)
x = seq(1,10,.01)
y = x * 7.95
points(x, y, col='gray', lty=3, type='l')
mtext(side=4, at=y[x==max(xlims)], line=0.25, las=2, text='LLQ', col='gray')
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, tck=-0.025, las=2, lwd=0, line=-0.25, labels=formatC(ybigs, format='d', big.mark=','), cex.axis=0.8)
mtext(side=2, line=2.5, text='β-syn (pg/mL)', cex=0.8)
for (smpl in unique(cv1$sample_name)) {
  
  subs = subset(bsqc, sample_name==smpl)
  color = subs$color[1]
  points(subs$dilution_factor, subs$bsyn, col=alpha(color, ci_alpha), pch=19)
  
  subs = subset(cv1, sample_name==smpl)
  points(subs$dilution_factor, subs$mn, type='l', col=color)
  points(subs$dilution_factor, subs$mn, pch='-', cex=2, col=color)
  
}
points(instudy_plasma_bsyn$dilution, instudy_plasma_bsyn$plasma_bsyn, pch=19, col=alpha(qcleg$color[qcleg$diagnosis=='study samples'],ci_alpha))
instudy_plasma_bsyn %>%
  group_by(dilution) %>%
  summarize(.groups='keep', n=n()) %>%
  ungroup() -> instudy_dilutions
text(instudy_dilutions$dilution, instudy_dilutions$dilution*bsyn_llq, pos=2, labels=paste0('N = ',instudy_dilutions$n), col=qcleg$color[qcleg$diagnosis=='study samples'])
csf_range = range(bsqc$bsyn[bsqc$biofluid=='CSF'])
plasma_range = range(c(bsqc$bsyn[bsqc$biofluid=='plasma'], instudy_plasma_bsyn$plasma_bsyn), na.rm=T)
axis(side=2, line=-1, tck=0.05, at=csf_range, labels=NA)
mtext(side=2, line=-1, at=sqrt(prod(csf_range)), text='CSF')
axis(side=2, line=-1, tck=0.05, at=plasma_range, labels=NA)
mtext(side=2, line=-1, at=sqrt(prod(plasma_range)), text='plasma')
par(xpd=T)
legend(x=max(xlims), y=max(ylims), qcleg$diagnosis, col=qcleg$color, pch=19, bty='n', cex=.6)
par(xpd=F)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5); panel = panel + 1

par(mar=c(3,3,3,1))
fluid_meta = tibble(biofluid=c('CSF','plasma'), color=c('#DCDCFF','#DAA520'), x_geom_offset=c(.95,1/.95))
xlims = c(1.5, 10)
ylims = c(0,0.1)
xats = 1:10
xbigs = c(2, 4, 8)
yats = seq(0,.1,.01)
ybigs = seq(0,.1,.05)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, xaxs='i', yaxs='i', log='x')
mean_cvs %>%
  inner_join(fluid_meta, by='biofluid') %>%
  mutate(xcent = dilution_factor * x_geom_offset) -> mean_cvs_anno
bar_geom_width = 1/.95
rect(xleft=mean_cvs_anno$xcent*bar_geom_width, xright=mean_cvs_anno$xcent/bar_geom_width,
     ybottom=rep(0, nrow(mean_cvs_anno)), ytop=mean_cvs_anno$mean_cv, col=mean_cvs_anno$color, border=NA)
axis(side=1, at=xats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
axis(side=1, at=xbigs, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-0.5)
mtext(side=1, line=2, text='dilution factor', cex=0.8)
axis(side=2, at=yats, tck=-0.025, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, labels=NA)
axis(side=2, at=ybigs, tck=-0.05, las=2, lwd=0, line=-0.5, labels=percent(ybigs), cex.axis=0.8)
mtext(side=2, line=2.5, text='mean %CV', cex=0.8)
par(xpd=T)
legend(x=5, y=.15, fluid_meta$biofluid, col=fluid_meta$color, pch=15, bty='n', cex=0.8)
par(xpd=F)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5); panel = panel + 1

par(mar=c(3,3,3,1))
xlims = c(0.5, 2.5)
ylims = c(0,600)
plot(NA, NA, xlim=xlims, ylim=ylims, axes=F, ann=F, yaxs='i', xaxs='i')
axis(side=1, at=xlims, lwd.ticks=0, labels=NA)
axis(side=1, at=c(1,2), labels=c('neat','CHAPS'), lwd.ticks=0, cex.axis=0.8)
for (i in 1:nrow(bsyn_chaps)) {
  points(1:2, c(bsyn_chaps$neat[i], bsyn_chaps$chaps[i]), pch=20)
  points(1:2, c(bsyn_chaps$neat[i], bsyn_chaps$chaps[i]), type='l')
}
axis(side=2, at=0:6*100, las=2)
mtext(side=2, line=2.5, text='CSF β-syn (pg/mL)', cex=0.8)
mtext(letters[panel], side=3, cex=2, adj = -0.1, line = 0.5); panel = panel + 1

tobj = t.test(bsyn_chaps$neat, bsyn_chaps$chaps, paired=T)
chaps_ratio = mean(bsyn_chaps$chaps / bsyn_chaps$neat)

write_stats('CSF beta-syn with vs. without CHAPS: mean ',percent(chaps_ratio-1,signed=T,digits=1),', T test P value = ',tobj$p.value)

unnecessary_message = dev.off()




plot_rtq_kinetic = function() {
  
  par(mar=c(3,4,3,1))
  
  xlims = c(0,24)
  xats = 0:24
  xbigs = seq(0,24,6)
  ylims = c(-0.025,1.05)
  yats = (0:4)/4
  
  
  plot(NA, NA, xlim=xlims, ylim=ylims, xaxs='i', yaxs='i', ann=FALSE, axes=FALSE)
  mtext(side=3, line=0.5, text='RT-QuIC')
  axis(side=1, at=xats, labels=NA, lwd=1, lwd.ticks=1, tck=-0.025)
  axis(side=1, at=xbigs, labels=NA, lwd=0, lwd.ticks=1, tck=-0.05)
  axis(side=1, at=xbigs, labels=xbigs, lwd=0, line=-0.5)
  axis(side=2, at=yats, labels=NA, lwd=1)
  axis(side=2, at=yats, labels=percent(yats), las=2, lwd=0, line=-0.5)
  
  mtext(side=1, line=1.5, text='reaction time (h)', cex=0.8)
  mtext(side=2, line=2.5, text='normalized fluorescence', cex=0.8)
  
  rtq_avkin_plot %>%
    mutate(conv_sort = case_when(converter & mut=='P102L' ~ 1, TRUE ~ 0)) %>%
    arrange(converter, conv_sort, indiv, iv) %>%
    distinct(iv) %>% 
    pull() -> ivs_in_order
  
  for (this_iv in ivs_in_order) {
    subs = subset(rtq_avkin_plot, iv == this_iv)
    points(x=subs$hours, y=subs$fluor, type='l', lwd=subs$lwd, col=subs$color)
  }
  
  
}


rtq_avkin_plot = read_tsv('../mgh_freeze2/nosync/analytic/rtq_avkin_plot.tsv', col_types=cols()) %>%
  select(-color, -mut, -cat) %>%
  mutate(indiv = as.integer(substr(iv, 1, 4))) %>%
  inner_join(gt_dem, by='indiv') %>%
  mutate(converter = indiv %in% conv_table$indiv) %>%
  mutate(group = case_when(converter & mut=='E200K' ~ 'E200K converter',
                           converter & mut=='P102L' ~ 'P102L converter',
                           !converter & mut=='none' ~ 'control',
                           !converter & mut!='none' ~ 'non-converting carrier')) %>%
  inner_join(leg, by=c('group'='disp'))


resx=300
png(paste0('for_slides/figure-rtq.png'),width=3.5*resx,height=3.5*resx,res=resx)
plot_rtq_kinetic()
unnecessary_message = dev.off()

####
# Figure 1
####

tell_user('done.\nCreating Figure 1...')

resx=600
png('display_items/figure-1.png', width=6.5*resx, height=8*resx, res=resx)
layout_matrix = matrix(c(rep(1:3, each=4),
                         rep(4:7, each=3),
                         rep(8:11, each=3),
                         rep(12:13,each=3), rep(14, each=6)),
                       nrow=4, 
                       byrow=T)
layout(layout_matrix)

panel = 1

plot_rtq_kinetic()
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5); panel = panel + 1
plot_rtq_endpoint()
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5); panel = panel + 1
plot_csf_prp_tr()
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5); panel = panel + 1
plot_biom_age('plasma_gfap')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5); panel = panel + 1
plot_biom_delta('plasma_gfap')
plot_biom_age('plasma_nfl')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5); panel = panel + 1
plot_biom_delta('plasma_nfl')
plot_biom_age('csf_nfl')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5); panel = panel + 1
plot_biom_delta('csf_nfl')
plot_biom_age('csf_tau')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5); panel = panel + 1
plot_biom_delta('csf_tau')
plot_biom_age('csf_bsyn')
mtext(letters[panel], side=3, cex=2, adj = -0.2, line = 0.5); panel = panel + 1
plot_biom_delta('csf_bsyn')

plot(NA, NA, xlim=0:1, ylim=0:1, xaxs='i', yaxs='i', axes=F, ann=F)
legend('center', legend=leg$disp, col=leg$color, pt.lwd=leg$ptlwd, lwd=leg$lwd, pch=leg$pch, bty='n', cex=1.25)

unnecessary_message = dev.off()



for (i in 1:nrow(bioms_meta)) {
  resx=600
  png(paste0('for_slides/biom_combined_',bioms_meta$varname[i],'.png'), width=6.5*resx, height=3.25*resx, res=resx)
  par(mfrow=c(1,2))
  plot_biom_age(bioms_meta$varname[i])
  plot_biom_delta(bioms_meta$varname[i])
  unnecessary_message = dev.off()
}





#####
# SUPPLEMENT
#####

tell_user('done.\nFinalizing supplementary tables...')

# write the supplement directory / table of contents
supplement_directory %>% rename(table_number = name, description=title) -> contents
addWorksheet(supplement,'contents')
bold_style = createStyle(textDecoration = "Bold")
writeData(supplement,'contents',contents,headerStyle=bold_style,withFilter=T)
freezePane(supplement,'contents',firstRow=T)
# move directory to the front
original_order = worksheetOrder(supplement)
n_sheets = length(original_order)
new_order = c(n_sheets, 1:(n_sheets-1))
worksheetOrder(supplement) = new_order
activeSheet(supplement) = 'contents'
# now save
saveWorkbook(supplement,supplement_path,overwrite = TRUE)

elapsed_time = Sys.time() - overall_start_time
cat(file=stderr(), paste0('done.\nAll tasks complete in ',round(as.numeric(elapsed_time),1),' ',units(elapsed_time),'.\n'))


