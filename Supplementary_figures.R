###############################################################################
### T cell quantification from DNA sequencing predicts immunotherapy        ###
###       response                                                          ###
###                                                                         ###
### Supplementary Information figures                                       ###
###                                                                         ###
### Author: Robert Bentham                                                  ###
### Date: 2021 July 30th                                                    ###
###############################################################################


# SF1 -----
# Toy examples
r_estimate1 <- function(tcell.fraction){
  return(log2(1-tcell.fraction))
}

r_estimate2 <- function(tcell.fraction, local.cN.non.t, local.cN){
  top <- local.cN.non.t*(1-tcell.fraction)
  bottom <- local.cN
  return(log2(top/bottom))
}


# Local tumour copy number = 1 to 6
# Real t-cell fraction = 0.05 to 0.4 (by 0.05)
# tumour purity = 0.75
ylines <- data.frame(intercept = seq(0.05, 0.4, by=0.05),
                     group = factor(seq(0.05, 0.4, by=0.05), levels = seq(0.05, 0.4, by=0.05)))


r.estimate.copy.example <- data.frame(tumour.ploidy = seq(1,6, by=0.01),
                                      tcell.fraction = rep(seq(0.05, 0.4, by=0.05), each = 501),
                                      tumour.purity = 0.75) %>%
  mutate(TCRA.local.copy = 2*(1-tumour.purity) + tumour.purity*tumour.ploidy) %>%
  mutate(adjusted.purity = tumour.purity/(1-tcell.fraction)) %>%
  mutate(TCRA.noT.local.copy = 2*(1-adjusted.purity) + adjusted.purity*tumour.ploidy) %>%
  mutate(copy.ratio = TCRA.local.copy/TCRA.noT.local.copy ) %>%
  mutate(r_original = r_estimate1(tcell.fraction)) %>%
  mutate(r_update = r_estimate2(tcell.fraction, TCRA.noT.local.copy, TCRA.local.copy)) %>%
  mutate(t_est = 1-(2^(r_update))) %>%
  mutate(t_update = ((1 - tumour.purity) + (tumour.purity*tumour.ploidy)/2) - 
           ((1 - tumour.purity) + (tumour.purity*tumour.ploidy)/2)*2^(r_update))


pdf('plots/supfigs/est_vs_cn_75.pdf', height = 6, width = 6)
r.estimate.copy.example %>%
  ggplot(aes(tumour.ploidy, t_est)) + 
  geom_hline(data = ylines, aes(yintercept = intercept, col = group)) +
  geom_point(aes(col = factor(tcell.fraction, levels = seq(0.05, 0.4, by=0.05)))) + 
  labs(col = "real tcell fraction") + ylab('estimated tcell fraction') + 
  xlab('local tumour copynumber') + ggtitle('Tumour purity = 0.75') +
  theme_bw()
dev.off()

r.estimate.copy.example <- data.frame(tumour.ploidy = seq(1,6, by=0.01),
                                      tcell.fraction = rep(seq(0.05, 0.4, by=0.05), each = 501),
                                      tumour.purity = 0.5) %>%
  mutate(TCRA.local.copy = 2*(1-tumour.purity) + tumour.purity*tumour.ploidy) %>%
  mutate(adjusted.purity = tumour.purity/(1-tcell.fraction)) %>%
  mutate(TCRA.noT.local.copy = 2*(1-adjusted.purity) + adjusted.purity*tumour.ploidy) %>%
  mutate(copy.ratio = TCRA.local.copy/TCRA.noT.local.copy ) %>%
  mutate(r_original = r_estimate1(tcell.fraction)) %>%
  mutate(r_update = r_estimate2(tcell.fraction, TCRA.noT.local.copy, TCRA.local.copy)) %>%
  mutate(t_est = 1-(2^(r_update))) %>%
  mutate(t_update = ((1 - tumour.purity) + (tumour.purity*tumour.ploidy)/2) - 
           ((1 - tumour.purity) + (tumour.purity*tumour.ploidy)/2)*2^(r_update))

pdf('plots/supfigs/est_vs_cn_5.pdf', height = 6, width = 6)
r.estimate.copy.example %>%
  ggplot(aes(tumour.ploidy, t_est)) + 
  geom_hline(data = ylines, aes(yintercept = intercept, col = group)) +
  geom_point(aes(col = factor(tcell.fraction, levels = seq(0.05, 0.4, by=0.05)))) + 
  labs(col = "real tcell fraction") + ylab('estimated tcell fraction') + 
  xlab('local tumour copynumber') + ggtitle('Tumour purity = 0.5') +
  theme_bw()
dev.off()

r.estimate.copy.example <- data.frame(tumour.ploidy = seq(1,6, by=0.01),
                                      tcell.fraction = rep(seq(0.05, 0.4, by=0.05), each = 501),
                                      tumour.purity = 0.25) %>%
  mutate(TCRA.local.copy = 2*(1-tumour.purity) + tumour.purity*tumour.ploidy) %>%
  mutate(adjusted.purity = tumour.purity/(1-tcell.fraction)) %>%
  mutate(TCRA.noT.local.copy = 2*(1-adjusted.purity) + adjusted.purity*tumour.ploidy) %>%
  mutate(copy.ratio = TCRA.local.copy/TCRA.noT.local.copy ) %>%
  mutate(r_original = r_estimate1(tcell.fraction)) %>%
  mutate(r_update = r_estimate2(tcell.fraction, TCRA.noT.local.copy, TCRA.local.copy)) %>%
  mutate(t_est = 1-(2^(r_update))) %>%
  mutate(t_update = ((1 - tumour.purity) + (tumour.purity*tumour.ploidy)/2) - 
           ((1 - tumour.purity) + (tumour.purity*tumour.ploidy)/2)*2^(r_update))

pdf('plots/supfigs/est_vs_cn_25.pdf', height = 6, width = 6)
r.estimate.copy.example %>%
  ggplot(aes(tumour.ploidy, t_est)) + 
  geom_hline(data = ylines, aes(yintercept = intercept, col = group)) +
  geom_point(aes(col = factor(tcell.fraction, levels = seq(0.05, 0.4, by=0.05)))) + 
  labs(col = "real tcell fraction") + ylab('estimated tcell fraction') + 
  xlab('local tumour copynumber') + ggtitle('Tumour purity = 0.25') +
  theme_bw()
dev.off()

# Get distributions of local tumour copy number and real tcell fraction

tx100.tcra.scores <- readRDS('data/Tx100_TCRA_scores_final_20201208.RDS')

pdf('plots/supfigs/tracerx100_copy.pdf', height = 4, width = 4)
tx100.tcra.scores %>%
  ggplot(aes(cnTotal)) + 
  geom_bar() + xlab('Local TCRA copy number') + ggtitle('TRACERx100')
dev.off()

pdf('plots/supfigs//tracerx100_purity.pdf', height = 4, width = 4)
tx100.tcra.scores %>% filter(purity!= 0) %>%
  ggplot(aes(purity)) + 
  geom_histogram() + xlab('Tumour purity') + ggtitle('TRACERx100')
dev.off()

pdf('plots/supfigs/tracerx100_naive_exact_copy.pdf', height = 6, width = 6)
tx100.tcra.scores %>%
  filter(cnTotal < 6 & cnTotal != 2) %>%
  ggplot(aes(TCRA.tcell.fraction, TCRA.naive.score)) +
  geom_point() +
  xlab('Exact t-cell fraction') + ylab('Naive t-cell fraction') +
  facet_wrap(~cnTotal) + 
  geom_abline(intercept = 0, slope = 1, colour = 'red', linetype = 'dashed') +
  ggtitle('TRACERx100 exact vs naive exact score by local copy number')
dev.off()


# Sup Figure 2 -----
library(tidyverse)

TCRA.opt.df <-readRDS('data/tracerx100_local_seg_opt.RData')


# Measures of comparison

pdf('plots/supfigs/scores_vs_segments.pdf', width = 16, height = 4)
TCRA.opt.df %>%
  gather('V_segments','Score',-sample) %>%
  mutate(V_segments = gsub('.cnCorrect.score','',V_segments)) %>%
  mutate(V_segments = gsub('TCRA.','',V_segments)) %>%
  mutate(V_segments = factor(V_segments, levels = paste0('V',seq(55)))) %>%
  ggplot(aes(V_segments, Score)) + #geom_point(aes(col = sample)) +
  geom_line(aes(col = sample, group = sample), alpha = 0.25) + 
  ggtitle('Scores vs segments used for normalisation') +
  xlab('') + ylab('Score') + theme_bw() + 
  theme(legend.position = 'none')
dev.off()

# 2. %above 0
pdf('plots/supfigs/nonzero_vs_segments.pdf', width = 16, height = 4)
data.frame(V_segments = factor(paste0('V',seq(55)), levels = paste0('V',seq(55))),
           percent_nonzero = as.numeric(colSums(TCRA.opt.df[,-1] > 0)/(dim(TCRA.opt.df)[1])),
           sample = 1) %>%
  ggplot(aes(V_segments, percent_nonzero)) + geom_line(aes(group = sample)) +
  xlab('') + ylab('Fraction') + ylim(0,1) + ggtitle('Fraction of non-zero samples') + 
  theme_bw() + 
  theme(legend.position = 'none') 
dev.off()
# 3. Correlation with tcell Danaher component
load('data/tracerx100.danaher.RData')

til_score_tx100_df_TCRA_new_plot <- list()
minus_log2_pvalue <- seq(55)

for(i in 1:55){
  col_to_select <- sym(colnames(TCRA.opt.df)[i+1])
  tx100.tcra.scores <- TCRA.opt.df %>%
    dplyr::select(sample, !!col_to_select)
  tx100.regions <- tx100.tcra.scores$sample
  
  
  til_score_tx100_df_TCRA <- til_score_tx100_df %>%
    rownames_to_column(var = 'sample') %>%
    left_join(tx100.tcra.scores) %>%
    filter(sample %in% tx100.regions) 
  
  til_score_tx100_df_TCRA_new_plot[[i]] <- til_score_tx100_df_TCRA %>%
    gather(Cell, Score, -sample, -!!col_to_select) %>%
    filter(Cell == 'tcells') %>%
    ggplot(aes(x=!!col_to_select, y=Score)) +
    facet_wrap(~Cell) +
    geom_point(aes(col = Cell), alpha = 0.5) +
    geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
    xlab('TCRA score') + 
    stat_cor(label.y = 12) +
    stat_regline_equation(label.y = 10) + 
    ylim(0,13) + theme_bw() 
  
  minus_log2_pvalue[i] <- -log2(ggplot_build(til_score_tx100_df_TCRA_new_plot[[i]])$data[[3]]$p)
}

pdf('plots/supfigs/danaher_vs_segments.pdf', width = 16, height = 4)
data.frame(V_segments = factor(paste0('V',seq(55)), levels = paste0('V',seq(55))),
           minus_log2 = minus_log2_pvalue,
           sample = 1) %>%
  ggplot(aes(V_segments, minus_log2)) + geom_line(aes(group = sample)) +
  xlab('') + ylab('-log2 pvalue')  + ggtitle('Correlation with Danaher tcells') + 
  theme_bw() + 
  theme(legend.position = 'none') 
dev.off()


# Figure S5A - Example of coverage bias ----


vdj.logR.df <- readRDS('data/tcga_exon_bias_example.RDS')

pdf('plots/supfigs/tcga_logR_example.pdf', height = 6, width = 12)
vdj.logR.df %>%
  ggplot(aes(pos, Ratio)) +
  geom_point() +
  geom_smooth() + 
  xlab('chr14 position') + ylab('log ratio') + ggtitle('TCGA LUAD example')
dev.off()


# Figure S5B - Identification of exons with coverage issues ----

tcga.luad.exon.ratio.df <- readRDS('data/TCGA_luad_exonRatio.RDS')
colnames(tcga.luad.exon.ratio.df) <- paste0( seq(192))

tcga.luad.exon.ratio.df.long <- gather(as.data.frame(tcga.luad.exon.ratio.df),
                                       key = 'exon', value = 'ratio')


median.tcga.luad.scores <- tcga.luad.exon.ratio.df.long  %>%
  dplyr::select(exon, ratio) %>%
  mutate(exon = factor(exon, levels = seq(192))) %>%
  group_by(exon) %>%
  summarise(TCGA.Median = median(ratio, na.rm = TRUE))


# Load Tx100 exons

tx100.exon.ratio.df <- readRDS('data/Tx100_meanExonRatio.RDS')

colnames(tx100.exon.ratio.df) <- paste0( seq(192))

tx100.exon.ratio.df.long <- gather(as.data.frame(tx100.exon.ratio.df),
                                   key = 'exon', value = 'ratio')


# Combine!
tx100.exon.ratio.df.long <- tx100.exon.ratio.df.long %>%
  mutate(samp1 = 'tx100')
tcga.luad.exon.ratio.df.long <-  tcga.luad.exon.ratio.df.long %>%
  mutate(samp1 = 'tcga.luad')

pdf('plots/supfigs/tcga_exon_ratio.pdf', height = 6, width = 15)
rbind(tx100.exon.ratio.df.long , tcga.luad.exon.ratio.df.long) %>% 
  filter(exon %in% as.character(c(121:192))) %>%
  mutate(exon = factor(exon, levels = seq(192))) %>%
  ggplot(aes(exon, ratio)) + geom_boxplot(aes(fill = samp1)) +
  ylab('log ratio') + xlab('TCRA exon')
dev.off()

# Get the locations to remove
median.tx100.scores <- tx100.exon.ratio.df.long %>%
  dplyr::select(exon, ratio) %>%
  mutate(exon = factor(exon, levels = seq(192))) %>%
  group_by(exon) %>%
  summarise(Tx421.Median = median(ratio, na.rm = TRUE))

median.exon.scores  <- median.tx100.scores  %>%
  left_join(median.tcga.luad.scores, by = 'exon') %>%
  mutate(medianDiff = TCGA.Median - Tx421.Median) %>%
  mutate(diffBig = (medianDiff) < -0.5)

remove.exons.0.5 <- as.numeric(as.character(median.exon.scores$exon[which(median.exon.scores$diffBig)]))


tx100.no.exons <- readRDS('data/tx100_TCRA_scores.exonsRemoved.RDS')
tx100.tcra.scores <- readRDS('data/Tx100_TCRA_scores_final_20201208.RDS')

# Join these two:
tx100.exon.score.df  <- tx100.tcra.scores %>%
  rename(All.exons = TCRA.orig.score) %>%
  select(sample, All.exons) %>%
  left_join(tx100.no.exons, 'sample') %>%
  rename(Remove.exons = TCRA.orig.score) %>%
  mutate(score.change = Remove.exons - All.exons)



# TCGA and Tx100 GC vs no GC
# load GC corrected scores:
# Older version of scores but these have not changed

TCGA.scores.final.gc <- readRDS('data/TCGA_TCRA_scores_final_gc_20210122.RDS')

Tx100.gc.final <-  readRDS('data/Tx100_TCRA_scores_gc_update_20210122.RDS')



load('data/tracerx100.clinical.RData')
tx100.hist.df <- clinical.data.tx100 %>%
  select(TRACERxID, histology = Histology) %>%
  mutate(histology = ifelse(histology == 'Invasive adenocarcinoma','LUAD',
                            ifelse(histology == 'Squamous cell carcinoma', 'LUSC', 'Other'))) 

# Use 

tcga.gc.df <- TCGA.scores.final.gc %>%
  select(sample = barcode, germline, TCRA.tcell.fraction.gc, TCRA.tcell.fraction,histology = tcga) %>%
  mutate(cohort = 'TCGA') %>% filter(!is.na(TCRA.tcell.fraction)) %>%
  mutate(histology = toupper(histology))

tx100.gc.df <- Tx100.gc.final %>%
  select(sample, germline, TCRA.tcell.fraction.gc, TCRA.tcell.fraction) %>%
  mutate(cohort = 'TRACERx') %>%
  mutate(TRACERxID = substr(sample, 3,8)) %>%
  left_join(tx100.hist.df, 'TRACERxID') %>%
  select(sample, germline, TCRA.tcell.fraction.gc, TCRA.tcell.fraction, histology, cohort)

all.gc.df <- rbind(tcga.gc.df, tx100.gc.df)


pre.gc.cor.df <- all.gc.df %>%
  filter(!germline) %>%
  filter(histology %in% c('LUAD','LUSC')) %>%
  filter(!is.na(TCRA.tcell.fraction.gc)) %>%   filter(!is.na(TCRA.tcell.fraction))

wilcox_test(pre.gc.cor.df[pre.gc.cor.df$histology=='LUAD',], TCRA.tcell.fraction ~ cohort)
wilcox_test(pre.gc.cor.df[pre.gc.cor.df$histology=='LUSC',], TCRA.tcell.fraction ~ cohort)

# GC correction:
pdf('plots/supfigs/TCGA_TX100_preGC_correct.pdf', width = 3, height = 2.5)
all.gc.df %>%
  filter(!germline) %>%
  filter(histology %in% c('LUAD','LUSC')) %>%
  filter(!is.na(TCRA.tcell.fraction.gc)) %>%   filter(!is.na(TCRA.tcell.fraction)) %>%
  ggplot(aes(cohort,  TCRA.tcell.fraction)) +
  geom_boxplot(outlier.size = 0.25) + 
  facet_wrap(~histology) + xlab('') + ylab('TCRA T cell fraction') +
  theme_bw() + stat_compare_means(size = 2) +
  theme(text = element_text(size = 7))
dev.off()

pdf('plots/supfigs/TCGA_TX100_postGC_correct.pdf', width = 3, height = 2.5)
all.gc.df %>%
  filter(!germline) %>%
  filter(histology %in% c('LUAD','LUSC')) %>%
  filter(!is.na(TCRA.tcell.fraction.gc)) %>%   filter(!is.na(TCRA.tcell.fraction)) %>%
  ggplot(aes(cohort, TCRA.tcell.fraction.gc)) +
  geom_boxplot(outlier.size = 0.25) + 
  facet_wrap(~histology) + xlab('') + ylab('TCRA T cell fraction') +
  theme_bw() + stat_compare_means(size = 2) +
  theme(text = element_text(size = 7))
dev.off()