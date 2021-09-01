###############################################################################
### T cell quantification from DNA sequencing predicts immunotherapy        ###
###       response                                                          ###
###                                                                         ###
### paper figure 1 + related extended figures                               ###
###                                                                         ###
### Author: Robert Bentham                                                  ###
### Date: 2021 July 30th                                                    ###
###############################################################################

library(tidyverse)
library(RColorBrewer)
library(readxl)
library(ggpubr)

# Extended Data Figure 1a ----

# TRACERx100 Breakpoints
load('data/tracerx100.ascat.chr14.RData')


known.end.positions <- tracerx.ascat.seg %>%
  group_by(sample) %>% summarise(last.pos = max(endpos)) %>%
  select(last.pos) %>% `[[`(1) %>% unique()

breakpoints.chr14 <- setdiff(tracerx.ascat.seg$endpos,  known.end.positions)

p.breakpoint1 <- data.frame(breakpoints = breakpoints.chr14) %>%
  mutate(TCRA = ifelse(breakpoints >= 22090057 & breakpoints <= 23221076, TRUE, FALSE)) %>%
  filter(breakpoints <= 1e8) %>%
  ggplot() + geom_histogram(aes(breakpoints, fill = TCRA), bins = 400) + 
  xlab('Genome position') +
  theme_classic() + theme(text = element_text(size=6),
                          legend.title = element_text( size = 6),
                          legend.text = element_text(size = 6)) +
  ggtitle('Breakpoints on chr14q and TCRA')

p.breakpoint2 <- data.frame(breakpoints = breakpoints.chr14) %>%
  filter(breakpoints > 22000000 & breakpoints < 23500000) %>%
  mutate(TCRA = ifelse(breakpoints >= 22090057 & breakpoints <= 23221076, TRUE, FALSE)) %>%
  ggplot() + geom_histogram(aes(breakpoints, fill = TCRA), bins = 80) + + 
  xlab('') + ylab('') +
  theme_classic() + theme(legend.position = 'none', text = element_text(size=6)) 

# Extended Figure 1a
pdf('plots/ext1/chr14q_breakpoints.pdf', width = 3.15, height =1.575)
p.breakpoint1 + annotation_custom(ggplotGrob(p.breakpoint2), ymin = 15, ymax = 55, 
                                  xmin = 3e7, xmax = 105e6)
dev.off()


# Extended figures 1b/c ----

TCRA_segments <- read_excel("data/TCRA_segments.xlsx")
TCRA.size <- 23021076 - 22090057

# To get coverage from bam files:
# samtools depth bam1 -r chr14:22090057-23021075 > CRUK0023_GL.txt
# samtools depth bam2 -r chr14:22090057-23021075 > CRUK0023_R1.txt
# and normalise for read depth

test.all2 <- readRDS('data/CRUK0023_example.RDS')

TRA_segment_find <- function(x){
  gene.loc <- setdiff(intersect(which(TCRA_segments$start_hg19 <= x),
                                which(x <= TCRA_segments$end_hg19)), c(1,57))
  return(ifelse(length(gene.loc) > 0, TCRA_segments$Gene[gene.loc], 'None'))
}

TRA_segment_find_v <- Vectorize(TRA_segment_find)

test.all3 <- test.all2 %>%
  mutate(TRA_segment = TRA_segment_find_v(X2))

test.all3$TRA_segment_short <- gsub('-','',gsub('[0-9]*','',gsub('TR[A-Z]','',test.all3$TRA_segment)))

test.all3a <- test.all3 %>%
  filter(TRA_segment_short != 'None') 
test.all3b <- test.all3 %>%
  filter(TRA_segment_short == 'None') 

pdf('plots/ext1/CRUK0023_R1_read_depth_ratio_with_leg.pdf', width = 3.5, height = 2.5)
test.all3a  %>%
  ggplot() + 
   geom_point(aes(x = X2, y=Ratio, col = TRA_segment_short), size = 0.25) + geom_smooth(aes(x = X2, y = Ratio)) +
  ggtitle('') + 
  geom_vline(aes(xintercept = 22090057), col = 'red', linetype = 'dashed') + 
  geom_vline(aes(xintercept = 23021075), col = 'red', linetype = 'dashed') +
  xlab('Chr 14 Location') +  labs(col = "TCRA gene segment") + theme_bw() +
  theme(text = element_text(size = 7)) 
dev.off()

# Reduced point version for easier display
pdf('plots/ext1/CRUK0023_R1_read_depth_ratio_reduced.pdf', width = 3.5, height = 2.5)
(test.all3a[sample(seq(30728 ),3000),]) %>%
  ggplot() + 
  # geom_point(data = test.all3b,
  #            aes(x = X2, y=Ratio, col = TRA_segment_short),alpha = 0.01) +
  geom_point(aes(x = X2, y=Ratio, col = TRA_segment_short), size = 0.25) + geom_smooth(aes(x = X2, y = Ratio)) +
  ggtitle('') + 
  geom_vline(aes(xintercept = 22090057), col = 'red', linetype = 'dashed') + 
  geom_vline(aes(xintercept = 23021075), col = 'red', linetype = 'dashed') +
  xlab('Chr 14') +  labs(col = "TCRA gene segment") + theme_bw() +
  theme(text = element_text(size = 7), legend.position = 'none') 
dev.off()

test.all2 <- readRDS('data/CRUK0039_example.RDS')

TRA_segment_find <- function(x){
  gene.loc <- setdiff(intersect(which(TCRA_segments$start_hg19 <= x),
                                which(x <= TCRA_segments$end_hg19)), c(1,57))
  return(ifelse(length(gene.loc) > 0, TCRA_segments$Gene[gene.loc], 'None'))
}

TRA_segment_find_v <- Vectorize(TRA_segment_find)

test.all3 <- test.all2 %>%
  mutate(TRA_segment = TRA_segment_find_v(X2))

test.all3$TRA_segment_short <- gsub('-','',gsub('[0-9]*','',gsub('TR[A-Z]','',test.all3$TRA_segment)))

test.all3a <- test.all3 %>%
  filter(TRA_segment_short != 'None') 
test.all3b <- test.all3 %>%
  filter(TRA_segment_short == 'None') 

pdf('plots/ext1/CRUK0039_R1_read_depth_ratio_reduced.pdf', width = 3.5, height = 2.5)
test.all3a[sample(seq(34100 ),3000),]  %>%
  ggplot() + 
  geom_point(aes(x = X2, y=Ratio, col = TRA_segment_short),size = 0.25) + geom_smooth(aes(x = X2, y = Ratio)) +
  ggtitle('') + 
  geom_vline(aes(xintercept = 22090057), col = 'red', linetype = 'dashed') + 
  geom_vline(aes(xintercept = 23021075), col = 'red', linetype = 'dashed') +
  xlab('Chr 14 Location') +  labs(col = "TCRA gene segment") + theme_bw()  +
  theme(text = element_text(size = 7),  legend.position = 'none') 
dev.off()


# Plot to extract legend for figure 
pdf('plots/supfig1/legend_extract.pdf', width = 3.5, height = 2.5)
test.all3a  %>%
  filter((X2 %% 50 == 1)) %>%
  ggplot() + 
  #  geom_point(data = test.all3b,
  #             aes(x = X2, y=Ratio, col = TRA_segment_short),alpha = 0.01) +
  geom_point(aes(x = X2, y=Ratio, col = TRA_segment_short),size = 0.25) + geom_smooth(aes(x = X2, y = Ratio)) +
  ggtitle('Read depth ratio of CRUK0039 R1 to matched germline within TCRA') + 
  geom_vline(aes(xintercept = 22090057), col = 'red', linetype = 'dashed') + 
  geom_vline(aes(xintercept = 23021075), col = 'red', linetype = 'dashed') +
  xlab('Chr 14 Location') +  labs(col = "TCRA gene segment") + theme_bw()  +
  theme(text = element_text(size = 7)) 
dev.off()


# Extended Figure 1e FFPE analysis ----

metrics_output <- readRDS('data/CPI_merged_data_orig_gc_20210517.RDS')
metrics_output <- metrics_output  %>% 
  rename(histology = Tum_Type)

ffpe.cpi.df <- metrics_output %>%
  filter(study %in% c('CRISTESCU_SCIENCE_2018','RIAZ_CELL_2017','VANALLEN_SCIENCE_2015',
                      'MARIATHASAN_NATURE_2018','MCDERMOT_NMED_2018','SNYDER_NEJM_2014')) %>%
  mutate(ffpe.status = ifelse(study %in% c('CRISTESCU_SCIENCE_2018','RIAZ_CELL_2017','VANALLEN_SCIENCE_2015'),
                              'FFPE','FreshFrozen')) %>%
  select(case, study, histology, ffpe.status, 
         tumour.tcell.fraction, normal.tcell.fraction) %>%
  filter(!is.na(tumour.tcell.fraction))

pdf('plots/ext1/ffpe_plot1.pdf', width = 2.5, height = 2.5)
ffpe.cpi.df %>%
  filter(histology %in% c('MELANOMA','BLADDER')) %>%
  ggplot(aes(ffpe.status, tumour.tcell.fraction)) +
  geom_boxplot(outlier.size = 0,outlier.alpha = 0) + 
  geom_sina(aes(colour = study), size = 0.25) +
  stat_compare_means(size = 3) + theme_bw() +
  ylab('TCRA T cell fraction') + xlab('') +
  theme(text = element_text(size = 7), legend.position = 'none') +
  facet_wrap(~histology)
dev.off()


# Extended Figure 1f Cell Lines ----
cellLine.df <- readRDS('data/cell_line_TCRA_scores_final_20210525.RDS')


blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

cellLine.pie <- cellLine.df  %>%
  select(sample, cellLineType, TCRA.TcellFraction = T.cell.fraction.raw.gc.smt) %>%
  mutate(nonTcellFraction = 1- TCRA.TcellFraction) 

cellLine.pie2 <- cellLine.pie %>%
  gather('Type', 'Fraction', -sample, -cellLineType) %>%
  mutate(Fraction = ifelse(Fraction < 0,0,Fraction))

cellLine.pie2$sample <- c(paste0('HCT116-R',seq(15)),
                          'JURKAT','PEER','HPB-ALL')

pdf('plots/ext1/cell_line.pdf', height = 2.5, width = 3.7)
cellLine.pie2 %>% 
  mutate(Type = ifelse(Type == 'TCRA.TcellFraction','T cell','non-T cell')) %>%
  ggplot() +
  geom_col(aes(x='', y = Fraction, fill = Type), position = 'fill') + 
  coord_polar("y", start=0) + facet_wrap(~sample) +  blank_theme +
  scale_fill_manual(values = c('#D89C70','#A6ACD7')) + 
  theme(axis.text.x=element_blank(), text = element_text(size = 6), legend.key.size = unit(0.2, 'cm'))
dev.off()

# Extended Figure 1h iDNA ----

load('data/tracerx100.clinical.RData')
load('data/tracerx100.danaher.RData')
tx100.tcra.scores <- readRDS('data/Tx100_TCRA_scores_gc_update_20210630.RDS')

tx100.tcra.scores <- tx100.tcra.scores %>%
  dplyr::select(sample, TCRA.tcell.fraction.gc)


iDNA.df2 <- readRDS('data/iDNA_summary_tx100.RDS')

pdf('plots/ext1/iDNA.pdf', height = 2.3, width = 2.3)
iDNA.df2 %>% 
  left_join(tx100.tcra.scores, 'sample') %>%
  ggplot(aes(iDNA.score, TCRA.tcell.fraction.gc)) +
  geom_point(size = 0.25) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size = 2) + theme_bw() + 
  xlab('CDR3 VDJ score') + ylab('TCRA T cell fraction') + 
  geom_smooth(method = 'lm') + 
  theme(text = element_text(size = 6)) + ggtitle('TRACERx100 cohort')
dev.off()

# Figure 1B TILs ----
load('data/tracerx.TIL.RData')

all.scores.tx100.til <- tx100.tcra.scores %>%
  left_join(Tx.til.long, by = 'sample')

hist.df <- clinical.data.tx100 %>%
  select(patient = TRACERxID, Histology) %>%
  mutate(Histology = ifelse(!Histology %in% c('Squamous cell carcinoma','Invasive adenocarcinoma'),
                            'Other',Histology)) 


til_score_tx100_df_TCRA <- til_score_tx100_df %>%
  rownames_to_column(var = 'sample') %>%
  left_join(tx100.tcra.scores)

cd8.danaher.df <- til_score_tx100_df_TCRA %>%
  dplyr::select(sample, cd8)

rna_seq_samps <- rownames(til_score_tx100_df)

til.tcra <- all.scores.tx100.til %>%
  filter(!is.na(Til)) %>%
  filter(sample %in% rna_seq_samps) %>% 
  left_join(hist.df,'patient') %>%
  ggplot(aes(Til, TCRA.tcell.fraction.gc)) +
  geom_point(aes(colour = Histology)) +
  geom_smooth(method = 'lm') +
  stat_cor(label.y.npc  = 1, method = 'spearman', cor.coef.name = 'rho', size = 2) +
  # stat_regline_equation(label.y.npc = 0.8) + 
  xlab('Pathological Til score') + ylab('TCRA t-cell score') +
  theme_bw()+ theme(legend.position = 'none', text = element_text(size=4))


til.danaher <- all.scores.tx100.til %>%
  filter(!is.na(Til)) %>%
  filter(sample %in% rna_seq_samps) %>% 
  left_join(hist.df,'patient') %>%
  left_join(cd8.danaher.df, 'sample') %>%
  ggplot(aes(Til, cd8)) +
  geom_point(aes(colour = Histology)) +
  geom_smooth(method = 'lm') +
  stat_cor(label.y.npc  = 1, method = 'spearman', cor.coef.name = 'rho', size = 2) +
  # stat_regline_equation(label.y.npc = 0.8) + 
  xlab('Pathological Til score') + ylab('Danaher CD8+') +
  theme_bw() + theme(legend.position = 'none', text = element_text(size=4))


# Others - davoli, xcell etc.
load('data/tracerx_immune_scores.RData')

# Make plots:
rnaseq.samps <- timer.cd8$sample

til.davoli <- Tx.til.long %>%
  filter(!is.na(Til)) %>%
  filter(sample %in% rnaseq.samps) %>%
  left_join(hist.df,'patient') %>%
  left_join(davoli.cd8, 'sample') %>%
  ggplot(aes(Til, cd8)) +
  geom_point(aes(colour = Histology)) +
  geom_smooth(method = 'lm') +
  stat_cor(label.y.npc  = 1, method = 'spearman', cor.coef.name = 'rho', size = 2) +
  xlab('Pathological Til score') + ylab('Davoli CD8+') +
  theme_bw() + theme(legend.position = 'none', text = element_text(size=4))


til.epic <- Tx.til.long %>%
  filter(!is.na(Til)) %>%
  filter(sample %in% rnaseq.samps) %>%
  left_join(hist.df,'patient') %>%
  left_join(epic.cd8, 'sample') %>%
  ggplot(aes(Til, cd8)) +
  geom_point(aes(colour = Histology)) +
  geom_smooth(method = 'lm') +
  stat_cor(label.y.npc  = 1, method = 'spearman', cor.coef.name = 'rho', size = 2) +
  xlab('Pathological Til score') + ylab('EPIC CD8+') +
  theme_bw() + theme(legend.position = 'none',  text = element_text(size=4))


til.xcell <- Tx.til.long %>%
  filter(!is.na(Til)) %>%
  filter(sample %in% rnaseq.samps) %>%
  left_join(hist.df,'patient') %>%
  left_join(xcell.cd8, 'sample') %>%
  ggplot(aes(Til, cd8)) +
  geom_point(aes(colour = Histology)) +
  geom_smooth(method = 'lm') +
  stat_cor(label.y.npc  = 1, method = 'spearman', cor.coef.name = 'rho',size = 2) +
  xlab('Pathological Til score') + ylab('xCell CD8+') +
  theme_bw() + theme(legend.position = 'none', text = element_text(size=4))

til.timer <- Tx.til.long %>%
  filter(!is.na(Til)) %>%
  filter(sample %in% rnaseq.samps) %>%
  left_join(hist.df,'patient') %>%
  left_join(timer.cd8, 'sample') %>%
  ggplot(aes(Til, cd8)) +
  geom_point(aes(colour = Histology)) +
  geom_smooth(method = 'lm') +
  stat_cor(label.y.npc  = 1, method = 'spearman', cor.coef.name = 'rho', size = 2) +
  # stat_regline_equation(label.y.npc = 0.8) + 
  xlab('Pathological Til score') + ylab('TIMER CD8+') +
  theme_bw() + theme(legend.position = 'none', text = element_text(size=4))

iDNA.df3 <- iDNA.df2



til.iDNA <- Tx.til.long %>%
  filter(!is.na(Til)) %>%
  filter(sample %in% rnaseq.samps) %>%
  left_join(hist.df,'patient') %>%
  left_join(iDNA.df3, 'sample') %>%
  ggplot(aes(Til, iDNA.score)) +
  geom_point(aes(colour = Histology)) +
  geom_smooth(method = 'lm') +
  stat_cor(label.y.npc  = 1, method = 'spearman', cor.coef.name = 'rho', size = 2) +
  # stat_regline_equation(label.y.npc = 0.8) + 
  xlab('Pathological Til score') + ylab('iDNA score') +
  theme_bw() + theme(legend.position = 'none', text = element_text(size=4))

# CIBERSORT

til.ciber <- Tx.til.long %>%
  filter(!is.na(Til)) %>%
  filter(sample %in% rnaseq.samps) %>%
  left_join(hist.df,'patient') %>%
  mutate(sample = gsub('_SU_T1-',':',sample)) %>%
  left_join(ciber.cd8, 'sample') %>%
  ggplot(aes(Til, CIBER_T_cells_CD8)) +
  geom_point(aes(colour = Histology)) +
  geom_smooth(method = 'lm') +
  stat_cor(label.y.npc  = 1, method = 'spearman', cor.coef.name = 'rho', size = 2) +
  xlab('Pathological Til score') + ylab('CIBER_T_cells_CD8') +
  theme_bw() + theme(legend.position = 'none', text = element_text(size=4))


til.list <- list()
til.list[['TCRA']] <- til.tcra
til.list[['danaher']] <- til.danaher
til.list[['davoli']] <- til.davoli
til.list[['epic']] <- til.epic
til.list[['xcell']] <- til.xcell
til.list[['timer']] <- til.timer
til.list[['iDNA']] <- til.iDNA
til.list[['ciber']] <- til.ciber

# Get rho values out
rho.out <- sapply(til.list, function(x) ggplot_build(x)$data[[3]]$r)
p.value.out  <- sapply(til.list, function(x) ggplot_build(x)$data[[3]]$p)

TIL.rho.df <- data.frame(tcell.name = c('TCRA','Danaher CD8+','Davoli CD8+','EPIC CD8+','xCell CD8+','TIMER CD8+', 'iDNA', 'CIBERSORT CD8+'),
                         rho = rho.out, p.value = p.value.out) %>%
  mutate(sig.level = ifelse(p.value > 0.05, 'ns',
                            ifelse( p.value > 0.01, '*',
                                    ifelse( p.value > 0.001, '**',
                                            ifelse(p.value > 0.0001,'***','****')))))
row.names(TIL.rho.df) <- NULL

library(scales)
TIL.rho.df <- TIL.rho.df %>%
  mutate(tcell.name = factor(tcell.name, levels = rev(c('Danaher CD8+','TCRA','Davoli CD8+','xCell CD8+', 'CIBERSORT CD8+', 'iDNA','TIMER CD8+','EPIC CD8+'))))

pdf('plots/fig1/Til_plot.pdf', width = 2.2, height = 3)
TIL.rho.df %>%
  ggplot() +
  geom_tile(aes(x=tcell.name, y=1, fill = rho)) + 
  xlab('') + ylab('') + theme_bw() + 
  scale_fill_gradient2(high = 'red',low = muted('blue')) +
  annotate("text", x = TIL.rho.df$tcell.name, y = 1, label = TIL.rho.df$sig.level, size = 2) +
  annotate("text", x = as.numeric(TIL.rho.df$tcell.name) - 0.15, y = 1,
           label = paste0('p = ',TIL.rho.df$p.value), size =1.75) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), 
        axis.ticks.y = element_blank(),
        panel.border = element_blank(),panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),text = element_text(size=7), legend.key.size = unit(0.25,'cm')) + coord_flip() +
  ggtitle('') + xlab('T cell estimation method')
dev.off()



# Figure 1C RNASeq ----
load('data/tracerx100.clinical.RData')
load('data/tracerx100.danaher.RData')
tx100.tcra.scores <- readRDS('data/Tx100_TCRA_scores_gc_update_20210630.RDS')

tx100.tcra.scores <- tx100.tcra.scores %>%
  dplyr::select(sample, TCRA.tcell.fraction.gc)

tx100.regions <- tx100.tcra.scores$sample

til_score_tx100_df_TCRA <- til_score_tx100_df %>%
  rownames_to_column(var = 'sample') %>%
  left_join(tx100.tcra.scores)

# Plots of three scores:
til_score_tx100_df_TCRA_new_plot <- til_score_tx100_df_TCRA %>%
  gather(Cell, Score, -sample, -TCRA.tcell.fraction.gc)

load('data/immune.subsets.combined.RData')

tx100.tcra.scores.immune <- tx100.tcra.scores %>%
  mutate(sample = gsub('_SU_T[0-9]-',':',sample)) %>%
  mutate(sample = gsub('_SU_',':',sample)) %>% 
  mutate(sample = gsub('_BS_',':',sample))

# Remove some of the less relevant columnes, e.g. residuals
cols.to.remove1 <- c('xCell_ImmuneScore', 'xCell_StromaScore','xCell_MicroenvironmentScore',
                     'cd8.to.treg.ratio','m1.to.m2.ratio','pro.to.anti.ratio',
                     'roh.immune.score','apc.score','hla.score',
                     'xCell_Smooth_muscle','xCell_Neurons','xCell_Skeletal_muscle',
                     'xCell_Astrocytes','xCell_Chondrocytes','xCell_Fibroblasts','xCell_GMP',
                     'xCell_Hepatocytes','xCell_HSC','xCell_Keratinocytes','xCell_ly_Endothelial_cells',
                     'xCell_Megakaryocytes','xCell_Melanocytes','xCell_MEP','xCell_Mesangial_cells',
                     'xCell_MPP','xCell_MSC', 'xCell_mv_Endothelial_cell','xCell_Myocytes','xCell_Osteoblast',
                     'xCell_Preadipocytes','xCell_Sebocytes')
cols.to.remove2 <- c(grep('resid',colnames(immune.subsets)),
                     which(colnames(immune.subsets) %in% cols.to.remove1))

immune.subsets.new <- immune.subsets[,-cols.to.remove2] %>%
  left_join(tx100.tcra.scores.immune, 'sample') %>%
  gather(Cell, Score, -sample, -TCRA.tcell.fraction.gc)

# Make plot with everything:
immune.test.plot <-immune.subsets.new  %>%
  filter(!is.na(TCRA.tcell.fraction.gc)) %>%
  mutate(TCRA.tcell.fraction.gc = ifelse(TCRA.tcell.fraction.gc < 0, 0, TCRA.tcell.fraction.gc)) %>%
  spread(Cell, Score) %>% 
  gather('Cell','Score', -sample, -TCRA.tcell.fraction.gc) %>%
  ggplot(aes(x=TCRA.tcell.fraction.gc, y=Score)) +
  facet_wrap(~`Cell`) +
  geom_point(aes(col = `Cell`), alpha = 0.5) +
  geom_smooth(method = "lm", se=FALSE, color="black", formula = y ~ x) +
  xlab('Tumour TCRA T cell fraction') + ylab('Danaher Score') +
  stat_cor(label.y = 8, method='spearman', cor.coef.name = 'rho',p.accuracy = 1e-100, p.digits = 3) +
  theme_bw() 


ggbuild1 <- ggplot_build(immune.test.plot)
danaher.associations <- ggplot_build(immune.test.plot)$data[[3]]
danaher.associations$danaher <- ggbuild1$layout$layout$Cell
danaher.associations$danaher[which(danaher.associations$danaher %in% colnames(TIMER.results$FMAT.all[[1]]))] <- paste0('TIMER_',danaher.associations$danaher[which(danaher.associations$danaher %in% colnames(TIMER.results$FMAT.all[[1]]))])

danaher.plot1.df <- danaher.associations %>%
  mutate(minus_logp = -log(p)) %>%
  arrange(desc(-r)) %>%
  mutate(danaher = factor(danaher, levels = danaher)) %>%
  mutate(danaher2 =  gsub('.danahaer','.d1',
                          gsub('.score','.d2',
                               gsub('TIMER_','', 
                                    gsub('.score.danaher','.d1',
                                         gsub('CIBER_','C_',
                                              gsub('EPIC_','E_',
                                                   gsub('xCell_','x_',danaher)))))))) %>%
  mutate(danaher2 = factor(danaher2, levels = danaher2))


metrics_plot <- data.frame(score = danaher.associations$danaher, 
                           minus_log_pvalue = -log2(danaher.associations$p.value)) %>%
  mutate(method = rep('Davoli',90)) %>%
  mutate(method = ifelse(grepl('xCell',score), 'xCell',
                         ifelse(grepl('danaher',score), 'Danaher',
                                ifelse(grepl('EPIC',score),'EPIC',
                                       ifelse(grepl('TIMER',score),'TIMER',
                                              ifelse(grepl('CIBER',score),'CIBERSORT',method)))))) %>%
  mutate(score2 = gsub('.danahaer','.d1',gsub('.score','.d2',gsub('TIMER_','',
                                                                  gsub('.score.danaher','.d1',
                                                                       gsub('CIBER_','C_',
                                                                            gsub('EPIC_','E_',
                                                                                 gsub('xCell_','x_',score)))))))) %>% 
  mutate(Tcell = ifelse(score2 %in% c('cd4.d2','cs4.d1','cd45.d1','cd8.exhausted.d1','cd8.d1','cd8.d2',
                                      'C_T_cells_CD8','E_CD4_Tcells','E_CD8_Tcells','T_cell.CD4','T_cell.CD8',
                                      'cyt.d2','cyto.d1',
                                      'tcells.d1','th1.d1','totat.til.d1','treg.d1','treg.d2',
                                      'x_CD4._memory_T.cells','x_CD4._naive_T.cells','x_CD4._T.cells',
                                      'x_CD4._Tcm','x_CD4._Tem','x_CD8._naive_T.cells','x_CD8._T.cells',
                                      'x_CD8._Tcm','x_CD8._Tem','x_Th1_cells','x_Th2_cells','x_Tregs'), TRUE, FALSE)) %>%
  select(score2, minus_log_pvalue, method, Tcell)

metrics_plot_join <- data.frame(danaher = danaher.associations$danaher,
                                RNAseqMethod = metrics_plot$method)

danaher.plot1.df2 <- danaher.associations %>%
  left_join(metrics_plot_join, 'danaher') %>%
  mutate(minus_logp = -log(p)) %>%
  arrange(RNAseqMethod, desc(-r)) %>%
  mutate(danaher = factor(danaher, levels = danaher)) %>%
  mutate(danaher2 =  gsub('.danahaer','.d1',
                          gsub('.score','.d2',
                               gsub('TIMER_','', 
                                    gsub('.score.danaher','.d1',
                                         gsub('CIBER_','C_',
                                              gsub('EPIC_','E_',
                                                   gsub('xCell_','x_',danaher)))))))) %>%
  mutate(danaher2 = factor(danaher2, levels = danaher2))

danaher.plot1 <- danaher.plot1.df2 %>%
  ggplot(aes(danaher2,r)) + geom_point(size = 0.5) +
  theme_bw()  + xlab('') + ylab("Spearman's rho") + 
  ggtitle('RNA-seq cell scores vs TCRA T cell fraction') + 
  theme(text = element_text(size = 7), axis.text.x = element_text(angle = 60,hjust = 1)) + 
  ylim(-0.5,0.7)

source('ggplot_add_tiles.R')

library(RColorBrewer)
library(gtable)
library(grid)
library(gridExtra)

plot.list <- list()
plot.list[[1]] <- danaher.plot1 
case.order <- as.character(danaher.plot1.df2$danaher2)

# Maybe put tracks above?
pdf('plots/fig1/rnaseq.pdf', height = 4, width = 9)
addTiles_ggplot(plot.list, metrics_plot, case.order,
                heights = c(4,0.3,0.3,0.3),
                widths = c(18,2), text_size = 6, legend_col_n = 1)
dev.off()


# Extended Figure 2a ----

# Simple implementation of T cell ExTRECT
vdj.logR.df <- readRDS('data/vdj_sim_example.RDS')

pdf('plots/supfig2/logR_example.pdf',
    height = 6, width = 12)
vdj.logR.df %>%
  ggplot(aes(pos, Ratio)) +
  geom_point() +
  geom_smooth() +
  ggtitle('Simulation proportions: 24% t-cells, 75% tumour (copy-number=1), 1% normal cells') +
  xlab('Chr14 position') + ylab('Log ratio')
dev.off()

#  Extended Figure 2b,c,d simulations  ----
TCRA.sim.final.df2 <- readRDS('data/sim_TCRA_scores_final_20201208.RDS')

# Extended figure 2b
TCRA.sim.final.df2 %>%
  filter(tumour.copynumber ==  2) %>%
  ggplot(aes(TCRA.naive.score, tcell)) +
  geom_point() +
  stat_cor(method = 'spearman', cor.coef.name = 'rho') +
  geom_abline(intercept = 0, slope = 1, colour = 'red', linetype = 'dashed') + 
  xlab('TCRA T cell fraction estimate') + ylab('Simulated T cell fraction')

# Extended figure 2c
pdf('plots/supfig2/naive_copynumber.pdf',
    height = 6, width = 6)
TCRA.sim.final.df2 %>%
  filter(tumour.copynumber !=  2) %>%
  filter(tcell!=0.25) %>%
  ggplot(aes(TCRA.naive.score, tcell)) +
  geom_point(aes(col = purity)) +
  facet_wrap(~tumour.copynumber) +
  geom_abline(intercept = 0, slope = 1, colour = 'red', linetype = 'dashed') + 
  xlab('Calculated naive score') + ylab('Simulated T ccell fraction') + 
  ggtitle('Simulation vs calculated naive score by local copy number')
dev.off()


# Extended figure 2d
pdf('plots/supfig2/exact_copynumber.pdf',
    height = 6, width = 6)
TCRA.sim.final.df2 %>%
  filter(tumour.copynumber !=  2) %>%
  filter(tcell!=0.25) %>%
  ggplot(aes(TCRA.t.cell.fraction, tcell)) +
  geom_point(aes(col = purity)) +
  facet_wrap(~tumour.copynumber)  +
  geom_abline(intercept = 0, slope = 1, colour = 'red', linetype = 'dashed') + 
  xlab('Calculated exact score') + ylab('Simulated t-cell fraction') + 
  ggtitle('Simulation vs calculated exact score by local copy number')
dev.off()


# Extended figure 2e/f - Simulated data depth ----

library(tidyverse)
library(ggpubr)

# Load output from downsampling + simulations!

Tx100.depth.df <- readRDS('data/TCRA.downasmp.CI.gc.20210615.RDS')
# Sim
sim.depth.df <- readRDS('data/TCRA.sim.CI.gc.20210615.RDS')


Tx100.depth.df$depth <- as.numeric(gsub('X','',
                                        sapply(Tx100.depth.df$sample, function(x) strsplit(as.character(x),'_')[[1]][5])))

Tx100.depth.df <- Tx100.depth.df %>%
  mutate(region = substr(sample,1,17))

# 

tx100.tcra.scores <- readRDS('data/Tx100_TCRA_scores_gc_update_20210630.RDS')

tx100.tcra.scores %>%
  dplyr::filter(sample %in% c("CRUK0040_SU_T1-R1", "CRUK0073_SU_T1-R2","CRUK0067_SU_T1-R1",
                              "CRUK0080_SU_T1-R3", "CRUK0068_SU_T1-R1")) %>%
  dplyr::select(sample, T.cell.fraction.raw.gc.smt)

tcra.coverage.df <- data.frame(sample = c("CRUK0040_SU_T1-R1", "CRUK0073_SU_T1-R2","CRUK0067_SU_T1-R1",
                                          "CRUK0080_SU_T1-R3", "CRUK0068_SU_T1-R1"),
                               patient = c('CRUK0040','CRUK0073', 'CRUK0067', 'CRUK0080', 'CRUK0068'),
                               region = c('R1','R2','R1','R3','R1'),
                               TCRA.tcell.fraction = c(0.29636720, 0.11610991, 0.05895609, 0.01105052, 0))

tmp1 <- tcra.coverage.df  %>%
  dplyr::select(region=sample, TCRA.actual.fraction=TCRA.tcell.fraction) 

Tx100.depth.df.final <-  Tx100.depth.df %>%
  left_join(tmp1, 'region')

# Extended data figure 2e
pdf('plots/ext2/depth_tx100.pdf', height = 6, width =6)
Tx100.depth.df.final %>%
  filter(depth!=5) %>%
  ggplot(aes(T.cell.fraction.raw.gc.smt, TCRA.actual.fraction)) +
  geom_point() +
  facet_wrap(~depth) +
  # geom_vline(xintercept = 0, linetype = 'dashed', colour = 'red') +
  geom_abline(slope = 1, intercept = 0) + 
  stat_cor(method = 'spearman', cor.coef.name = 'rho') + 
  xlab('downsampled TCRAT cell fraction') + ylab('TCRA T cell fraction') +
  theme_bw() + ggtitle('TCRA t-cell fraction at different depths (TRACERx)')
dev.off()


sim.depth.df2 <- sim.depth.df %>%
  mutate(sample = gsub('_TCRA.txt','',sample)) %>%
  mutate(tcell.fraction = as.numeric(gsub('X[0-9]*','',gsub('tcell.','',sample)))) %>%
  mutate(depth = as.numeric(gsub('tcell.0.[0-9]*X','',sample)))

# Extended data figure 2f
pdf('plots/ext2/depth_sim.pdf', height = 6, width =6)
sim.depth.df2 %>%
  filter(depth!=5) %>%
  ggplot(aes(T.cell.fraction.raw.smt, tcell.fraction)) +
  geom_point() + facet_wrap(~depth,nrow = 2) +
  xlab('TCRA T cell fraction') + ylab('simulated T cell fraction') +
  geom_abline(slope = 1, intercept = 0) +   stat_cor(method = 'spearman', cor.coef.name = 'rho')  + 
  theme_bw() + ggtitle('TCRA T cell fraction at different depths (Simulation)')
dev.off()
