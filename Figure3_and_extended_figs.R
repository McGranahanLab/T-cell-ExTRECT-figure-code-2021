###############################################################################
### T cell quantification from DNA sequencing predicts immunotherapy        ###
###       response                                                          ###
###                                                                         ###
### paper figure 3 + related extended figures                               ###
###                                                                         ###
### Author: Robert Bentham                                                  ###
### Date: 2021 July 30th                                                    ###
###############################################################################
library(tidyverse)
library(ggpubr)
library(ggforce)
library(rstatix)
library(survival)
library(survminer)
library(survplot)

# Figure 3 KM plots Tx100 mean ----

load('data/tracerx100.clinical.RData')
tx100.tcra.scores <- readRDS('data/Tx100_TCRA_scores_gc_update_20210630.RDS')
load('data/tracerx100.danaher.RData')
load('data/tracerx.survival.15052018.RData')


tx100.mean <- tx100.tcra.scores %>%
  mutate(germline = grepl('GL',sample)) %>%
  filter(!germline) %>%
  select(TCRA.tcell.fraction.gc) %>% `[[`(1) %>% mean(na.rm = TRUE)

thershold.TCRA.score <- tx100.mean

region.below.df <- tx100.tcra.scores %>%
  mutate(germline = grepl('GL',sample)) %>%
  filter(!germline) %>%
  dplyr::select(sample, TCRA.tcell.fraction.gc,TCRA.tcell.fraction.gc.upr.95,TCRA.tcell.fraction.gc.lwr.95) %>%
  mutate(REGTrialNo = substr(sample, 1,8)) %>%
  mutate(TCRA.tcell.fraction.gc = ifelse(TCRA.tcell.fraction.gc < 0, 0, TCRA.tcell.fraction.gc)) %>%
  mutate(belowThreshold = TCRA.tcell.fraction.gc <= thershold.TCRA.score )


region.risk.df1 <- region.below.df %>%
  group_by(REGTrialNo) %>%
  filter(!is.na(TCRA.tcell.fraction.gc)) %>%
  summarise(coldRegions = sum(belowThreshold , na.rm = TRUE),
            total.regions = n()) %>%
  mutate(risk = ifelse(coldRegions > 0, '>= 1 cold regions','< 1 cold regions'))

region.risk.df2 <- region.below.df %>%
  group_by(REGTrialNo) %>%
  summarise(coldRegions = sum(belowThreshold , na.rm = TRUE),
            total.regions = n()) %>%
  mutate(risk = ifelse(coldRegions > 1, '>= 2 cold regions','< 2 cold regions'))

region.risk.df3 <- region.below.df %>%
  group_by(REGTrialNo) %>%
  summarise(coldRegions = sum(belowThreshold , na.rm = TRUE),
            total.regions = n()) %>%
  mutate(risk = ifelse(coldRegions > 2, '>= 3 cold regions','< 3 cold regions'))


vdj.surv.df <- surv.df %>%
  left_join(region.risk.df1 , 'REGTrialNo')

vdj.surv.df2 <- surv.df %>%
  left_join(region.risk.df2 , 'REGTrialNo') %>%
  filter(total.regions >= 2)

vdj.surv.df3 <- surv.df %>%
  left_join(region.risk.df3 , 'REGTrialNo') %>%
  filter(total.regions >= 3)


vdj.surv.df.luad <- vdj.surv.df %>% filter(histology_group == "Adenocarcinoma")
vdj.surv.df.lusc <- vdj.surv.df %>% filter(histology_group == "Squamous cell carcinoma")

pdf('plots/fig3/tx100_survival_luad_1cold.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.luad, 
         main = "TRACERx100 - LUAD", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()


pdf('plots/fig3/tx100_survival_lusc_1cold.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.lusc, 
         main = "TRACERx100 - LUSC", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()


vdj.surv.df.luad2 <- vdj.surv.df2 %>% filter(histology_group == "Adenocarcinoma")
vdj.surv.df.lusc2 <- vdj.surv.df2 %>% filter(histology_group == "Squamous cell carcinoma")

pdf('plots/fig3/tx100_survival_luad_2cold.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.luad2, 
         main = "TRACERx100 - LUAD", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()


pdf('plots/fig3/tx100_survival_lusc_2cold.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.lusc2, 
         main = "TRACERx100 - LUSC", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()


vdj.surv.df.luad3 <- vdj.surv.df3 %>% filter(histology_group == "Adenocarcinoma")
vdj.surv.df.lusc3 <- vdj.surv.df3 %>% filter(histology_group == "Squamous cell carcinoma")

pdf('plots/fig3/tx100_survival_luad_3cold.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.luad3, 
         main = "TRACERx100 - LUAD", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()


pdf('plots/fig3/tx100_survival_lusc_3cold.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.lusc3, 
         main = "TRACERx100 - LUSC", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()

# Extended Figure 5 Tx100 median----
tx100.median <- tx100.tcra.scores %>%
  mutate(germline = grepl('GL',sample)) %>%
  filter(!germline) %>%
  select(TCRA.tcell.fraction.gc) %>% `[[`(1) %>% median(na.rm = TRUE)

thershold.TCRA.score <- tx100.median 

region.below.df <- tx100.tcra.scores %>%
  mutate(germline = grepl('GL',sample)) %>%
  filter(!germline) %>%
  dplyr::select(sample, TCRA.tcell.fraction.gc,TCRA.tcell.fraction.gc.upr.95,TCRA.tcell.fraction.gc.lwr.95) %>%
  mutate(REGTrialNo = substr(sample, 1,8)) %>%
  mutate(TCRA.tcell.fraction.gc = ifelse(TCRA.tcell.fraction.gc < 0, 0, TCRA.tcell.fraction.gc)) %>%
  mutate(belowThreshold = TCRA.tcell.fraction.gc <= thershold.TCRA.score )


region.risk.df1 <- region.below.df %>%
  group_by(REGTrialNo) %>%
  filter(!is.na(TCRA.tcell.fraction.gc)) %>%
  summarise(coldRegions = sum(belowThreshold , na.rm = TRUE),
            total.regions = n()) %>%
  mutate(risk = ifelse(coldRegions > 0, '>= 1 cold regions','< 1 cold regions'))

region.risk.df2 <- region.below.df %>%
  group_by(REGTrialNo) %>%
  summarise(coldRegions = sum(belowThreshold , na.rm = TRUE),
            total.regions = n()) %>%
  mutate(risk = ifelse(coldRegions > 1, '>= 2 cold regions','< 2 cold regions'))

region.risk.df3 <- region.below.df %>%
  group_by(REGTrialNo) %>%
  summarise(coldRegions = sum(belowThreshold , na.rm = TRUE),
            total.regions = n()) %>%
  mutate(risk = ifelse(coldRegions > 2, '>= 3 cold regions','< 3 cold regions'))


vdj.surv.df <- surv.df %>%
  left_join(region.risk.df1 , 'REGTrialNo')

vdj.surv.df2 <- surv.df %>%
  left_join(region.risk.df2 , 'REGTrialNo') %>%
  filter(total.regions >= 2)

vdj.surv.df3 <- surv.df %>%
  left_join(region.risk.df3 , 'REGTrialNo') %>%
  filter(total.regions >= 3)


vdj.surv.df.luad <- vdj.surv.df %>% filter(histology_group == "Adenocarcinoma")
vdj.surv.df.lusc <- vdj.surv.df %>% filter(histology_group == "Squamous cell carcinoma")

pdf('plots/ext5/tx100_survival_luad_1cold_median.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.luad, 
         main = "TRACERx100 - LUAD", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()


pdf('plots/ext5/tx100_survival_lusc_1cold_median.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.lusc, 
         main = "TRACERx100 - LUSC", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()


vdj.surv.df.luad2 <- vdj.surv.df2 %>% filter(histology_group == "Adenocarcinoma")
vdj.surv.df.lusc2 <- vdj.surv.df2 %>% filter(histology_group == "Squamous cell carcinoma")

pdf('plots/ext5/tx100_survival_luad_2cold_median.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.luad2, 
         main = "TRACERx100 - LUAD", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()


pdf('plots/ext5tx100_survival_lusc_2cold_median.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.lusc2, 
         main = "TRACERx100 - LUSC", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()


vdj.surv.df.luad3 <- vdj.surv.df3 %>% filter(histology_group == "Adenocarcinoma")
vdj.surv.df.lusc3 <- vdj.surv.df3 %>% filter(histology_group == "Squamous cell carcinoma")

pdf('plots/ext5/tx100_survival_luad_3cold_median.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.luad3, 
         main = "TRACERx100 - LUAD", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()


pdf('plots/ext5/tx100_survival_lusc_3cold_median.pdf', height = 6, width = 6)
survplot(Surv(DFS_time_days, DFS_censor_variable) ~ risk, 
         data = vdj.surv.df.lusc3, 
         main = "TRACERx100 - LUSC", 
         xlab = 'Time (days)', 
         ylab = "Disease free survival",
         col = c('darkred', 'darkblue'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()

# Extended Figure 5 TCGA ----
TCGA.scores.final <- readRDS('data/TCGA_GC_CI_final_20210520.RDS')

TCGA.scores.final <- TCGA.scores.final %>%
  mutate(sample_type = as.character(sample_type)) %>%
  mutate(sample_type = ifelse(sample_type == 'Primary Tumor','primary.tumour',
                              ifelse(sample_type == 'Blood Derived Normal','blood.normal',
                                     ifelse(sample_type == 'Solid Tissue Normal','tissue.normal',sample_type))))

luad.lusc.pfs <- readRDS('data/luad_lusc_pfs.RDS')


TCGA.scores.final <-  TCGA.scores.final %>%
  left_join(luad.lusc.pfs, 'case_id')



library(TCGAbiolinks)

clin.luad <- GDCquery_clinic("TCGA-LUAD", "clinical")
clin.lusc <- GDCquery_clinic("TCGA-LUSC", "clinical")


# TCGA mean/median threshold
tcga.threshold.median  <- TCGA.scores.final %>%
  filter(sample_type == 'primary.tumour')  %>% 
  filter(tcga == 'luad') %>%
  select(TCRA.tcell.fraction.gc.smt) %>% `[[`(1) %>% median(na.rm = TRUE)

tcga.threshold.mean  <- TCGA.scores.final %>%
  filter(sample_type == 'primary.tumour')  %>% 
  filter(tcga == 'luad') %>%
  select(TCRA.tcell.fraction.gc.smt) %>% `[[`(1) %>% mean(na.rm = TRUE)


TCGA.groups.luad <- TCGA.scores.final %>%
  filter(sample_type == 'primary.tumour') %>%
  filter(tcga == 'luad') %>%
  dplyr::select(bcr_patient_barcode = case_id, TCRA.tcell.fraction.gc.smt, pfs.days, pfs.censor) %>%
  filter(bcr_patient_barcode %in% clin.luad$bcr_patient_barcode) %>%
  mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > tcga.threshold.mean, 'Immune hot','Immune cold'))

TCGA.groups.luad.median <- TCGA.scores.final %>%
  filter(sample_type == 'primary.tumour') %>%
  filter(tcga == 'luad') %>%
  dplyr::select(bcr_patient_barcode = case_id, TCRA.tcell.fraction.gc.smt, pfs.days, pfs.censor) %>%
  filter(bcr_patient_barcode %in% clin.luad$bcr_patient_barcode) %>%
  mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > tcga.threshold.median, 'Immune hot','Immune cold'))


TCGA.groups.lusc <- TCGA.scores.final %>%
  filter(sample_type == 'primary.tumour') %>%
  filter(tcga == 'lusc') %>%
  dplyr::select(bcr_patient_barcode = case_id, TCRA.tcell.fraction.gc.smt,pfs.days, pfs.censor) %>%
  filter(bcr_patient_barcode %in% clin.lusc$bcr_patient_barcode) %>%
  # mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > lusc.median  , 'Immune hot','Immune cold'))
  mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > tcga.threshold.mean  , 'Immune hot','Immune cold'))

TCGA.groups.all <- TCGA.scores.final %>%
  filter(sample_type %in% c('primary.tumour')) %>%
  dplyr::select(bcr_patient_barcode = case_id, TCRA.tcell.fraction.gc.smt, pfs.days, pfs.censor) %>%
  filter(bcr_patient_barcode %in% c(clin.luad$bcr_patient_barcode,clin.lusc$bcr_patient_barcode)) %>%
  # mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > lusc.median  , 'Immune hot','Immune cold'))
  mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > tcga.threshold.mean , 'Immune hot','Immune cold'))

clin.luad.upd <- clin.luad %>%
  mutate(death.censor = ifelse(vital_status == 'Dead',1,0)) %>%
  mutate(days = ifelse(!is.na(days_to_death),days_to_death,days_to_last_follow_up)) %>%
  left_join(TCGA.groups.luad, 'bcr_patient_barcode') %>%
  filter(!is.na(TCRA.tcell.group))

clin.luad.upd.median <- clin.luad %>%
  mutate(death.censor = ifelse(vital_status == 'Dead',1,0)) %>%
  mutate(days = ifelse(!is.na(days_to_death),days_to_death,days_to_last_follow_up)) %>%
  left_join(TCGA.groups.luad.median, 'bcr_patient_barcode') %>%
  filter(!is.na(TCRA.tcell.group))


clin.lusc.upd <- clin.lusc %>%
  mutate(death.censor = ifelse(vital_status == 'Dead',1,0)) %>%
  mutate(days = ifelse(!is.na(days_to_death),days_to_death,days_to_last_follow_up)) %>%
  left_join(TCGA.groups.lusc, 'bcr_patient_barcode') %>%
  filter(!is.na(TCRA.tcell.group))

clin.all.upd <- rbind(clin.luad,clin.lusc) %>%
  mutate(death.censor = ifelse(vital_status == 'Dead',1,0)) %>%
  mutate(days = ifelse(!is.na(days_to_death),days_to_death,days_to_last_follow_up)) %>%
  left_join(TCGA.groups.all, 'bcr_patient_barcode') %>%
  filter(!is.na(TCRA.tcell.group))

clin.all.upd_scaled <- clin.all.upd %>%
  mutate( TCRA.tcell.fraction.gc.smt = scale( TCRA.tcell.fraction.gc.smt))

clin.luad.upd_scaled <- clin.luad.upd %>%
  mutate( TCRA.tcell.fraction.gc.smt = scale( TCRA.tcell.fraction.gc.smt))

clin.lusc.upd_scaled <- clin.lusc.upd %>%
  mutate( TCRA.tcell.fraction.gc.smt = scale( TCRA.tcell.fraction.gc.smt))

# Using TCGA LUAD mean 
pdf('plots/ext5/tcga_luad_os_surv_mean.pdf',
    width = 6, height = 6)
survplot(Surv(days, death.censor) ~ TCRA.tcell.group, 
         data = clin.luad.upd, 
         main = "TCGA LUAD", 
         xlab = 'Time (days)', 
         ylab = "Overall survival",
         col = c('darkblue','darkred'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()

pdf('plots/ext5/tcga_luad_pfs_surv_mean.pdf',
    width = 6, height = 6)
survplot(Surv(pfs.days, pfs.censor) ~ TCRA.tcell.group, 
         data = clin.luad.upd, 
         main = "TCGA LUAD", 
         xlab = 'Time (days)', 
         ylab = "PFS",
         col = c('darkblue','darkred'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()

# Sup Figure 4C  TCGA KM plots ----
pdf('plots/ext5/tcga_lusc_os_surv_mean.pdf',
    width = 6, height = 6)
survplot(Surv(days, death.censor) ~ TCRA.tcell.group, 
         data = clin.lusc.upd, 
         main = "TCGA LUSC", 
         xlab = 'Time (days)', 
         ylab = "Overall survival",
         col = c('darkblue','darkred'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()

pdf('plots/ext5/tcga_lusc_pcf_surv_mean.pdf',
    width = 6, height = 6)
survplot(Surv(pfs.days, pfs.censor) ~ TCRA.tcell.group, 
         data = clin.lusc.upd, 
         main = "TCGA LUSC", 
         xlab = 'Time (days)', 
         ylab = "PFS",
         col = c('darkblue','darkred'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()

pdf('plots/ext5/tcga_all_os_surv_mean.pdf',
    width = 6, height = 6)
survplot(Surv(days, death.censor) ~ TCRA.tcell.group, 
         data = clin.all.upd, 
         main = "TCGA LUAD & LUSC", 
         xlab = 'Time (days)', 
         ylab = "Overall survival",
         col = c('darkblue','darkred'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()

pdf('plots/ext5/tcga_all_pcf_surv_mean.pdf',
    width = 6, height = 6)
survplot(Surv(pfs.days, pfs.censor) ~ TCRA.tcell.group, 
         data = clin.all.upd, 
         main = "TCGA LUAD & LUSC", 
         xlab = 'Time (days)', 
         ylab = "PFS",
         col = c('darkblue','darkred'),
         lwd = 1, las = 1, cex.axis = 0.8)
dev.off()



# Sup Figure 4B  TCGA ranges ----
thershold.TCRA.score.vec <- seq(0,0.16,0.0025)
tcga_1 <- list()
tcga_luad_1 <- list()
tcga_lusc_1 <- list()

TCGA.scores.final <- readRDS('data/TCGA_GC_CI_final_20210520.RDS')



TCGA.scores.final <- TCGA.scores.final %>%
  mutate(sample_type = as.character(sample_type)) %>%
  mutate(sample_type = ifelse(sample_type == 'Primary Tumor','primary.tumour',
                              ifelse(sample_type == 'Blood Derived Normal','blood.normal',
                                     ifelse(sample_type == 'Solid Tissue Normal','tissue.normal',sample_type)))) %>%
  mutate(TCRA.tcell.fraction.gc.smt = ifelse(TCRA.tcell.fraction.gc.smt.lwr > maxPossible,
                                             T.cell.fraction.raw.gc.smt, TCRA.tcell.fraction.gc.smt)) %>%
  mutate(TCRA.tcell.fraction.gc.smt.lwr = ifelse(TCRA.tcell.fraction.gc.smt.lwr > maxPossible,
                                                 T.cell.fraction.raw.gc.smt.lwr, TCRA.tcell.fraction.gc.smt.lwr)) %>%
  mutate(TCRA.tcell.fraction.gc.smt.upr = ifelse(TCRA.tcell.fraction.gc.smt.lwr > maxPossible,
                                                 T.cell.fraction.raw.gc.smt.upr, TCRA.tcell.fraction.gc.smt.upr)) %>%
  mutate(TCRA.tcell.fraction.gc = TCRA.tcell.fraction.gc.smt) %>%
  left_join(luad.lusc.pfs, 'case_id')


TCGA.groups.luad <- TCGA.scores.final %>%
  filter(sample_type == 'primary.tumour') %>%
  filter(tcga == 'luad') %>%
  dplyr::select(bcr_patient_barcode = case_id, TCRA.tcell.fraction.gc.smt, pfs.days, pfs.censor) %>%
  filter(bcr_patient_barcode %in% clin.luad$bcr_patient_barcode) %>%
  # mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > lusc.median  , 'Immune hot','Immune cold'))
  mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > 0.1  , 'Immune hot','Immune cold'))

TCGA.groups.lusc <- TCGA.scores.final %>%
  #filter(smt.plausible) %>%
  # filter(sample_type == 'primary.tumour') %>%
  filter(sample_type == 'primary.tumour') %>%
  filter(tcga == 'lusc') %>%
  dplyr::select(bcr_patient_barcode = case_id, TCRA.tcell.fraction.gc.smt,pfs.days, pfs.censor) %>%
  filter(bcr_patient_barcode %in% clin.lusc$bcr_patient_barcode) %>%
  # mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > lusc.median  , 'Immune hot','Immune cold'))
  mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > 0.1  , 'Immune hot','Immune cold'))

TCGA.groups.all <- TCGA.scores.final %>%
  # filter(sample_type == 'primary.tumour') %>%
  filter(sample_type %in% c('primary.tumour')) %>%
  dplyr::select(bcr_patient_barcode = case_id, TCRA.tcell.fraction.gc.smt, pfs.days, pfs.censor) %>%
  filter(bcr_patient_barcode %in% c(clin.luad$bcr_patient_barcode,clin.lusc$bcr_patient_barcode)) %>%
  # mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > lusc.median  , 'Immune hot','Immune cold'))
  mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc.smt > 0.1  , 'Immune hot','Immune cold'))


clin.luad.upd <- clin.luad %>%
  mutate(death.censor = ifelse(vital_status == 'Dead',1,0)) %>%
  mutate(days = ifelse(!is.na(days_to_death),days_to_death,days_to_last_follow_up)) %>%
  left_join(TCGA.groups.luad, 'bcr_patient_barcode') %>%
  filter(!is.na(TCRA.tcell.group))

clin.lusc.upd <- clin.lusc %>%
  mutate(death.censor = ifelse(vital_status == 'Dead',1,0)) %>%
  mutate(days = ifelse(!is.na(days_to_death),days_to_death,days_to_last_follow_up)) %>%
  left_join(TCGA.groups.lusc, 'bcr_patient_barcode') %>%
  filter(!is.na(TCRA.tcell.group))

clin.all.upd <- rbind(clin.luad,clin.lusc) %>%
  mutate(death.censor = ifelse(vital_status == 'Dead',1,0)) %>%
  mutate(days = ifelse(!is.na(days_to_death),days_to_death,days_to_last_follow_up)) %>%
  left_join(TCGA.groups.all, 'bcr_patient_barcode') %>%
  filter(!is.na(TCRA.tcell.group))


for(i in 1:length(thershold.TCRA.score.vec)){
  threshold.value <- thershold.TCRA.score.vec[i]
  
  TCGA.groups.luad <- TCGA.scores.final %>%
    filter(sample_type == 'primary.tumour') %>%
    filter(tcga == 'luad') %>%
    dplyr::select(bcr_patient_barcode = case_id, TCRA.tcell.fraction.gc, pfs.days, pfs.censor, TCRA.tcell.fraction.gc.smt.lwr) %>%
    filter(bcr_patient_barcode %in% clin.luad$bcr_patient_barcode) %>%
    mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc > threshold.value  , 'Immune hot','Immune cold'))
  
  
  TCGA.groups.lusc <- TCGA.scores.final %>%
    filter(sample_type == 'primary.tumour') %>%
    filter(tcga == 'lusc') %>%
    dplyr::select(bcr_patient_barcode = case_id, TCRA.tcell.fraction.gc, pfs.days, pfs.censor, TCRA.tcell.fraction.gc.smt.lwr) %>%
    filter(bcr_patient_barcode %in% clin.lusc$bcr_patient_barcode) %>%
    mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc > threshold.value  , 'Immune hot','Immune cold'))
  
  TCGA.groups.all <- TCGA.scores.final %>%
    filter(sample_type == 'primary.tumour') %>%
    dplyr::select(bcr_patient_barcode = case_id, TCRA.tcell.fraction.gc, pfs.days, pfs.censor,TCRA.tcell.fraction.gc.smt.lwr) %>%
    filter(bcr_patient_barcode %in% c(clin.luad$bcr_patient_barcode,clin.lusc$bcr_patient_barcode)) %>%
    mutate(TCRA.tcell.group = ifelse(TCRA.tcell.fraction.gc > threshold.value  , 'Immune hot','Immune cold'))
  
  clin.luad.upd <- clin.luad %>%
    mutate(death.censor = ifelse(vital_status == 'Dead',1,0)) %>%
    mutate(days = ifelse(!is.na(days_to_death),days_to_death,days_to_last_follow_up)) %>%
    left_join(TCGA.groups.luad, 'bcr_patient_barcode') %>%
    filter(!is.na(TCRA.tcell.group))
  
  clin.lusc.upd <- clin.lusc %>%
    mutate(death.censor = ifelse(vital_status == 'Dead',1,0)) %>%
    mutate(days = ifelse(!is.na(days_to_death),days_to_death,days_to_last_follow_up)) %>%
    left_join(TCGA.groups.lusc, 'bcr_patient_barcode') %>%
    filter(!is.na(TCRA.tcell.group))
  
  clin.all.upd <- rbind(clin.luad,clin.lusc) %>%
    mutate(death.censor = ifelse(vital_status == 'Dead',1,0)) %>%
    mutate(days = ifelse(!is.na(days_to_death),days_to_death,days_to_last_follow_up)) %>%
    left_join(TCGA.groups.all, 'bcr_patient_barcode') %>%
    filter(!is.na(TCRA.tcell.group))
  
  tcga_luad_1[[i]] <- survplot(Surv(days, death.censor) ~ TCRA.tcell.group, 
                               data = clin.luad.upd, 
                               main = "TCGA LUAD", 
                               xlab = 'Time (days)', 
                               ylab = "Overall survival",
                               col = c('darkblue','darkred'),
                               lwd = 1, las = 1, cex.axis = 0.8)
  
  tcga_lusc_1[[i]] <- survplot(Surv(days, death.censor) ~ TCRA.tcell.group, 
                               data = clin.lusc.upd, 
                               main = "TCGA LUSC", 
                               xlab = 'Time (days)', 
                               ylab = "Overall survival",
                               col = c('darkblue','darkred'),
                               lwd = 1, las = 1, cex.axis = 0.8)
  
  
  tcga_1[[i]] <- survplot(Surv(days, death.censor) ~ TCRA.tcell.group, 
                          data = clin.all.upd, 
                          main = "TCGA LUAD & LUSC", 
                          xlab = 'Time (days)', 
                          ylab = "Overall survival",
                          col = c('darkblue','darkred'),
                          lwd = 1, las = 1, cex.axis = 0.8)
  
}


tcga_hr <- gsub(')','',gsub('\\([0-9]*.[0-9]*[e+00]*-[0-9]*.[0.9]*[nf]*','',gsub(' ','',gsub('HR = ','',sapply(tcga_1 , `[`,1))))) %>% as.numeric()
tcga_luad_hr <- gsub(')','',gsub('\\([0-9]*.[0-9]*[e+00]*-[0-9]*.[0.9]*[nf]*','',gsub(' ','',gsub('HR = ','',sapply(tcga_luad_1, `[`,1))))) %>% as.numeric()
tcga_lusc_hr <- gsub(')','',gsub('\\([0-9]*.[0-9]*[e+00]*-[0-9]*.[0.9]*[nf]*','',gsub(' ','',gsub('HR = ','',sapply(tcga_lusc_1, `[`,1))))) %>% as.numeric()

tcga_p1 <- as.numeric(gsub('logrank P = ','',sapply(tcga_1, `[`,2)))
tcga_luad_p2 <- as.numeric(gsub('logrank P = ','',sapply(tcga_luad_1, `[`,2)))
tcga_lusc_p3 <- as.numeric(gsub('logrank P = ','',sapply(tcga_lusc_1, `[`,2)))

tcga_survival_df_results <- data.frame(threshold = rep(thershold.TCRA.score.vec,3),
                                       HR = c(tcga_hr, tcga_luad_hr, tcga_lusc_hr),
                                       pvalue = c(tcga_p1,tcga_luad_p2,tcga_lusc_p3),
                                       regions = c(rep('tcga',65),rep('tcga-luad',65),rep('tcga-lusc',65)))

tcga_survival_df_results$HR[which(tcga_survival_df_results$HR > 1000)] <- NA

tcga_survival_df_results <- tcga_survival_df_results %>%
  mutate(sig.level = ifelse(pvalue < 0.0001,'****',
                            ifelse(pvalue < 0.001,'***',
                                   ifelse(pvalue < 0.01, '**',
                                          ifelse(pvalue < 0.05, '*','')))))

library(scales)

# Do OS and PFS separately with code above
pdf('plots/ext5/tcga_os_threshold.pdf',
    width = 2.5, height = 3.5)
tcga_survival_df_results  %>%
  ggplot() +
  geom_tile(aes(x=regions, y=threshold, fill = log2(HR))) + 
  xlab('') + ylab('') + theme_bw() + 
  scale_fill_gradientn(limits = c(-2,2),
                       colours=c(muted('blue'), "white", "red")) + 
  annotate("text", x = tcga_survival_df_results$regions,
           y = tcga_survival_df_results$threshold,
           label = tcga_survival_df_results$sig.level) + 
  ggtitle('Hazard ratio of immune hot vs cold groups') + ylab('Threshold') + xlab('Cohort') + 
  theme(text = element_text(size = 6), legend.key.size = unit(0.25,'cm'))
dev.off()


# Extended data 5  Cox models ----

# Tx100:
region.summary.df1 <- region.below.df %>%
  group_by(REGTrialNo) %>%
  filter(!is.na(TCRA.tcell.fraction.gc)) %>%
  summarise(coldRegions = sum(belowThreshold , na.rm = TRUE),
            mean.TCRA.tcell.fraction = mean(TCRA.tcell.fraction.gc, na.rm = TRUE),
            sd.TCRA.tcell.fraction = sd(TCRA.tcell.fraction.gc, na.rm = TRUE),
            median.TCRA.tcell.fraction = median(TCRA.tcell.fraction.gc, na.rm = TRUE),
            max.TCRA.tcell.fraction = max(TCRA.tcell.fraction.gc, na.rm = TRUE),
            min.TCRA.tcell.fraction = min(TCRA.tcell.fraction.gc, na.rm = TRUE),
            min.TCRA.tcell.fraction2 = min(TCRA.tcell.fraction.gc.upr.95, na.rm = TRUE),
            relative.region.immune.escape = max(TCRA.tcell.fraction.gc, na.rm = TRUE)/min(TCRA.tcell.fraction.gc.upr.95, na.rm = TRUE),
            total.regions = n()) %>%
  mutate(risk = ifelse(coldRegions > 1, '>= 1 cold regions','< 1 cold regions'))

# Cox model first:
forest.surv <- surv.df %>%
  left_join(region.summary.df1 , 'REGTrialNo') %>%
  mutate(TRACERxID = gsub('^[A-Z]_LTX','LTX0', REGTrialNo)) %>%
  filter(TRACERxID %in% clinical.data.tx100$TRACERxID)

forest.surv_scaled <- forest.surv %>%
  mutate(mean.TCRA.tcell.fraction = scale(mean.TCRA.tcell.fraction),
         min.TCRA.tcell.fraction2 = scale(min.TCRA.tcell.fraction2),
         max.TCRA.tcell.fraction = scale(max.TCRA.tcell.fraction),
         relative.region.immune.escape = scale(relative.region.immune.escape)) 

# Following set for TRACERx100, TRACERx100 LUAD and TRACERx100 LUSC 
# For combining:
ggforest_2 <- function (model, data = NULL, main = "Hazard ratio", cpositions = c(0.02, 
                                                                                  0.22, 0.4), 
                        fontsize = 2, refLabel = "reference", noDigits = 2, min_x, max_x, breaks) 
{
  conf.high <- conf.low <- estimate <- NULL
  stopifnot(inherits(model, "coxph"))
  data <- survminer:::.get_data(model, data = data)
  terms <- attr(model$terms, "dataClasses")[-1]
  coef <- as.data.frame(tidy(model, conf.int = TRUE))
  gmodel <- broom.mixed::glance(model)
  allTerms <- lapply(seq_along(terms), function(i) {
    var <- names(terms)[i]
    if (terms[i] %in% c("factor", "character")) {
      adf <- as.data.frame(table(data[, var]))
      cbind(var = var, adf, pos = 1:nrow(adf))
    }
    else if (terms[i] == "numeric") {
      data.frame(var = var, Var1 = "", Freq = nrow(data), 
                 pos = 1)
    }
    else {
      vars = grep(paste0("^", var, "*."), coef$term, value = TRUE)
      data.frame(var = vars, Var1 = "", Freq = nrow(data), 
                 pos = seq_along(vars))
    }
  })
  allTermsDF <- do.call(rbind, allTerms)
  colnames(allTermsDF) <- c("var", "level", "N", "pos")
  inds <- apply(allTermsDF[, 1:2], 1, paste0, collapse = "")
  rownames(coef) <- gsub(coef$term, pattern = "`", replacement = "")
  toShow <- cbind(allTermsDF, coef[inds, ])[, c("var", "level", 
                                                "N", "p.value", "estimate", "conf.low", "conf.high", 
                                                "pos")]
  toShowExp <- toShow[, 5:7]
  toShowExp[is.na(toShowExp)] <- 0
  toShowExp <- format(exp(toShowExp), digits = noDigits)
  toShowExpClean <- data.frame(toShow, pvalue = signif(toShow[, 
                                                              4], noDigits + 1), toShowExp)
  toShowExpClean$stars <- paste0(round(toShowExpClean$p.value, 
                                       noDigits + 1), " ", ifelse(toShowExpClean$p.value < 0.05, 
                                                                  "*", ""), ifelse(toShowExpClean$p.value < 0.01, "*", 
                                                                                   ""), ifelse(toShowExpClean$p.value < 0.001, "*", ""))
  toShowExpClean$ci <- paste0("(", toShowExpClean[, "conf.low.1"], 
                              " - ", toShowExpClean[, "conf.high.1"], ")")
  toShowExpClean$estimate.1[is.na(toShowExpClean$estimate)] = refLabel
  toShowExpClean$stars[which(toShowExpClean$p.value < 0.001)] = "<0.001 ***"
  toShowExpClean$stars[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$ci[is.na(toShowExpClean$estimate)] = ""
  toShowExpClean$estimate[is.na(toShowExpClean$estimate)] = 0
  toShowExpClean$var = as.character(toShowExpClean$var)
  toShowExpClean$var[duplicated(toShowExpClean$var)] = ""
  toShowExpClean$N <- paste0("(N=", toShowExpClean$N, ")")
  toShowExpClean <- toShowExpClean[nrow(toShowExpClean):1, 
  ]
  rangeb <- range(toShowExpClean$conf.low, toShowExpClean$conf.high, 
                  na.rm = TRUE)
  # breaks <- axisTicks(rangeb/2, log = TRUE, nint = 7)
  # breaks <- c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,1.9,2.3)
  rangeplot <- rangeb
  rangeplot[1] <- rangeplot[1] - diff(rangeb)
  rangeplot[2] <- rangeplot[2] + 0.15 * diff(rangeb)
  width <- diff(rangeplot)
  y_variable <- rangeplot[1] + cpositions[1] * width
  y_nlevel <- rangeplot[1] + cpositions[2] * width
  y_cistring <- rangeplot[1] + cpositions[3] * width
  y_stars <- rangeb[2]
  x_annotate <- seq_len(nrow(toShowExpClean))
  annot_size_mm <- fontsize 
  
  p <- ggplot(toShowExpClean, aes(seq_along(var), exp(estimate))) + 
    geom_rect(aes(xmin = seq_along(var) - 0.5, xmax = seq_along(var) +0.5,
                  ymin = min_x, ymax = max_x, 
                  fill = ordered(seq_along(var)%%2 + 1))) + 
    scale_fill_manual(values = c("#FFFFFF33","#00000033"), guide = "none") +
    geom_point(pch = 15,size = 1) + 
    geom_errorbar(aes(ymin = exp(conf.low), ymax = exp(conf.high)),width = 0.15) + 
    geom_hline(yintercept = 1, linetype = 3) + 
    coord_flip(ylim = c(min_x,max_x), clip = 'off') + 
    ggtitle(main) + 
    scale_y_log10(name = "",labels = sprintf("%g", breaks), expand = c(0.02, 0.02),
                  breaks = breaks) + theme_light() + 
    theme(panel.grid.minor.y = element_blank(), panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_blank(), legend.position = "none", panel.border = element_blank(),
          axis.title.y = element_blank(), axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), plot.title = element_text(hjust = 0.5),
          axis.text.x = element_text(size = 6),
          plot.margin = unit(c(-0.5, 0.5, -0.5, 0.5), "cm")) +
    xlab("") + 
    # name
    annotate(geom = "text", x = x_annotate,
             #y = exp(y_variable),
             y = min_x, 
             label = toShowExpClean$var, fontface = "bold", hjust = 0,
             size = annot_size_mm) + 
    # estimate
    annotate(geom = "text", x = x_annotate,y = exp(y_nlevel),
             hjust = 0, label = toShowExpClean$level,
             vjust = -0.1, size = annot_size_mm) + 
    #annotate(geom = "text", x = x_annotate, y = exp(y_nlevel),
    #         label = toShowExpClean$N,fontface = "italic", hjust = 0,
    #         vjust = ifelse(toShowExpClean$level == "", 0.5, 1.1), size = annot_size_mm) + 
    # CI
    annotate(geom = "text",x = x_annotate, y = exp(y_cistring),
             label = toShowExpClean$estimate.1,  size = annot_size_mm,
             vjust = ifelse(toShowExpClean$estimate.1 == "reference", 0.5, -0.1)) +
    annotate(geom = "text", x = x_annotate, y = exp(y_cistring), label = toShowExpClean$ci,
             size = annot_size_mm, vjust = 1.1, fontface = "italic") + 
    # Sig
    annotate(geom = "text", x = x_annotate, y = max_x, 
             label = toShowExpClean$stars, size = annot_size_mm, 
             hjust = -0.2, fontface = "italic") 
  # annotate(geom = "text",  x = 0.5, y = exp(y_variable),
  #         label = paste0("# Events: ", gmodel$nevent, "; Global p-value (Log-Rank): ",
  #                          format.pval(gmodel$p.value.log, eps = ".001"), " \nAIC: ",
  #                         round(gmodel$AIC,  2), "; Concordance Index: ",
  #                        round(gmodel$concordance, 2)),
  #        size = annot_size_mm, hjust = 0, vjust = 1.2, fontface = "italic")
  
  gt <- ggplot_gtable(ggplot_build(p))
  gt$layout$clip[gt$layout$name == "panel"] <- "off"
  ggpubr::as_ggplot(gt)
  
  
}


forest_plot_list <- list()

forest_plot_list[[1]] <- ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                             mean.TCRA.tcell.fraction,
                                           data = forest.surv_scaled), main = '',min_x = 0.1, max_x = 2.5,
                                    breaks = c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,1.9,2.3))


forest_plot_list[[2]] <- ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                             min.TCRA.tcell.fraction2,
                                           data = forest.surv_scaled), main = '',min_x = 0.1, max_x = 2.5,
                                    breaks = c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,1.9,2.3))

forest_plot_list[[3]] <-ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                            max.TCRA.tcell.fraction,
                                          data = forest.surv_scaled), main = '',min_x = 0.1, max_x = 2.5,
                                   breaks = c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,1.9,2.3))

forest_plot_list[[4]] <-ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                            relative.region.immune.escape,
                                          data = forest.surv_scaled), main = '',min_x = 0.1, max_x = 2.5,
                                   breaks = c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,1.9,2.3))

forest_plot_list[[5]] <- ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                             min.TCRA.tcell.fraction2 + max.TCRA.tcell.fraction,
                                           data = forest.surv_scaled), main = '',min_x = 0.1, max_x = 2.5,
                                    breaks = c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,1.9,2.3))

# Join these plots together:

source('ggplot_add_tiles.R')

library(RColorBrewer)
library(gtable)
library(grid)
library(gridExtra)

ggplot_combine_list <- function(plot.list, heights){
  plot_list_ggplotGrop <- lapply(1:length(plot.list), function(x){
    p <- plot.list[[x]]
    if(x == 1){
      p <- ggplotGrob(p + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")))
    } else if(x == length(plot.list)) {
      p <- ggplotGrob(p + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")))
    } else{
      p <- ggplotGrob(p + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")))
    }
    return(p)
  })
  
  all_widths <- lapply(plot_list_ggplotGrop, function(x) {x$widths})
  plot_list_alignedWiths <- lapply(plot_list_ggplotGrop, function(x){
    x$widths <- do.call(unit.pmax, all_widths)
    return(x)
  })
  
  
  #set heights
  g <- do.call(gtable_rbind, plot_list_alignedWiths)
  id_panels_h <- unique(g$layout[g$layout$name=="panel", "t"])
  g$heights[id_panels_h] <- grid::unit(heights, "null")
  
  grid.arrange(g)
  
}

pdf('plots/ext5/forest_continuous_tx100.pdf',
    width = 3, height = 3)
ggplot_combine_list(forest_plot_list, rep(1,5))
dev.off()

# LUAD:
forest_plot_list2 <- list()

forest_plot_list2[[1]] <- ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                              mean.TCRA.tcell.fraction,
                                            data = forest.surv_scaled[forest.surv_scaled$histology_group == 'Adenocarcinoma',]),
                                     main = '',min_x = 0.1, max_x = 4.5,
                                     breaks = c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,2,3,4))

forest_plot_list2[[2]] <- ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                              min.TCRA.tcell.fraction2,
                                            data = forest.surv_scaled[forest.surv_scaled$histology_group == 'Adenocarcinoma',]),
                                     main = '',min_x = 0.1, max_x = 4.5,
                                     breaks = c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,2,3,4))

forest_plot_list2[[3]] <-ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                             max.TCRA.tcell.fraction,
                                           data = forest.surv_scaled[forest.surv_scaled$histology_group == 'Adenocarcinoma',]),
                                    main = '',min_x = 0.1, max_x = 4.5,
                                    breaks = c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,2,3,4))

forest_plot_list2[[4]] <-ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                             relative.region.immune.escape,
                                           data = forest.surv_scaled[forest.surv_scaled$histology_group == 'Adenocarcinoma',]),
                                    main = '',min_x = 0.1, max_x = 4.5,
                                    breaks = c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,2,3,4))

forest_plot_list2[[5]] <- ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                              min.TCRA.tcell.fraction2 + max.TCRA.tcell.fraction,
                                            data = forest.surv_scaled[forest.surv_scaled$histology_group == 'Adenocarcinoma',]),
                                     main = '',min_x = 0.1, max_x = 4.5,
                                     breaks = c(0.1,0.3,0.5,0.7,0.9,1.1,1.5,2,3,4))

pdf('plots/ext5/forest_continuous_tx100_luad.pdf',
    width = 3, height = 3)
ggplot_combine_list(forest_plot_list2, rep(1,5))
dev.off()

# LUSC:
forest_plot_list3 <- list()


forest_plot_list3[[1]] <- ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                              mean.TCRA.tcell.fraction,
                                            data = forest.surv_scaled[forest.surv_scaled$histology_group == 'Adenocarcinoma',]),
                                     main = '',min_x = 0.1, max_x = 2.5,
                                     breaks = c(0.1,0.3,0.5,0.7,1.5,1.9,2.3))

forest_plot_list3[[2]] <- ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                              min.TCRA.tcell.fraction2,
                                            data = forest.surv_scaled[forest.surv_scaled$histology_group == 'Squamous cell carcinoma',]),
                                     main = '',min_x = 0.1, max_x = 2.5,
                                     breaks = c(0.1,0.3,0.5,0.7,1.5,1.9,2.3))

forest_plot_list3[[3]] <-ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                             max.TCRA.tcell.fraction,
                                           data = forest.surv_scaled[forest.surv_scaled$histology_group == 'Squamous cell carcinoma',]),
                                    main = '',min_x = 0.1, max_x = 2.5,
                                    breaks = c(0.1,0.3,0.5,0.7,1.5,1.9,2.3))

forest_plot_list3[[4]] <-ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                             relative.region.immune.escape,
                                           data = forest.surv_scaled[forest.surv_scaled$histology_group == 'Squamous cell carcinoma',]),
                                    main = '',min_x = 0.1, max_x = 2.5,
                                    breaks = c(0.1,0.3,0.5,0.7,1.5,1.9,2.3))

forest_plot_list3[[5]] <- ggforest_2(coxph( Surv(DFS_time_days, DFS_censor_variable) ~ 
                                              min.TCRA.tcell.fraction2 + max.TCRA.tcell.fraction,
                                            data = forest.surv_scaled[forest.surv_scaled$histology_group == 'Squamous cell carcinoma',]),
                                     main = '',min_x = 0.1, max_x = 2.5,
                                     breaks = c(0.1,0.3,0.5,0.7,1.5,1.9,2.3))

pdf('plots/ext5/forest_continuous_tx100_lusc.pdf',
    width = 2, height = 3)
ggplot_combine_list(forest_plot_list3, rep(1,5))
dev.off()

# TCGA: 

forest_plot_list4 <- list()

forest_plot_list4[[1]] <- ggforest_2(coxph( Surv(days, death.censor) ~  TCRA.tcell.fraction.gc.smt,
                                            data = clin.all.upd_scaled),min_x = 0.5, max_x = 1.5, breaks = c(0.5,0.75,1,1.2),main = '')
forest_plot_list4[[2]] <- ggforest_2(coxph( Surv(pfs.days, pfs.censor) ~  TCRA.tcell.fraction.gc.smt,
                                            data = clin.all.upd_scaled),min_x = 0.5, max_x = 1.5, breaks = c(0.5,0.75,1,1.2),main = '')

forest_plot_list4[[3]] <-ggforest_2(coxph( Surv(days, death.censor) ~  TCRA.tcell.fraction.gc.smt,
                                           data = clin.luad.upd_scaled),min_x = 0.5, max_x = 1.5, breaks = c(0.5,0.75,1,1.2),main = '')
forest_plot_list4[[4]] <-ggforest_2(coxph( Surv(pfs.days, pfs.censor) ~  TCRA.tcell.fraction.gc.smt,
                                           data = clin.luad.upd_scaled),min_x = 0.5, max_x = 1.5, breaks = c(0.5,0.75,1,1.2),main = '')

forest_plot_list4[[5]] <- ggforest_2(coxph( Surv(days, death.censor) ~  TCRA.tcell.fraction.gc.smt,
                                            data = clin.lusc.upd_scaled),min_x = 0.5, max_x = 1.5, breaks = c(0.5,0.75,1,1.2),main = '')
forest_plot_list4[[6]] <- ggforest_2(coxph( Surv(pfs.days, pfs.censor) ~  TCRA.tcell.fraction.gc.smt,
                                            data = clin.lusc.upd_scaled),min_x = 0.5, max_x = 1.5, breaks = c(0.5,0.75,1,1.2),main = '')

pdf('plots/ext5/forest_continuous_tcga.pdf',
    width = 2, height = 3)
ggplot_combine_list(forest_plot_list4, rep(1,5))
dev.off()

