###############################################################################
### T cell quantification from DNA sequencing predicts immunotherapy        ###
###       response                                                          ###
###                                                                         ###
### paper figure 2 + related extended figures                               ###
###                                                                         ###
### Author: Robert Bentham                                                  ###
### Date: 2021 July 30th                                                    ###
###############################################################################

library(tidyverse)
library(ggforce)
library(ggpubr)
library(rstatix)

# Extended Figure 3a ----
load('data/tracerx100.clinical.RData')
tx100.tcra.scores <- readRDS('data/Tx100_TCRA_scores_gc_update_20210630.RDS')

tx100.tcra.scores <- tx100.tcra.scores  %>%
  mutate(TCRA.tcell.fraction = TCRA.tcell.fraction.gc)

age.df <- clinical.data.tx100 %>%
  dplyr::select(patient = TRACERxID, Age, Gender, Histology) 


response.plot.data <- tx100.tcra.scores %>%
  dplyr::select(sample, TCRA.tcell.fraction) %>%
  mutate(Sample_type = ifelse(grepl('GL',sample), 'blood',
                              ifelse(grepl('LN',sample),'lymph node','tumour'))) %>%
  filter(Sample_type != 'lymph node') %>%
  mutate(patient = substr(sample, 1,8)) %>% 
  left_join(age.df, 'patient')

# LUAD/LUSC
blood.plot.data2 <- response.plot.data %>%
  filter(Sample_type == 'blood') %>%
  filter(Histology %in% c('Squamous cell carcinoma', 'Invasive adenocarcinoma')) %>%
  mutate(Histology = ifelse(Histology == 'Invasive adenocarcinoma', 'LUAD','LUSC')) 

wilcox_test(blood.plot.data2, TCRA.tcell.fraction ~ Histology)
ES.plot2 <- wilcox_effsize(blood.plot.data2, TCRA.tcell.fraction ~ Histology)$effsize[[1]]

pdf('plots/ext3/blood_hist.pdf', height = 2, width = 1.75)
blood.plot.data2  %>%
  ggplot(aes(Histology, TCRA.tcell.fraction)) +
  xlab('') + ylab('Blood TCRA T cell fraction') +
  geom_boxplot() + geom_sina(size = 0.5) +
  stat_compare_means(size = 2) + 
  annotate('text', x =1.5 , y=0.3 ,
           label = paste0('ES = ',format(ES.plot2,digits = 2)),size = 2) +
  theme_bw() + theme(text = element_text(size = 7)) + ggtitle('')
dev.off()

# Figure 2a
# Sex
blood.plot.data3 <-response.plot.data %>%
  filter(Sample_type == 'blood')
ES.plot3 <- wilcox_effsize(blood.plot.data3, TCRA.tcell.fraction ~ Gender)$effsize[[1]]

pdf('plots/fig2/blood_sex.pdf', height = 2, width = 1.75)
blood.plot.data3 %>%
  ggplot(aes(Gender, TCRA.tcell.fraction)) +
  xlab('') + ylab('Blood TCRA T cell fraction') +
  geom_boxplot() + geom_sina(size = 0.5) +
  stat_compare_means(size = 2) +
  annotate('text', x =1.5 , y=0.3 ,
           label = paste0('ES = ',format(ES.plot3,digits = 2)),size = 2) +
  theme_bw() + theme(text = element_text(size = 7)) + ggtitle('')
dev.off()

# Tumour infiltration
tumour.infiltration.df <- tx100.tcra.scores %>%
  mutate(patient = substr(sample, 1,8)) %>%
  filter(!germline) %>%
  mutate(LN = grepl('LN', sample)) %>%
  filter(!LN) %>%
  group_by(patient) %>%
  summarise(mean.infiltration = mean(TCRA.tcell.fraction, na.rm = TRUE))

blood.plot.data4 <- response.plot.data %>%
  filter(Sample_type == 'blood') %>%
  left_join(tumour.infiltration.df, 'patient') 

pdf('plots/fig2/blood_infilt.pdf', height = 2, width = 2)
blood.plot.data4 %>%
  ggplot(aes(mean.infiltration, TCRA.tcell.fraction)) +
  xlab('Mean tumour TCRA fraction') + ylab('Blood TCRA T cell fraction ') +
  geom_point(size = 0.5) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size = 2) + 
  geom_smooth(method = 'lm') +
  theme_bw() + theme(text = element_text(size = 7)) + ggtitle('')
dev.off()

# Extended Figure 3b TCGA blood ----
TCGA.scores <- readRDS('data/TCGA_GC_CI_final_20210520.RDS')
TCGA.scores2 <- TCGA.scores  %>%
  mutate(TCRA.tcell.fraction = ifelse(TCRA.tcell.fraction.gc.smt < 0, 0,
                                      ifelse(TCRA.tcell.fraction.gc.smt.lwr > maxPossible, maxPossible, TCRA.tcell.fraction.gc.smt))) %>%
  select(case_id, sample_type, TCRA.tcell.fraction) %>%
  mutate(sample_type = as.character(sample_type)) %>%
  mutate(sample_type = ifelse(sample_type == 'Solid Tissue Normal','tissue.normal',
                              ifelse(sample_type == 'Blood Derived Normal','blood.normal',
                                     ifelse(sample_type == 'Primary Tumor','primary.tumour', sample_type)))) %>%
  filter(sample_type %in% c('primary.tumour','tissue.normal','blood.normal')) %>%
  group_by(case_id, sample_type) %>%
  summarise(TCRA.tcell.fraction = mean(TCRA.tcell.fraction, na.rm = TRUE)) %>%
  spread(sample_type, TCRA.tcell.fraction)

TCGA.scores.clin <- TCGA.scores %>%
  select(case_id, ascat_sex, age, tcga) %>% distinct()

TCGA.scores2 <- TCGA.scores2 %>%
  left_join(TCGA.scores.clin, 'case_id')

# LUAD/LUSC
blood.tcga.data1 <- TCGA.scores2 %>%
  filter(!is.na(tcga)) %>%
  mutate(tcga = toupper(tcga)) %>%
  filter(!is.na(blood.normal)) %>%
  as.data.frame()
ES.plot4<- wilcox_effsize(blood.tcga.data1 , blood.normal ~ tcga)$effsize[[1]]

pdf('plots/ext3/blood_hist_tcga.pdf', height = 2, width = 1.75)
blood.tcga.data1 %>%
  ggplot(aes(tcga, blood.normal)) +
  xlab('') + ylab('Blood TCRA T cell fraction') +
  geom_boxplot(outlier.size = 0, outlier.alpha = 0) + geom_sina(size = 0.25, alpha = 0.5) +
  stat_compare_means(size = 2) +
  annotate('text', x =1.5 , y=0.5 ,
           label = paste0('ES = ',format(ES.plot4,digits = 2)),size = 2) +
  theme_bw() + theme(text = element_text(size = 7)) + ggtitle('')
dev.off()

# SEX
blood.tcga.data2 <- TCGA.scores2 %>%
  filter(!is.na(ascat_sex)) %>%
  mutate(ascat_sex = ifelse(ascat_sex == 'XY','Male','Female')) %>%
  as.data.frame()
ES.plot5<- wilcox_effsize(blood.tcga.data1 , blood.normal ~ ascat_sex)$effsize[[1]]

pdf('plots/ext3/blood_sex_tcga.pdf', height = 2, width = 1.75)
blood.tcga.data2 %>%
  ggplot(aes(ascat_sex, blood.normal)) + 
  xlab('') + ylab('Blood TCRA T cell fraction') +
  geom_boxplot(outlier.size = 0, outlier.alpha = 0) + geom_sina(size = 0.25, alpha = 0.5) +
  stat_compare_means(size = 2) +
  annotate('text', x =1.5 , y=0.5 ,
           label = paste0('ES = ',format(ES.plot5,digits = 2)),size = 2) +
  theme_bw() + theme(text = element_text(size = 7)) + ggtitle('')
dev.off()

# Tumour infiltration
pdf('plots/ext3/blood_tcga_infilt.pdf', height = 2, width = 2)
TCGA.scores2 %>%
  ggplot(aes(primary.tumour, blood.normal)) + 
  xlab('Mean tumour TCRA fraction') + ylab('Blood TCRA T cell fraction ') +
  geom_point(size = 0.5) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size = 2) + 
  geom_smooth(method = 'lm') +
  theme_bw() + theme(text = element_text(size = 7)) + ggtitle('')
dev.off()


# Extended figure 3d PNE summary ----
normal.oesophagus.df  <- readRDS('data/normal_oesophagus_20210525.RDS')

normal.oesophagus.df.cat <- normal.oesophagus.df %>%
  filter(phenotype == 'nondysplasia') %>%
  group_by(subject_id) %>%
  summarise(meanT = mean(T.cell.fraction.raw.gc.smt)) %>%
  mutate(InfiltrationCategory = ifelse(meanT < 0.001,'< 0.001',
                                       ifelse(meanT >= 0.001 & meanT < 0.01, '0.001 - 0.01',
                                              ifelse(meanT >= 0.01 & meanT < 0.05, '0.01 - 0.05','>= 0.05')))) %>%
  mutate(InfiltrationCategory = factor(InfiltrationCategory, levels = rev(c('< 0.001','0.001 - 0.01','0.01 - 0.05','>= 0.05'))))
# 0.001 is threshold here

pdf('plots/ext3d/pne_counts.pdf', height = 1.75, width = 2)
normal.oesophagus.df.cat %>% 
  ggplot(aes(x= InfiltrationCategory)) + geom_bar() +
  theme_bw() + coord_flip() + xlab('mean TCRA T cell fraction') + 
  ggtitle('PNE TCRA T cell fraction') + 
  theme(text = element_text(size = 7))
dev.off()

# Extended figure 3e PNE summary----

pdf('plots/ext3/pne_overview.pdf', height = 3, width = 6)
normal.oesophagus.df %>%
  filter(phenotype %in% c('nondysplasia','blood')) %>%
  filter(subject_id %in% sample.order) %>%
  mutate(subject_id = factor(subject_id, levels = sample.order)) %>%
  mutate(phenotype = ifelse(phenotype == 'nondysplasia','PNE',phenotype)) %>%
  ggplot(aes(subject_id, T.cell.fraction.raw.gc.smt)) +
  geom_point(aes(col = phenotype), size = 0.5) + 
  #coord_flip() + 
  theme_bw() + xlab('') + ylab('TCRA T cell fraction') +
  ggtitle('') + theme(text = element_text(size = 7),legend.position = 'top',
                      axis.text.x = element_text(angle = 60, vjust = 0.5))
dev.off()

# Figure 2c PNE blood influence ----

samps1 <- sample.order[seq(30)]
samps2 <- sample.order[-seq(30)]

pne.data1 <- normal.oesophagus.df %>%
  filter(phenotype %in% c('nondysplasia','blood')) %>%
  filter(subject_id %in% sample.order) %>%
  mutate(zeroInfiltrationStatus = ifelse(subject_id %in% samps1,'TRUE','FALSE')) %>%
  group_by(subject_id,zeroInfiltrationStatus) %>%
  filter(phenotype == 'blood') %>%
  summarise(mean_infiltration = mean(T.cell.fraction.raw.gc.smt)) 

ES.plot8<- wilcox_effsize(as.data.frame(pne.data1) , mean_infiltration~ zeroInfiltrationStatus)$effsize[[1]]

pdf('plots/fig2/pne_inf.pdf', height = 2, width = 1.2)
pne.data1 %>%
  ggplot(aes(zeroInfiltrationStatus,mean_infiltration)) + 
  ylab('Blood TCRA T cell fraction') + xlab('T cell infiltration in PNE') +
  geom_boxplot(outlier.size = 0.5) + geom_sina(size = 0.5) +
  stat_compare_means(size = 2) +
  annotate('text', x =1.5 , y=0.5 ,
           label = paste0('ES = ',format(ES.plot8,digits = 2)),size = 2) +
  theme_bw() + theme(text = element_text(size = 7)) + ggtitle('')
dev.off()

# Extended data figure 3f,g ----

pne.data1 %>%
  ggplot(aes(zeroInfiltrationStatus,mean_infiltration)) +
  geom_boxplot() + stat_compare_means() +
  ylab('Blood TCRA T cell fraction') + xlab('T cell infiltration in PNE') + theme_bw()

linear.model.oesophagus.df <- normal.oesophagus.df %>%
  filter(phenotype %in% c('nondysplasia','blood')) %>%
  filter(subject_id %in% sample.order) %>%
  group_by(subject_id, phenotype) %>%
  summarise(TcellFraction = mean(T.cell.fraction.raw.gc.smt),
            age_years = unique(age_years), sex = unique(sex),
            alcohol = unique(alcohol), smoking = unique(smoking),
            ESCC_risk = unique(ESCC_risk), ESCC_HGD_status = unique(ESCC_HGD_status)) %>%
  spread(phenotype, TcellFraction) %>%
  mutate(Infiltration = nondysplasia > 0.001) %>%
  mutate(age_cat = ifelse(age_years < 60, '< 60','>= 60')) %>%
  mutate(ESCC_HGD_status = gsub(' \\(after ER\\)','',ESCC_HGD_status)) %>%
  filter(!is.na(blood)) %>%
  ungroup() %>%
  mutate(blood = scale(blood), age_years = scale(age_years)) %>%
  rename(Age = age_years, Smoke = smoking, PNE.infiltration = Infiltration, Sex = sex) %>%
  mutate(Sex = ifelse(Sex == 'M','Male','Female')) %>%
  mutate(Smoke = factor(Smoke, levels = c('Non or light smoker','Heavy smoker')))

summary(lm(blood ~ age_cat + sex + smoking + Infiltration + ESCC_HGD_status,
           data =linear.model.oesophagus.df))

# linear.model.oesophagus.df2 <- linear.model.oesophagus.df %>%
pdf('plots/ext3/pne_lm.pdf', height = 1.7, width = 2)
plot_model(lm(blood ~ Age + Sex + Smoke + PNE.infiltration  + ESCC_HGD_status,
              data =linear.model.oesophagus.df),dot.size =0.5, line.size = 0.5) + theme(text = element_text(size = 6)) + 
  ggtitle('Linear model for Blood TCRA fraction') 
dev.off()


nondysplaysia.lm.df.scaled <-
  nondysplaysia.lm.df %>%
  mutate(T.cell.fraction.raw.gc.smt = scale(T.cell.fraction.raw.gc.smt),
         Size  = scale(Size), mut_num = scale(mut_num), mut_num_dr = scale(mut_num_dr),
         TCF_max = scale(TCF_max), TCF_max_dr = scale(TCF_max_dr),  mut_num_dr_density = scale( mut_num_dr_density),
         mut_num_density = scale(mut_num_density), `Number of NOTCH1 mutation` = scale(`Number of NOTCH1 mutation`),
         `Number of TP53 mutations` = scale(`Number of TP53 mutations`), `Number of PPM1D mutations`  = scale(`Number of PPM1D mutations` ),
         TCF_max_TP53 = scale(TCF_max_TP53), TCF_max_NOTCH1 = scale(TCF_max_NOTCH1), TCF_max_PPM1D = scale(TCF_max_PPM1D),
         SignatureA_number= scale(SignatureA_number ), SignatureB_number = scale(SignatureB_number), 
         SignatureC_number = scale(SignatureC_number), SignatureD_number = scale(SignatureD_number))

genomic.lm.pne <- lm(T.cell.fraction.raw.gc.smt ~ Size + mut_num + mut_num_dr + TCF_max + 
                       TCF_max_dr + mut_num_dr_density +mut_num_density + `Number of NOTCH1 mutation` + 
                       `Number of TP53 mutations` + `Number of PPM1D mutations` + 
                       TCF_max_TP53 + TCF_max_NOTCH1 + TCF_max_PPM1D + SignatureA_number + SignatureB_number + 
                       SignatureC_number + SignatureD_number,
                     nondysplaysia.lm.df.scaled )

pdf('~/plots/ext3/pne_lm_genomic.pdf', height = 2.5, width = 2)
plot_model(lm(T.cell.fraction.raw.gc.smt ~ Size + mut_num + mut_num_dr + TCF_max + 
                TCF_max_dr + mut_num_dr_density +mut_num_density + `Number of NOTCH1 mutation` + 
                `Number of TP53 mutations` + `Number of PPM1D mutations` + 
                TCF_max_TP53 + TCF_max_NOTCH1 + TCF_max_PPM1D + SignatureA_number + SignatureB_number + 
                SignatureC_number + SignatureD_number,
              nondysplaysia.lm.df.scaled ),dot.size =0.5, line.size = 0.5) + theme(text = element_text(size = 6)) + 
  ggtitle('Linear model for TCRA T cell fraction for PNE from genomics') 
dev.off()


# Figure 2c Kraken ----

blood.kraken2 <- readRDS('data/tcga_kraken_blood_lm.RDS')

blood.kraken3 <- blood.kraken2 %>%
  ungroup() %>%
  mutate(age_at_diagnosis = scale(age_at_diagnosis)) %>%
  mutate(TCRA.tcell.fraction = scale(TCRA.tcell.fraction)) 

blood.kraken.lm <- lm(TCRA.tcell.fraction~investigation + gender + age_at_diagnosis + kraken.mean.cpm.norm.cat,
                      blood.kraken3) 
plot_model(blood.kraken.lm )
summary(blood.kraken.lm )

ES.plot6<- wilcox_effsize(as.data.frame(blood.kraken2) , TCRA.tcell.fraction ~ kraken.mean.cpm.norm.cat)$effsize[[1]]

pdf('plots/fig2/blood_kraken_tcga.pdf', height = 2, width = 1.75)
blood.kraken2 %>%
  ggplot(aes(kraken.mean.cpm.norm.cat, TCRA.tcell.fraction)) + 
  xlab('KRAKEN microbiome category') + ylab('Blood TCRA T cell fraction') +
  geom_boxplot() + geom_sina(size = 0.5) +
  stat_compare_means(size = 2) +
  annotate('text', x =1.5 , y=0.5 ,
           label = paste0('ES = ',format(ES.plot6,digits = 2)),size = 2) +
  theme_bw() + theme(text = element_text(size = 7)) + ggtitle('')
dev.off()

# Extended data figure 3h,i,j,k Kraken ----
kraken.df.summary2 <- readRDS('data/kraken_summary.RDS')

kraken.tumour.median <- kraken.df.summary2 %>%
  dplyr::filter(sample_type_key %in% c('tumour')) %>% ungroup() %>%
  dplyr::select(kraken.total.counts) %>% `[[`(1) %>% median(na.rm = TRUE)


tumour.kraken3  <- kraken.df.summary2 %>%
  dplyr::filter(sample_type_key %in% c('tumour')) %>%
  dplyr::filter(!is.na(kraken.total.counts)) %>%
  mutate(kraken.mean.cpm.norm.cat = ifelse(kraken.total.counts < kraken.tumour.median,'low','high')) 
ES.plot7<- wilcox_effsize(as.data.frame(tumour.kraken3) ,
                          TCRA.tcell.fraction ~ kraken.mean.cpm.norm.cat)$effsize[[1]]

pdf('plots/ext3/tumour_kraken_tcga.pdf', height = 2, width = 1.75)
tumour.kraken3 %>%
  ggplot(aes(kraken.mean.cpm.norm.cat, TCRA.tcell.fraction)) + 
  xlab('Kraken microbiome category') + ylab('Tumour TCRA T cell fraction') +
  geom_boxplot(outlier.size = 0, outlier.alpha = 0) + geom_sina(size = 0.25, alpha = 0.5) +
  stat_compare_means(size = 2) +
  annotate('text', x =1.5 , y=0.5 ,
           label = paste0('ES = ',format(ES.plot7,digits = 2)),size = 2) +
  theme_bw() + theme(text = element_text(size = 7)) + ggtitle('')
dev.off()


all.kraken.data <- readRDS('data/tcga_all_kraken_data.RDS')
# Uses Kraken-TCGA-Voom-SNM-Most-Stringent-Filtering-Data.csv and must have > 1000 raw reads across all the lung samples - 59 species remaining

# Download this file from: ftp://ftp.microbio.me/pub/cancer_microbiome_analysis/.
kraken.data.raw <- read_csv('data/Kraken-TCGA-Raw-Data-17625-Samples.csv') %>%
  dplyr::filter(X1 %in% all.kraken.data$X1)

all.kraken.blood <- all.kraken.data %>% filter(sample_type_key == 'blood')
all.kraken.normal <- all.kraken.data %>% filter(sample_type_key == 'normal')
all.kraken.tumour <- all.kraken.data %>% filter(sample_type_key == 'tumour')
all.kraken.lusc.tumour <- all.kraken.data %>% filter(investigation == 'TCGA-LUSC') %>% 
  filter(sample_type_key == 'tumour')
all.kraken.lusc.normal <- all.kraken.data %>% filter(investigation == 'TCGA-LUSC') %>%
  filter(sample_type_key == 'normal')
all.kraken.lusc.blood <- all.kraken.data %>% filter(investigation == 'TCGA-LUSC') %>% 
  filter(sample_type_key == 'blood')
all.kraken.luad.tumour <- all.kraken.data %>% filter(investigation == 'TCGA-LUAD') %>% 
  filter(sample_type_key == 'tumour')
all.kraken.luad.normal <- all.kraken.data %>% filter(investigation == 'TCGA-LUAD') %>%
  filter(sample_type_key == 'normal')
all.kraken.luad.blood <- all.kraken.data %>% filter(investigation == 'TCGA-LUAD') %>% 
  filter(sample_type_key == 'blood')

blood.cor.vector<- seq(59)
tumour.cor.vector<- seq(59)
normal.cor.vector<- seq(59)
blood.cor.vector<- seq(59)
luad.tumour.cor.vector<- seq(59)
luad.normal.cor.vector<- seq(59)
luad.blood.cor.vector<- seq(59)
lusc.tumour.cor.vector<- seq(59)
lusc.normal.cor.vector<- seq(59)
lusc.blood.cor.vector<- seq(59)

for(i in 2:60){
  index <- i-1
  blood.cor.vector[index] <- cor.test(all.kraken.blood[[i]], all.kraken.blood$TCRA.tcell.fraction, method = 'spearman')$p.value
  tumour.cor.vector[index] <- cor.test(all.kraken.tumour[[i]], all.kraken.tumour$TCRA.tcell.fraction, method = 'spearman')$p.value
  normal.cor.vector[index] <- cor.test(all.kraken.normal[[i]], all.kraken.normal$TCRA.tcell.fraction, method = 'spearman')$p.value
  
  luad.blood.cor.vector[index] <- cor.test(all.kraken.luad.blood[[i]], all.kraken.luad.blood$TCRA.tcell.fraction, method = 'spearman')$p.value
  luad.tumour.cor.vector[index] <- cor.test(all.kraken.luad.tumour[[i]], all.kraken.luad.tumour$TCRA.tcell.fraction, method = 'spearman')$p.value
  luad.normal.cor.vector[index] <- cor.test(all.kraken.luad.normal[[i]], all.kraken.luad.normal$TCRA.tcell.fraction, method = 'spearman')$p.value
  
  lusc.blood.cor.vector[index] <- cor.test(all.kraken.lusc.blood[[i]], all.kraken.lusc.blood$TCRA.tcell.fraction, method = 'spearman')$p.value
  lusc.tumour.cor.vector[index] <- cor.test(all.kraken.lusc.tumour[[i]], all.kraken.lusc.tumour$TCRA.tcell.fraction, method = 'spearman')$p.value
  lusc.normal.cor.vector[index] <- cor.test(all.kraken.lusc.normal[[i]], all.kraken.lusc.normal$TCRA.tcell.fraction, method = 'spearman')$p.value
  
}

all.logp.df <- data.frame(species = colnames(all.kraken.data)[c(2:60)],
                          blood.logp = -log10(blood.cor.vector),
                          normal.logp = -log10(normal.cor.vector),
                          tumour.logp = -log10(tumour.cor.vector))

luad.logp.df <- data.frame(species = colnames(all.kraken.data)[c(2:60)],
                           luad.blood.logp = -log10(luad.blood.cor.vector),
                           luad.normal.logp = -log10(luad.normal.cor.vector),
                           luad.tumour.logp = -log10(luad.tumour.cor.vector))

lusc.logp.df <- data.frame(species = colnames(all.kraken.data)[c(2:60)],
                           lusc.blood.logp = -log10(lusc.blood.cor.vector),
                           lusc.normal.logp = -log10(lusc.normal.cor.vector),
                           lusc.tumour.logp = -log10(lusc.tumour.cor.vector))

pdf('plots/ext3/luad_all_hits.pdf', height = 1.5, width = 3.5)
luad.logp.df %>%
  gather('comparison','minus_logp',-species) %>%
  filter(comparison != 'luad.normal.logp') %>%
  mutate(comparison = ifelse(comparison == 'luad.blood.logp','blood','tumour')) %>%
  ggplot(aes(species, minus_logp)) + 
  geom_point(aes(col = comparison), size = 0.5) + 
  geom_hline(yintercept = -log10(0.05/(length(blood.cor.vector)*2)),
             colour = 'red') + 
  theme_bw() + ggtitle('LUAD') + theme(axis.text.x = element_blank(), legend.key.size = unit(0.25,'cm')) +
  theme(text = element_text(size = 7))
dev.off()

# Only 18 lusc blood samples so remove as there is not any power
pdf('plots/ext3/lusc_all_hits.pdf', height = 1.5, width = 3.5)
lusc.logp.df %>%
  dplyr::select(species, lusc.blood.logp, lusc.tumour.logp) %>%
  gather('comparison','minus_logp',-species) %>%
  filter(comparison != 'luad.normal.logp') %>%
  mutate(comparison = ifelse(comparison == 'lusc.blood.logp','blood','tumour')) %>%
  ggplot(aes(species, minus_logp)) + 
  geom_point(aes(col = comparison), size = 0.5) + 
  geom_hline(yintercept = -log10(0.05/(length(blood.cor.vector)*2)),
             colour = 'red') + 
  theme_bw() + ggtitle('LUSC') + theme(axis.text.x = element_blank(), legend.key.size = unit(0.25,'cm')) +
  theme(text = element_text(size = 7))
dev.off()

# Load raw.counts.df....
raw.counts.df <- kraken.data.raw %>%
  select(X1, 
         Williamsia.raw.count = k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Corynebacteriales.f__Williamsiaceae.g__Williamsia,
         Paeniclostridium.raw.count = k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Peptostreptococcaceae.g__Paeniclostridium)


pdf('plots/ext3/luad_williamsia_tcga.pdf', height = 2, width = 2.5)
all.kraken.luad.tumour %>%
  left_join(raw.counts.df, 'X1') %>%
  mutate(Reads = ifelse(Williamsia.raw.count != 0, '> 0','0')) %>%
  ggplot(aes(k__Bacteria.p__Actinobacteria.c__Actinobacteria.o__Corynebacteriales.f__Williamsiaceae.g__Williamsia, TCRA.tcell.fraction)) +
  geom_point(aes(col = Reads), size = 0.5) + 
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size =2) + 
  xlab('Williamsia (log-cpm post SNM)')  + ylab('TCRA T cell fraction') +
  theme_bw() + ggtitle('TCGA LUAD tumours')  + theme(text = element_text(size = 7))
dev.off()

pdf('plots/ext3/lusc_Paeniclostridium_tcga.pdf', height = 2, width = 1.75)
all.kraken.lusc.tumour %>%
  left_join(raw.counts.df, 'X1') %>%
  ggplot(aes(k__Bacteria.p__Firmicutes.c__Clostridia.o__Clostridiales.f__Peptostreptococcaceae.g__Paeniclostridium, TCRA.tcell.fraction)) +
  geom_point(size = 0.5) + stat_cor(method = 'spearman', cor.coef.name = 'rho', size = 2) + 
  xlab('Paeniclostridium (log-cpm post SNM)')  + 
  ylab('TCRA T cell fraction') + theme_bw() + ggtitle('TCGA LUSC tumours')  + theme(text = element_text(size = 7))
dev.off()


# Pan cancer setup ----
load('data/pairwise_immune_copy_20210208.RData')
tx100.tcra.scores <- readRDS('data/Tx100_TCRA_scores_gc_update_20210630.RDS')
mets.tcell.output.final <- readRDS('data/MetsPrimary_TCRA_scores.CI.final.gc.20210517.RDS')
cytoband.df <- readRDS('data/20210214_gain_loss_tot_raw_cpn_by_cytoband.rds')


# Fix scores - if implausible -> NA 
tx100.tcra.scores <- tx100.tcra.scores %>%
  mutate(TCRA.tcell.fraction.gc = ifelse(TCRA.tcell.fraction.gc < 0, 0, TCRA.tcell.fraction.gc)) %>% 
  mutate(TCRA.tcell.fraction.gc = ifelse(TCRA.tcell.fraction.gc.lwr.95 > maxPossible, 
                                         T.cell.fraction.raw.gc.smt, TCRA.tcell.fraction.gc)) 

mets.tcell.output.final <- mets.tcell.output.final %>%
  mutate(TCRA.tcell.fraction.gc.smt = ifelse(TCRA.tcell.fraction.gc.smt < 0, 0,
                                             TCRA.tcell.fraction.gc.smt)) %>% 
  mutate(TCRA.tcell.fraction.gc.smt = ifelse(TCRA.tcell.fraction.gc.smt.lwr > maxPossible,
                                             T.cell.fraction.raw.gc.smt, TCRA.tcell.fraction.gc.smt))


group1.tcell.tx100 <- tx100.tcra.scores %>%
  select(group1 = sample, group1.tcell.fraction = TCRA.tcell.fraction.gc)
group2.tcell.tx100 <- tx100.tcra.scores %>%
  select(group2 = sample, group2.tcell.fraction = TCRA.tcell.fraction.gc)
group1.tcell <- mets.tcell.output.final %>%
  select(group1 = sample, group1.tcell.fraction = TCRA.tcell.fraction.gc.smt)
group2.tcell <- mets.tcell.output.final %>%
  select(group2 = sample, group2.tcell.fraction = TCRA.tcell.fraction.gc.smt)


a.df1 <- mets.pairwise.het %>%
  select(group_id, comparisonstatus, histology, data_set_name, range)
a.df2 <- tx100.pairwise.het %>%
  mutate(comparisonstatus = 'PRI-PRI') %>%
  select(group_id, comparisonstatus, histology, data_set_name, range)

all.ranges <- rbind(a.df2, a.df1) 

group1.tcell.all <- rbind(group1.tcell.tx100, group1.tcell)
group2.tcell.all <- rbind(group2.tcell.tx100, group2.tcell)

pairwise.everrything.df <- readRDS('data/pan_cancer_pairwise_all.RDS')

all.patients <- mets.tcell.output.final %>%
  filter(!germline) %>%
  filter(!is.na(TCRA.tcell.fraction.gc.smt)) %>%
  filter(sample_timepoint %in% c('PRI','PRI_POST')) %>% 
  select(patient_id) %>% `[[`(1)

multiple.primary.samples <- names(which(summary(factor(all.patients), maxsum = 1000) > 1))

load('data/tracerx100.clinical.RData')

tx100.hist.df <- clinical.data.tx100 %>%
  select(TRACERxID, Histology) %>%
  mutate(Histology = ifelse(Histology == 'Invasive adenocarcinoma','LUAD',
                            ifelse(Histology == 'Squamous cell carcinoma','LUSC','LUNG-OTHER'))) %>%
  mutate(data_set_name = 'TX100') %>%
  mutate(TRACERxID = gsub('LTX0','LTX', TRACERxID))

tx100.tcra.scores <- tx100.tcra.scores %>%
  mutate(germline = grepl('GL',sample)) %>%
  mutate(TRACERxID = substr(sample,3,8)) %>%
  filter(!germline) %>%
  select(TRACERxID, sample, TCRA.tcell.fraction.gc, 
         TCRA.tcell.fraction.gc.lwr.95, TCRA.tcell.fraction.gc.upr.95)

tx100.tcra.scores.hetero <- tx100.tcra.scores %>%
  left_join(tx100.hist.df, 'TRACERxID') %>%
  filter(!is.na(TCRA.tcell.fraction.gc)) %>%
  mutate(LN = grepl('LN', sample)) %>%
  filter(!LN)

# Extended Figure 4a ----

grid.plot.df <- mets.tcell.output.final %>% 
  filter(!is.na(TCRA.tcell.fraction.gc.smt)) %>%
  filter(patient_id %in% multiple.primary.samples) %>%
  filter(!germline) %>%
  dplyr::select(sample, patient_id, data_set_name, histology, TCRA.tcell.fraction.gc.smt,
                TCRA.tcell.fraction.gc.smt.lwr, TCRA.tcell.fraction.gc.smt.upr) %>%
  arrange(desc(TCRA.tcell.fraction.gc.smt)) 




# Add in Tx100
tx100.df.add <- tx100.tcra.scores.hetero %>%
  mutate(data_set_name = 'Tx100') %>%
  select(sample, patient_id = TRACERxID, data_set_name,
         histology = Histology, TCRA.tcell.fraction.gc.smt = TCRA.tcell.fraction.gc,
         TCRA.tcell.fraction.gc.smt.lwr = TCRA.tcell.fraction.gc.lwr.95,
         TCRA.tcell.fraction.gc.smt.upr = TCRA.tcell.fraction.gc.upr.95)

grid.plot.df <-rbind(grid.plot.df, tx100.df.add) %>%
  arrange(desc(TCRA.tcell.fraction.gc.smt)) 

TCRA.mean <- mean(grid.plot.df$TCRA.tcell.fraction.gc.smt)

mult.primary.heterogeneity.summary <- mets.tcell.output.final %>%
  filter(!germline) %>%
  filter(patient_id %in% multiple.primary.samples) %>%
  filter(!data_set_name %in% c('findlay','thomsen','murugaesu')) %>%
  group_by(patient_id) %>%
  summarise(histology = unique(histology), data_set_name = unique(data_set_name),
            min = min(TCRA.tcell.fraction.gc.smt, na.rm = TRUE), max = max(TCRA.tcell.fraction.gc.smt, na.rm = TRUE),
            sd = sd(TCRA.tcell.fraction.gc.smt, na.rm = TRUE),
            range = max(TCRA.tcell.fraction.gc.smt, na.rm = TRUE) - min(TCRA.tcell.fraction.gc.smt, na.rm = TRUE),
            hotRegions = sum(TCRA.tcell.fraction.gc.smt> TCRA.mean, na.rm = TRUE),
            coldRegions = sum(TCRA.tcell.fraction.gc.smt <= TCRA.mean, na.rm = TRUE),
            total = n())

mult.primary.heterogeneity.summary <- mult.primary.heterogeneity.summary  %>%
  mutate(category = ifelse(coldRegions == 0, 'hot',
                           ifelse(hotRegions == 0, 'cold', 'heterogeneous')))

grid.plot.df$x.loc <- c(seq_len(sum(grid.plot.df$TCRA.tcell.fraction.gc.smt !=0)),
                        rep(sum(grid.plot.df$TCRA.tcell.fraction.gc.smt !=0) + 1,
                            sum(grid.plot.df$TCRA.tcell.fraction.gc.smt ==0)))

category.df <- mult.primary.heterogeneity.summary %>%
  dplyr::select(patient_id, category)

grid.sample.order <- grid.plot.df %>%
  group_by(patient_id) %>%
  summarise(Score = median(TCRA.tcell.fraction.gc.smt),
            data_set_name = unique(data_set_name), histology = unique(histology)) %>%
  left_join(category.df, 'patient_id') %>%
  mutate(category = factor(category, levels = c('hot','heterogeneous','cold'))) %>%
  arrange(histology, category, desc(Score))

grid.sample.order$y.loc <- rev(seq(dim(grid.sample.order)[1]))

grid.sample.order2 <- grid.sample.order %>% dplyr::select(patient_id, y.loc)

grid.plot.df <- grid.plot.df %>%
  left_join(grid.sample.order2, 'patient_id')

# Change GMBLGG to GBM (mistake in annotation)
grid.plot.df <- grid.plot.df %>%
  mutate(histology = gsub(' $','',histology)) %>%
  mutate(histology = ifelse(histology == 'GBMLGG','GMB',histology)) %>%
  mutate(histology = factor(histology))

# Create histology colour scale
library(RColorBrewer)
# Define the number of colors you want
nb.cols <- 13
hist.colours <- colorRampPalette(brewer.pal(12, "Paired"))(nb.cols)
names(hist.colours) <- (levels(grid.plot.df$histology))

grid1.plot <- grid.plot.df  %>%
  ggplot(aes(x.loc, y.loc)) +
  geom_line(aes(group = patient_id), size = 0.1) +
  geom_point(aes(colour = histology), size = 0.25) + theme_bw() +
  scale_color_manual(values = hist.colours) +
  xlab('Region') + ylab('Patient')  + 
  theme(axis.text.x = element_blank(), axis.text.y = element_blank(),
        axis.ticks.x = element_blank(), axis.ticks.y = element_blank(),
        text = element_text(size=6),legend.position = 'none') +
  geom_vline(xintercept = which(grid.plot.df$TCRA.tcell.fraction.gc.smt < TCRA.mean)[1] - 0.5,
             col = 'red', linetype = 'dashed')

grid2.plot <- grid.plot.df %>%
  ggplot(aes(x.loc, TCRA.tcell.fraction.gc.smt)) +
  geom_col(aes(fill = histology)) + theme_bw() + 
  scale_fill_manual(values = hist.colours) +
  xlab('') + ylab('T cell fraction') + 
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = 'none',
        text = element_text(size=6))

# COMBINE
require(grid)
require(gridExtra)

plot_list_withoutLegends <- list(grid2.plot, grid1.plot)

#convert into ggplotGrob object
plot_list_ggplotGrop <- lapply(1:length(plot_list_withoutLegends), function(x){
  p <- plot_list_withoutLegends[[x]]
  if(x == 1){
    p <- ggplotGrob(p + theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")))
  } else if(x == length(plot_list_withoutLegends)) {
    p <- ggplotGrob(p + theme(plot.margin = unit(c(0, 0.5, 1, 0.5), "cm")))
  } else{
    p <- ggplotGrob(p + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")))
  }
  return(p)
})

#align widths 
all_widths <- lapply(plot_list_ggplotGrop, function(x) {x$widths})
plot_list_alignedWiths <- lapply(plot_list_ggplotGrop, function(x){
  x$widths <- do.call(unit.pmax, all_widths)
  return(x)
})

#set heights
g <- do.call(gtable_rbind, plot_list_alignedWiths)
id_panels_h <- unique(g$layout[g$layout$name=="panel", "t"])
g$heights[id_panels_h] <- grid::unit(c(1,3), "null")


pdf('plots/ext4/multi_primary_overview_without_legend_20210312_update.pdf', width = 3.1, height =4)
grid.arrange(g)
dev.off()

# Figure 2d Pan cancer heterogeneity ----
tx100.heterogeneity.summary <- tx100.tcra.scores.hetero %>%
  group_by(TRACERxID) %>%
  summarise(histology = unique(Histology), data_set_name = unique(data_set_name),
            min = min(TCRA.tcell.fraction.gc), max = max(TCRA.tcell.fraction.gc),
            sd = sd(TCRA.tcell.fraction.gc),
            range = max(TCRA.tcell.fraction.gc) - min(TCRA.tcell.fraction.gc),
            hotRegions = sum(TCRA.tcell.fraction.gc> TCRA.mean), 
            coldRegions = sum(TCRA.tcell.fraction.gc <= TCRA.mean),
            total = n())

tx100.heterogeneity.summary <- tx100.heterogeneity.summary %>%
  mutate(category = ifelse(coldRegions == 0, 'hot',
                           ifelse(hotRegions == 0, 'cold', 'heterogeneous')))

all.heterogeneity.summary <- tx100.heterogeneity.summary %>%
  dplyr::rename(patient_id = TRACERxID) %>%
  rbind(mult.primary.heterogeneity.summary)

all.heterogeneity.summary$histology <- gsub(' $','',all.heterogeneity.summary$histology)

total.hist <- all.heterogeneity.summary %>% group_by(histology) %>% summarise(total = n())

histology.order <- all.heterogeneity.summary %>%
  filter(histology %in% c('LUAD','LUSC','LUNG-OTHER','BLCA','BRCA ER+',
                          'BRCA HER2+','GBMLGG','KIRC')) %>%
  group_by(category, histology) %>% summarise(num=n()) %>%
  filter(category == 'heterogeneous') %>%
  left_join(total.hist, 'histology') %>% 
  mutate(percent = num/total) %>% arrange(desc(percent))  %>% 
  ungroup() %>%
  select(histology) %>% `[[`(1)

all.heterogeneity.summary %>%
  filter(histology %in% c('LUAD','LUSC','LUNG-OTHER','BLCA','BRCA ER+',
                          'BRCA HER2+','GBMLGG','KIRC')) %>%
  filter(!is.na(category)) %>% select(histology) %>% `[[`(1) %>% factor() %>%summary()

pdf('plots/fig2/multi_primary_category_no_legend.pdf', width = 3.5, height = 1.75)
all.heterogeneity.summary %>%
  filter(histology %in% c('LUAD','LUSC','LUNG-OTHER','BLCA','BRCA ER+',
                          'BRCA HER2+','GBMLGG','KIRC')) %>%
  filter(!is.na(category)) %>%
  mutate(histology = factor(histology, levels = histology.order)) %>%
  mutate(category = factor(category, levels = c('hot','cold','heterogeneous'))) %>%
  mutate(count = 1) %>%
  ggplot(aes(histology,count)) +
  geom_col(aes(fill = category), position = 'fill')  + 
  scale_fill_manual(values=c('#BD2131', '#30ABDF', '#df9b2f')) +
  theme_bw() + xlab('') + ylab('proportion') +
  ggtitle('Chi-squared test: p-value = 1.616e-07') +
  theme(text = element_text(size=7),
        legend.key.size = unit(0.3, 'cm'),
        axis.text.x = element_text(angle=45, hjust = 1),
        legend.position = 'none') +
  annotate("text", x = 1:8, y = rep(1.05,8), label = paste0('n = ',c(6,5,7,24,59,22,12,32)), size = 2)
dev.off()

all.heterogeneity.table <- all.heterogeneity.summary  %>%
  filter(histology %in% c('LUAD','LUSC','LUNG-OTHER','BLCA','BRCA ER+',
                          'BRCA HER2+','GBMLGG','KIRC')) %>%
  filter(!is.na(category)) %>%
  select(histology, category) %>% table()

all.heterogeneity.table/rowSums(all.heterogeneity.table)

chisq.test(all.heterogeneity.table)

# Figure 2e Pan cancer SCNA pairwise ----
multi.primary.samples <- unique(grid.plot.df$sample)
multi.primary.patients <- unique(grid.plot.df$patient_id)
tx100.patients <- unique(pairwise.everrything.df$patient_id)[grepl('LTX',unique(pairwise.everrything.df$patient_id))]

immuneDiff.df <- pairwise.everrything.df %>%
  filter(patient_id %in% c(multi.primary.patients, tx100.patients)) %>%
  mutate(immuneDiff = (abs(group1.tcell.fraction - group2.tcell.fraction)))

pair.mean <- mean(immuneDiff.df$immuneDiff, na.rm = TRUE)

median.tcell.diff <- pair.mean

pairwise.everything.df <- pairwise.everrything.df %>%
  filter(patient_id %in% c(multi.primary.patients, tx100.patients)) %>%
  mutate(immuneComparison = ifelse(group1.tcell.fraction < pair.mean & group2.tcell.fraction < pair.mean,
                                   'Cold-Cold',
                                   ifelse(group1.tcell.fraction >= pair.mean& group2.tcell.fraction >= pair.mean,
                                          'Hot-Hot','Cold-Hot'))) %>%
  mutate(immuneDiff = ifelse(abs(group1.tcell.fraction - group2.tcell.fraction) > median.tcell.diff,
                             'large','small'))

paired.immuneDiff <- pairwise.everything.df %>%
  filter(!is.na(immuneDiff)) %>%
  group_by(event.type, patient_id, immuneDiff) %>%
  summarise(prop_ab_in_one_diff_abs = mean(prop_ab_in_one_diff_abs, na.rm = TRUE),
            prop_ab_in_one = mean(prop_ab_in_one, na.rm = TRUE),
            tcell.diff = mean(abs(tcell.diff)))

hetero.patients <- paired.immuneDiff %>% filter(event.type == 'all') %>% ungroup()  %>% 
  select(patient_id)  %>% `[[`(1) %>% factor() %>% summary(maxsum=400)  


paired.plot.df <- paired.immuneDiff %>%
  filter(event.type %in% c('all','gain','loh_loss')) %>%
  filter(patient_id %in% names(which(hetero.patients == 2))) %>%
  mutate(immuneDiff = ifelse(immuneDiff == 'large', '>= 0.065', '< 0.065')) %>%
  mutate(event.type = ifelse(event.type == 'all','All events',
                             ifelse(event.type == 'gain','Gain events','Loss or LOH events')))


wilcox_effsize(ungroup(paired.plot.df[paired.plot.df$event.type == 'All events',]),
               prop_ab_in_one ~immuneDiff, paired = TRUE) # 0.347
wilcox_effsize(ungroup(paired.plot.df[paired.plot.df$event.type == 'Gain events',]),
               prop_ab_in_one ~immuneDiff, paired = TRUE)  # 0.318 
wilcox_effsize(ungroup(paired.plot.df[paired.plot.df$event.type == 'Loss or LOH events',]),
               prop_ab_in_one ~immuneDiff, paired = TRUE)  # 0.253

wilcox.es.df <- data.frame(event.type = c('All events','Gain events','Loss or LOH events'),
                           x = 1, y = 0.45, lab = paste0('Effect size = ', c(0.347, 0.318,0.253)))

pdf('plots/fig2/pairwise_scna_tcell.pdf', width = 3.5, height = 1.75)
paired.plot.df  %>%
  ggplot(aes(immuneDiff, prop_ab_in_one)) +
  geom_boxplot(outlier.size = 0.2) + geom_line(aes(group = patient_id), alpha = 0.2, size = 0.2) + 
  stat_compare_means(size = 1, paired = TRUE, label.y = c(0.6,0.5,0.4)) +
  # geom_sina() + 
  theme_bw() + ylab('Pairwise SCNA heterogeniety') + 
  facet_wrap(~event.type) + xlab('TCRA T cell fraction difference') + 
  theme(text = element_text(size = 5)) +
  geom_text(data = wilcox.es.df, mapping =  aes(x , y , label = lab), size = 1) +
  ggtitle('Multi-region primary tumours with both small and large pairwise TCRA T cell differences between regions (n = 76)')
dev.off()



# Extended Figure 4c and Figure 2f 12q24.31-32 ----

all.tcell.df <- group1.tcell.all %>%
  rename(region = group1, tcell.fraction = group1.tcell.fraction)
all.hist.df <- pairwise.everything.df %>%
  select(patient_id, histology) %>% distinct()


all.cytobands <- unique(cytoband.df$cytoband_name)
all.primary.samples <- unique(c(pairwise.everything.df$group2, pairwise.everything.df$group1))

cytoband.subclonality.df <- cytoband.df %>%
  filter(!is.na(cpn_call_vs_ploidy)) %>%
  filter(sample %in% all.primary.samples) %>% # take out all met data... 
  group_by(cytoband_name, patient_id) %>%
  summarise(num.regions = length(unique(sample)),
            loss.regions = sum(cpn_call_vs_ploidy %in% c('loss','deep_loss')),
            neutral.regions = sum(cpn_call_vs_ploidy == 'neutral'),
            gain.regions = sum(cpn_call_vs_ploidy %in% c('gain','amp'))) %>%
  mutate(loss.subclonal = ifelse(loss.regions < num.regions & loss.regions != 0, TRUE, FALSE)) %>%
  mutate(gain.subclonal = ifelse(gain.regions < num.regions & gain.regions != 0, TRUE, FALSE))

cytoband.subclonality.summary.df <- cytoband.subclonality.df %>%
  summarise(total.loss.subclonal = sum(loss.subclonal),
            total.gain.subclonal = sum(gain.subclonal)) 

# Need to get cytoband locations and make a genome plot - sort of like the CN plots I have made before
# wget -qO- "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" | gunzip -c > cytoband.hg19.txt
cytoband.loc.df <- read_tsv('data/cytoBand.hg19.txt',col_names = FALSE)

colnames(cytoband.loc.df) <- c('chr','start','stop','cytoband','X5')

cytoband.loc.df2 <- cytoband.loc.df %>%
  mutate(chr = gsub('chr','',chr)) %>%
  mutate(cytoband_name = paste0(chr, cytoband)) %>%
  select(chr, cytoband_name, start,stop)
# Add chromosome lengths to make
cytoband.max <- cytoband.loc.df2 %>%
  group_by(chr) %>% summarise(max.pos = max(stop))

cytoband.max <- cytoband.max %>%
  mutate(chr = ifelse(chr =='X','23',ifelse(chr == 'Y','24',chr))) %>%
  mutate(chr = as.numeric(chr)) %>%
  arrange(chr) %>%
  mutate(cumulative.len = cumsum(max.pos)) %>%
  mutate(chr = chr + 1) 

cytoband.max <-rbind(data.frame(chr = 1, max.pos = 0, cumulative.len = 0),
                     cytoband.max)
cytoband.max <- cytoband.max %>%
  filter(chr <=22)

cytoband.loc.df3 <- cytoband.loc.df2 %>%
  mutate(chr = ifelse(chr =='X','23',ifelse(chr == 'Y','24',chr))) %>%
  mutate(chr = as.numeric(chr)) %>%
  left_join(cytoband.max, 'chr') %>%
  mutate(start.adj = start + cumulative.len,
         stop.adj = stop + cumulative.len) %>%
  select(cytoband_name,chr, start, stop, start.adj, stop.adj )

# Combine + plot 
cytoband.subclonality.plot.df <- cytoband.subclonality.summary.df %>%
  left_join(cytoband.loc.df3, 'cytoband_name') %>%
  mutate(total.loss.subclonal = -total.loss.subclonal)

# Do test for every cytoband region above a certain threshold:
loss.cytobands <- cytoband.subclonality.summary.df %>% 
  filter(total.loss.subclonal >= 30) %>%
  select(cytoband_name) %>% `[[`(1)

gain.cytobands <- cytoband.subclonality.summary.df %>% 
  filter(total.gain.subclonal >= 30) %>%
  select(cytoband_name) %>% `[[`(1)

# group together by region
unique.loss.cytobands <- sapply(loss.cytobands,
                                function(x) strsplit(as.character(x),'\\.')[[1]][1]) %>%
  unique()
unique.gain.cytobands <- sapply(gain.cytobands,
                                function(x) strsplit(as.character(x),'\\.')[[1]][1]) %>%
  unique()

loss.cytobands.groups <- lapply(unique.loss.cytobands ,
                                function(x) loss.cytobands[grepl(paste0('^',x), loss.cytobands)])
gain.cytobands.groups <- lapply(unique.gain.cytobands ,
                                function(x) gain.cytobands[grepl(paste0('^',x), gain.cytobands)])


gain.loss.plot.list <- function(cytoband, alteration){
  
  alterations.1 <- cytoband.df %>%
    filter(sample %in% all.primary.samples) %>%
    filter(cytoband_name %in% cytoband) %>%
    filter(cpn_call_vs_ploidy %in% alteration)
  
  # Identify patients with subclonal gains
  alterations.1.patients <- unique(alterations.1$patient_id)
  alterations.1.regions <- unique(alterations.1$sample)
  
  alterations.1.subclonal.patients <- cytoband.df %>%
    filter(sample %in% all.primary.samples) %>%
    filter(patient_id %in% alterations.1.patients) %>%
    filter(!sample %in% alterations.1.regions) %>% ungroup() %>%
    select(patient_id) %>% `[[`(1) %>% unique()
  
  alterations.1.subclonal.regions <-  cytoband.df %>%
    filter(sample %in% all.primary.samples) %>%
    filter(patient_id %in% alterations.1.subclonal.patients) %>%
    filter(cytoband_name %in% cytoband) %>%
    filter(cpn_call_vs_ploidy %in% alteration) %>% ungroup() %>%
    select(sample) %>% `[[`(1) %>% unique()
  
  alterations.1.no.subclonal.regions <-  cytoband.df %>%
    filter(sample %in% all.primary.samples) %>%
    filter(patient_id %in% alterations.1.subclonal.patients) %>%
    filter(!sample %in% alterations.1.subclonal.regions) %>%ungroup() %>%
    select(sample) %>% `[[`(1) %>% unique()
  
  # Make data frame:
  alterations.1.df <- data.frame(region = c(alterations.1.no.subclonal.regions,
                                            alterations.1.subclonal.regions),
                                 alterations.1 = c(rep(FALSE, length(alterations.1.no.subclonal.regions)),
                                                   rep(TRUE, length(alterations.1.subclonal.regions))))
  
  alterations.1.df.summary <- alterations.1.df %>%
    mutate(patient_id = substr(region,1,8)) %>%
    left_join(all.tcell.df, 'region') %>%
    filter(!is.na(tcell.fraction)) %>% group_by(patient_id, alterations.1)  %>%
    summarise(tcell.fraction = mean(tcell.fraction, na.rm = TRUE))  
  
  alterations.1.df.summary %>%
    filter(patient_id %in% names(which(summary(factor(alterations.1.df.summary$patient_id)) == 2))) %>%
    left_join(all.hist.df, 'patient_id') %>%
    ggplot(aes(alterations.1, tcell.fraction)) +
    geom_boxplot() + geom_line(aes(colour = histology, group = patient_id)) +
    stat_compare_means(paired = TRUE, size = 2) + xlab(paste(paste0(cytoband, collapse = ','),'-',
                                                             paste0(alteration, collapse = ','))) %>%
    return()
}

cpn.plot.gains.list <- list()
for(i in seq_len(length(gain.cytobands.groups))){
  cpn.plot.gains.list[[i]] <- gain.loss.plot.list(gain.cytobands.groups[[i]],c('gain','amp'))
}
cpn.plot.gains.pvalues <- (sapply(cpn.plot.gains.list, function(x) ggplot_build(x)$data[[3]]$p.adj))
#
cpn.plot.loss.list <- list()
for(i in seq_len(length(loss.cytobands.groups))){
  cpn.plot.loss.list[[i]] <- gain.loss.plot.list(loss.cytobands.groups[[i]],c('loss','deep_loss'))
}

cpn.plot.loss.pvalues <- (sapply(cpn.plot.loss.list, function(x) ggplot_build(x)$data[[3]]$p.adj))

# Build data frame to plot log p values - also get location of each
get_cytoband_pos <- function(cytoband_regions){
  tmp.df <- cytoband.loc.df3 %>%
    filter(cytoband_name %in% cytoband_regions)
  return((min(tmp.df$start.adj) + max(tmp.df$stop.adj))/2)
}

chr.position <- sapply(seq(22), function(x){
  tmp.df <- cytoband.loc.df3 %>%
    filter(chr %in% x)
  return((min(tmp.df$start.adj) + max(tmp.df$stop.adj))/2)
})


logp.df <- rbind(data.frame(genome.pos = sapply(loss.cytobands.groups , get_cytoband_pos),
                            minus_logp = -log10(cpn.plot.loss.pvalues),
                            status = 'loss'),
                 data.frame(genome.pos = sapply(gain.cytobands.groups , get_cytoband_pos),
                            minus_logp = -log10(cpn.plot.gains.pvalues),
                            status = 'gain'))

logp.plot <- logp.df %>%
  ggplot(aes(genome.pos, minus_logp)) +
  geom_point(aes(colour = status), size = 0.25) +
  geom_hline(yintercept = -log10(0.05/(length(cpn.plot.loss.pvalues) + length(cpn.plot.gains.pvalues))),
             size = 0.25, colour = 'red') + 
  scale_color_manual(values = c(loss = 'darkblue', gain = 'darkred')) +
  geom_vline(data = cytoband.max, aes(xintercept = cumulative.len),
             linetype = 'dashed', size = 0.25) +
  theme_bw() + xlab('') +
  theme(text = element_text(size = 7), legend.position = 'none',
        axis.text.x = element_blank(), axis.ticks.x = element_blank()) 

cytoband.plot <- cytoband.subclonality.plot.df %>% 
  ggplot(aes(xmin = start.adj, xmax = stop.adj)) +
  geom_hline(yintercept = c(-30,30), size = 0.25) +
  geom_rect(aes(ymin = 0, ymax = total.gain.subclonal),fill = 'darkred') +
  geom_rect(aes(ymin = total.loss.subclonal, ymax = 0),fill = 'darkblue') +
  theme_bw() +xlab('') + ylab('Number of patients with subclonal gains/losses') +
  geom_vline(data = cytoband.max, aes(xintercept = cumulative.len), linetype = 'dashed', size = 0.25) +
  xlab('Genome position') +
  annotate('text', x = chr.position,
           y = rep(c(-50,50),11), label = cytoband.max$chr, size =2) +
  theme(text = element_text(size = 7),
        axis.text.x = element_blank(), axis.ticks.x = element_blank())

plot_list_withoutLegends <- list(logp.plot, cytoband.plot)

#convert into ggplotGrob object
plot_list_ggplotGrop <- lapply(1:length(plot_list_withoutLegends), function(x){
  p <- plot_list_withoutLegends[[x]]
  if(x == 1){
    p <- ggplotGrob(p + theme(plot.margin = unit(c(1, 0.5, 0, 0.5), "cm")))
  } else if(x == length(plot_list_withoutLegends)) {
    p <- ggplotGrob(p + theme(plot.margin = unit(c(0, 0.5, 1, 0.5), "cm")))
  } else{
    p <- ggplotGrob(p + theme(plot.margin = unit(c(0, 0.5, 0, 0.5), "cm")))
  }
  return(p)
})

#align widths 
all_widths <- lapply(plot_list_ggplotGrop, function(x) {x$widths})
plot_list_alignedWiths <- lapply(plot_list_ggplotGrop, function(x){
  x$widths <- do.call(unit.pmax, all_widths)
  return(x)
})

#set heights
g <- do.call(gtable_rbind, plot_list_alignedWiths)
id_panels_h <- unique(g$layout[g$layout$name=="panel", "t"])
g$heights[id_panels_h] <- grid::unit(c(1.5,3), "null")

pdf('plots/ext4/subclonal_cytoband_20210312_update.pdf', width = 5, height = 3)
grid.arrange(g)
dev.off()

gain.loss.df.make <- function(cytoband, alteration){
  
  alterations.1 <- cytoband.df %>%
    filter(sample %in% all.primary.samples) %>%
    filter(cytoband_name %in% cytoband) %>%
    filter(cpn_call_vs_ploidy %in% alteration)
  
  # Identify patients with subclonal gains
  alterations.1.patients <- unique(alterations.1$patient_id)
  alterations.1.regions <- unique(alterations.1$sample)
  
  alterations.1.subclonal.patients <- cytoband.df %>%
    filter(patient_id %in% alterations.1.patients) %>%
    filter(!sample %in% alterations.1.regions) %>% ungroup() %>%
    select(patient_id) %>% `[[`(1) %>% unique()
  
  alterations.1.subclonal.regions <-  cytoband.df %>%
    filter(patient_id %in% alterations.1.subclonal.patients) %>%
    filter(cytoband_name %in% cytoband) %>%
    filter(cpn_call_vs_ploidy %in% alteration) %>% ungroup() %>%
    select(sample) %>% `[[`(1) %>% unique()
  
  alterations.1.no.subclonal.regions <-  cytoband.df %>%
    filter(patient_id %in% alterations.1.subclonal.patients) %>%
    filter(!sample %in% alterations.1.subclonal.regions) %>%ungroup() %>%
    select(sample) %>% `[[`(1) %>% unique()
  
  # Make data frame:
  alterations.1.df <- data.frame(region = c(alterations.1.no.subclonal.regions,
                                            alterations.1.subclonal.regions),
                                 alterations.1 = c(rep(FALSE, length(alterations.1.no.subclonal.regions)),
                                                   rep(TRUE, length(alterations.1.subclonal.regions))))
  
  alterations.1.df.summary <- alterations.1.df %>%
    mutate(patient_id = substr(region,1,8)) %>%
    left_join(all.tcell.df, 'region') %>%
    filter(!is.na(tcell.fraction)) %>% group_by(patient_id, alterations.1)  %>%
    summarise(tcell.fraction = mean(tcell.fraction, na.rm = TRUE))  
  
  alterations.1.df.summary %>%
    filter(patient_id %in% names(which(summary(factor(alterations.1.df.summary$patient_id)) == 2))) %>%
    left_join(all.hist.df, 'patient_id') %>%
    ungroup() %>%
    return()
}


loss_12q24.df <- gain.loss.df.make(c('12q24.31','12q24.32'), c('loss','deep_loss'))


wilcox_test(loss_12q24.df, formula = tcell.fraction~alterations.1, paired = TRUE)
wilcox_effsize(loss_12q24.df, formula = tcell.fraction~alterations.1, paired = TRUE)
loss.cytobands.groups[[69]]

pdf('plots/fig2/loss_12q24.pdf', width = 2, height = 3)
cpn.plot.loss.list[[69]] +
  scale_color_manual(values = hist.colours) +
  theme_bw() +
  theme(text = element_text(size = 7)) +
  annotate('text',x=1, y =0.25, label = 'Effect size = 0.75', size = 2) +
  xlab('Subclonal 12q24.31-32 loss') +
  ylab('TCRA T cell fraction') + theme(legend.position = 'none')
dev.off()

# Extended Figure 4d SPPL3 volcano ----

tx100.subclonal.12q.volcano.df2 <- read_csv('data/_subclonal_subtypelimma_voom_results.csv')


library(EnhancedVolcano)

png('plots/supfig3volcano_loss_12q24.png', width = 4, height = 2.5, units = 'in',res = 150)
EnhancedVolcano(toptable = as.data.frame(tx100.subclonal.12q.volcano.df2),
                lab = as.data.frame(tx100.subclonal.12q.volcano.df2)$gene_id,
                x = 'logFC',
                transcriptLabhjust = c(1.2,1,1.6),
                transcriptLabvjust = c(1.2,-1,1.6),
                y = 'adj.P.Val',
                pCutoff = 0.05,
                FCcutoff = 0.25,
                title = NULL,
                subtitle = NULL,
                ylim = c(0,2.5),
                selectLab = c('SPPL3','ABCB9','OGFOD2'),
                legendPosition = 'topt',
                legendLabSize = 6,
                legendIconSize = 0.5,
                colConnectors = 'black',
                drawConnectors = TRUE,
                widthConnectors = 0.25, titleLabSize = 8,
                subtitleLabSize = 6,axisLabSize = 8,captionLabSize = 6,
                transcriptLabSize = 2)
dev.off()
