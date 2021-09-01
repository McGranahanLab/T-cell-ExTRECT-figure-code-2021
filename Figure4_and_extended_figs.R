###############################################################################
### T cell quantification from DNA sequencing predicts immunotherapy        ###
###       response                                                          ###
###                                                                         ###
### paper figure 4 + related extended figures                               ###
###                                                                         ###
### Author: Robert Bentham                                                  ###
### Date: 2021 July 30th                                                    ###
###############################################################################

library(tidyverse)
library(ggpubr)
library(meta)
library(ggforce)
library(xgboost)
library(ROCR)
library(pROC)
library(rstatix)
library(scales)

###############################################################################
### Functions                                                               ###
###############################################################################


create_mva_scaled <- function(columns_to_select, metrics_input, studies){
  mva_ind_set<- metrics_input[,colnames(metrics_input) %in% columns_to_select]
  mva_ind_set<-cbind(mva_ind_set,
                     metrics_input$case,
                     metrics_input$histology,
                     metrics_input$study,
                     metrics_input$response)
  names(mva_ind_set)[ncol(mva_ind_set)-3]<-"case"
  names(mva_ind_set)[ncol(mva_ind_set)-2]<-"histology"
  names(mva_ind_set)[ncol(mva_ind_set)-1]<-"study"
  names(mva_ind_set)[ncol(mva_ind_set)]<-"response"
  mva_ind_set$hist_study<-paste0(mva_ind_set$histology,"-",mva_ind_set$study)
  mva_ind_set<-mva_ind_set[!is.na(mva_ind_set$Clonal_TMB), ]
  mva_ind_set_all <- mva_ind_set
  mva_ind_set <- mva_ind_set[complete.cases(mva_ind_set),]
  mva_ind_set_scaled <- mva_ind_set[0,]
  
  mva_set<-mva_ind_set_all
  mva_set_scaled<-mva_set[0,]
  
  # Do not scale either TMB or gender
  selected.cols <- which(colnames(mva_set) %in% columns_to_select)
  do.not.scale <- which(colnames(mva_set) %in% c('TMB','Gender', 'tumour.tcell.fraction'))
  cols_to_scale <- setdiff(selected.cols, do.not.scale)
  #scale variables within each study
  for(j in 1:length(studies)){
    print(studies[j])
    #opt_tmb_auc<-c()
    #tmb_auc<-c()
    this_study<-mva_set[mva_set$hist_study==as.character(studies[j]),]
    for(k in 1:length(cols_to_scale)){
      this_study[cols_to_scale[k]]<-as.numeric(scale(this_study[cols_to_scale[k]]))}	
    mva_tmp<-mva_set_scaled
    mva_set_scaled<-rbind(mva_tmp,this_study)
  }
  mva_set_scaled$Response_binary <- 1
  mva_set_scaled$Response_binary[mva_set_scaled$response=="no_response"]<-0
  return(mva_set_scaled)
}

tuneplot <- function(x, probs = 1) {
  ggplot(x) +
    coord_cartesian(ylim = c(min(x$results$ROC), quantile(x$results$ROC, probs = probs))) + 
    theme_bw()
}

auc.pred.fun <- function(model, testData, cols_to_select, glm = FALSE){
  selected.cols <- which(colnames(testData) %in% cols_to_select)
  do.not.scale <- which(colnames(testData) %in% c('TMB'))
  col1 <- setdiff(selected.cols, do.not.scale)
  if(glm){
    testData1 <- testData[,col1]
  }else{
    testData1 <- as.matrix(testData[,col1])
  }
  if(glm){
    pred <- predict(model,
                    testData1, type = 'response')
    # pred <- ifelse(as.numeric(pred) > 0.5, 1,0)
    no.na.loc <- !is.na(pred)
    pr <- prediction(as.numeric(pred)[no.na.loc], testData$Response_binary[no.na.loc])
  }else{
    pred <- predict(model,
                    testData1)
    pr <- prediction(as.numeric(pred), testData$Response_binary)
  }
  
  auc <- performance(pr, measure = 'auc')
  auc_val <- auc@y.values[[1]]
  auc_tpr <- performance(pr, "tpr","fpr")
  return(list(auc, auc_val, auc_tpr, pred))
}

auc.pred.tmb.fun <- function(testData){
  pred <- ifelse(testData[,1]>=199,1,0)
  pr2 <- prediction(pred,
                    testData$Response_binary)
  auc <- performance(pr2, measure = "auc")
  auc_val<-auc@y.values[[1]]
  auc_tpr <- performance(pr2,"tpr","fpr")
  return(list(auc, auc_val, auc_tpr, pred))
}

###############################################################################
### Figures                                                                 ###
###############################################################################

metrics_output <- readRDS('data/CPI_merged_data_orig_gc_20210517.RDS')
metrics_output <- metrics_output  %>% 
  rename(histology = Tum_Type)

# Remove samples that have had previous CPI treatment
library(readxl)
CPI.treatment <- read_excel('data/CPI1000_treatment.xlsx')

treatment.patients <- CPI.treatment$Patient[which(CPI.treatment$`Biopsy timepoint` %in% c('on_treatment','prior_CTLA-4'))]

metrics_output <- metrics_output %>%
  filter(!case %in% treatment.patients)

###############################################################################
### Statistics + numbers                                                    ###
############################################################################### 

# Numbers 
# 1. Remove SNYDER PLOSMED 2017
metrics_final <- metrics_output %>%
  filter(study != 'SNYDER_PLOSMED_2017') 



summary(metrics_final$study)


# How many with RNA-seq
summary(!is.na(metrics_final$CD8A))


# How many with WES
summary(metrics_final$WXS)


summary(factor(metrics_final$histology[metrics_final$WXS]))
summary(factor(metrics_final$histology[!is.na(metrics_final$CD8A)]))
summary(factor(metrics_final$histology[!is.na(metrics_final$CD8A) & metrics_final$WXS]))

summary(factor(metrics_final$histology[!is.na(metrics_final$CD8A) & !is.na(metrics_final$tumour.tcell.fraction)]))

# How many with TCRA score 
summary(!is.na(metrics_final$tumour.tcell.fraction))


# Responders
responders.table <- metrics_final %>%
  filter(!is.na(tumour.tcell.fraction.gc)) %>%
  mutate(immune.cold = tumour.tcell.fraction.gc < median(metrics_final$tumour.tcell.fraction.gc, na.rm = TRUE)) %>%
  dplyr::select(Response_responder, immune.cold) %>%
  table()
responders.table


fisher.test(responders.table)


metrics_final <- metrics_final %>%
  filter(WXS)

case.order <- metrics_final %>%
  arrange(histology, study, !is.na(CD8A), Response_responder,desc(tumour.tcell.fraction.gc)) %>%
  select(case) %>% `[[`(1)

# Plot of tumour.tcell.fraction

tcell.fract.plot <- metrics_final %>%
  dplyr::select(case, tumour.tcell.fraction.gc, purity) %>%
  mutate(tumour.tcell.fraction.gc = ifelse(tumour.tcell.fraction.gc < 0 ,0, tumour.tcell.fraction.gc)) %>%
  mutate(case = factor(case, levels = case.order)) %>%
  ggplot(aes(case, tumour.tcell.fraction.gc)) +
  geom_col() + ylab('Tumour T cell fraction estimate') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

tcell.fract.norm.plot <- metrics_final %>%
  dplyr::select(case, normal.tcell.fraction.gc) %>%
  mutate(case = factor(case, levels = case.order)) %>%
  ggplot(aes(case, normal.tcell.fraction.gc)) +
  geom_col() + ylab('Blood T cell fraction estimate') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank())

tmb.plot <- metrics_final %>%
  dplyr::select(case, Clonal_TMB) %>%
  mutate(Clonal_TMB = log(ifelse(is.na(Clonal_TMB),0,Clonal_TMB) +1)) %>%
  mutate(case = factor(case, levels = case.order)) %>%
  ggplot(aes(case, Clonal_TMB)) +
  geom_col() + ylab('log(clonal TMB +1)') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank())



metrics_plot <- metrics_final %>%
  mutate(RNAseq = !is.na(CD8A)) %>%
  mutate(Clonal_TMB = log(ifelse(is.na(Clonal_TMB),0,Clonal_TMB) +1)) %>%
  dplyr::select(sample = case,study,histology, 
                RNAseq,Clonal_TMB, 
                response = Response_responder)

# Tiles:
# WES - TRUE/FALSE
# RNA-seq - TRUE/FALSE
# Study - categorical
# Tumour type - categorical
# TCRA tcell-fraction (tumour)
# TCRA tcell-fraction (blood)
# Response_responder

# Order by:
# Study, Tumour type, RNA-seq present,TMB

source('ggplot_add_tiles.R')

library(RColorBrewer)
library(gtable)
library(grid)
library(gridExtra)

plot.list <- list()
plot.list[[1]] <- tcell.fract.plot 
plot.list[[2]] <- tcell.fract.norm.plot

pdf('plots/ext6/CPI_overview_gc.pdf', height = 6.6, width = 10.38)
addTiles_ggplot(plot.list, metrics_plot, case.order,
                heights = c(3,3,0.5,0.5,0.5,0.5,0.5),
                widths = c(10,2,2), text_size = 10, legend_col_n = 2,
                c('Paired', 'Set2', 'Set3', 'Accent'))
dev.off()

# Figure 4b -----

tmp1 <- metrics_final %>%
  filter(Clonal_TMB > 0)



pdf('plots/fig4/CPI_TMBvsTCRA_v3_gc.pdf', height = 1.65, width = 2.5)
tmp1 %>%
  ggplot(aes(Clonal_TMB, tumour.tcell.fraction.gc)) +
  geom_point(aes(col = response), size = 0.25) + 
  # stat_cor(method = 'spearman', cor.coef.name = 'rho', size =2) +
  theme_bw() + theme(text = element_text(size=6), legend.position = 'none')  +
  xlab('Clonal TMB') + ylab('Tumour TCRA T cell') +
  scale_x_continuous(trans = 'log2',
                     breaks = c(1,10,20,50,100,200,500,1000,2000,5000,10000)) +
  scale_y_continuous(trans = 'log2',
                     breaks = c(0.01,0.025,0.05,0.1,0.25,0.5, 1),limits = c(0.001,3)) +
  geom_hline(yintercept = median(tmp1$tumour.tcell.fraction.gc, na.rm = TRUE), linetype = 'dashed', col = 'black') +
  geom_vline(xintercept = median(tmp1$Clonal_TMB, na.rm = TRUE), linetype = 'dashed', col = 'black')
dev.off()


# Make pieplots
# Make pieplots
blank_theme <- theme_minimal()+
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.border = element_blank(),
    panel.grid=element_blank(),
    axis.ticks = element_blank(),
    plot.title=element_text(size=14, face="bold")
  )

TMB.threshold <- median(tmp1$Clonal_TMB, na.rm = TRUE)
Tcell.threshold <- median(tmp1$tumour.tcell.fraction.gc, na.rm = TRUE)

tmp1 %>%
  filter(Clonal_TMB >= TMB.threshold  & 
           tumour.tcell.fraction.gc >= Tcell.threshold ) %>%
  group_by(response) %>%
  summarise(n())

pdf('plots/fig4/piechart_qtr_gc.pdf', height = 1, width = 1)
tmp1 %>%
  filter(Clonal_TMB >= TMB.threshold  & tumour.tcell.fraction.gc >= Tcell.threshold ) %>%
  ggplot() +
  geom_bar(aes(x='',fill = response), position = 'fill') +
  coord_polar("y", start=0) + 
  blank_theme +
  theme(axis.text.x=element_blank(), legend.position = 'none') +
  # geom_text(aes(x = 0,y = 0, label = '')) + 
  annotate('text',x = 1,y = 0.25, label = percent(103/(103+126)), size=2) +
  annotate('text',x = 1,y = 0.75, label = percent(126/(103+126)), size=2)
dev.off()

tmp1 %>%
  filter(Clonal_TMB >= TMB.threshold  & 
           tumour.tcell.fraction.gc < Tcell.threshold ) %>%
  group_by(response) %>%
  summarise(n())


pdf('plots/fig4/piechart_qbr_gc.pdf', height = 1, width = 1)
tmp1 %>%
  filter(Clonal_TMB >= TMB.threshold  & tumour.tcell.fraction.gc < Tcell.threshold ) %>%
  ggplot() +
  geom_bar(aes(x ='',fill = response), position = 'fill') +
  coord_polar("y", start=0) +  blank_theme +
  theme(axis.text.x=element_blank(), legend.position = 'none') +
  annotate('text',y = 0.25, x =1,  label = percent(65/(65 + 155)), size=2) + 
  annotate('text',y = 0.75, x =1,  label = percent(155/(65  + 155)), size=2)
dev.off()

tmp1 %>%
  filter(Clonal_TMB < TMB.threshold  & tumour.tcell.fraction.gc >=  Tcell.threshold ) %>%
  group_by(response) %>%
  summarise(n())

pdf('plots/fig4/piechart_qtl.pdf', height = 1, width = 1)
tmp1 %>%
  filter(Clonal_TMB < TMB.threshold  & tumour.tcell.fraction.gc  >= Tcell.threshold) %>%
  ggplot() +
  geom_bar(aes(x='',fill = response), position = 'fill') +
  coord_polar("y", start=0) + 
  blank_theme +
  theme(axis.text.x=element_blank(), legend.position = 'none') +
  annotate('text',y = 0.125, x =1,  label = percent(44/(44 + 165)), size = 2) + 
  annotate('text',y = 0.625, x =1,  label = percent(165/(44 + 165)), size = 2)
dev.off()

tmp1 %>%
  filter(Clonal_TMB < TMB.threshold  & 
           tumour.tcell.fraction.gc <  Tcell.threshold ) %>%
  group_by(response) %>%
  summarise(n())

pdf('plots/fig4/piechart_qbl.pdf', height = 1, width = 1)
tmp1 %>%
  filter(Clonal_TMB < TMB.threshold & tumour.tcell.fraction.gc < Tcell.threshold ) %>%
  ggplot() +
  geom_bar(aes(x='',fill = response), position = 'fill') +
  coord_polar("y", start=0) + 
  blank_theme +
  theme(axis.text.x=element_blank(), legend.position = 'none') +
  annotate('text', y = 0.09, x =1,  label = percent(17/(200 + 17)), size = 2) + 
  annotate('text',y = 0.625, x =1,  label = percent(200/(200 + 17)), size = 2)
dev.off()

pdf('plots/fig4/CPI_TMBvsTCRA.pdf', height = 2, width = 2)
metrics_final %>%
  filter(Clonal_TMB > 0) %>%
  ggplot(aes(Clonal_TMB, tumour.tcell.fraction)) +
  geom_point(aes(col = response), size = 0.25) + 
  stat_cor(method = 'spearman', cor.coef.name = 'rho', size = 2) +
  theme_bw() + theme(text = element_text(size=6))  +
  xlab('Clonal TMB') + ylab('Tumour T cell fraction estimate') +
  scale_x_continuous(trans = 'log2',
                     breaks = c(1,10,20,50,100,250,500,1000,2500,5000,10000)) +
  geom_hline(yintercept = 0.1, linetype = 'dashed', col = 'black')
dev.off()

pdf('plots/fig4/CPI_TMBvsTCRAnorm.pdf', height = 2, width = 2)
metrics_final %>%
  filter(Clonal_TMB > 0) %>%
  ggplot(aes(log(Clonal_TMB + 1), normal.tcell.fraction)) +
  geom_point(size = 0.25) + stat_cor(method = 'spearman', cor.coef.name = 'rho', size = 2) +
  theme_bw() + theme(text = element_text(size=6))  +
  xlab('log(Clonal TMB + 1)') + ylab('Blood T cell fraction estimate')
dev.off()




# Figure 4a - responders vs non-responders ----

pdf('plots/fig4/CPI_response_cncorrect_gc.pdf', width = 2.5, height = 1.65)
metrics_output %>%
  filter(response %in% c('response', 'no_response')) %>%
  filter(study != 'SNYDER_PLOSMED_2017') %>%
  # filter(tumour.original < 0.5) %>%
  ggplot(aes(response, tumour.tcell.fraction.gc)) +
  geom_violin(aes(fill = response)) +  geom_sina(col = rgb(0,0,0,0.2), size = 0.5) +
  geom_hline(yintercept = mean(metrics_final$tumour.tcell.fraction.gc, na.rm = TRUE), linetype = 'dashed') +
  # ylim(0,0.25) + 
  stat_compare_means(size = 2) + ylab('Tumour TCRA T cell fraction') + theme_bw() + 
  theme(legend.position = 'none', text = element_text(size = 7))
dev.off()


response.plot.data <- metrics_output %>%
  filter(response %in% c('response', 'no_response')) %>%
  filter(study != 'SNYDER_PLOSMED_2017') 
wilcox_test(response.plot.data, tumour.tcell.fraction.gc ~response)
wilcox_effsize(response.plot.data, tumour.tcell.fraction.gc ~response)
# ES = 0.174

mean(response.plot.data$tumour.tcell.fraction.gc[which(response.plot.data$response == 'response')], na.rm = TRUE) # 0.097
mean(response.plot.data$tumour.tcell.fraction.gc[which(response.plot.data$response == 'no_response')], na.rm = TRUE) # 0.058


# Figure 4C - Univariate analysis  ----

metrics_output <- metrics_output %>%
  mutate(tcell.fraction.change.gc = tumour.tcell.fraction.gc - normal.tcell.fraction.gc)

cols_for_test <- which(colnames(metrics_output) %in% 
                         c('CD8A','Clonal_TMB','purity',
                           'tumour.tcell.fraction.gc','normal.tcell.fraction.gc',
                           'tcell.fraction.change.gc'))

studies<-c("MELANOMA-SNYDER_NEJM_2014", #57
           "MELANOMA-VANALLEN_SCIENCE_2015", #100
           "MELANOMA-HUGO_CELL_2016", #38
           "MELANOMA-RIAZ_CELL_2017", #80
           "MELANOMA-CRISTESCU_SCIENCE_2018", #89
           "BLADDER-MARIATHASAN_NATURE_2018", #346
           "BLADDER-CRISTESCU_SCIENCE_2018", #17
           "RENAL-MCDERMOT_NMED_2018", # 105
           "HEAD AND NECK-CRISTESCU_SCIENCE_2018", # 107
           "BREAST-CRISTESCU_SCIENCE_2018") # 14


ORs_out <- data.frame(matrix(0, ncol = length(studies), nrow = length(cols_for_test)))
Pvals_out <- data.frame(matrix(1, ncol = length(studies), nrow = length(cols_for_test)))
Estimates_out <- data.frame(matrix(0, ncol = length(studies), nrow = length(cols_for_test)))
SEs_out <- data.frame(matrix(0, ncol = length(studies), nrow = length(cols_for_test)))

metrics_output <- metrics_output %>%
  filter(!is.na(CD8A)) %>% filter(!is.na(tumour.tcell.fraction.gc)) %>%
  filter(!is.na(Clonal_TMB))

cols.loc <- which(colnames(metrics_output) %in% c('histology','study','response'))

for(i in 1:length(cols_for_test)){
  this_tab<-metrics_output[,c(1,cols_for_test[i],cols.loc)]
  names(this_tab)[2]<-"metric"
  for(j in 1:length(studies)){
    this_hist<-this_tab[this_tab$histology==as.character(data.frame(strsplit(studies[j],"-"))[1,1]),]
    this_study<-this_hist[this_hist$study==as.character(data.frame(strsplit(studies[j],"-"))[2,1]),]
    this_study<-this_study[complete.cases(this_study),]
    print(studies[j]);print(table(this_study$response))
    
    if(nrow(this_study)>=10 && sum(this_study$metric)){
      this_study<-this_study[!(is.na(this_study$metric)),]
      this_study$metric_scaled<-scale(this_study$metric)
      #this_study$metric_scaled<-this_study$metric
      
      #use logistic regression
      this_study$Response_binary<-1
      this_study$Response_binary[this_study$response=="no_response"]<-0
      b<-summary(glm(Response_binary~metric_scaled,data=this_study,family=binomial))
      Pvals_out[(i),(j)]<-b$coefficients[2,4]
      if(b$coefficients[2,1]>log(10)){Estimates_out[(i),(j)]<-log(10)}
      if(b$coefficients[2,1]<log(0.1)){Estimates_out[(i),(j)]<-log(0.1)}
      if(b$coefficients[2,1]>=log(0.1) && b$coefficients[2,1]<=log(10)){Estimates_out[(i),(j)]<-b$coefficients[2,1]}
      SEs_out[(i),(j)]<-b$coefficients[2,2]
    }               
  }
}

names(Pvals_out)<-studies
names(Estimates_out)<-studies
names(SEs_out)<-studies
row_names<-colnames(metrics_output)[cols_for_test]
rownames(Estimates_out)<-row_names
rownames(Pvals_out)<-row_names
rownames(SEs_out)<-row_names
Estimates_ORs_out_set1<-exp(Estimates_out)


meta_set1<-data.frame(matrix(1, ncol = 4, nrow = nrow(Estimates_out)))
rownames(meta_set1)<-row_names
colnames(meta_set1)<-c("OR","lower_OR","upper_OR","Meta_Pval")

study.log2.values <- list()

for(k in 1:nrow(Estimates_out)){
  #print(rownames(Estimates_out)[k])
  ######need to do a look here over all the measures, do meta to get the ORs per study, per metirc, this then becomes what gets plotted (on log2 scale) in heatmap
  #####can then also do forest with ORs from meta-analysis for each measure, with p-value <0.001 etc
  #meta_dat<-data.frame(t(Estimates_out[(k-15),]),t(SEs_out[(k-15),]),studies)
  meta_dat<-data.frame(t(Estimates_out[(k),]),t(SEs_out[(k),]),studies)
  names(meta_dat)[1]<-"Effect"
  names(meta_dat)[2]<-"SE"
  #print(metagen(Effect,SE, studlab = studies,sm = "OR",data = meta_dat))
  a<-metagen(Effect,SE, studlab = studies,sm = "OR",data = meta_dat)
  meta_set1[k,1]<-exp(a$TE.fixed)
  meta_set1[k,2]<-exp(a$lower.fixed)
  meta_set1[k,3]<-exp(a$upper.fixed)
  meta_set1[k,4]<-a$pval.fixe
  
  study.log2.values[[k]] <- data.frame(variable = rownames(Estimates_out)[k],
                                       log2OR = a$data$.TE,
                                       study = a$data$.studlab)
}
print(meta_set1)

study.log2.values <- Reduce(rbind, study.log2.values)

pdf('plots/fig4/CPI_forest_gc.pdf', width = 1.5, height = 2.65)
meta_set1 %>%
  rownames_to_column('Variable') %>%
  ggplot(aes(Variable,OR, ymin=lower_OR, ymax = upper_OR)) +
  geom_pointrange(size = 0.25) +
  geom_hline(aes(yintercept =1), linetype = 2) +
  # ylim(0.75,2.5)+
  geom_text(mapping = aes(x = Variable, y=2, label = paste0('p=',formatC(Meta_Pval))), size = 2) +
  coord_flip(ylim = c(0.75, 2), clip = 'off') + xlab('') + theme_bw() +
  theme( text = element_text(size=6)) + 
  ggtitle('Meta analysis across all cohorts')
dev.off()

# Study heatmap - note only included studies with everything!
pdf('plots/fig4/CPI_study_heat_univariate_gc.pdf', width = 2, height = 2.65)
study.log2.values %>%
  mutate(log2OR = ifelse(log2OR == 0, NA, log2OR)) %>%
  filter(!study %in% c('CRC-Diaz_MSI_CRC','LUNG-HELLMAN_ALL_2018','LUNG-RIZVI_SCIENCE_2015')) %>%
  mutate(variable = factor(variable, levels = rev(c('tumour.tcell.fraction.gc', 'tcell.fraction.change.gc','purity',
                                                    'normal.tcell.fraction.gc','Clonal_TMB','CD8A')))) %>%
  ggplot(aes(study, variable)) +
  geom_tile(aes(fill = log2OR)) +
  scale_fill_gradient2(low = ("blue"),mid = "white",
                       high = ("red"),midpoint = 0, na.value = 'grey') +
  theme_bw() + xlab('') + ylab('') +
  labs(fill = "log2 OR") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.title = element_text(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=6))
dev.off()

# Extended data figure 6 - Multivariate analysis  ----

studies<-c("MELANOMA-HUGO_CELL_2016","MELANOMA-RIAZ_CELL_2017",
           "MELANOMA-CRISTESCU_SCIENCE_2018","LUNG-RIZVI_SCIENCE_2015",
           "LUNG-HELLMAN_ALL_2018", #"OTHER-CRISTESCU_SCIENCE_2018",
           "BLADDER-MARIATHASAN_NATURE_2018","BLADDER-CRISTESCU_SCIENCE_2018",
           "RENAL-MCDERMOT_NMED_2018","CRC-Diaz_MSI_CRC",
           "HEAD AND NECK-CRISTESCU_SCIENCE_2018","BREAST-CRISTESCU_SCIENCE_2018")

cols_uni_sig_orig <- c("Clonal_TMB","TMB","CD8A",
                       "tumour.tcell.fraction.gc")


mva_set_scaled_orig <- create_mva_scaled(cols_uni_sig_orig, metrics_output,
                                         studies)

mva_set_scaled_orig <- mva_set_scaled_orig %>%
  filter(!is.na(CD8A)) %>% filter(!is.na(tumour.tcell.fraction.gc)) %>%
  filter(!is.na(Clonal_TMB))

mva_set_scaled_orig %>%
  ggplot(aes(CD8A, tumour.tcell.fraction.gc)) +
  geom_point(aes(col = response)) +
  stat_cor(method = 'spearman', cor.coef.name = 'rho') +
  xlab('scaled CD8A') + ylab('scaled TCRA T cell fraction') +theme_bw()

# 578 samples with all
# Select samples with all data
mva_mod_glm1 <- glm(Response_binary ~ Clonal_TMB,
                    data=mva_set_scaled_orig,
                    family=binomial())
mva_mod_glm2 <- glm(Response_binary ~ Clonal_TMB + CD8A,
                    data=mva_set_scaled_orig,
                    family=binomial())
mva_mod_glm3 <- glm(Response_binary ~ Clonal_TMB + tumour.tcell.fraction.gc,
                    data=mva_set_scaled_orig,
                    family=binomial())

mva_mod_glm4 <- glm(Response_binary ~ Clonal_TMB + tumour.tcell.fraction.gc + CD8A,
                    data=mva_set_scaled_orig,
                    family=binomial())

summary(mva_mod_glm1)
summary(mva_mod_glm2)
summary(mva_mod_glm3)
summary(mva_mod_glm4)

model.df <- data.frame(model.name = c('clonal.TMB + CD8A', 'clonal.TMB + TCRA',
                                      'clonal.TMB + TCRA + CD8A', 'clonal.TMB + TCRA + CD8A'),
                       AIC.change = -c(mva_mod_glm2$aic - mva_mod_glm1$aic,
                                       mva_mod_glm3$aic - mva_mod_glm1$aic,
                                       mva_mod_glm4$aic - mva_mod_glm1$aic,
                                       mva_mod_glm4$aic - mva_mod_glm1$aic),
                       logPvalue =c(-log2(summary(mva_mod_glm2)$coefficients[3,4]), 
                                    -log2(summary(mva_mod_glm3)$coefficients[3,4]),
                                    -log2(summary(mva_mod_glm4)$coefficients[3,4]),
                                    -log2(summary(mva_mod_glm4)$coefficients[4,4])
                       ),
                       variable = c('CD8A', 'TCRA','TCRA','CD8A'))



auc_orig1 <- auc.pred.fun(mva_mod_glm1, mva_set_scaled_orig,
                          cols_uni_sig_orig, glm = TRUE)
auc_orig2 <- auc.pred.fun(mva_mod_glm2, mva_set_scaled_orig,
                          cols_uni_sig_orig, glm = TRUE)
auc_orig3 <- auc.pred.fun(mva_mod_glm3, mva_set_scaled_orig,
                          cols_uni_sig_orig, glm = TRUE)
auc_orig4 <- auc.pred.fun(mva_mod_glm4, mva_set_scaled_orig,
                          cols_uni_sig_orig, glm = TRUE)


pdf('plots/ext6/glm_model_rocplot_gc.pdf', width = 5, height = 5)
plot(auc_orig1[[3]],col="darkblue",lwd=3,xaxt="n") + 
  mtext(paste0("Clonal TMB: ",
               round(auc_orig1[[2]],2)),
        col="darkblue",line = -15, at = 0.75,cex=0.8) + 
  abline(coef = c(0,1)) 

plot(auc_orig2[[3]],col="darkgreen",add=T, lwd=3,xaxt="n") + 
  mtext(paste0("Clonal TMB + CD8A: ",
               round(auc_orig2[[2]],2)),
        col="darkgreen",line = -16, at = 0.75,cex=0.8) + 
  abline(coef = c(0,1)) 

plot(auc_orig3[[3]],col="darkred",add=T, lwd=3,xaxt="n") + 
  mtext(paste0("Clonal TMB + TCRA: ",
               round(auc_orig3[[2]],2)),
        col="darkred",line = -17, at = 0.75,cex=0.8) + 
  abline(coef = c(0,1)) 

plot(auc_orig4[[3]],col="red",add=T, lwd=3,xaxt="n") + 
  mtext(paste0("Clonal TMB + TCRA + CD8A: ",
               round(auc_orig4[[2]],2)),
        col="darkred",line = -18, at = 0.75,cex=0.8) + 
  abline(coef = c(0,1)) 

dev.off()

roc1 <- roc(mva_set_scaled_orig$Response_binary,
            as.numeric(auc_orig1[[4]]))
roc2 <- roc(mva_set_scaled_orig$Response_binary,
            as.numeric(auc_orig2[[4]]))
roc3 <- roc(mva_set_scaled_orig$Response_binary,
            as.numeric(auc_orig3[[4]]))
roc4 <- roc(mva_set_scaled_orig$Response_binary,
            as.numeric(auc_orig4[[4]]))
print(roc.test(roc1, roc2)) # p.value = 0.112
print(roc.test(roc1, roc3)) # p-value = 0.01121
print(roc.test(roc2, roc3)) # p-value = 0.1939
print(roc.test(roc4, roc3)) # 0.91


# Extended figure 6  Lung cohort overview  ----

metrics_output <- readRDS('data/CPI_merged_data_withShim_gc_CI_20210428.RDS')

metrics_output <- metrics_output %>%
  mutate(tcell.fraction.change.gc = tumour.tcell.fraction.gc - normal.tcell.fraction.gc)

metrics_final <- metrics_output %>%
  filter(study %in% c("RIZVI_SCIENCE_2015","HELLMAN_ALL_2018",
                      "SHIM_AOO_2020")) 

metrics_final$histology = 'LUNG'
metrics_final$WXS = TRUE

case.order <- metrics_final %>%
  arrange(response, desc(tumour.tcell.fraction.gc), study) %>%
  select(case) %>% `[[`(1)

# Add median lines for responders/non-responders
median(metrics_final$tumour.tcell.fraction.gc[!metrics_final$response == 'response'], na.rm = TRUE) # 0.08355291
median(metrics_final$tumour.tcell.fraction.gc[metrics_final$response == 'response'], na.rm = TRUE) # 0.1161409

tcell.fract.plot <- metrics_final %>%
  dplyr::select(case, tumour.tcell.fraction.gc) %>%
  mutate(case = factor(case, levels = case.order)) %>%
  ggplot(aes(case, tumour.tcell.fraction.gc)) +
  geom_col() + ylab('Tumour TCRA T cell fraction') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank()) +
  annotate("segment", x = 0, xend = 177, y = 0.08355291, yend = 0.08355291,
           colour = "red") + 
  annotate("segment", x = 177, xend = 267, y = 0.1161409, yend = 0.1161409,
           colour = "red")

tcell.fract.blood.plot <- metrics_final %>%
  dplyr::select(case, normal.tcell.fraction.gc) %>%
  mutate(case = factor(case, levels = case.order)) %>%
  ggplot(aes(case, normal.tcell.fraction.gc)) +
  geom_col() + ylab('Normal T cell fraction estimate') + 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        axis.title.x = element_blank())


metrics_plot <- metrics_final %>%
  mutate(RNAseq = !is.na(CD8A)) %>%
  dplyr::select(sample = case,study,
                response)

# Tiles:
# WES - TRUE/FALSE
# RNA-seq - TRUE/FALSE
# Study - categorical
# Tumour type - categorical
# TCRA tcell-fraction (tumour)
# TCRA tcell-fraction (blood)
# Response_responder

# Order by:
# Study, Tumour type, RNA-seq present,TMB

source('ggplot_add_tiles.R')

library(RColorBrewer)
library(gtable)
library(grid)
library(gridExtra)

plot.list <- list()
plot.list[[1]] <- tcell.fract.plot
plot.list[[2]] <-tcell.fract.blood.plot

metrics_plot$response <- factor(metrics_plot$response, levels = c('response', 'no_response'))

pdf('plots/ext6/CPI_lung_overview.pdf', height = 4, width = 10)
addTiles_ggplot(plot.list, metrics_plot, case.order,
                heights = c(3,3, 0.5,0.5),
                widths = c(10,2), text_size = 10, legend_col_n = 1)
dev.off()

#Parameter information:
# *ggplot_list = List with ggplots to which the tiles should be added
# *tile_data = dataframe with samples in first column and all tiles that should be included as columns (column names is used as label)
# *order_samples = samples in order they are supposed to be plotted
# *heights and weights = heights and weights for grid plot
# *text_size = size for text in plots
# *categorical and numerical = color palletes for tile plots
# *legend_col_n = value from 1 to 3 to classify in how many columns the legend should be split

# Figure 4d - Lung cohort univariate analysis  ----

cols_for_test <- which(colnames(metrics_output) %in% 
                         c(
                           'purity',
                           #'tumour.tcell.fraction','normal.tcell.fraction',
                           #'tcell.fraction.change',
                           'tumour.tcell.fraction.gc','normal.tcell.fraction.gc',
                           'tcell.fraction.change.gc'
                         ))

studies<-c("LUNG-RIZVI_SCIENCE_2015","LUNG-HELLMAN_ALL_2018",
           "LUNG-SHIM_AOO_2020")


ORs_out <- data.frame(matrix(0, ncol = length(studies), nrow = length(cols_for_test)))
Pvals_out <- data.frame(matrix(1, ncol = length(studies), nrow = length(cols_for_test)))
Estimates_out <- data.frame(matrix(0, ncol = length(studies), nrow = length(cols_for_test)))
SEs_out <- data.frame(matrix(0, ncol = length(studies), nrow = length(cols_for_test)))


metrics_output <- metrics_output %>%
  filter(study!='SNYDER_PLOSMED_2017') %>%
  filter(!is.na(tumour.tcell.fraction)) %>% filter(is.na(CD8A))  %>%
  filter(study %in% c( 'HELLMAN_ALL_2018', 'RIZVI_SCIENCE_2015', 'SHIM_AOO_2020'))
# filter(study %in% c('Diaz_MSI_CRC', 'HELLMAN_ALL_2018', 'RIZVI_SCIENCE_2015'))
#filter(!study %in% c('RIZVI_SCIENCE_2015','HUGO_CELL_2016'))
# 1. Use everything


cols.loc <- which(colnames(metrics_output) %in% c('histology','study','response'))

for(i in 1:length(cols_for_test)){
  this_tab<-metrics_output[,c(1,cols_for_test[i],cols.loc)]
  names(this_tab)[2]<-"metric"
  for(j in 1:length(studies)){
    this_hist<-this_tab[this_tab$histology==as.character(data.frame(strsplit(studies[j],"-"))[1,1]),]
    this_study<-this_hist[this_hist$study==as.character(data.frame(strsplit(studies[j],"-"))[2,1]),]
    this_study<-this_study[complete.cases(this_study),]
    print(studies[j]);print(table(this_study$response))
    
    if(nrow(this_study)>=10 && sum(this_study$metric)){
      this_study<-this_study[!(is.na(this_study$metric)),]
      this_study$metric_scaled<-scale(this_study$metric)
      #this_study$metric_scaled<-this_study$metric
      
      #use logistic regression
      this_study$Response_binary<-1
      this_study$Response_binary[this_study$response=="no_response"]<-0
      b<-summary(glm(Response_binary~metric_scaled,data=this_study,family=binomial))
      Pvals_out[(i),(j)]<-b$coefficients[2,4]
      if(b$coefficients[2,1]>log(10)){Estimates_out[(i),(j)]<-log(10)}
      if(b$coefficients[2,1]<log(0.1)){Estimates_out[(i),(j)]<-log(0.1)}
      if(b$coefficients[2,1]>=log(0.1) && b$coefficients[2,1]<=log(10)){Estimates_out[(i),(j)]<-b$coefficients[2,1]}
      SEs_out[(i),(j)]<-b$coefficients[2,2]
    }               
  }
}

names(Pvals_out)<-studies
names(Estimates_out)<-studies
names(SEs_out)<-studies
row_names<-colnames(metrics_output)[cols_for_test]
rownames(Estimates_out)<-row_names
rownames(Pvals_out)<-row_names
rownames(SEs_out)<-row_names
Estimates_ORs_out_set1<-exp(Estimates_out)


meta_set1<-data.frame(matrix(1, ncol = 4, nrow = nrow(Estimates_out)))
rownames(meta_set1)<-row_names
colnames(meta_set1)<-c("OR","lower_OR","upper_OR","Meta_Pval")

study.log2.values <- list()

for(k in 1:nrow(Estimates_out)){
  #print(rownames(Estimates_out)[k])
  ######need to do a look here over all the measures, do meta to get the ORs per study, per metirc, this then becomes what gets plotted (on log2 scale) in heatmap
  #####can then also do forest with ORs from meta-analysis for each measure, with p-value <0.001 etc
  #meta_dat<-data.frame(t(Estimates_out[(k-15),]),t(SEs_out[(k-15),]),studies)
  meta_dat<-data.frame(t(Estimates_out[(k),]),t(SEs_out[(k),]),studies)
  names(meta_dat)[1]<-"Effect"
  names(meta_dat)[2]<-"SE"
  #print(metagen(Effect,SE, studlab = studies,sm = "OR",data = meta_dat))
  a<-metagen(Effect,SE, studlab = studies,sm = "OR",data = meta_dat)
  meta_set1[k,1]<-exp(a$TE.fixed)
  meta_set1[k,2]<-exp(a$lower.fixed)
  meta_set1[k,3]<-exp(a$upper.fixed)
  meta_set1[k,4]<-a$pval.fixe
  
  study.log2.values[[k]] <- data.frame(variable = rownames(Estimates_out)[k],
                                       log2OR = a$data$.TE,
                                       study = a$data$.studlab)
}
print(meta_set1)

study.log2.values <- Reduce(rbind, study.log2.values)

pdf('plots/fig4/CPI_forest_lung.pdf', width = 1.5, height = 2)
meta_set1 %>%
  rownames_to_column('Variable') %>%
  ggplot(aes(Variable,OR, ymin=lower_OR, ymax = upper_OR)) +
  geom_pointrange(size = 0.25) +
  geom_hline(aes(yintercept =1), linetype = 2) +
  # ylim(0.75,2.5)+
  geom_text(mapping = aes(x = Variable , y=2, label = paste0('p=',formatC(Meta_Pval))), nudge_x = 0.2, size = 2) +
  coord_flip(clip = 'off') + xlab('') + theme_bw() +
  theme(text = element_text(size=6)) + 
  ggtitle('Meta analysis across lung cohorts')
dev.off()

pdf('plots/fig4/CPI_study_heat_univariate_lung.pdf', width = 1.5, height = 3)
study.log2.values %>%
  mutate(log2OR = ifelse(log2OR == 0, NA, log2OR)) %>%
  # filter(study %in% c('CRC-Diaz_MSI_CRC', 'LUNG-HELLMAN_ALL_2018', 'LUNG-RIZVI_SCIENCE_2015')) %>%
  mutate(variable = factor(variable, levels = rev(c('tumour.tcell.fraction.gc', 'tcell.fraction.change.gc','purity',
                                                    'normal.tcell.fraction.gc','Clonal_TMB')))) %>%
  ggplot(aes(study, variable)) +
  geom_tile(aes(fill = log2OR)) +
  scale_fill_gradient2(low = ("blue"),mid = "white",
                       high = ("red"),midpoint = 0, na.value = 'grey', limits = c(-2,2)) +
  theme_bw() + xlab('') + ylab('') +
  labs(fill = "log2 OR") +
  theme(axis.text.x = element_text(angle = 45,hjust = 1),legend.title = element_text(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), text = element_text(size=6))
dev.off()
