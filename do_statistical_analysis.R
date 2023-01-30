library(tidyverse)
library(ggplot2)
####################
# Load data and build dataset

# Clinical Metadata table containing demographic and medical history info, SARS-CoV-2 variant (e.g. Alpha,Delta,BA.1,etc.), SampleDate, identified by 'sample_id' and 'ssi_id'
M <- read.table(file='samples_plus_metadata.tsv',sep='\t',header=T)

# WGS metadata table containing information on coverage and SARS-CoV-2 variant (e.g. Alpha,Delta,BA.1,etc.) with the columns 'ssi_id','variant','sample_id','cov'
N <- read.table(file='all_WGS_metadata.tsv',sep='\t',header=T)

# Subgenomic RNA read occurrence table with the columns 'sample' 'pos' 'count' 'lineage' 'mapped'
O <- read.table(file='results/sites_PE_occ_all_allR1_allsites_w_Nmapped_expanded_pseudocount_N_tilAUG.tsv',sep='\t',header=T)


library(data.table)

# WGS linelist containing the columns 'ssi_id' 'date_sequencing'
linelist <- fread(file='linelist.tsv',sep='\t',header=T)

####################
Building the dataset

M <- merge(M,linelist[,c('ssi_id','date_sequencing')],by='ssi_id',all.x=T)
M$sequencing_delay <- as.Date(M$date_sequencing) - as.Date(M$date_sampling)

library(lubridate)
M$yearmonth <- paste0(year(M$SampleDate),'M',ifelse(month(M$SampleDate)<10,paste0('0',month(M$SampleDate)),month(M$SampleDate))) 
M <- merge(M,N3,by='ssi_id', all.x=T)
M <- M[M$sample_id %in% O$sample,]

M$variant <- as.character(M$variant)
# BA.5 is split into BA.5 and late BA.5 as samples collected in April-May 2022 have a suboptimal amplicon1 coverage. in the main study we keep samples later than June 12.
M$variant[M$variant == 'BA.5'] <- ifelse(as.Date(M$SampleDate[M$variant == 'BA.5']) > as.Date('2022-06-12'),'Late BA.5','BA.5' )

# Late Delta corresponds to the delta samples collected synchonously with the first Omicron samples in fall 2021
M$variant[M$variant == 'Delta2'] <- 'Late Delta'

M$variant <- factor(M$variant,levels=c('Alpha','Delta','Late Delta','BA.1','BA.2','BA.2.12.1','BA.4','BA.5','Late BA.5'))

M$RUKS_Astma <-as.logical(M$RUKS_Astma)
M$RUKS_Diab1 <-as.logical(M$RUKS_Diab1)
M$RUKS_Diab2 <-as.logical(M$RUKS_Diab2)
M$RUKS_KOL <-as.logical(M$RUKS_KOL)
M$newlyAdm[is.na(M$newlyAdm)==T] <- 0
M$AdmResp[is.na(M$AdmResp)==T] <- 0

M$Death30Days_final[is.na(M$Death30Days_final)==T] <-0
M$Death60Days_final[is.na(M$Death60Days_final)==T] <-0

# list of diseases to include
dis_list <- c('RUKS_Astma','RUKS_Diab1','RUKS_Diab2','RUKS_KOL','Diabet',
              'Adipos','Cancer','Neuro','Nyre','Haem_c','Card_dis','Resp_dis','Immu_dis')

# NA means 'super healthy therefore can be set as FALSE= healthy, TRUE means having had the named condition within the past 5 years
for (d in dis_list){
  M[which(is.na(M[,d])==T),d] <- FALSE 
}


P <- merge(O,M[,setdiff(colnames(M),colnames(O))],by.x='sample',by.y='sample_id')
P$lineage <-factor(P$lineage, levels=c("Alpha",'Delta','Delta2','BA.1','BA.2','BA.2.12.1','BA.4','BA.5'))

# Pool positions together for each subgenomic RNA and reorder
P$site[P$pos %in% c(21552:21562)] <- 'S'
P$site[P$pos %in% c(28250:29903)] <- 'ORF9'
P$site[P$pos %in% c(26237:26244)] <- 'E'
P$site[P$pos %in% c(26468:26522)] <- 'M'
P$site[P$pos %in% c(27041:27202)] <- 'ORF6'
P$site[P$pos %in% c(27385:27393)] <- 'ORF7a'
P$site[P$pos %in% c(27884:27893)] <- 'ORF8'
P$site[P$pos %in% c(25382:25392)] <- 'ORF3a'
P$site <- factor(P$site, levels=c('S','ORF3a','E','M','ORF6','ORF7a','ORF8','ORF9'))

allres_sites_agg <- aggregate(x = P$count, by = list(P$sample,P$site), FUN= "sum")
colnames(allres_sites_agg) <- c('sample','site','count')
allres_sites_agg <- merge(allres_sites_agg,P[,c('sample','lineage','variant','mapped','yearmonth')],by='sample',all.x=T)
allres_sites_agg <- allres_sites_agg[-which(duplicated(allres_sites_agg)),]

allres_sites_agg$PE151 <- ifelse(allres_sites_agg$sample %in% N[N$readlen>76,c('sample_id')],TRUE,FALSE)

allres_sites_agg$site <- factor(allres_sites_agg$site, levels=c('S','ORF3a','E','M','ORF6','ORF7a','ORF8','ORF9'))
allres_sites_agg$lineage <- factor(allres_sites_agg$lineage, levels=c('Alpha','Delta','Delta2','BA.1','BA.2','BA.2.12.1','BA.4','BA.5'))
allres_sites_agg$variant <- factor(allres_sites_agg$variant,levels=c('Alpha','Delta','Late Delta','BA.1','BA.2','BA.2.12.1','BA.4','BA.5','Late BA.5'))
allres_sites_agg_meta <- merge(allres_sites_agg,P[,c('sample',setdiff(colnames(P),colnames(allres_sites_agg)))],by='sample')
allres_sites_agg_meta <- allres_sites_agg_meta[,-which(colnames(allres_sites_agg_meta)=='pos')]

# the merge created duplicates, remove them here
allres_sites_agg_meta <- allres_sites_agg_meta[-which(duplicated(allres_sites_agg_meta[,c('sample','site','count','variant')])),]

# count with different normalizations
allres_sites_agg_meta$count_per_1000_genomic <- allres_sites_agg_meta$count/(allres_sites_agg_meta$leader_depth)*1000
allres_sites_agg_meta$count_per_million <- allres_sites_agg_meta$count/(allres_sites_agg_meta$mapped)*1000000
allres_sites_agg_meta$count_per_2million <- allres_sites_agg_meta$count/(allres_sites_agg_meta$mapped)*2000000

# take out the samples with PE151 if any

allres_sites_agg <- allres_sites_agg[allres_sites_agg$PE151 == FALSE,]
allres_sites_agg_meta <- allres_sites_agg_meta[allres_sites_agg_meta$PE151 == FALSE,]

##########################################
# exclude BA.2.12.1, BA.4 and BA.5 groups from main analysis
allres_sites_agg_meta_excl <- allres_sites_agg_meta[-which(allres_sites_agg_meta$variant %in% c('BA.2.12.1','BA.4','BA.5')),]
allres_sites_agg_meta_excl$variant <- as.character(allres_sites_agg_meta_excl$variant)

allres_sites_agg_meta_excl$variant[allres_sites_agg_meta_excl$variant == 'Late BA.5'] <- 'BA.5'
allres_sites_agg_meta_excl$variant <- factor(allres_sites_agg_meta_excl$variant,levels=c('Alpha','Delta','Late Delta','BA.1','BA.2','BA.5'))

allres_sites_agg_meta_excl <- allres_sites_agg_meta_excl[-which(allres_sites_agg_meta_excl$yearmonth %in% c('2021M03','2021M04')),]
allres_sites_agg_meta_excl <- allres_sites_agg_meta_excl[-which(allres_sites_agg_meta_excl$yearmonth %in% c('2022M04','2022M05')),]


# check
table(unique(allres_sites_agg_meta[,c('sample','variant')])$variant)
table(unique(allres_sites_agg_meta_excl[,c('sample','variant')])$variant)



###########################################

# dummy variables
allres_sites_agg_meta_dummy <- allres_sites_agg_meta_excl
allres_sites_agg_meta_dummy$SexM <- ifelse(allres_sites_agg_meta_dummy$Sex == 'M',1,0)
allres_sites_agg_meta_dummy$lineageAlpha <- ifelse(allres_sites_agg_meta_dummy$lineage == 'Alpha',1,0)
allres_sites_agg_meta_dummy$lineageDelta <- ifelse(allres_sites_agg_meta_dummy$variant == 'Delta',1,0)
allres_sites_agg_meta_dummy$lineageLateDelta <- ifelse(allres_sites_agg_meta_dummy$variant == 'Late Delta',1,0)
allres_sites_agg_meta_dummy$lineageBA.1 <- ifelse(allres_sites_agg_meta_dummy$variant == 'BA.1',1,0)
allres_sites_agg_meta_dummy$lineageBA.2 <- ifelse(allres_sites_agg_meta_dummy$variant == 'BA.2',1,0)
allres_sites_agg_meta_dummy$lineageBA.5 <- ifelse(allres_sites_agg_meta_dummy$variant == 'BA.5',1,0)

allres_sites_agg_meta_dummy$SiteS <- ifelse(allres_sites_agg_meta_dummy$site == 'S',1,0)
allres_sites_agg_meta_dummy$SiteORF3a <- ifelse(allres_sites_agg_meta_dummy$site == 'ORF3a',1,0)
allres_sites_agg_meta_dummy$SiteE <- ifelse(allres_sites_agg_meta_dummy$site == 'E',1,0)
allres_sites_agg_meta_dummy$SiteM <- ifelse(allres_sites_agg_meta_dummy$site == 'M',1,0)
allres_sites_agg_meta_dummy$SiteORF6 <- ifelse(allres_sites_agg_meta_dummy$site == 'ORF6',1,0)
allres_sites_agg_meta_dummy$SiteORF7a <- ifelse(allres_sites_agg_meta_dummy$site == 'ORF7a',1,0)
allres_sites_agg_meta_dummy$SiteORF8 <- ifelse(allres_sites_agg_meta_dummy$site == 'ORF8',1,0)
allres_sites_agg_meta_dummy$SiteORF9 <- ifelse(allres_sites_agg_meta_dummy$site == 'ORF9',1,0)

allres_sites_agg_meta_dummy$age_group0_5 <- ifelse(allres_sites_agg_meta_dummy$age_group == '0-5',1,0)
allres_sites_agg_meta_dummy$age_group10_14 <- ifelse(allres_sites_agg_meta_dummy$age_group == '10-14',1,0)
allres_sites_agg_meta_dummy$age_group15_19 <- ifelse(allres_sites_agg_meta_dummy$age_group == '15-19',1,0)
allres_sites_agg_meta_dummy$age_group20_29 <- ifelse(allres_sites_agg_meta_dummy$age_group == '20-29',1,0)
allres_sites_agg_meta_dummy$age_group30_39 <- ifelse(allres_sites_agg_meta_dummy$age_group == '30-39',1,0)
allres_sites_agg_meta_dummy$age_group40_49 <- ifelse(allres_sites_agg_meta_dummy$age_group == '40-49',1,0)
allres_sites_agg_meta_dummy$age_group50_64 <- ifelse(allres_sites_agg_meta_dummy$age_group == '50-64',1,0)
allres_sites_agg_meta_dummy$age_group6_9 <- ifelse(allres_sites_agg_meta_dummy$age_group == '6-9',1,0)
allres_sites_agg_meta_dummy$age_group65_79 <- ifelse(allres_sites_agg_meta_dummy$age_group == '65-79',1,0)
allres_sites_agg_meta_dummy$age_group80 <- ifelse(allres_sites_agg_meta_dummy$age_group == '80+',1,0)

allres_sites_agg_meta_dummy$Vaccination_status_Booster <- ifelse(allres_sites_agg_meta_dummy$Vaccination_status == 'Booster vaccinated full effect',1,0)
allres_sites_agg_meta_dummy$Vaccination_status_no <- ifelse(allres_sites_agg_meta_dummy$Vaccination_status == 'Not vaccinated',1,0)
allres_sites_agg_meta_dummy$Vaccination_status_primary_full <- ifelse(allres_sites_agg_meta_dummy$Vaccination_status == 'Primary program full effect',1,0)
allres_sites_agg_meta_dummy$Vaccination_status_primary_start <- ifelse(allres_sites_agg_meta_dummy$Vaccination_status == 'Started primary program',1,0)

allres_sites_agg_meta_dummy$Region_Hovedstaden <- ifelse(allres_sites_agg_meta_dummy$Region == 'Hovedstaden',1,0)
allres_sites_agg_meta_dummy$Region_Midtjylland <- ifelse(allres_sites_agg_meta_dummy$Region == 'Midtjylland',1,0)
allres_sites_agg_meta_dummy$Region_Nordjylland <- ifelse(allres_sites_agg_meta_dummy$Region == 'Nordjylland',1,0)
allres_sites_agg_meta_dummy$Region_Sjælland <- ifelse(allres_sites_agg_meta_dummy$Region == 'Sjælland',1,0)
allres_sites_agg_meta_dummy$Region_Syddanmark <- ifelse(allres_sites_agg_meta_dummy$Region == 'Syddanmark',1,0)

# convert to numeric
allres_sites_agg_meta_dummy[,dis_list] <- apply(allres_sites_agg_meta_dummy[,dis_list],2,as.numeric )

# remove all cases without symptoms and without ct
allres_sites_agg_meta_dummy <- allres_sites_agg_meta_dummy[-which(is.na(allres_sites_agg_meta_dummy$HarSymptomer_STPS)==T),]
allres_sites_agg_meta_dummy <- allres_sites_agg_meta_dummy[-which(is.na(allres_sites_agg_meta_dummy$ct)==T),]
allres_sites_agg_meta_dummy$HarSymptomer_STPS <- as.numeric(allres_sites_agg_meta_dummy$HarSymptomer_STPS)
allres_sites_agg_meta_dummy$count_per_1000_genomic <- allres_sites_agg_meta_dummy$count/(allres_sites_agg_meta_dummy$leader_depth- allres_sites_agg_meta_dummy$count)*1000


# check

table(unique(allres_sites_agg_meta_dummy[,c('sample','variant')])$variant)


#########################################
# PCA biplot

dfPC <-reshape(allres_sites_agg_meta_excl[,c('sample','variant','count_per_2million','site')],direction = "wide",timevar = "site",idvar = c('sample','variant'))
dfPC[is.na(dfPC)] <- 0

#colnames(dfPC) <- c('sample','lineage','ORF6','E','M','S','N','ORF7a','ORF3a','ORF8')
colnames(dfPC) <- c('sample','lineage','S','ORF3a','ORF6','ORF9','E','M','ORF8','ORF7a')
pc <- prcomp(dfPC[,-c(1,2)],scale. = T,center = T)
library(ggfortify)
library(cluster)

pdf(file = 'results/PCAbiplot_count_vs_lineage_figure2_v2.pdf',width = 16,height = 9)
autoplot(pam(dfPC[,-c(1,2)],6),data=dfPC,colour='lineage',loadings=T,loadings.label=T,loadings.label.size = 6,loadings.label.colour = 'black',loadings.color = 'black',size=2, frame = TRUE, frame.type = 'norm')+ theme_bw() +theme(text=element_text(size=21)) + scale_color_manual(values = c("#F8766DFF","#B79F00FF", "#A0BA00FF","#00BFC4FF","#619CFFFF","#7C64F5FF")) #+labs(colour='lineage')
dev.off()

##########################################
# PCA biplot for disease history

tmp <- as.data.frame(allres_sites_agg_meta_excl[setdiff(colnames(allres_sites_agg_meta_excl), c("pos", "count_per_2million","count_per_million",'count','leader_depth','mapped','ssi_id','CPR','date_sampling','SampleDate','cpr_episode','date_diff',
                                                                                                                      'Episodekey','First_VaccineDate','First_VaccineName','Second_VaccineDate','Second_VaccineName','Re_vaccineDate','Re_vaccineName',
                                                                                                                      'NewlyAdmDate','LastDsc','ICUAdm','ICUDsc','AdmResp','EpilprResp_start', "Death30Days_final" , "Death60Days_final",  "DateOfDeath_final",
                                                                                                                      "countgrp",'lineage','Region','age_group','Haem_c','site','yearmonth','PE151','count_per_1000_genomic', 'date_sequencing', 'sequencing_delay','ct' ))])

tmp <- unique(tmp[complete.cases(tmp),])

tmp2 <- tmp[,-c(1,2,20)]

tmp2[,2:17] <-apply(tmp2[2:17],2,function(x){as.numeric(factor(x))})

pc <- prcomp(tmp2,scale. = T)
library(ggfortify)
 pdf(file = 'results/pca_medical_history_supplFig1.pdf',width = 16,height = 9)
 autoplot(pc,data=tmp,colour='variant',loadings=TRUE,loadings.label=T,loadings.label.size = 6,size=4)+ theme_bw() +theme(text=element_text(size=21))+labs(colour='lineage')
 dev.off()

####################################################################
# linear regression analysis and elastic net and variable importance plot

library(glmnet)
library(doParallel)
library(foreach)
library(pROC)

registerDoParallel(cores = 4)

n <- length(unique(allres_sites_agg_meta_dummy$sample))

set.seed(1234)
it <- 50
varimportance <- c()
for (it in 1:50){
sample <- sample(seq(n), size = n * 0.5, replace = FALSE)

train <- allres_sites_agg_meta_dummy[which(allres_sites_agg_meta_dummy$sample %in% unique(allres_sites_agg_meta_dummy$sample)[sample]), -1]  
test <- allres_sites_agg_meta_dummy[which(allres_sites_agg_meta_dummy$sample %in% unique(allres_sites_agg_meta_dummy$sample)[sample]), -1]


mdlY <- data.matrix(train["count_per_2million"])
mdlX <- data.matrix(train[setdiff(colnames(allres_sites_agg_meta_dummy), c("pos","sample", "count_per_million",'count_per_2million','count','leader_depth','mapped','ssi_id','CPR','date_sampling','SampleDate','cpr_episode','date_diff',
                                               'Episodekey','SampleAge','First_VaccineDate','First_VaccineName','Second_VaccineDate','Second_VaccineName','Re_vaccineDate','Re_vaccineName',
                                               'NewlyAdmDate','LastDsc','ICUAdm','ICUDsc','AdmResp','EpilprResp_start', "Death30Days_final" , "Death60Days_final",  "DateOfDeath_final",
                                               "countgrp",'lineage','Region','age_group','Sex','variant','Vaccination_status','site','yearmonth','cov','PE151','count_per_1000_genomic', 'date_sequencing', 'sequencing_delay' ))])


a <- seq(0.1, 0.9, 0.05)
search <- foreach(i = a, .combine = rbind) %dopar% {
  cv <- cv.glmnet(mdlX, mdlY, family = "gaussian", nfold = 10, type.measure = "deviance", paralle = TRUE, alpha = i,standardize = T)
  data.frame(cvm = cv$cvm[cv$lambda == cv$lambda.1se], lambda.1se = cv$lambda.1se, alpha = i)
}
cv3 <- search[search$cvm == min(search$cvm), ]
md3 <- glmnet(mdlX, mdlY, family = "gaussian", lambda = cv3$lambda.1se, alpha = cv3$alpha)
#coef(md3)

library(caret)
varImp <- function(object, lambda = NULL, ...) {
  
  ## skipping a few lines
  
  beta <- predict(object, s = lambda, type = "coef")
  if(is.list(beta)) {
    out <- do.call("cbind", lapply(beta, function(x) x[,1]))
    out <- as.data.frame(out, stringsAsFactors = TRUE)
  } else out <- data.frame(Overall = beta[,1])
  out <- abs(out[rownames(out) != "(Intercept)",,drop = FALSE])
  out
}
varimportance <- bind_cols(varimportance,varImp(md3, lambda = md3$lambda.min))
}
varimportance$variable <- row.names(varimportance)

varimportance_melt <- melt(varimportance)
colnames(varimportance_melt) <- c('variable','iteration','value')

varimportance_agg <- aggregate(varimportance_melt[,3],by=list(varimportance_melt$variable),FUN = mean)
colnames(varimportance_agg) <- c('variable','value')

pdf(file='results/variableimportance_plot_10CV_50rep_Figure4_square.pdf',width = 16,height=16)
ggplot(varimportance_agg,aes(x=fct_reorder(variable,value),y=value)) + geom_bar(stat='identity', fill = 'deepskyblue4') + coord_flip() + theme_bw() + theme(text=element_text(size=18)) + ylab('Importance') + xlab('')

dev.off()

# fit the selected variables in a linear regression

fitselect <- glm(allres_sites_agg_meta_dummy, formula=count_per_2million ~ lineageAlpha + SiteORF9 + SiteM + SiteE+ SiteORF3a +SiteORF8 + SiteORF7a+ SiteORF6 + ct+ lineageLateDelta + lineageDelta + Vaccination_status_no + Vaccination_status_primary_full  ,family = 'gaussian')
summary(fitselect)


####################################################
# Boxplots
# 

dfplot <- aggregate(count ~ sample + site + variant,data =allres_sites_agg_meta_excl,FUN=function(x){sum(which(x==1))})
dfplot2 <- aggregate(count ~ site + variant,data =dfplot,FUN=sum)
sampsize <- table(unique(allres_sites_agg_meta_excl[,c('sample','variant')])$variant)
sampsize <- t(rbind(names(sampsize),sampsize))
colnames(sampsize) <- c('variant','N')
dfplot2 <- merge(dfplot2,sampsize,by='variant',all.x=T)
dfplot2$freqmissing <- dfplot2$count/dfplot2$N

pdf(file='results/missingness_pattern_supplfig4.pdf',width = 16, height = 9)
ggplot(dfplot2,aes(x=site,y=count/as.numeric(as.character(N)),fill=variant))+geom_bar(stat = 'identity',position = 'dodge') +ggtitle('')  + theme_bw() + ylab('Missingness')+ xlab('') +theme(text = element_text(size=18))+labs(fill='lineage')
dev.off()

# to AUG, hlines at integer counts

pdf(file='results/boxplot_countper2million_AUG_for_figure1_horiz.pdf', width=16,height=10)
plotdata <- allres_sites_agg_meta_excl
ggplot(plotdata,aes(x=reorder(site,desc(site)),y=log2(count_per_2million),fill=variant))+geom_boxplot(outlier.shape = NA)+ geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3),alpha=0.1,size=0.03) +geom_hline(yintercept = log2(1),linetype='dotted') +geom_hline(yintercept = log2(2),linetype='dotted') +geom_hline(yintercept = log2(3),linetype='dotted') +geom_hline(yintercept = log2(4),linetype='dotted') +geom_hline(yintercept = log2(5),linetype='dotted') +geom_hline(yintercept = log2(6),linetype='dotted')+xlab('') +geom_hline(yintercept = log2(7),linetype='dotted') +geom_hline(yintercept = log2(8),linetype='dotted')+geom_hline(yintercept = log2(9),linetype='dotted')+geom_hline(yintercept = log2(10),linetype='dotted')  + theme_classic() + ylab('Count per 2 million reads') +theme(text = element_text(size=18))+labs(fill='lineage')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank() )+ annotate("text",                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 x = length(table(plotdata$site))+0.5, y = log2(1:10), label = as.character(1:10), col = "black", size=4.5, vjust =0,5) + coord_flip()
dev.off()



pdf(file='results/boxplot_countper2million_AUG_for_supplfigure3_alllineages_horiz.pdf',width = 16,height = 10)
plotdata <- allres_sites_agg_meta
ggplot(plotdata,aes(x=reorder(site,desc(site)),y=log2(count_per_2million),fill=variant))+geom_boxplot(outlier.shape = NA)+ geom_point(position = position_jitterdodge(dodge.width = 0.75, jitter.width = 0.3),alpha=0.1,size=0.03) +geom_hline(yintercept = log2(1),linetype='dotted') +geom_hline(yintercept = log2(2),linetype='dotted') +geom_hline(yintercept = log2(3),linetype='dotted') +geom_hline(yintercept = log2(4),linetype='dotted') +geom_hline(yintercept = log2(5),linetype='dotted') +geom_hline(yintercept = log2(6),linetype='dotted')+xlab('') +geom_hline(yintercept = log2(7),linetype='dotted') +geom_hline(yintercept = log2(8),linetype='dotted')+geom_hline(yintercept = log2(9),linetype='dotted')+geom_hline(yintercept = log2(10),linetype='dotted')  + theme_classic() + ylab('Count per 2 million reads') +theme(text = element_text(size=18))+labs(fill='lineage')+theme(axis.text.x=element_blank(),axis.ticks.x=element_blank() )+ annotate("text", x = length(table(plotdata$site))+0.5, y = log2(1:10), label = as.character(1:10), col = "black", size=3.5, vjust =0,5) + coord_flip()
dev.off()

######

#####################
# peak and bump picking algorithm from https://rdrr.io/cran/peakPick/src/R/peakpicking.R
# Weber, C.M., Ramachandran, S., and Henikoff, S. (2014). Nucleosomes are context-specific, H2A.Z-modulated barriers to RNA polymerase. Molecular Cell 53, 819-830.
source('peakpicking.R')

smtroub <-density(log2(allres_sites_agg_meta_excl_alpha_only$count_per_2million))
alphathr <- smtroub$x[which(peakpick(-smtroub$y,neighlim = 1))][1]

smtroub <-density(log2(allres_sites_agg_meta_excl_BA1_only$count_per_2million))
ba1thr <- smtroub$x[which(peakpick(-smtroub$y,neighlim = 1))][1]

smtroub <-density(log2(allres_sites_agg_meta_excl_BA2_only$count_per_2million))
ba2thr <- smtroub$x[which(peakpick(-smtroub$y,neighlim = 1))][1]

smtroub <-density(log2(allres_sites_agg_meta_excl_BA5_only$count_per_2million))
ba5thr <- smtroub$x[which(peakpick(-smtroub$y,neighlim = 1))][1]

low_abundance_N_alpha <-allres_sites_agg_meta_excl[which(log2(allres_sites_agg_meta_excl$count_per_2million) < alphathr & allres_sites_agg_meta_excl$variant=='Alpha' & allres_sites_agg_meta_excl$site == 'ORF9'),'sample']
low_abundance_M_BA1 <-allres_sites_agg_meta_excl[which(log2(allres_sites_agg_meta_excl$count_per_2million) < ba1thr & allres_sites_agg_meta_excl$variant=='BA.1' & allres_sites_agg_meta_excl$site == 'M'),'sample']
low_abundance_M_BA2 <-allres_sites_agg_meta_excl[which(log2(allres_sites_agg_meta_excl$count_per_2million) < ba2thr & allres_sites_agg_meta_excl$variant=='BA.2' & allres_sites_agg_meta_excl$site == 'M'),'sample']
low_abundance_M_BA5 <-allres_sites_agg_meta_excl[which(log2(allres_sites_agg_meta_excl$count_per_2million) < ba5thr & allres_sites_agg_meta_excl$variant=='BA.5' & allres_sites_agg_meta_excl$site == 'M'),'sample']

allres_sites_agg_meta_excl$low_abundance_N_alpha <- ifelse(allres_sites_agg_meta_excl$sample %in% low_abundance_N_alpha,1,0)
allres_sites_agg_meta_excl$low_abundance_N_alpha <- factor(allres_sites_agg_meta_excl$low_abundance_N_alpha,levels=c("1","0"))
allres_sites_agg_meta_excl$low_abundance_M_BA1 <- ifelse(allres_sites_agg_meta_excl$sample %in% low_abundance_M_BA1,1,0)
allres_sites_agg_meta_excl$low_abundance_M_BA1 <- factor(allres_sites_agg_meta_excl$low_abundance_M_BA1,levels=c("1","0"))
allres_sites_agg_meta_excl$low_abundance_M_BA2 <- ifelse(allres_sites_agg_meta_excl$sample %in% low_abundance_M_BA2,1,0)
allres_sites_agg_meta_excl$low_abundance_M_BA2 <- factor(allres_sites_agg_meta_excl$low_abundance_M_BA2,levels=c("1","0"))

allres_sites_agg_meta_excl$low_abundance_M_BA5 <- ifelse(allres_sites_agg_meta_excl$sample %in% low_abundance_M_BA5,1,0)
allres_sites_agg_meta_excl$low_abundance_M_BA5 <- factor(allres_sites_agg_meta_excl$low_abundance_M_BA5,levels=c("1","0"))

allres_sites_agg_meta_excl_BA1_only <- allres_sites_agg_meta_excl[which(allres_sites_agg_meta_excl$variant == 'BA.1' & allres_sites_agg_meta_excl$site =='M'),]
allres_sites_agg_meta_excl_BA2_only <- allres_sites_agg_meta_excl[which(allres_sites_agg_meta_excl$variant == 'BA.2' & allres_sites_agg_meta_excl$site =='M'),]
allres_sites_agg_meta_excl_BA5_only <- allres_sites_agg_meta_excl[which(allres_sites_agg_meta_excl$variant == 'BA.5' & allres_sites_agg_meta_excl$site =='M'),]
allres_sites_agg_meta_excl_alpha_only <- allres_sites_agg_meta_excl[which(allres_sites_agg_meta_excl$variant == 'Alpha' & allres_sites_agg_meta_excl$site =='ORF9'),]

library(ggpubr)
library(cowplot)

xplot <- ggboxplot(unique(allres_sites_agg_meta_excl_alpha_only), x="low_abundance_N_alpha", y="ct", 
                   color = "low_abundance_N_alpha", fill = "low_abundance_N_alpha", palette = c("#E7B800","#00AFBB"),
                   alpha = 0.5, ggtheme = theme_bw())+theme(text = element_text(size=16),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.y = element_text(size=11)) + rotate() + xlab("Subgenomic RNA level")+ ylab('')+ 
  scale_x_discrete(labels = c('Low','high'))+ stat_compare_means(label =  "p.signif", label.x = 2)

sp <- ggplot(unique(allres_sites_agg_meta_excl_alpha_only[which(allres_sites_agg_meta_excl_alpha_only$variant=='Alpha'),c('ct','variant','sample','count_per_2million')]),aes(x=ct,y=log2(count_per_2million)))+stat_density2d_filled()+theme_bw()+labs(color='Lineage')+theme(text = element_text(size=20))+ggtitle("A) Alpha lineage, subgenomic Orf9") + xlab('Ct')+ylab("log2(count per 2M)") + geom_hline(yintercept = alphathr,linetype='dotted')
sp <- sp + rremove("legend")


xplot <- xplot + rremove("legend")+ rremove("legend")


xplot2 <- ggboxplot(unique(allres_sites_agg_meta_excl_BA5_only), x="low_abundance_M_BA5", y="ct", 
                   color = "low_abundance_M_BA5", fill = "low_abundance_M_BA5", palette = c("#E7B800","#00AFBB"),
                   alpha = 0.5, ggtheme = theme_bw())+theme(text = element_text(size=18),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + rotate() + xlab("")+ ylab('')+ 
  scale_x_discrete(labels = c('Low','High')) + stat_compare_means(label =  "p.signif", label.x = 2)

sp2 <- ggplot(unique(allres_sites_agg_meta_excl_BA5_only[which(allres_sites_agg_meta_excl_BA5_only$variant=='BA.5'),c('ct','variant','sample','count_per_2million')]),aes(x=ct,y=log2(count_per_2million)))+stat_density2d_filled()+theme_bw()+labs(color='Lineage')+theme(text = element_text(size=20))+ggtitle("B) BA.5 lineage, subgenomic M") + xlab('Ct')+ylab("") + geom_hline(yintercept = ba5thr,linetype='dotted')
sp2 <- sp2 + rremove("legend")



xplot2 <- xplot2 + rremove("legend")+ rremove("legend")



xplot3 <- ggboxplot(unique(allres_sites_agg_meta_excl_BA2_only), x="low_abundance_M_BA2", y="ct", 
                    color = "low_abundance_M_BA2", fill = "low_abundance_M_BA2", palette = c("#E7B800","#00AFBB"),
                    alpha = 0.5, ggtheme = theme_bw())+theme(text = element_text(size=18),axis.ticks.x=element_blank(),axis.text.x=element_blank(),axis.title.y = element_text(size=11)) + rotate() + xlab("Subgenomic RNA level")+ ylab('')+ 
  scale_x_discrete(labels = c('Low','High')) + stat_compare_means(label =  "p.signif", label.x = 2)

sp3 <- ggplot(unique(allres_sites_agg_meta_excl_BA2_only[which(allres_sites_agg_meta_excl_BA2_only$variant=='BA.2'),c('ct','variant','sample','count_per_2million')]),aes(x=ct,y=log2(count_per_2million)))+stat_density2d_filled()+theme_bw()+labs(color='Lineage')+theme(text = element_text(size=20))+ggtitle("C) BA.2 lineage, subgenomic M") + xlab('Ct')+ylab("log2(count per 2M)") + geom_hline(yintercept = ba2thr,linetype='dotted')
sp3 <- sp3 + rremove("legend")



xplot3 <- xplot3 + rremove("legend")+ rremove("legend")


xplot4 <- ggboxplot(unique(allres_sites_agg_meta_excl_BA1_only), x="low_abundance_M_BA1", y="ct", 
                    color = "low_abundance_M_BA1", fill = "low_abundance_M_BA1", palette = c("#E7B800","#00AFBB"),
                    alpha = 0.5, ggtheme = theme_bw())+theme(text = element_text(size=16),axis.ticks.x=element_blank(),axis.text.x=element_blank()) + rotate() + xlab("")+ ylab('')+ 
  scale_x_discrete(labels = c('Low','High')) + stat_compare_means(label =  "p.signif", label.x = 2)

sp4 <- ggplot(unique(allres_sites_agg_meta_excl_BA1_only[which(allres_sites_agg_meta_excl_BA1_only$variant=='BA.1'),c('ct','variant','sample','count_per_2million')]),aes(x=ct,y=log2(count_per_2million)))+stat_density2d_filled()+theme_bw()+labs(color='Lineage')+theme(text = element_text(size=20))+ggtitle("D) BA.1 lineage, subgenomic M") + xlab('Ct')+ylab("") + geom_hline(yintercept = ba1thr,linetype='dotted')
sp4 <- sp4 + rremove("legend")



#xplot <- xplot + clean_theme() + rremove("legend")
xplot4 <- xplot4 + rremove("legend")+ rremove("legend")



pdf(file='results/Ct_vs_count_alpha_omicron_Figure5_v3.pdf',width = 16,height = 16)

plot_grid(sp,sp2, NULL,xplot, xplot2, NULL,sp3,sp4,NULL,xplot3,xplot4,NULL, ncol = 3, align = "hv", rel_widths = c(10, 10,1), rel_heights = c(4, 1,4,1))

dev.off()


## output summary statistics


sgsite <- c('S','ORF3a','E','M','ORF6','ORF7a','ORF8','ORF9')
lineagegrp <- c('Alpha','Delta','Late Delta','BA.1','BA.2','BA.5')

outtable <- c()
for (s in sgsite){
  for (l in lineagegrp){
    
    
    sumtmp <- summary(as.numeric(allres_sites_agg_meta_excl$count_per_2million[allres_sites_agg_meta_excl$site == s & allres_sites_agg_meta_excl$variant == l]))
    
    
    sdtmp <- sd(as.numeric(allres_sites_agg_meta_excl$count_per_2million[allres_sites_agg_meta_excl$site == s & allres_sites_agg_meta_excl$variant == l]))
    
    outtable <- rbind(outtable,c(as.character(s),as.character(l),sprintf("%.3f", sumtmp[3]),sprintf("%.3f", sdtmp),sprintf("%.3f", sumtmp[1]),sprintf("%.3f", sumtmp[6])))
    
    
  }
}
colnames(outtable) <- c('Subgenomic RNA','Lineage','Median(count per 2 million)','Standard Deviation','Min','Max')
View(outtable)

write.table(file='results/boxplot_summarystats_suppltable4.tsv',outtable,col.names = T,row.names = F,sep='\t',quote = F)
 
