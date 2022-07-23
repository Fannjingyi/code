rm(list = ls())
library(stringr)
library(tidyverse)
met <- read.table('GDF1011.txt',sep = '\t',row.names = 1,header = T)
met <- t(met)
met1 <- as.data.frame(met)
row.names(met1) <-  str_replace_all(row.names(met),"\\.","-")


met <- data.frame(patient_ID = str_sub(row.names(met1), 1, 12),met1)  %>%
       dplyr::select(patient_ID, everything())%>%
dplyr::filter(str_detect(row.names(met1), "-01"))
#####################
all <- met
  
  colnames(all)[1] <- 'Symbol'
#平均数
all[,2:ncol(all)] <- apply(all[,2:ncol(all)],2,function(x){as.numeric(x)})
all.unique <- aggregate(.~Symbol,all,mean)
all.final <- all.unique[,-1] 
rownames(all.final) <- all.unique$Symbol
met <- all.final
##########################################
sur <- read.table('TCGA-UCEC.survival.tsv',sep = '\t',header = T)%>%
  dplyr::filter(str_detect(sample, "-01")) %>%
  dplyr::select(X_PATIENT, everything())%>%select(-2)
all <- sur
#转换ID,只保留symbol
colnames(all)[1] <- 'Symbol'
#平均数
all[,2:ncol(all)] <- apply(all[,2:ncol(all)],2,function(x){as.numeric(x)})
all.unique <- aggregate(.~Symbol,all,mean)
all.final <- all.unique[,-1] 
rownames(all.final) <- all.unique$Symbol
sur <- all.final
sur$OS.time <- sur$OS.time/365
##提取患病样本
samSample=intersect(row.names(sur),row.names(met))
sur=sur[samSample,]
met <- met[samSample,]
surl <- cbind(sur,met )
################批量画生存曲线++
#单个基因
library(survival)
library(survminer)
gs=names(surl)[3:18]
  splots <- lapply(gs, function(g){
  surl$gene = ifelse(surl[,g] > median(surl[,g]),'high','low')
  sfit1=survfit(Surv(OS.time,OS)~gene, data=surl)
  
  p <- ggsurvplot(sfit1,pval =TRUE, data = surl, risk.table = TRUE
                  ,conf.int = TRUE
                  ,title=g
                  
  )
  pdf(paste0(g, "_surv.pdf"),width = 8, height =10)
  print(p, newpage = FALSE)
  dev.off()

  }) 












