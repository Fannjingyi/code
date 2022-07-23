
#install.packages("glmnet")
#install.packages("survival")

load("mRNA.data")
library("glmnet")
library("survival")
library(stringr)
library(tidyverse)
set.seed(2020)
gene <-read.table('all.txt',sep = '\t')[,1] 
rt=mRNA_FPKM    #读取文件
sur <- read.table('TCGA-UCEC.survival.tsv',header = T,row.names = 1)
samSample=intersect(row.names(rt),row.names(sur))
sur=sur[samSample,]
rt=rt[samSample,]
#####################
rt1 <- cbind(rt,sur)
rt1 <- rt1%>%
dplyr::filter(str_detect(row.names(rt), "-01")) %>%
  dplyr::select(X_PATIENT,OS,OS.time, everything())

rt1 <- rt1[,gene]
all <- rt1 

#转换ID,只保留symbol
colnames(all)[1] <- 'Symbol'
#取基因的平均数
all[,2:ncol(all)] <- apply(all[,2:ncol(all)],2,function(x){as.numeric(x)})
all.unique <- aggregate(.~Symbol,all,mean)
all.final <- all.unique[,-1] 
rownames(all.final) <- all.unique$Symbol
rt1 <- all.final
rt <- rt1

##################
rt$OS.time=rt$OS.time/365

set.seed(2020)

x=as.matrix(rt[,c(3:ncol(rt))])
y=data.matrix(Surv(rt$OS.time,rt$OS))

fit <- glmnet(x, y, family = "cox", maxit = 1000)
pdf("lambda.pdf")
plot(fit, xvar = "lambda", label = TRUE)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

cvfit <- cv.glmnet(x, y, family="cox", maxit = 1000)
pdf("cvfit.pdf")
plot(cvfit)
abline(v=log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef <- coef(fit, s = cvfit$lambda.min)
index <- which(coef != 0)
actCoef <- coef[index]
lassoGene=row.names(coef)[index]
lassoGene=c("OS.time","OS",lassoGene)
lassoSigExp=rt[,lassoGene]
lassoSigExp=cbind(id=row.names(lassoSigExp),lassoSigExp)
write.table(lassoSigExp,file="lassoSigExp.txt",sep="\t",row.names=F,quote=F)


