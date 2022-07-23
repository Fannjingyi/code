library(tidyverse)
library(survival)

rt <- read.table('lassoSigExp.txt',sep = '\t',header = T)
set.seed(2020)
library(caret)
inTrain<-createDataPartition(y=rt[,3],p=0.75,list=F)
train<-rt[inTrain,]

write.table(train,file="train.txt",sep="\t",quote=F,row.names=F)
write.table(test,file="test.txt",sep="\t",quote=F,row.names=F)

rt <- read.table('lassoSigExp.txt',sep = '\t',header = T)
set.seed(2020)
library(caret)
inTrain<-createDataPartition(y=rt[,3],p=0.75,list=F)
train<-rt[inTrain,]
test<-rt[-inTrain,]
write.table(train,file="train.txt",sep="\t",quote=F,row.names=F)
write.table(test,file="test.txt",sep="\t",quote=F,row.names=F)
rt <- read.table('lassoSigExp.txt',sep = '\t',header = T)
set.seed(2020)
library(caret)
inTrain<-createDataPartition(y=rt[,3],p=0.75,list=F)
train<-rt[inTrain,]
test<-rt[-inTrain,]
write.table(train,file="train.txt",sep="\t",quote=F,row.names=F)
write.table(test,file="test.txt",sep="\t",quote=F,row.names=F)

multiCox=coxph(Surv(OS.time,OS) ~ ., data = train)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)

multiCox=coxph(Surv(OS.time,OS) ~ ., data = train)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
train <-train[,-1]
multiCox=coxph(Surv(OS.time,OS) ~ ., data = train)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
outTab=data.frame()
outTab=cbind(
coef=multiCoxSum$coefficients[,"coef"],
HR=multiCoxSum$conf.int[,"exp(coef)"],
HR.95L=multiCoxSum$conf.int[,"lower .95"],
HR.95H=multiCoxSum$conf.int[,"upper .95"],
pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)
ggforest(multiCox,
main = "Hazard ratio",
cpositions = c(0.02,0.22, 0.4),
fontsize = 0.7,
refLabel = "reference",
noDigits = 2)
outTab=data.frame()
outTab=cbind(
coef=multiCoxSum$coefficients[,"coef"],
HR=multiCoxSum$conf.int[,"exp(coef)"],
HR.95L=multiCoxSum$conf.int[,"lower .95"],
HR.95H=multiCoxSum$conf.int[,"upper .95"],
pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)