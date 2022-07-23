rm(list = ls())
library(tidyverse)
gene <- read.table('comm.xls',sep = '\t')
gene1 <- as.data.frame(t(gene))
colnames(gene1) <-gene1[1,]
gene1 <- gene1[-1,]
set.seed(2020)
table(str_sub(gene1[,1],14,15))
group_list = ifelse(as.numeric(str_sub(gene1[,1],14,15)) < 10,'tumor','normal')
group_list = factor(group_list,levels = c("normal","tumor"))
table(group_list)
write_tsv(gene1,'gene1.xls')


tumor <- read.table('tumor.txt',header = T,row.names  = 1,sep = '\t')
###取交集有生存信息的
sur <-read.table('TCGA-UCEC.survival.tsv',sep = '\t',row.names = 1,header = T) 
###
samSample=intersect(row.names(tumor),row.names(sur))
tumor =tumor[samSample,]
sur=sur[samSample,]
###
mul <- cbind(sur,tumor)
#####
mul[,"OS.time"]=mul[,"OS.time"]/365
mul <- data.frame(row.names(mul),mul)
##############取平均值
mul <- mul %>%mutate(patient_ID = str_sub(mul$row.names.mul., 1, 12)) %>%
  dplyr::select(patient_ID, everything())
####################取平均值


library(glmnet)
library(survival)
v1<-as.matrix(mul[,c(3:ncol(mul))])
v2 <- as.matrix(Surv(mul$OS.time,mul$OS))
myfit <- glmnet(v1, v2, family = "cox")
plot(myfit, xvar = "lambda", label = TRUE)
myfit2 <- cv.glmnet(v1, v2, family="cox")
plot(myfit2)
abline(v=log(c(myfit2$lambda.min,myfit2$lambda.1se)),lty="dashed")
myfit2$lambda.min
coe <- coef(myfit, s = myfit2$lambda.min)
act_index <- which(coe != 0)
act_coe <- coe[act_index]
row.names(coe)[act_index]
##rm(list = ls())




###############################################cox多因素



a <- row.names(coe)[act_index]
library(survival)
mul1 <- mul[,a]
mul1 <- cbind(sur,mul1)
mul1[,"OS.time"]=mul1[,"OS.time"]/365
rt <- mul1
#划分训练集和测试集
set.seed(2020)
library(caret)
inTrain<-createDataPartition(y=rt[,3],p=0.75,list=F)
train<-rt[inTrain,]
test<-rt[-inTrain,]
write.table(train,file="train.txt",sep="\t",quote=F,row.names=F)
write.table(test,file="test.txt",sep="\t",quote=F,row.names=F)













##################################
###########################################################多因素cox

multiCox=coxph(Surv(OS.time,OS) ~ ., data = train)
multiCox=step(multiCox,direction = "both")
multiCoxSum=summary(multiCox)
###########################################################



#输出模型相关信息
outTab=data.frame()
outTab=cbind(
  coef=multiCoxSum$coefficients[,"coef"],
  HR=multiCoxSum$conf.int[,"exp(coef)"],
  HR.95L=multiCoxSum$conf.int[,"lower .95"],
  HR.95H=multiCoxSum$conf.int[,"upper .95"],
  pvalue=multiCoxSum$coefficients[,"Pr(>|z|)"])
outTab=cbind(id=row.names(outTab),outTab)
write.table(outTab,file="multiCox.xls",sep="\t",row.names=F,quote=F)
#绘制森林图
pdf(file="forest.pdf",
    width = 8,             #图片的宽度
    height = 5,            #图片的高度
)
ggforest(multiCox,
         main = "Hazard ratio",
         cpositions = c(0.02,0.22, 0.4), 
         fontsize = 0.7, 
         refLabel = "reference", 
         noDigits = 2)
dev.off()


################################################
#输出train组风险文件
riskScore=predict(multiCox,type="risk",newdata=train)           #利用train得到模型预测train样品风险
coxGene=rownames(multiCoxSum$coefficients)
coxGene=gsub("`","",coxGene)
outCol=c("OS.time","OS",coxGene)
medianTrainRisk=median(riskScore)
risk=as.vector(ifelse(riskScore>medianTrainRisk,"high","low"))
write.table(cbind(id=rownames(cbind(train[,outCol],riskScore,risk)),cbind(train[,outCol],riskScore,risk)),
            file="riskTrain.txt",
            sep="\t",
            quote=F,
            row.names=F)

#输出test组风险文件
rtTest=test   #读取test输入文件

riskScoreTest=predict(multiCox,type="risk",newdata=rtTest)      #利用train得到模型预测test样品风险
riskTest=as.vector(ifelse(riskScoreTest>medianTrainRisk,"high","low"))
write.table(cbind(id=rownames(cbind(rtTest[,outCol],riskScoreTest,riskTest)),cbind(rtTest[,outCol],riskScore=riskScoreTest,risk=riskTest)),
            file="riskTest.txt",
            sep="\t",
            quote=F,
            row.names=F)









###############################################cox

library(survival)

#绘制train组生存曲线
rt=read.table("riskTrain.txt",header=T,sep="\t")
diff=survdiff(Surv(OS.time,OS) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(OS.time,OS) ~ risk, data = rt)
summary(fit)    #查看五年生存率
pdf(file="survivalTrain.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()

#绘制test组生存曲线
rt=read.table("riskTest.txt",header=T,sep="\t")
diff=survdiff(Surv(OS.time,OS) ~risk,data = rt)
pValue=1-pchisq(diff$chisq,df=1)
pValue=signif(pValue,4)
pValue=format(pValue, scientific = TRUE)
fit <- survfit(Surv(OS.time,OS) ~ risk, data = rt)
summary(fit)    #查看五年生存率
pdf(file="survivalTest.pdf",width=5.5,height=5)
plot(fit, 
     lwd=2,
     col=c("red","blue"),
     xlab="Time (year)",
     ylab="Survival rate",
     main=paste("Survival curve (p=", pValue ,")",sep=""),
     mark.time=T)
legend("topright", 
       c("high risk", "low risk"),
       lwd=2,
       col=c("red","blue"))
dev.off()
################################################################ROC曲线
###install.packages("survivalROC")
library(survivalROC)
rt=read.table("riskTrain.txt",header=T,sep="\t",check.names=F,row.names=1)
pdf(file="rocTrain.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()

rt=read.table("riskTest.txt",header=T,sep="\t",check.names=F,row.names=1)
pdf(file="rocTest.pdf",width=6,height=6)
par(oma=c(0.5,1,0,1),font.lab=1.5,font.axis=1.5)
roc=survivalROC(Stime=rt$OS.time, status=rt$OS, marker = rt$riskScore, 
                predict.time =1, method="KM")
plot(roc$FP, roc$TP, type="l", xlim=c(0,1), ylim=c(0,1),col='red', 
     xlab="False positive rate", ylab="True positive rate",
     main=paste("ROC curve (", "AUC = ",round(roc$AUC,3),")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)
abline(0,1)
dev.off()
###################################################

library(pheatmap)

rt=read.table("riskTrain.txt",sep="\t",header=T,row.names=1,check.names=F)      #读取train输入文件
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)

annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmapTrain.pdf",width = 12,height = 5)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         fontsize_col=3,
         show_colnames= F,
         color = colorRampPalette(c("green", "black", "red"))(50) )
dev.off()

#绘制test组风险热图
rt=read.table("riskTest.txt",sep="\t",header=T,row.names=1,check.names=F)      #读取test输入文件
rt=rt[order(rt$riskScore),]
rt1=rt[c(3:(ncol(rt)-2))]
rt1=t(rt1)
annotation=data.frame(type=rt[,ncol(rt)])
rownames(annotation)=rownames(rt)
pdf(file="heatmapTest.pdf",width = 12,height = 5)
pheatmap(rt1, 
         annotation=annotation, 
         cluster_cols = FALSE,
         fontsize_row=11,
         fontsize_col=3,
         show_colnames= F,
         color = colorRampPalette(c("green", "black", "red"))(50) )

dev.off()
###################################################################################
#绘制train组风险图
rt=read.table("riskTrain.txt",header=T,sep="\t",check.names=F,row.names=1)     #读取train输入文件
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file="riskScoreTrain.pdf",width = 12,height = 5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
trainMedianScore=median(rt$riskScore)
abline(h=trainMedianScore,v=lowLength,lty=2)
dev.off()

#绘制test组风险图
rt=read.table("riskTest.txt",header=T,sep="\t",check.names=F,row.names=1)       #读取test输入文件
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
line=rt[,"riskScore"]
line[line>10]=10
pdf(file="riskScoreTest.pdf",width = 12,height = 5)
plot(line,
     type="p",
     pch=20,
     xlab="Patients (increasing risk socre)",
     ylab="Risk score",
     col=c(rep("green",lowLength),
           rep("red",highLength)))
abline(h=trainMedianScore,v=lowLength,lty=2)
dev.off()
##################################################
#绘制train组生存状态图
rt=read.table("riskTrain.txt",header=T,sep="\t",check.names=F,row.names=1)        #读取train输入文件
rt=rt[order(rt$riskScore),]



color=as.vector(rt$OS)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStatTrain.pdf",width = 12,height = 5)
plot(rt$OS.time,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()

#绘制test组生存状态图
rt=read.table("riskTest.txt",header=T,sep="\t",check.names=F,row.names=1)         #读取test输入文件
rt=rt[order(rt$riskScore),]
riskClass=rt[,"risk"]
lowLength=length(riskClass[riskClass=="low"])
highLength=length(riskClass[riskClass=="high"])
color=as.vector(rt$OS)
color[color==1]="red"
color[color==0]="green"
pdf(file="survStatTest.pdf",width = 12,height = 5)
plot(rt$OS.time,
     pch=19,
     xlab="Patients (increasing risk socre)",
     ylab="Survival time (years)",
     col=color)
abline(v=lowLength,lty=2)
dev.off()












