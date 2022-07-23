##加载包
library('tidyverse')
options(stringsAsFactors = F)
##读取信息
train <- read.table('riskTrain.txt',sep = '\t',row.names = 1,header = T)
###表型
fileEncoding='UCS-2LE'
clinical<- read.table('TCGA-UCEC.GDC_phenotype.tsv',sep = '\t',header = T,quote = '',fill = T)
colnames(clin_data)
clin_data <- clinical %>% 
  dplyr::select(patient_ID = submitter_id,
                submitter_id.samples,
                age = age_at_diagnosis.diagnoses,
                Stage = clinical_stage,
                Grade= neoplasm_histologic_grade
                           ) 
clin_data$age <- (as.numeric(clin_data$age))/365
#####
clin_data <- clin_data%>% 
dplyr::filter(str_detect(submitter_id.samples, "-01A")) %>% 
  column_to_rownames('patient_ID')
#########取交集
samSample=intersect(row.names(clin_data),row.names(train ))
train =train [samSample,]
clin_data=clin_data[samSample,]
input <- cbind(train,clin_data)
tab <- data.frame(row.names(input),input)
write.table(tab ,file="tab.xls",sep="\t",row.names=F,quote=F)
#临床相关性分析，输出图形结果
library(limma)
library(ggpubr)
rt <- read.table('iput.txt', header=T, sep="\t", check.names=F, row.names=1)
for(clinical in colnames(rt[,2:ncol(rt)])){
  data=rt[c("riskScore", clinical)]
  colnames(data)=c("riskScore", "clinical")
  data=data[(data[,"clinical"]!="unknow"),]
  #设置比较组
  group=levels(factor(data$clinical))
  data$clinical=factor(data$clinical, levels=group)
  comp=combn(group,2)
  my_comparisons=list()
  for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
  #绘制箱线图
  boxplot=ggboxplot(data, x="clinical", y="riskScore", color="clinical",
                    xlab=clinical,
                    ylab="Risk score",
                    legend.title=clinical,
                    add = "jitter")+ 
    stat_compare_means(comparisons = my_comparisons)
  #输出图片
  pdf(file=paste0(clinical, ".pdf"), width=5.5, height=5)
  print(boxplot)
  dev.off()
}