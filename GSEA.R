# 清空变量
remove(list = ls())

# 加载所需R包
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(tidyverse)
library(limma)
rm(list=ls())  
options(stringsAsFactors = F)

library(clusterProfiler)
library(org.Hs.eg.db)
library(tidyverse)

##deg
rt <- read.table('riskTrain.txt',header = T)[,c(1,8)]
exp <- read.table('uniq.symbol.txt',header = T,row.names = 1)
group_list <- rt$risk
group_list = factor(group_list,
                    levels = c("low","high"))
group_list
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
colnames(deg)
write.table(deg, file='results.txt',
            sep = "\t",row.names = T,col.names = NA,quote = F)


input <- read.table("results.txt",sep="\t",header=T,check.names=F)[,c(1,2)]
names(input)[1] <- 'id'
head(input)
genelist <- bitr(input$id , fromType="SYMBOL",
                 toType="ENTREZID", OrgDb='org.Hs.eg.db')
DEG <- inner_join(input,genelist,by=c("id"="SYMBOL"))



dim(input)

# 按照logFC值对基因进行排序
## 1: 提取logFC值，并储存在一个向量中
geneList =DEG[,2]
## 2: 对geneList进行命名
names(geneList) = as.character(DEG[,3])
head(geneList)
## 基因名：ENTREZID
# 4312     8318    10874    55143    55388      991 
## 基因对应的logFC值
# 4.572613 4.514594 4.418218 4.144075 3.876258 3.677857 

## 3: 根据logFC值降序排列
geneList = sort(geneList, decreasing = TRUE)

########## GO的GSEA富集分析：gseGO ##########
gsego <- gseGO(geneList     = geneList,#排序后的基因列表，一般根据logFC进行排序
              OrgDb        = org.Hs.eg.db,
              ont          = "ALL",#可选择bp.MF,CC,ALL
           #   nPerm        = 1000,#置换检验的次数，默认1000，保持默认即可
            #  minGSSize    = 100,#最小基因集的基因数
             # maxGSSize    = 500,#最大基因集的基因数
              pvalueCutoff = 0.05,#p值的阈值
              verbose      = FALSE)#是否输出提示信息，默认为false
head(gsego)

library(clusterProfiler)
library(enrichplot)
library(ReactomePA)



gseaplot2(gsego,
          1:5, #绘制前5个
          ) # 






# 保存GO的GSEA分析结果到文件
write.table(gsego,file="GSEA_GO_result.txt",sep="\t",
            quote=F,row.names = F)

########## KEGG的GSEA富集分析：gseKEGG ##########
gsekk <- gseKEGG(geneList     = geneList,
               organism     = 'hsa',
             #  nPerm        = 1000,
               #minGSSize    = 50,
           pvalueCutoff =0.05,
               verbose      = T)
head(gsekk)

# 保存KEGG的GSEA分析结果到文件
write.table(gsekk,file="GSEA_KEGG_result.txt",sep="\t",
            quote=F,row.names = F)
gseaplot2(gsekk,
          1:5, #绘制前3个
          ) 
