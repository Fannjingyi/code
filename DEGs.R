
rm(list = ls())
set.seed(2020)
library(stringr)
###################
gene_exp <- read.table('gene_exp.xls',sep = '\t',header = T,row.names = 1)

colnames(gene_exp) <- str_replace_all(colnames(gene_exp),"\\.","-")
exp=gene_exp
exp <- as.numeric(unlist(exp))
range(exp)
exp=gene_exp



#write.table(gene_exp, file="gene_exp.xls",row.names=T, col.names=T,quote=FALSE,sep="\t")
sample_info <- read.table('type.txt',sep = '\t',header = T,row.names = 1)

##########差异分析
library(limma)

design <- model.matrix(~ 0 + sample_info$TN)
colnames(design) <- levels(factor(sample_info$TN))
rownames(design) <- rownames(sample_info)

contrasts <- makeContrasts(
 TN=tumor-normal,
  levels = design)



fit <- lmFit(exp,design)
fit <- contrasts.fit(fit, contrasts)
fit <- eBayes(fit)
de_result <- topTable(fit,
                      coef = 'TN',
                      number = Inf) # 保留所有基因
de_result[1:5, c('logFC', 'P.Value', 'adj.P.Val')]

library(tidyverse)
de_result <- rownames_to_column(
  as.data.frame(de_result), var = 'sym')
de_result <-  mutate(de_result,direction = if_else(
      P.Value>=0.05, 'ns', if_else(
      abs(logFC) <= 1, 'ns', if_else(
        logFC > 0.5, 'up', 'down')
    )))%>%
  arrange(desc(abs(logFC)))
###排序
  
table(de_result$direction)
de <-  dplyr::select(de_result,-t, -B) 
#ddown    ns    up 
#1260 17461   991  
write.table(de,'diff.xls',row.names = F,quote = F,sep = '\t',col.names = T)



#################画图
de_exp_top <- dplyr::slice(de_result, 1:50)[,1]
#差异基因表达热图
library(pheatmap)

group_list = factor(sample_info$TN,
                    levels = c("normal","tumor"))
annotation_col=data.frame(group=group_list)
rownames(annotation_col)=row.names(sample_info)
a <- exp[de_exp_top,]
############
s <- read.table('type.txt',sep = '\t',header = T)
t <- subset(s,type==1)[,1]
n <-  subset(s,type==0)[,1]
rt=cbind(exp[,t], exp[,n])
pheatmap(rt[de_exp_top,],
         annotation_col=annotation_col,
         scale = "row",
         show_rownames = T,
         show_colnames = F,
         color = colorRampPalette(c("green","white","red"))(100),
         cluster_cols =F,
         fontsize = 10,
         fontsize_row=10,
         fontsize_col=3)

#######绘制火山图
de_result$logP <- -log10(de_result$P.Value)
library(ggplot2)
p <- ggplot(data = de_result, 
            aes(x = logFC, y = logP)) +
  geom_point(alpha=0.4, size=3.5, 
             aes(color=direction)) +
  ylab("-log10(pvalue)")+
  scale_color_manual(values=c("navy", "grey", "red"))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.8) +
  theme_bw()
