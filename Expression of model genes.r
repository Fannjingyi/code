#######
rm(list=ls())
library(stringr)
library(parallel)
library(tidyverse)
library(ggpubr)      
options(stringsAsFactors = FALSE)
gene <- c('B4GALNT3',	'GREB1'	,'DNAJC22')
##
exp <- read.table('gene_exp.xls',header = T)


tp <- read.table('type.txt',header = T)
a <- exp[gene,]
a <- t(a)
 row.names(a) <- str_replace_all(row.names(a),"\\.","-")

###
a <- data.frame(row.names(a),a)
as.data.frame(a)

names(a)[1] <- names(tp)[1] 
 
rt<- inner_join(a,tp)
###
group=levels(factor(rt$TN))
rt$TN=factor(rt$TN, levels=group)
comp=combn(group,2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
B4GALNT3=ggboxplot(rt, x="TN", y="B4GALNT3", color="TN",
                  xlab="Type",
               
                  legend.title="TN",
                  palette = c("blue","red"),
                  add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons
                ,   method="wilcox.test"
                     )
dev.off()


GREB1=ggboxplot(rt, x="TN", y="GREB1", color="TN",
                   xlab="Type",
                 
                   legend.title="TN",
                   palette = c("blue","red"),
                   add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons
                     ,   method="wilcox.test"
  )





DNAJC22<- ggboxplot(rt, x="TN", y="DNAJC22", color="TN",
         xlab="Type",
        
         legend.title="TN",
         palette = c("blue","red"),
         add = "jitter")+ 
  stat_compare_means(comparisons = my_comparisons
                     ,   method="wilcox.test"
  )

  ##±£´æÍ¼Æ¬
pdf(file='B4GALNT3.pdf', width=6, height=5)
print(B4GALNT3)
dev.off()

pdf(file='DNAJC22.pdf', width=6, height=5)
print(DNAJC22)
dev.off()
pdf(file='GREB1.pdf', width=6, height=5)
print(GREB1)
dev.off()



