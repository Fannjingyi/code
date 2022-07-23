rm(list=ls())   #清空环境变量
options(stringsAsFactors = F)

#install.packages('e1071')
#install.packages('parallel')

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("preprocessCore", version = "3.13")

#函数下载：https://content.cruk.cam.ac.uk/fmlab/sivakumar2016/Cibersort.R
source("source.R")   #注释文件

#整理基因表达数据
load("CESC_FPKM_tumor.Rda")
ciber_input <- log2(CESC_FPKM_tumor_final + 1)
write.table(ciber_input, file = "cibersort_input.txt",
            sep = "\t", row.names = T,col.names = NA,quote = F)

#CIBERSORT计算
sig_matrix <- "LM22.txt"   #注释文件名
mixture_file = 'uniq.symbol.txt'   #表达数据文件名
res_cibersort <- CIBERSORT(sig_matrix, mixture_file, perm=100, QN=TRUE)
save(res_cibersort,file = "res_cibersort.Rdata")   #保存结果

#可视化展示
rm(list=ls())
load("res_cibersort.Rdata")
res_cibersort <- res_cibersort[,1:22]   #取前22列为细胞丰度数据
ciber.res <- res_cibersort[,colSums(res_cibersort) > 0]   #去除丰度全为0的细胞

#barplot图
mycol <- ggplot2::alpha(rainbow(ncol(ciber.res)), 0.7) #创建彩虹色板（带70%透明度）
par(bty="o", mgp = c(2.5,0.3,0), mar = c(2.1,4.1,2.1,10.1),tcl=-.25,las = 1,xpd = F)
barplot(as.matrix(t(ciber.res)),
        border = NA, # 柱子无边框
        names.arg = rep("",nrow(ciber.res)), # 无横坐标样本名
        yaxt = "n", # 先不绘制y轴
        ylab = "Relative percentage", # 修改y轴名称
        col = mycol) # 采用彩虹色板
axis(side = 2, at = c(0,0.2,0.4,0.6,0.8,1), # 补齐y轴添加百分号
     labels = c("0%","20%","40%","60%","80%","100%"))
legend(par("usr")[2]-20, # 这里-20要根据实际出图的图例位置情况调整
       par("usr")[4], 
       legend = colnames(ciber.res), 
       xpd = T,
       fill = mycol,
       cex = 0.8, 
       border = NA, 
       y.intersp = 1,
       x.intersp = 0.2,
       bty = "n")
dev.off()   #关闭画板

#相关性热图
M <- round(cor(ciber.res),2) # 计算相关性矩阵并保留两位小数

library(corrplot)
corrplot.mixed(M,
               lower.col = "black", #左下方字体颜色为黑色
               tl.pos = "lt",  #标签出现在左侧和顶部
               number.cex = 0.5, #左下方字号为0.5
               tl.cex = 0.5) #标签字号为0.5
dev.off()   #关闭画板




