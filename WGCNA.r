####wgcna
####6个预后相关的
rm(list = ls())
trait <- read.table('GDF1011.txt',sep="\t",header=T,check.names=F,row.names = 1)
trait <-as.data.frame( t(trait))
##特异性位点
c <- c(
'cg02974931',
'cg04110601',
'cg07773116',
'cg14720763',
'cg20186445',
'cg24636477')

trait <- trait[,c]
###############
library(tidyverse)
gene_info <- read_tsv("gencode.gene.info.v22.tsv")

gene_info_TCGA <- gene_info %>%
  dplyr::select("gene_id", "gene_name", "gene_type") 
FPKM <- read.table('TCGA-UCEC.htseq_fpkm.tsv',row.names = 1,header = T)
FPKM1 <- FPKM%>%
  rownames_to_column("gene_id")%>%
  inner_join(gene_info_TCGA, "gene_id")%>%
  dplyr::select("gene_id",  "gene_name", "gene_type", everything())
table(FPKM1$gene_type)
####################################################
mRNA_FPKM <- FPKM1 %>%
  subset(gene_type == "protein_coding") %>%
  dplyr::select(-c(1,3))%>%
  aggregate(. ~ gene_name, . , mean) 
######################################################################################################
 row.names(mRNA_FPKM) <- mRNA_FPKM[,1]

######################################################################################################
mRNA_FPKM <-mRNA_FPKM[,-1]
mRNA_FPKM <-t(mRNA_FPKM)
mRNA_FPKM <- as.data.frame(mRNA_FPKM)
row.names(mRNA_FPKM) <- str_replace_all(row.names(mRNA_FPKM),"\\.","-")
save(mRNA_FPKM,file = 'mRNA.data')
##############################################################################################
##取交集转录组
samSample=intersect(row.names(mRNA_FPKM),row.names(trait))
trait=trait[samSample,]
mRNA_FPKM=mRNA_FPKM[samSample,]
#########################################
library(WGCNA)

datExpr0 <- mRNA_FPKM
##check missing value
gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

###############################Sample cluster##########样品聚类#################### 
sampleTree = hclust(dist(datExpr0), method = "average")
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.


#sizeGrWindow(12,9)

pdf(file = "1_sampleClustering.pdf", width =30, height = 9)
par(cex = 0.6)
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

### Plot a line to show the cut
#abline(h = 160, col = "red")##剪切高度不确定，故无红线

dev.off()
##############所以这一步不一定能够做，剪切高度问题,这个根据实际设置后可用


### Determine cluster under the line
#clust = cutreeStatic(sampleTree, cutHeight =160, minSize = 10)
#table(clust)


### clust 1 contains the samples we want to keep.
#keepSamples = (clust==1)
#datExpr0 = datExpr0[keepSamples, ]


############### 载入性状数据## input trait data###############
#Loading clinical trait data
traitData = trait########trait file name can be changed######性状数据文件名，根据实际修改，如果工作路径不是实际性状数据路径，需要添加正确的数据路径
dim(traitData)
names(traitData)
# remove columns that hold information we do not need.
allTraits = traitData
dim(allTraits)
names(allTraits)

# Form a data frame analogous to expression data that will hold the clinical traits.
fpkmSamples = rownames(datExpr0)
traitSamples =rownames(allTraits)
traitRows = match(fpkmSamples, traitSamples)
datTraits = allTraits[traitRows,]
rownames(datTraits) 
collectGarbage()

# Re-cluster samples
sampleTree2 = hclust(dist(datExpr0), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(datTraits, signed = FALSE)
# Plot the sample dendrogram and the colors underneath.

#sizeGrWindow(12,12)
pdf(file="2_Sample dendrogram and trait heatmap.pdf",width=30,height=12)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")
dev.off()


##################################
#save(datExpr0, file = "fpkm_forAnalysis.RData")
#save(datTraits, file="trait_forAnalysis.RData")





#############################network constr########################################

# Allow multi-threading within WGCNA. At present this call is necessary.
# Any error here may be ignored but you may want to update WGCNA if you see one.
# Caution: skip this line if you run RStudio or other third-party R environments.
# See note above.
enableWGCNAThreads(48)


# Choose a set of soft-thresholding powers
powers = c(1:20)

# Call the network topology analysis function
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)

# Plot the results:
#sizeGrWindow(9, 5)
pdf(file="3_Scale independence.pdf",width=9,height=5)
par(mfrow = c(1,2))
cex1 = 0.9
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.80,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()


######chose the softPower


softPower =8  ##大于0.8
adjacency = adjacency(datExpr0, power = softPower)

##### Turn adjacency into topological overlap
TOM = TOMsimilarity(adjacency);
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)

#sizeGrWindow(12,9)
pdf(file="4_Gene clustering on TOM-based dissimilarity.pdf",width=12,height=9)
plot(geneTree, xlab="", sub="", main = "Gene clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)
dev.off()


# We like large modules, so we set the minimum module size relatively high:
minModuleSize = 30
# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 3, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
pdf(file="5_Dynamic Tree Cut.pdf",width=8,height=6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()


# Calculate eigengenes
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average")
# Plot the result
#sizeGrWindow(7, 6)
pdf(file="6_Clustering of module eigengenes.pdf",width=7,height=6)
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
MEDissThres = 0.1######剪切高度可修改
# Plot the cut line into the dendrogram
abline(h=MEDissThres, col = "red")
dev.off()


# Call an automatic merging function
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# The merged module colors
mergedColors = merge$colors
# Eigengenes of the new merged modules:
mergedMEs = merge$newMEs

#sizeGrWindow(12, 9)
pdf(file="7_merged dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
######################################
pdf(file="7.1_merged dynamic.pdf", width = 9, height = 6)
plotDendroAndColors(geneTree,  mergedColors,
                     "Colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()




#######################################
# Rename to moduleColors
moduleColors = mergedColors
# Construct numerical labels corresponding to the colors
colorOrder = c("grey", standardColors(50))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs

# Save module colors and labels for use in subsequent parts
#save(MEs, TOM, dissTOM,  moduleLabels, moduleColors, geneTree, sft, file = "networkConstruction-stepByStep.RData")




##############################relate modules to external clinical triats######################################
# Define numbers of genes and samples
nGenes = ncol(datExpr0)
nSamples = nrow(datExpr0)

moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

#sizeGrWindow(10,6)
pdf(file="8_Module-trait relationships.pdf",width=15,height=30)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")

dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3))

# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


######## Define variable weight containing all column of datTraits

###MM and GS


# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="")

#names of those trait
traitNames=names(datTraits)

geneTraitSignificance = as.data.frame(cor(datExpr0, datTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))

names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")


####plot MM vs GS for each trait vs each module


##########example:royalblue and CK
#module="royalblue"
#column = match(module, modNames)
#moduleGenes = moduleColors==module

#trait="CK"
#traitColumn=match(trait,traitNames)

#sizeGrWindow(7, 7)

#par(mfrow = c(1,1))
#verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
abs(geneTraitSignificance[moduleGenes, traitColumn]),
xlab = paste("Module Membership in", module, "module"),
ylab = paste("Gene significance for ",trait),
main = paste("Module membership vs. gene significance\n"),
cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
######

for (trait in traitNames){
  traitColumn=match(trait,traitNames)
  
  for (module in modNames){
    column = match(module, modNames)
    moduleGenes = moduleColors==module
    
    if (nrow(geneModuleMembership[moduleGenes,]) > 1){####进行这部分计算必须每个模块内基因数量大于2，由于前面设置了最小数量是30，这里可以不做这个判断，但是grey有可能会出现1个gene,它会导致代码运行的时候中断，故设置这一步
      
      #sizeGrWindow(7, 7)
      pdf(file=paste("9_", trait, "_", module,"_Module membership vs gene significance.pdf",sep=""),width=7,height=7)
      par(mfrow = c(1,1))
      verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                         abs(geneTraitSignificance[moduleGenes, traitColumn]),
                         xlab = paste("Module Membership in", module, "module"),
                         ylab = paste("Gene significance for ",trait),
                         main = paste("Module membership vs. gene significance\n"),
                         cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
      dev.off()
    }
  }
}

#####
names(datExpr0)
probes = names(datExpr0)


#################export GS and MM############### 

geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)

for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]

write.table(geneInfo, file = "10_GS_and_MM.xls",sep="\t",row.names=F)

