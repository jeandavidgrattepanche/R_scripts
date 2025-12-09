setwd('/home/tuk61790/Network/')
library('ape')
library('phytools')
library('phyloseqCompanion')
library("ggplot2")
library("plyr")
library(stringr)
#library(factoextra)
library(WGCNA)
source('/home/tuk61790/Network/multiplot.R', chdir = TRUE)
source('/home/tuk61790/Network/get_top_taxa.R', chdir = TRUE)

options(stringsAsFactors = FALSE)

data.all <- read.table('OTUtable_ingroup_100.txt',sep="\t",header=TRUE,row.names=1)
#data.all <- read.table('SPtable.txt')
env.all <- read.table('NBP1910_envdata_v4.5.txt',sep="\t",header=TRUE,row.names=1)

datExpr0 = as.data.frame(t(cbind(data.all[c(0)],data.all[c(10:105)])))
datExpr0[] <- lapply(datExpr0, as.numeric)
annot = as.data.frame((cbind(data.all[c(0:10)])))

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#if TRUE continue

sampleTree = hclust(dist(datExpr0), method = "average")
pdf("results/sample_clust.pdf")
plot(sampleTree)
dev.off()

# if outliner remove it

env.sub <- cbind(env.all[c(0)], env.all[c(7:41)], env.all[c(5)], env.all[c(7)])
traitRows = match(rownames(datExpr0), rownames(env.sub))
env.red = env.sub[traitRows, -1]
env.red[] <- lapply(env.red, as.numeric)
env.reds <- scale(env.red)
names(env.reds) <- names(env.red)
traitColors = numbers2colors(env.reds, signed = TRUE)
pdf("results/sample_clust2.pdf", width =15, height=9)
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(env.reds))
dev.off()

enableWGCNAThreads(20)

powers = c(c(1:10), seq(from = 12, to=30, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5)
pdf("results/test_WGCNA_clust_1.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)", ylab = "Scale Free topology Model Fit, signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], labels=powers,cex=0.9,col="red")
abline(h=0.90,col="red")
dev.off()

pdf("results/test_WGCNA_clust_2.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)", ylab = "Main Connectivity",type="n",main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,cex=0.9,col="red")
dev.off()

datExpr0[] <- lapply(datExpr0, as.numeric)
net2 = blockwiseModules(datExpr0, power=5,TOMType ="unsigned", deepSplit= 0 ,minModuleSize=30, maxBlockSize = 30000, reassignThreshold = 1e-6 , mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "results/test_TOM", verbose= 3, nThreads = 10)


pdf("results/test_WGCNA_clust_3.pdf", width = 15, height= 9)
mergedColors = labels2colors(net2$colors)
plotDendroAndColors(net2$dendrograms[[1]], mergedColors[net2$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang =0.05)
dev.off()


moduleLabels = net2$colors
moduleColors = labels2colors(net2$colors)
MEs = net2$MEs
geneTree = net2$dendrograms[[1]]
nSamples = nrow(datExpr0)
nTaxa = ncol(datExpr0)
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
nMEs2 = as.data.frame(t(cbind(table(moduleEigengenes(datExpr0, moduleColors)$validColors))))
MEs2 = orderMEs(MEs0) 
moduleTraitCor = cor(MEs2, env.red, use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("results/test_WGCNA_clust_4.pdf", width = 20, height= 36)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(env.red), yLabels = names(MEs2), ySymbols = names(MEs2), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins= FALSE, cex.text=0.5, zlim= c(-1,1), main = paste("Module-Environment relationships"))
dev.off()

SST = as.data.frame(env.red$ice)
names(SST) = "icecoverage"
modNames = substring(names(MEs2), 3) 
geneModuleMembership = as.data.frame(cor(datExpr0, MEs2, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
names(geneModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="") 
geneTraitSignificance = as.data.frame(cor(datExpr0, SST, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", names(SST), sep="") 
names(GSPvalue) = paste("p.GS.", names(SST), sep="") 

module = "blue"
column = match(module, modNames)
moduleGenes = moduleColors==module
pdf("results/test_WGCNA_clust_5.pdf")
par(mfrow = c(1,1))
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),abs(geneTraitSignificance[moduleGenes, 1]), xlab = paste("Module Membership in", module, "module"), ylab = "Gene significance for Ice Coverage", main = paste("Module membership vs. taxa significance\n"), cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module) 
dev.off()

print("heat map done. Create a summary table\n")

#probes2annot = match(names(datExpr0), annot$Bsp) # for SPtable
probes2annot = match(names(datExpr0), rownames(annot)) # work for bot OTUtable and SPtable

# if okay
sum(is.na(probes2annot)) 
# is equal to 0

print("merge Taxo and module")

TaxaInfo0 = data.frame(Taxon = names(datExpr0), TaxaSymbol = rownames(annot)[probes2annot], LinkID = annot$T2[probes2annot], moduleColor = moduleColors, geneTraitSignificance, GSPvalue)

modOrder = order(-abs(cor(MEs2, SST, use = "p")))
print("rank Taxo and module")

for (mod in 1:ncol(geneModuleMembership))
{
	oldNames = names(TaxaInfo0)
	TaxaInfo0 = data.frame(TaxaInfo0, geneModuleMembership[, modOrder[mod]], MMPvalue[,modOrder[mod]])
	names(TaxaInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], spe=""), paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

print("ranking done! saving the file")
TaxaOrder = order(TaxaInfo0$moduleColor, -abs(TaxaInfo0$GS.icecoverage))
TaxaInfo = TaxaInfo0[TaxaOrder,]
write.csv(TaxaInfo, file = "results/TaxaInfo.csv")
















