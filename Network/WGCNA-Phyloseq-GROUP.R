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
env.all <- read.table('NBP1910_envdata_v4.5.txt',sep="\t",header=TRUE,row.names=1)

mytable = otu_table(cbind(data.all[c(0)],data.all[c(10:105)]), taxa_are_rows=TRUE,errorIfNULL=TRUE)
envdata = sample_data(env.all)
testb <- as.matrix(data.all[c(0:10)])
TAX <- tax_table(testb)
physeq <- phyloseq(mytable, envdata, TAX)

#merge 3 size fractions for same sample:
station = as.character(get_variable(physeq, "station"))
layer = as.character(get_variable(physeq,"layer"))
sample_data(physeq)$StLay <- mapply(paste0,station,layer,collapse =" _")
physeq.group <- merge_samples(physeq, "StLay")

#need to be repeat for each size fraction
#datatoreduce <- physeq.group
#norare <- prune_taxa(taxa_sums(datatoreduce) > 0.05, datatoreduce)
#norare <- filter_taxa(physeq, function(x) sum(x>5) > (0.1*length(x)), TRUE)
#Hellinger transformation
physeq.group.Hell = transform_sample_counts(physeq.group, function(x) sqrt(x /sum(x)))   

#norare <- filter_taxa(physeq.group.Hell, function(x) sum(x>0.005) > (0.25*length(x)), TRUE)
MB <-subset_samples(physeq.group.Hell, group=="4")
MBc <- prune_taxa(taxa_sums(MB) > 0.0, MB)

#datExpr0 = as.data.frame(t(otu_table(norare)))
datExpr0 = as.data.frame(t(otu_table(MBc)))
datExpr0[] <- lapply(datExpr0, as.numeric)
annot = as.data.frame((cbind(data.all[c(0:10)])))
env.r = as.data.frame(sample_data(MBc))

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

#if TRUE continue

sampleTree = hclust(dist(datExpr0), method = "ward.D2")
pdf("results_2/sample_clust.pdf")
plot(sampleTree)
dev.off()

# if outliner remove it

#env.sub <- cbind(env.all[c(0)], env.all[c(7:35)], env.all[c(5)], env.all[c(7)])
#env.sub <- as.data.frame(cbind(env.r$size_num, env.r$group, env.r$latitude, env.r$longitude, env.r$time, env.r$bottom, env.r$depth, env.r$ZML_TS, env.r$Ze, env.r$ZCM, env.r$ice, env.r$air_temp, env.r$water_temp, env.r$conductivity, env.r$salinity, env.r$PAR, env.r$surf_PAR, env.r$oxygen, env.r$O_saturation, env.r$beam_trans, env.r$NH4, env.r$NO2NO3, env.r$PO4, env.r$fluorescence, env.r$Chla, env.r$Pprod_PAR, env.r$Babun, env.r$Bprod, env.r$Observed, env.r$Chao1, env.r$se.chao1, env.r$Shannon, env.r$picoDiversity, env.r$nanoDiversity, env.r$microDiversity))

#env.sub <- as.data.frame(cbind(env.r$sample, env.r$latitude, env.r$longitude, env.r$time, env.r$bottom, env.r$depth, env.r$ZML_TS, env.r$Ze, env.r$ZCM, env.r$ice, env.r$air_temp, env.r$water_temp, env.r$conductivity, env.r$salinity, env.r$PAR, env.r$surf_PAR, env.r$oxygen, env.r$O_saturation, env.r$beam_trans, env.r$NH4, env.r$NO2NO3, env.r$PO4, env.r$fluorescence, env.r$Chla, env.r$Pprod_PAR, env.r$Babun, env.r$Bprod, env.r$Observed, env.r$Chao1, env.r$se.chao1, env.r$Shannon, env.r$picoDiversity, env.r$nanoDiversity, env.r$microDiversity))

#env.sub <- as.data.frame(cbind(env.r$StLay,env.r$latitude, env.r$longitude, env.r$time, env.r$bottom, env.r$depth, env.r$ZML_TS, env.r$Ze, env.r$ZCM, env.r$ice, env.r$air_temp, env.r$water_temp, env.r$conductivity, env.r$salinity, env.r$PAR, env.r$surf_PAR, env.r$oxygen, env.r$O_saturation, env.r$beam_trans, env.r$NH4, env.r$NO2NO3, env.r$PO4, env.r$fluorescence, env.r$Chla, env.r$Pprod_PAR, env.r$Babun, env.r$Bprod, env.r$Observed, env.r$Chao1, env.r$se.chao1, env.r$Shannon, env.r$picoDiversity, env.r$nanoDiversity, env.r$microDiversity))

rownames(env.sub) <- rownames(env.r)
# readd size if needed
names(env.sub) <- c("Sample","latitude","longitude","time", "bottom", "depth", "ZML_TS", "Ze", "ZCM", "ice", "air_t", "water_t", "conductivity", "salinity", "PAR", "surf_PAR", "oxygen", "O_saturation", "beam_trans", "NH4", "NO2NO3", "PO4", "fluorescence", "Chla", "Pprod_PAR", "Babun", "Bprod", "Div", "Chao1", "Ch1se", "Shannon", "picoDiv", "nanoDiv", "microDiv")
#names(env.sub) <- c("StLayer","latitude","longitude","time", "bottom", "depth", "ZML_TS", "Ze", "ZCM", "ice", "air_t", "water_t", "conductivity", "salinity", "PAR", "surf_PAR", "oxygen", "O_saturation", "beam_trans", "NH4", "NO2NO3", "PO4", "fluorescence", "Chla", "Pprod_PAR", "Babun", "Bprod", "Div", "Chao1", "Ch1se", "Shannon", "picoDiv", "nanoDiv", "microDiv")
traitRows = match(rownames(datExpr0), rownames(env.sub)) # check if t table or table
#traitRows = match(names(datExpr0), rownames(env.sub))
env.red = env.sub[traitRows, -1]
env.red[] <- lapply(env.red, as.numeric)
env.reds <- scale(env.red)
names(env.reds) <- names(env.red)
traitColors = numbers2colors(env.reds, signed = TRUE)
pdf("results_2/sample_clust2.pdf", width =15, height=9)
plotDendroAndColors(sampleTree, traitColors, groupLabels = names(env.reds))
dev.off()

enableWGCNAThreads(20)

powers = c(c(1:10), seq(from = 5, to=50, by=2))
sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, networkType= "signed")
pdf("results/test_WGCNA_clust_1.pdf")
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)", ylab = "Scale Free topology Model Fit, signed R^2",type="n",main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2], labels=powers,cex=0.9,col="red")
abline(h=0.90,col="red")
dev.off()

pdf("results/test_WGCNA_clust_2.pdf")
plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)", ylab = "Main Connectivity",type="n",main = paste("Mean Connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,cex=0.9,col="red")
dev.off()

##### new method #### based on Wilson et al, 2018

softpower = 25
adjacency = adjacency(datExpr0, power = softpower, type = "signed")
TOM = TOMsimilarity(adjacency, TOMType= "signed")
dissTOM= 1-TOM
TaxaTree = hclust(as.dist(dissTOM), method= "average")

pdf("results/test_WGCNA_14clust_3_v2.pdf", width = 15, height= 9)
#?mergedColors = labels2colors(net2$colors)
plot(TaxaTree, xlab="", sub="", main="Tax clustering on TOM-based dissimilarity", labels = FALSE, hang = 0.04)
dev.off()

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = TaxaTree, distM = dissTOM, deepSplit = 2, pamRespectsDendro = FALSE, minClusterSize = minModuleSize)
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

pdf("results/test_WGCNA_14clust_3_v2b.pdf", width = 15, height= 9)
plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang =0.05, main = "Taxa dendogram and module colors")
dev.off()

MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = hclust(as.dist(MEDiss), method = "average")
pdf("results/test_WGCNA_14clust_4_v2c.pdf", width = 15, height= 9)
plot(METree, main=" Clustering of module eigengenes", xlab="", sub="")
dev.off()

MEDissThres = 0.0 #0.3
merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight=MEDissThres, verbose =3)
mergedColors = merge$colors
mergedMEs = merge$newMEs
pdf("results/test_WGCNA_14clust_3_v2bMergedc.pdf", width = 15, height= 9)
plotDendroAndColors(TaxaTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang =0.05, main = "Taxa dendogram and module colors")
dev.off()


colorOrder = c("grey", standardColors(50))
moduleLabels = match(mergedColors, colorOrder)-1
MEs = mergedMEs

save(MEs, moduleLabels, mergedColors, TaxaTree, file= "NBP1910_networkConstruction_stepbystepb.RData")

nSamples = nrow(datExpr0)
nTaxa = ncol(datExpr0)
MEs0 = moduleEigengenes(datExpr0, mergedColors)$eigengenes
MEs =orderMEs(MEs0)
moduleTraitCor = cor(MEs, env.red, use="p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
pdf("results/test_WGCNA_14clust_4c.pdf", width = 20, height= 8)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(env.red), yLabels = names(MEs), ySymbols = names(MEs), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins= FALSE, cex.text=0.5, zlim= c(-1,1), main = paste("Module-Environment relationships"))
dev.off()

MEs2 = MEs
moduleColors = mergedColors


###### Old #######
powertouse = 9 
print(c("soft threshold (power) used =", powertouse))
datExpr0[] <- lapply(datExpr0, as.numeric)
net2 = blockwiseModules(datExpr0, power=powertouse ,TOMType ="signed", deepSplit= 2 ,minModuleSize=10, maxBlockSize = 30000, reassignThreshold = 1e-12 , mergeCutHeight = 0.25, numericLabels = TRUE, pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "results/test_TOM", verbose= 3, nThreads = 10)


pdf("results/test_WGCNA_14clust_3.pdf", width = 15, height= 9)
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
pdf("results/test_WGCNA_14clust_4.pdf", width = 20, height= 8)
par(mar = c(6, 8.5, 3, 3))
labeledHeatmap(Matrix = moduleTraitCor, xLabels = names(env.red), yLabels = names(MEs2), ySymbols = names(MEs2), colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix, setStdMargins= FALSE, cex.text=0.5, zlim= c(-1,1), main = paste("Module-Environment relationships"))
dev.off()

###### same for both methods ######


SST = as.data.frame(env.red$water_t)
names(SST) = "water_t"
modNames = substring(names(MEs2), 3) 
TaxaModuleMembership = as.data.frame(cor(datExpr0, MEs2, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaModuleMembership), nSamples))
names(TaxaModuleMembership) = paste("MM", modNames, sep="")
names(MMPvalue) = paste("p.MM", modNames, sep="") 
TaxaTraitSignificance = as.data.frame(cor(datExpr0, SST, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaTraitSignificance), nSamples))
names(TaxaTraitSignificance) = paste("GS.", names(SST), sep="") 
names(GSPvalue) = paste("p.GS.", names(SST), sep="") 

module = "darkturquoise"
column = match(module, modNames)
#if NA => colors not in the module list
moduleTaxa = moduleColors==module
pdf("testWGCNA/test_WGCNA_clust_5c.pdf")
par(mfrow = c(1,1))
verboseScatterplot(abs(TaxaModuleMembership[moduleTaxa, column]),abs(TaxaTraitSignificance[moduleTaxa, 1]), xlab = paste("Module Membership in", module, "module"), ylab = "Taxa significance for water Temp", main = paste("Module membership vs. Taxa significance\n"), cex.main = 1.2 , cex.lab = 1.2, cex.axis = 1.2, col = module) 
dev.off()

print("heat map done. Create a summary table\n")

#probes2annot = match(names(datExpr0), annot$Bsp) # for SPtable
probes2annot = match(names(datExpr0), rownames(annot)) # work for bot OTUtable and SPtable

# if okay
sum(is.na(probes2annot)) 
# is equal to 0

print("merge Taxo and module")

TaxaInfo0 = data.frame(Taxon = names(datExpr0), TaxaSymbol = rownames(annot)[probes2annot], LinkID = annot$Btaxo_rank3[probes2annot], moduleColor = moduleColors, TaxaTraitSignificance, GSPvalue)

modOrder = order(-abs(cor(MEs2, SST, use = "p")))
print("rank Taxo and module")

for (mod in 1:ncol(TaxaModuleMembership))
{
	oldNames = names(TaxaInfo0)
	TaxaInfo0 = data.frame(TaxaInfo0, TaxaModuleMembership[, modOrder[mod]], MMPvalue[,modOrder[mod]])
	names(TaxaInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""), paste("p.MM.", modNames[modOrder[mod]], sep=""))
}

print("ranking done! saving the file")
TaxaOrder = order(TaxaInfo0$moduleColor, -abs(TaxaInfo0$GS.water_t))
TaxaInfo = TaxaInfo0[TaxaOrder,]
write.csv(TaxaInfo, file = "testWGCNA/TaxaInfoc.csv")



##network
library(cooccur)
library(visNetwork)
library(igraph)
datExpr1 <- t((datExpr0>0) *1L)
co <- cooccur(datExpr1, spp_names = TRUE)
cox <- print(co)
cox[,"sp1_name"] == rownames(datExpr1)[cox$sp1]
cox[,"sp2_name"] == rownames(datExpr1)[cox$sp2]
nodes <- data.frame(id = 1:nrow(datExpr1), label = rownames(datExpr1), color = "#606482", shadow= FALSE, size=rowSums(datExpr1))
edges <- data.frame(from=cox$sp1, to = cox$sp2, color = ifelse(cox$p_lt <= 0.05, "red", "#3C3F51"), dashes = ifelse(cox$p_lt <= 0.05, TRUE, FALSE), weight = cox$p_lt)

#p <- visNetwork(nodes = nodes, edges= edges) %>% 
#visIgraphLayout(layout = "layout_with_kk")


G <- graph_from_data_frame(d=edges, vertices = nodes, directed= TRUE)
#deg <- degree(G, mode="all")
#V(G)$size <- deg*3

pdf("testWGCNA/test_WGCNA_network_v2.pdf" , width=20, height=20)
plot(G, vertex.size=1, vertex.label.family="Helvetica",vertex.label.cex=0.2, edge.arrow.size=0.1)
dev.off()

Isolated = which(degree(G)==0)
G2 = delete.vertices(G,Isolated)
LO2 = layout_with_fr(G)[-Isolated,]
pdf("testWGCNA/test_WGCNA_network2.pdf", width=20, height=20)
plot(G2, vertex.label.family="Helvetica",vertex.label.cex=0.2, edge.arrow.size=0.1, vertex.size=1, layout=LO2)
dev.off()

toremove = which(E(G2)$color == "poor")
G3 = delete.edges(G2,toremove)
LO3 = layout_with_fr(G2)[-toremove,]
pdf("testWGCNA/test_WGCNA_network3.pdf")
plot(G3, vertex.label.family="Helvetica",vertex.label.cex=0.2, edge.arrow.size=0)
dev.off()




## correlation network:
library(Hmisc)
library(Matrix)
library(igraph)
#datExpr0 = as.matrix(cbind(data.all[c(0)],data.all[c(16:111)]))
#datExpr0[] <- lapply(datExpr0, as.numeric)
datExpr3 = as.matrix((otu_table(physeq.group)))
datExpr3[] <- lapply(datExpr3, as.numeric)
datExpr0 = t(datExpr3)
#shortlist = datExpr0[0:30,0:10]
shortlist=datExpr0[ rowSums(datExpr0) >= 100, ]
renamesd <- c(gsub("[.]", "-", names(shortlist)))
names(shortlist) <- renamesd

corr <- rcorr(t(shortlist), type="spearman")
corr.pval <- forceSymmetric(corr$P)
#tax <- tax_table(norare)
testbc <- data.all[c(0:17)]
sel.taxa <- testbc[rownames(corr.pval),, drop=FALSE]
all.equal(rownames(sel.taxa), rownames(corr.pval))
## should be TRUE to continue

# pvalue cut off 0.001
p.ok <- corr.pval < 0.001
r.val = corr$r
p.ok.rval <- r.val*p.ok

#r cut off 0.75
p.ok.r.sel <- abs(p.ok.rval)>0.85
p.ok.rval.sel <- p.ok.r.sel*r.val
###try with only + and oly -
HERE
pp.ok.r.sel <- (p.ok.rval)>0.85
np.ok.r.sel <- (p.ok.rval)<(-0.85)
pp.ok.rval.sel <- pp.ok.r.sel*r.val
np.ok.rval.sel <- np.ok.r.sel*r.val

#creat adjacency matrix
adjm <- as.matrix(p.ok.rval.sel)
colnames(adjm) <- as.vector(sel.taxa$Bsp)
rownames(adjm) <- as.vector(sel.taxa$Bsp)

#net grap with igraph
net_grph=graph.adjacency(adjm, mode="undirected",weighted=TRUE,diag=FALSE)
edgew <- E(net_grph)$weight
edgew2 <- (1+edgew)/2
bad.vs<-V(net_grph)[degree(net_grph) == 0]
net_grph <- delete.vertices(net_grph,bad.vs)
pdf("testWGCNA/correlationnetwork_SL100r_0.75c.pdf", width=20, height=20)
#plot(net_grph, vertex.size=2, vertex.frame.color="black", edge.curved=F, edge.width=1.5, layout=layout.fruchterman.reingold, edge.color=ifelse(edgew <0, "red", "blue"), vertex.label.cex=0.5)

#dev.off()

plot(net_grph, vertex.size=2, vertex.frame.color="black", edge.curved=F, edge.width=0.5, layout=layout.reingold.tilford, edge.color=ifelse(edgew <(-0.9), "red", ifelse ( (-0.9) <  edgew & edgew<0, "darkorange", ifelse( 0<edgew & edgew<0.9, "forestgreen", ifelse ( edgew> 0.9, "blue","grey")))), vertex.label.cex=0.5, vertex.label.angle=0.5)
dev.off()

c_scale <- colorRamp(c('red', 'forestgreen', 'blue'))

E(net_grph)$color = apply(c_scale(edgew2), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))  
pdf("testWGCNA/correlationnetwork_SL100r_0.75c.pdf", width=20, height=20)
plot(net_grph, vertex.size=2, vertex.frame.color="black", edge.curved=F, edge.width=abs(edgew), layout=layout.circle, vertex.label.cex=1)
dev.off()





