Setwd('/home/tuk61790/Network/')
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
sample_data(physeq)$StLay <- mapply(paste0,station,layer,collapse ="-")
physeq.group <- merge_samples(physeq, "StLay")

#Hellinger transformation
physeq.group.Hell = transform_sample_counts(physeq.group, function(x) sqrt(x /sum(x)))   

#norare <- filter_taxa(physeq.group.Hell, function(x) sum(x>0.005) > (0.25*length(x)), TRUE)
#all 0.45, MB 0.15, StLay 0.32

#network after Hellinger transformation (better than just % or relative abundance) for samples and taxa using all samples
network1 <- make_network(physeq.group.Hell, type="taxa",distance="jaccard",max.dist= 0.32)
netplot1 <- plot_network(network1, physeq.group.Hell, type="taxa",color="Trophic")
pdf("results_2/nettaxa_StLay.pdf",width=13, height=8)
plot(netplot1)
dev.off()

network2 <- make_network(physeq.group.Hell, type="samples", max.dist= 0.6)
netplot2 <- plot_network(network2, physeq.group.Hell, color="station")
pdf("results_2/netsamples_StLay.pdf",width=13, height=8)
plot(netplot2)
dev.off()

MB <-subset_samples(physeq.group.Hell, group=="4")
MBc <- prune_taxa(taxa_sums(MB) > 0.05, MB)
MB2 <-subset_samples(physeq.group, group=="4")
MBred <- prune_taxa(taxa_sums(MB2) > 100.0, MB2)
MBred.Hell = transform_sample_counts(MBred, function(x) sqrt(x /sum(x))) # not sure transformation appropriate here

GS2 <-subset_samples(physeq.group, group=="1")
GSred <- prune_taxa(taxa_sums(GS2) > 100.0, GS2)
GSred.Hell = transform_sample_counts(GSred, function(x) sqrt(x /sum(x))) # not sure transformation appropriate here



network3h <- make_network(MBc, type="taxa",distance="jaccard",max.dist= 0.15)
netplot3h <- plot_network(network3h, MBc, type="taxa",color="Trophic")
pdf("results_2/nettaxa_MB.Hell.pdf",width=13, height=8)
plot(netplot3h)
dev.off()


network3 <- make_network(MBred, type="taxa",distance="jaccard",max.dist= 0.3)
netplot3 <- plot_network(network3, MBred, type="taxa",color="Trophic")
pdf("results_2/nettaxa_MB2.pdf",width=13, height=8)
plot(netplot3)
dev.off()
network4 <- make_network(GSred, type="taxa",distance="jaccard",max.dist= 0.4)
netplot4 <- plot_network(network4, MBc, type="taxa",color="Trophic")
pdf("results_2/nettaxa_GS2.pdf",width=13, height=8)
plot(netplot4)
dev.off()


### other network
library(cooccur)
library(visNetwork)
library(igraph)
datExpr0 = as.data.frame(t(otu_table(MBc)))
datExpr0[] <- lapply(datExpr0, as.numeric)
annot = as.data.frame((cbind(data.all[c(0:10)])))
env.r = as.data.frame(sample_data(MBc))

gsg = goodSamplesGenes(datExpr0, verbose = 3)
gsg$allOK

datExpr1 <- ((datExpr0>0) *1L)
co <- cooccur(datExpr1, spp_names = TRUE) #can take awhile if correctly done on OTUs even if only a few OTUs
cox <- print(co)
cox[,"sp1_name"] == rownames(datExpr1)[cox$sp1]
cox[,"sp2_name"] == rownames(datExpr1)[cox$sp2]
nodes <- data.frame(id = 1:nrow(datExpr1), label = rownames(datExpr1), color = "#606482", shadow= FALSE, size=rowSums(datExpr1))
edges <- data.frame(from=cox$sp1, to = cox$sp2, color = ifelse(cox$p_lt <= 0.5, "red", "#3C3F51"), dashes = ifelse(cox$p_lt <= 0.05, TRUE, FALSE), weight = cox$p_lt)

G <- graph_from_data_frame(d=edges, vertices = nodes, directed= TRUE)
#deg <- degree(G, mode="all")
#V(G)$size <- deg*3

pdf("results_2/test_network_v2.pdf" , width=20, height=20)
plot(G, vertex.size=1, vertex.label.family="Helvetica",vertex.label.cex=0.2, edge.arrow.size=0.1)
dev.off()

Isolated = which(degree(G)==0)
G2 = delete.vertices(G,Isolated)
LO2 = layout_with_fr(G)[-Isolated,]
pdf("results_2/test_network_v2-sp.2.pdf", width=10, height=10)
plot(G2, vertex.label.family="Helvetica",vertex.label.cex=0.5, edge.arrow.size=0.1, vertex.size=2, layout=LO2)
dev.off()

toremove = which(E(G2)$color == "poor")
G3 = delete.edges(G2,toremove)
LO3 = layout_with_fr(G2)[-toremove,]
pdf("results_2/test_network_v2.3.pdf")
plot(G3, vertex.label.family="Helvetica",vertex.label.cex=0.2, edge.arrow.size=0)
dev.off()


## correlation network:
library(Hmisc)
library(Matrix)
library(igraph)
#datExpr0 = as.matrix(cbind(data.all[c(0)],data.all[c(16:111)]))
#datExpr0[] <- lapply(datExpr0, as.numeric)
MB <-subset_samples(physeq.group.Hell, group=="4")
MBc <- prune_taxa(taxa_sums(MB) > 0.1, MB)
#0.05 -> 882
#0.1 -> 481
#0.25 -> 105
datExpr3 = as.matrix((otu_table(MBc)))

corr <- rcorr((datExpr3), type="spearman")
corr.pval <- forceSymmetric(corr$P)
#tax <- tax_table(norare)

# pvalue cut off 0.001
p.ok <- corr.pval < 0.001 #identify OTUs couple with an ok pvalue
r.val = corr$r
p.ok.rval <- r.val*p.ok #extract rvalue for identify OTUs couple with an ok pvalue

#r cut off 0.75
p.ok.r.sel <- abs(p.ok.rval)>0.85 #identify OTUs couple with an ok rvalue
p.ok.rval.sel <- p.ok.r.sel*r.val #extract rvalue for identify OTUs couple with an ok rvalue

#creat adjacency matrix
adjm <- as.matrix(p.ok.rval.sel)

net_grph=graph.adjacency(adjm, mode="undirected",weighted=TRUE,diag=FALSE)
edgedt <- as_data_frame(net_grph, what="edges")
testbc <- data.all[c(0,2:9)]
sel.taxa2 <- testbc[rownames(p.ok.rval.sel),, drop=FALSE]
sel.taxa2['name'] = rownames(sel.taxa2)
sel.taxa2 <- sel.taxa2[,c(9,1,2,3,4,6,5,7,8)] #really important to have the rowname and the first column (also the name matching the edge dataframe) with the same OTU name to build the network
net2 <- graph_from_data_frame(edgedt, sel.taxa2, directed=F)
#c("Klepto?", "Mixo", "mixo?", "Parasite", "Phagotroph", "Phototroph")
#c("Gray50","red",  "orange","purple","blue","forestgreen")
V(net2)$color <- ifelse(V(net2)$Trophic == "Klepto?","Gray50", ifelse(V(net2)$Trophic == "Mixo","red",ifelse(V(net2)$Trophic == "mixo?","orange",ifelse(V(net2)$Trophic == "Parasite","purple", ifelse(V(net2)$Trophic == "Phagotroph","blue", ifelse(V(net2)$Trophic == "Phototroph","forestgreen","black"))))))
bad.vs<-V(net2)[degree(net2) == 0]
net2 <- delete.vertices(net2,bad.vs)
edgew <- E(net2)$weight
pdf("results_2/correlationnetwork_MBc_v2.pdf", width=20, height=20)
plot(net2, vertex.size=log(V(net2)$readnumber)/2, edge.curved=F, edge.width=5, layout=layout.graphopt, edge.color=ifelse(edgew <(-0.9), "red", ifelse ( (-0.9) <=  edgew & edgew<=0, "darkorange", ifelse( 0<edgew & edgew<=0.9, "forestgreen", ifelse ( edgew> 0.9, "blue","grey")))), vertex.label.cex=0.5, vertex.label.angle=0.5)
dev.off()
#colnames(adjm) <- as.vector(sel.taxa$Genus.species)
#rownames(adjm) <- as.vector(sel.taxa$Genus.species)


#net grap with igraph
net_grph=graph.adjacency(adjm, vertices=sel.taxa, mode="undirected",weighted=TRUE,diag=FALSE)
edgew <- E(net_grph)$weight
V(net_grph)$trophic <- as.vector(sel.taxa$Trophic)
edgew2 <- (1+edgew)/2
bad.vs<-V(net_grph)[degree(net_grph) == 0]
net_grph <- delete.vertices(net_grph,bad.vs)
pdf("results_2/correlationnetwork_MBc_0.85c.pdf", width=15, height=15)
#layout=layout.reingold.tilford
#layout=layout.fruchterman.reingold
#vertex.size=2,
plot(net_grph, vertex.size=2, vertex.frame.color="black", edge.curved=F, edge.width=0.5, layout=layout.graphopt, edge.color=ifelse(edgew <(-0.95), "red", ifelse ( (-0.95) <  edgew & edgew<0, "darkorange", ifelse( 0<edgew & edgew<0.95, "forestgreen", ifelse ( edgew> 0.95, "blue","grey")))), vertex.label.cex=0.5, vertex.label.angle=0.5)
dev.off()

c_scale <- colorRamp(c('red', 'forestgreen', 'blue'))

E(net_grph)$color = apply(c_scale(edgew2), 1, function(x) rgb(x[1]/255,x[2]/255,x[3]/255))
pdf("testWGCNA/correlationnetwork_SL100r_0.75c.pdf", width=20, height=20)
plot(net_grph, vertex.size=2, vertex.frame.color="black", edge.curved=F, edge.width=abs(edgew), layout=layout.circle, vertex.label.cex=1)
dev.off()
