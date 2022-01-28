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
co <- cooccur(datExpr1, spp_names = TRUE)
cox <- print(co)
cox[,"sp1_name"] == rownames(datExpr1)[cox$sp1]
cox[,"sp2_name"] == rownames(datExpr1)[cox$sp2]
nodes <- data.frame(id = 1:nrow(datExpr1), label = rownames(datExpr1), color = "#606482", shadow= FALSE, size=rowSums(datExpr1))
edges <- data.frame(from=cox$sp1, to = cox$sp2, color = ifelse(cox$p_lt <= 0.05, "red", "#3C3F51"), dashes = ifelse(cox$p_lt <= 0.05, TRUE, FALSE), weight = cox$p_lt)

G <- graph_from_data_frame(d=edges, vertices = nodes, directed= TRUE)
#deg <- degree(G, mode="all")
#V(G)$size <- deg*3

pdf("results_2/test_network_v2.pdf" , width=20, height=20)
plot(G, vertex.size=1, vertex.label.family="Helvetica",vertex.label.cex=0.2, edge.arrow.size=0.1)
dev.off()
