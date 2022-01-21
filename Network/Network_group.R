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
MB <-subset_samples(physeq.group.Hell, group=="4")
MBc <- prune_taxa(taxa_sums(MB) > 0.0, MB)

network <- make_network(MBc, type="samples")

network2 <- make_network(MBc, type="taxa",distance="jaccard",max.dist= 0.15)
netplot <- plot_network(network2, MBc, type="taxa",color="Trophic")
pdf("results_2/nettaxa_MB.pdf",width=13, height=8)
plot(netplot)
dev.off()
#all 0.45, MB 0.15, StLay 0.32
network2 <- make_network(physeq.group.Hell, type="taxa",distance="jaccard",max.dist= 0.32)
netplot <- plot_network(network2, physeq.group.Hell, type="taxa",color="Trophic")
pdf("results_2/nettaxa_StLay.pdf",width=13, height=8)
plot(netplot)
dev.off()
network2 <- make_network(physeq.group.Hell, type="samples", max.dist= 0.6)
netplot <- plot_network(network2, physeq.group.Hell, color="station")
pdf("results_2/netsamples_StLay.pdf",width=13, height=8)
plot(netplot)
dev.off()
