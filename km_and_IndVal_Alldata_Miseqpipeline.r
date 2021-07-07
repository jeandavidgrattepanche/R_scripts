setwd('/home/tuk61790/HPCr/')
library('phyloseqCompanion')
library("ggplot2")
library("plyr")
library(stringr)
library(stats)
library(indicspecies)
source('/home/tuk61790/HPCr/multiplot.R', chdir = TRUE)
source('/home/tuk61790/HPCr/get_top_taxa.R', chdir = TRUE)

#create phyloseq/physeq objects
data.all <- read.table('OTUtable_ingroup_noTtree.txt')
mytable = otu_table(cbind(data.all[c(0)],data.all[c(16:111)]), taxa_are_rows=TRUE,errorIfNULL=TRUE)
env.all <- read.table('NBP1910_envdata_v3b_size_noNaN.txt')
envdata = sample_data(env.all)
testb <- as.matrix(data.all[c(0:15)])
TAX <- tax_table(testb)
physeq <- phyloseq(mytable, envdata, TAX)


#convert abundance in percentage:
physeq.perc = transform_sample_counts(physeq, function(x) 100 * x/sum(x))

#reduce the number of OTU based on abundance (1e-05 here)
physeq.red <- filter_taxa(physeq.perc, function(x) var(x) > 1e-05, TRUE)

datatorun <- physeq.red
# S or size, G or Group

km = kmeans(as(t(otu_table(physeq.red)),"matrix"), centers= 5)

pt = multipatt(as(t(otu_table(physeq.red)),"matrix"), km$cluster, control = how(nperm=999))

cat("indval_km",capture.output(summary(pt, indvalcomp=TRUE)),file="results_km_indval.txt",sep="\t")

// #create network plot for each size fraction. size can be replace by another varaible
// sizes <- c("micro","nano","pico")
// Groups <- c("other","northern","southern","offshore")
// ind <- "jaccard" 
// list <- list()
// i = 0
// for(Gp in Groups){
// 	for(sizedata in sizes){
// 		i = i + 1
// 		named <- paste(Gp,sizedata,sep = "-")
// 		print(named)
// 		network <- make_network(subset_samples(subset_samples(datatorun, size == sizedata), Group == Gp), type="taxa", distance = ind, max.dist = 0.3, keep.isolates=FALSE)
// 		ordplotPm <- plot_network(network, subset_samples(subset_samples(datatorun, size == sizedata), Group == Gp), type ="taxa", color="Btaxo_rank4", label=NULL ,title=named)
// 		assign(paste("p",i,sep=""), ordplotPm)
// 		list[[length(list)+1]] <- paste("p",i,sep="")
// 		plot(get(paste("p",i,sep="")))
// 	}
// }
// pdf('Network_size_group.pdf')
// multiplot(p1 + theme(legend.position="none"),p2+ theme(legend.position="none"),p3+ theme(legend.position="none"),p4+ theme(legend.position="none"),p5+ theme(legend.position="none"),p6+ theme(legend.position="none"),p7+ theme(legend.position="none"),p8+ theme(legend.position="none"),p9+ theme(legend.position="none"),p10+ theme(legend.position="none"),p11+ theme(legend.position="none"),p12+ theme(legend.position="none"), p1,cols=6)
// #multiplot(p1 ,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,cols=4)
// dev.off()
