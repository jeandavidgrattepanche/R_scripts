setwd('/home/tuk61790/HPCr/')
library('ape')
library('phytools')
library('phyloseqCompanion')
library("ggplot2")
library("plyr")
library(stringr)
#library(factoextra)
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
physeq.bc = physeq # backup in case 

#read distribution OTUs and Samples
readnumb = data.frame(readn = sort(taxa_sums(physeq), TRUE), sorted = 1:ntaxa(physeq), type = "OTUS")
readnumb = rbind(readnumb, data.frame(readn = sort(sample_sums(physeq), TRUE), sorted = 1:nsamples(physeq), type = "Samples"))
titlenb = "read number distribution"
p = ggplot(readnumb, aes(x= sorted, y = readn)) + geom_bar(stat="identity")
plotp = p + ggtitle(titlenb) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
jpeg(filename='barplot_read.jpeg', width=60, height=15, unit="cm", res=100)
print(plotp)
dev.off()

#check if empty sample - can be useful after next step
any(sample_sums(physeq) == 0)
physeq = prune_samples(sample_sums(physeq) >0, physeq)

#remove singleton ! be careful if samples are already rarefy
physeq = prune_taxa(taxa_sums(physeq) > 1, physeq)

#read distribution OTUs and Samples no singelton
readnumb = data.frame(readn = sort(taxa_sums(physeq), TRUE), sorted = 1:ntaxa(physeq), type = "OTUS")
readnumb = rbind(readnumb, data.frame(readn = sort(sample_sums(physeq), TRUE), sorted = 1:nsamples(physeq), type = "Samples"))
titlenbns = "read number distribution no singleton"
p = ggplot(readnumb, aes(x= sorted, y = readn)) + geom_bar(stat="identity")
plotp = p + ggtitle(titlenbns) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")
jpeg(filename='barplot_read_nosing.jpeg', width=60, height=15, unit="cm", res=100)
print(plotp)
dev.off()