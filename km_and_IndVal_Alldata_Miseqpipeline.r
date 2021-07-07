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

#reduce the number of OTU based on variance (1e-05 here)
physeq.red <- filter_taxa(physeq.perc, function(x) var(x) > 1e-05, TRUE)

scaled_data <- as(t(otu_table(physeq.red)),"matrix")
# S or size, G or Group

km = kmeans(scaled_data, centers= 5)

pt = multipatt(scaled_data, km$cluster, control = how(nperm=999))

cat("indval_km",capture.output(summary(pt, indvalcomp=TRUE)),file="results_km_indval.txt",sep="\t")

library(factoextra)
library(cluster)
library(tidyverse)

tot_withinss <- map_dbl(1:20, function(k) {
	model <- kmeans(as(t(otu_table(physeq.red)),"matrix"), centers = k)
	model$tot.withinss
})

bothkw <- data.frame( k = 1:20, tot_withinss = tot_withinss)

pdf("kmeans_Elbow.pdf")
ggplot(bothkw, aes(x = k, y = tot_withinss)) + 
	geom_line() + geom_point() + scale_x_continuous(breaks= 1:20)
dev.off()


set.seed(123)
gap_stat <- clusGap(as(t(otu_table(physeq.red)),"matrix"), FUN = kmeans, nstart = 25, K.max = 20, B = 50)

pdf("kmeans_gap.pdf")
fviz_gap_stat(gap_stat)
dev.off()

library(factoextra)
library(NbClust)
# Elbow method
plEl <- fviz_nbclust(scaled_data, kmeans, method = "wss") +
  geom_vline(xintercept = 2, linetype = 2)+
  labs(subtitle = "Elbow method")
# Silhouette method
plSil <- fviz_nbclust(scaled_data, kmeans, method = "silhouette")+
  labs(subtitle = "Silhouette method")
# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
plGap <- fviz_nbclust(scaled_data, kmeans, nstart = 25,  method = "gap_stat", nboot = 50)+
  labs(subtitle = "Gap statistic method")

pdf('kmeans_selection.pdf')
multiplot(plEl,plSil,plGap,cols=1)
dev.off()
