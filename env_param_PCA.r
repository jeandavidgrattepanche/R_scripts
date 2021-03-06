setwd('/Users/jaydiii/Desktop/JD_desktop/R_MiSeq/R_all')
library('ape')
library('phytools')
library('phyloseqCompanion')
library("ggplot2")
library("plyr")
library(factoextra)
env.all <- read.table('NBP1910_envdataCTD_2021d.txt')
envdata = sample_data(data.frame(cbind(data.frame(scale(env.all[c(5:22)])),env.all[c(1:4)])))
envdata2 = sample_data(data.frame(scale(env.all[c(5:22)])))
row.names(envdata2) <- env.all$sample
env.pca <- prcomp(envdata2, center = TRUE, scale = TRUE)
summary(env.pca)
varexp <- env.pca$sdev^2/sum(env.pca$sdev^2)
env.pca$x %>%
as.data.frame %>%
ggplot(aes(x=PC1,y=PC2)) + geom_point(size=4) + theme_bw(base_size=32) + labs(x=paste0("PC1: ", round(varexp[1]*100,1),"%"), y=paste0("PC2: ",round(varexp[2]*100,1), "%")) + theme(legend.position="top")
