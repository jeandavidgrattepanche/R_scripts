#if issue with phyloseq and ape
#install.packages("devtools")
#devtools::install_github("joey711/phyloseq")
#library("phyloseq")


setwd("~/Desktop/JD_OG_paper")
library('ape')
library("phyloseq")
library("foreach")
library("iterators")
library("parallel")
library("doParallel")
library("ggplot2")
library("plyr")
library("vegan")
data.all <- read.table('Table_Rubisco.txt')
mytable = otu_table(data.all, taxa_are_rows=TRUE,errorIfNULL=TRUE)
env.all <- read.table('env_Memi_august2018_pooled.txt')
envdata = sample_data(env.all)
tree <- read.tree('Rubisco_memi.tree')
phy_tree(tree, errorIfNULL=TRUE)
physeq <- phyloseq(mytable, envdata, tree)


treedraw <- plot_tree(physeq, color="all", size="Abundance",label.tips="taxa_names", ladderize = TRUE, justify="left", sizebase=5, base.spacing=0.03, title="RuBisCO", text.size=2)
quartz()
plot(treedraw)
