setwd('/home/tuk61790/Network/')
library('ape')
library('phytools')
library('phyloseqCompanion')
library("ggplot2")
library("plyr")
library(gtools)
library(stringr)
#library(factoextra)

options(stringsAsFactors = FALSE)

data.all <- read.table('OTUtable_ingroup_100.txt',sep="\t",header=TRUE,row.names=1)
env.all <- read.table('NBP1910_envdata_v4.5.txt',sep="\t",header=TRUE,row.names=1)
mytable = otu_table(cbind(data.all[c(0)],data.all[c(10:105)]), taxa_are_rows=TRUE,errorIfNULL=TRUE)
envdata = sample_data(env.all)
testb <- as.matrix(data.all[c(0:10)])
TAX <- tax_table(testb)
physeq <- phyloseq(mytable, envdata, TAX)
pico <- subset_samples(physeq, size == "pico")
nano <- subset_samples(physeq, size == "nano")
micro <- subset_samples(physeq, size == "micro")

Test AICc script from https://github.com/kdyson/R_Scripts/blob/master/AICc_PERMANOVA.R

source('/home/tuk61790/Network/AICc_table_generation_edited.R',chdir =TRUE)
source('AICc_PERMANOVA.R',chdir=T)
#testvar <-c("group","size","layer","latitude")
#testvar2 <- permutations(4,4,testvar) => not working properly yet
#AICc.table.Nvar(testvar2, matrix.char=distance(physeq,"jaccard"),perm=99,n.var=4,method="jaccard", df=as(sample_data(physeq),"data.frame"))

env.var <- c("group","size","longitude","latitude","time","bottom","depth","ZML_TS","Ze","ZCM","ice","air_temp","water_temp","conductivity","salinity","PAR","surf_PAR","oxygen","O_saturation","beam_trans","fluorescence")
#env.var2 <- permutations(length(env.var),length(env.var),env.var) #no role of order here so skipped
#AICresults <- AICc.table.Nvar(env.var, matrix.char=distance(physeq,"jaccard"),perm=9999,n.var=21,method="jaccard", df=as(sample_data(physeq),"data.frame"))
#write.table(AICresults, file="AIC_Table.txt", append=TRUE, sep="\t")
AICresults <- AICc.table.all(env.var, matrix.char=distance(physeq,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(physeq),"data.frame"),comb.incl=c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21))
write.table(AICresults, file="AIC_Table.txt", append=FALSE, sep="\t")
#select order based on AIC value (1=lowest and last= highest)
#use version below to add 1 parameter at the time => need to edit script to add most significant paramters each time.=> need to close R and restart
AICresults2 <- AICc.table.all(env.var, matrix.char=distance(physeq,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(physeq),"data.frame"),comb.incl=c(1), control.var.char = c("size + latitude"))
write.table(AICresults2, file="AIC_Table2-3.txt", append=FALSE, sep="\t")
#use this order to run adonis
Adonis2.results = adonis2(formula = distance(physeq, "jaccard") ~ group + layer + latitude + size  , data = as(sample_data(physeq), "data.frame"), permutations = 9999, by = "Term") 
cat("Adonis_results", capture.output(Adonis2.results), file="summary_of_Adonis2bmargin2_results021722.txt", sep="\n", append=TRUE)

#should be repeat for each size fraction
