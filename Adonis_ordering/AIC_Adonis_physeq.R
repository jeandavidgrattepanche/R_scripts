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
source('/home/tuk61790/Network/AICc_table_generation_edited.R',chdir =TRUE)
source('AICc_PERMANOVA.R',chdir=T)

data.all <- read.table('OTUtable_ingroup_100.txt',sep="\t",header=TRUE,row.names=1)
#env.all <- read.table('NBP1910_envdata_v4.5.txt',sep="\t",header=TRUE,row.names=1)
env.all <- read.table('NBP1910_envdata_v4.5_edited.txt',sep="\t",header=TRUE,row.names=1)                                                              
env.all.sc <- data.frame(cbind(env.all[c(0:7)],scale(env.all[c(8:41)]))) 
envdata = sample_data(env.all.sc)
mytable = otu_table(cbind(data.all[c(0)],data.all[c(10:105)]), taxa_are_rows=TRUE,errorIfNULL=TRUE)
envdata = sample_data(env.all)
testb <- as.matrix(data.all[c(0:10)])
TAX <- tax_table(testb)
physeq <- phyloseq(mytable, envdata, TAX)
physeq = transform_sample_counts(physeq, function(x) 100 * x/sum(x)) 
# rerun without B, D and bucket => missing nuts  
#N because not used
physeqb = subset_samples(physeq, station != "B")
physeqb = subset_samples(physeqb, station != "D")
physeqb = subset_samples(physeqb, layer != "Bucket")

torun <- subset_samples(physeqb, size == "pico")
env.var.pico <- c( "latitude", "salinity","NO2NO3","Chla","conductivity", "ZML_TS","beam_trans", "fluorescence", "Ze","oxygen" ,"water_temp",  "Bprod", "ice","NH4","O_saturation","depth","PO4", "bottom")
AICresults.pico <- AICc.table.all(env.var.pico, matrix.char=distance(torun,"jaccard"),perm=999,method="jaccard", df=as(sample_data(torun),"data.frame"),comb.incl=c(1,2,3,4,5,6,7))
write.table(AICresults.pico, file="AIC_Table_pico_scaled.txt", append=F, sep="\t")

torun <- subset_samples(physeqb, size == "nano")
env.var.pico <- c( "latitude","NO2NO3", "ice", "ZML_TS", "Chla", "oxygen", "NH4", "Ze", "fluorescence", "beam_trans", "PO4", "salinity", "bottom")
AICresults.pico <- AICc.table.all(env.var.pico, matrix.char=distance(torun,"jaccard"),perm=999,method="jaccard", df=as(sample_data(torun),"data.frame"),comb.incl=c(1,2,3,4,5,6,7))
write.table(AICresults.pico, file="AIC_Table_nano_scaled.txt", append=F, sep="\t")

torun <- subset_samples(physeqb, size == "micro")
env.var.pico <- c( "latitude","Ze", "NO2NO3","oxygen","ZML_TS","Chla", "fluorescence","beam_trans","Bprod","salinity","ice","conductivity", "O_saturation", "water_temp", "PO4", "NH4") 
AICresults.pico <- AICc.table.all(env.var.pico, matrix.char=distance(torun,"jaccard"),perm=999,method="jaccard", df=as(sample_data(torun),"data.frame"),comb.incl=c(1,2,3,4,5,6,7))
write.table(AICresults.pico, file="AIC_Table_micro_scaled.txt", append=F, sep="\t")

torun <- subset_samples(physeqb, size == "pico")
Adonis2.results = adonis2(formula = distance(torun, "jaccard") ~ latitude + NO2NO3  , data = as(sample_data(torun), "data.frame"), permutations = 9999, by = "terms")
cat("Adonis_results_pico", capture.output(Adonis2.results), file="Adonis2_results022522.txt", sep="\n", append=TRUE)