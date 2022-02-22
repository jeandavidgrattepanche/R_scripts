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
#scale env data
#mytable = otu_table(cbind(data.all[c(0)],data.frame(scale(data.all[c(10:105)]))), taxa_are_rows=TRUE,errorIfNULL=TRUE)
# env.all.sc <- data.frame(scale(env.all))
# env.all.sc["sample"] <- row.names(env.all)
# envdata = sample_data(env.all.sc)
envdata = sample_data(env.all)
testb <- as.matrix(data.all[c(0:10)])
TAX <- tax_table(testb)
physeq <- phyloseq(mytable, envdata, TAX)
physeq = transform_sample_counts(physeq, function(x) 100 * x/sum(x)) # better than scaling (no negative abundance)
pico <- subset_samples(physeq, size == "pico")
nano <- subset_samples(physeq, size == "nano")
micro <- subset_samples(physeq, size == "micro")

#Test AICc script from https://github.com/kdyson/R_Scripts/blob/master/AICc_PERMANOVA.R

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
Adonis2.results = adonis2(formula = distance(physeq, "jaccard") ~ group + layer + latitude + size  , data = as(sample_data(physeq), "data.frame"), permutations = 9999, by = "terms")
cat("Adonis_results", capture.output(Adonis2.results), file="summary_of_Adonis2bmargin2_results021722.txt", sep="\n", append=TRUE)

#should be repeat for each size fraction
data.all <- read.table('OTUtable_ingroup_100.txt',sep="\t",header=TRUE,row.names=1)
env.all <- read.table('NBP1910_envdata_v4.5_edited.txt',sep="\t",header=TRUE,row.names=1)
# use average to replace missing data
mytable = otu_table(cbind(data.all[c(0)],data.all[c(10:105)]), taxa_are_rows=TRUE,errorIfNULL=TRUE)
envdata = sample_data(env.all)
testb <- as.matrix(data.all[c(0:10)])
TAX <- tax_table(testb)
physeq <- phyloseq(mytable, envdata, TAX)
pico <- subset_samples(physeq, size == "pico")
nano <- subset_samples(physeq, size == "nano")
micro <- subset_samples(physeq, size == "micro")
env.var.micro <- c("time" , "bottom" , "depth" , "ZML_TS" , "Ze" , "ZCM" , "ice" , "air_temp" , "water_temp" , "conductivity" , "salinity" , "PAR" , "surf_PAR" , "oxygen" , "O_saturation" , "beam_trans" , "NH4" , "NO2NO3" , "PO4" , "fluorescence" , "Chla" , "Pprod_Sun" , "Pprod_PAR" , "Babun" , "Bprod" , "picoDiversity" , "nanoDiversity")
AICresults.micro <- AICc.table.all(env.var.micro, matrix.char=distance(micro,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(micro),"data.frame"),comb.incl=c(1,2))
write.table(AICresults.micro, file="AIC_Table_micro.txt", append=FALSE, sep="\t")
AICresults.micro <- AICc.table.all(env.var.micro, matrix.char=distance(micro,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(micro),"data.frame"),comb.incl=c(1), control.var.char = c("NO2NO3 + ZML_TS + time + Babun + salinity + Ze + O_saturation"))
write.table(AICresults.micro, file="AIC_Table_micro2.txt", append=T, sep="\t")


# micro: Three parameters played a role in shaping the microplankton community: the nitrate/nitrite, then the depth of the mized layer and the time of sampling were important regarding the compostion of the microplankton. Additional parameters are not adding explanation to the variance, suggesting negligiable role.

env.var.nano <- c("time" , "bottom" , "depth" , "ZML_TS" , "Ze" , "ZCM" , "ice" , "air_temp" , "water_temp" , "conductivity" , "salinity" , "PAR" , "surf_PAR" , "oxygen" , "O_saturation" , "beam_trans" , "NH4" , "NO2NO3" , "PO4" , "fluorescence" , "Chla" , "Pprod_Sun" , "Pprod_PAR" , "Babun" , "Bprod" , "picoDiversity" ,  "microDiversity")
AICresults.nano <- AICc.table.all(env.var.nano, matrix.char=distance(nano,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(nano),"data.frame"),comb.incl=c(1))
write.table(AICresults.nano, file="AIC_Table_nano.txt", append=FALSE, sep="\t")

env.var.nano <- c("bottom" , "Ze" , "ZCM" , "ice" , "air_temp" , "water_temp" , "conductivity" , "salinity" , "PAR" , "surf_PAR" , "oxygen" , "O_saturation" , "NH4" , "PO4" , "fluorescence" , "Chla" , "Pprod_Sun" , "Pprod_PAR" , "Babun" , "Bprod" , "picoDiversity" ,  "microDiversity")
AICresults.nano <- AICc.table.all(env.var.nano, matrix.char=distance(nano,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(nano),"data.frame"),comb.incl=c(1), control.var.char = c("NO2NO3 + time + ZML_TS + depth + beam_trans"))
write.table(AICresults.nano, file="AIC_Table_nano.txt", append=T, sep="\t")

# nano: Four parameters were important regarding the nanoplankton community: as for the microplankton, the nitrate/nitrite, then the time of sampling and the depth of the mixed layer have also played for the nanoplankton. In additon, the depth of sampling was also responsible of the nanoplankton composition. Additional parameters are not adding explanation to the variance. suggesting negligiable role.

env.var.pico <- c("latitude","longitude","time" , "bottom" , "depth" , "ZML_TS" , "Ze" , "ZCM" , "ice" , "air_temp" , "water_temp" , "conductivity" , "salinity" , "PAR" , "surf_PAR" , "oxygen" , "O_saturation" , "beam_trans" , "NH4" , "NO2NO3" , "PO4" , "fluorescence" , "Chla" , "Pprod_Sun" , "Pprod_PAR" , "Babun" , "Bprod" , "nanoDiversity" ,  "microDiversity")
AICresults.pico <- AICc.table.all(env.var.pico, matrix.char=distance(pico,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(pico),"data.frame"),comb.incl=c(1))
write.table(AICresults.pico, file="AIC_Table_pico2.txt", append=T, sep="\t")

env.var.pico <- c( "water_temp", "time", "conductivity")
AICresults.pico <- AICc.table.all(env.var.pico, matrix.char=distance(pico,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(pico),"data.frame"),comb.incl=c(3), control.var.char = c("NO2NO3" ))
write.table(AICresults.pico, file="AIC_Table_pico3.txt", append=T, sep="\t")

#pico: Four parameters were important regarding the nanoplankton community: as for the microplankton and nanoplankton, the nitrate/nitrite, then the time of sampling have also played for the picoplankton. In additon, the water temperature and conductivity were impacting the picoplankton composition. Additional parameters are not adding explanation to the variance. suggesting negligiable/unclear role.
#removing parameters with missing data (nut and prod) lat, Bprod, time and O sat are the most important factor regarding the pico plankton
#adonis2(formula = distance(pico, "jaccard") ~ NO2NO3 + water_temp + time + conductivity + ice + depth + bottom + ZML_TS + beam_trans, data = as(sample_data(pico), "data.frame"), permutations = 9999)
Df SumOfSqs      R2      F Pr(>F)
NO2NO3        1   1.4991 0.15810 8.9175 0.0001 ***
water_temp    1   0.8768 0.09247 5.2157 0.0001 ***
time          1   0.7888 0.08319 4.6922 0.0001 ***
conductivity  1   0.6646 0.07009 3.9533 0.0002 ***
ice           1   0.4248 0.04480 2.5271 0.0040 **
depth         1   0.3826 0.04035 2.2760 0.0094 **
bottom        1   0.3798 0.04005 2.2592 0.0101 *
ZML_TS        1   0.2956 0.03117 1.7584 0.0478 *
beam_trans    1   0.3033 0.03198 1.8039 0.0371 *
Residual     23   3.8666 0.40778
Total        32   9.4821 1.00000

Permutation test for adonis under reduced model
Terms added sequentially (first to last)
Permutation: free
Number of permutations: 9999

adonis2(formula = distance(pico, "jaccard") ~ latitude + Bprod + time + O_saturation + water_temp + ice + conductivity, data = as(sample_data(pico), "data.frame"), permutations = 9999)
             Df SumOfSqs      R2       F Pr(>F)
latitude      1   1.7311 0.18678 10.3476 0.0001 ***
Bprod         1   0.8796 0.09491  5.2580 0.0001 ***
time          1   0.8805 0.09500  5.2632 0.0001 ***
O_saturation  1   0.6026 0.06502  3.6021 0.0003 ***
water_temp    1   0.4359 0.04703  2.6056 0.0035 **
ice           1   0.3468 0.03742  2.0730 0.0223 *
conductivity  1   0.2094 0.02259  1.2514 0.2135
Residual     25   4.1824 0.45126
Total        32   9.2684 1.00000
---
# rerun without B, D and bucket => missing nuts  
N because not used
physeqb = subset_samples(physeq, station != "B")
physeqb = subset_samples(physeqb, station != "D")
physeqb = subset_samples(physeqb, layer != "Bucket")
torun <- subset_samples(physeqb, size == "pico")
env.var <- c("bottom" , "depth" , "ZML_TS" , "Ze" , "ZCM" , "ice" , "air_temp" , "water_temp" , "conductivity" , "salinity" , "PAR" , "surf_PAR" , "oxygen" , "O_saturation" , "beam_trans" , "NH4" , "NO2NO3" , "PO4" , "fluorescence" , "Chla" , "Pprod_Sun" , "Pprod_PAR" , "Babun" , "Bprod" , "picoDiversity" , "nanoDiversity")
source('/home/tuk61790/Network/AICc_table_generation_edited.R',chdir =TRUE)
source('AICc_PERMANOVA.R',chdir=T)
AICresults.pico <- AICc.table.all(env.var, matrix.char=distance(torun,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(torun),"data.frame"),comb.incl=c(1))
write.table(AICresults.pico, file="AIC_Table_pico_BD.txt", append=T, sep="\t")
env.var.pico <- c("bottom" , "depth" , "ZML_TS" , "Ze" , "ZCM" , "ice" , "air_temp" , "conductivity" , "salinity" , "PAR" , "surf_PAR" , "oxygen" , "O_saturation" , "beam_trans" , "NH4" , "PO4" , "fluorescence" ,  "Bprod" )
AICresults.pico <- AICc.table.all(env.var.pico, matrix.char=distance(torun,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(torun),"data.frame"),comb.incl=c(1), control.var.char = c("NO2NO3 +  water_temp" ))
write.table(AICresults.pico, file="AIC_Table_pico_BD.txt", append=T, sep="\t")
env.var.pico <- c( "NO2NO3", "water_temp", "time", "conductivity" , "Bprod", "ice", "depth", "ZML_TS", "oxygen")
AICresults.pico <- AICc.table.all(env.var.pico, matrix.char=distance(torun,"jaccard"),perm=9999,method="jaccard", df=as(sample_data(torun),"data.frame"),comb.incl=c(1,2,3,4,5,6,7,8,9))
write.table(AICresults.pico, file="AIC_Table_pico_BD2.txt", append=T, sep="\t")

env.var.pico <- c( "NO2NO3", "water_temp", "conductivity" , "Bprod", "ice", "depth", "ZML_TS", "oxygen", "beam_trans", "bottom")
AICresults.pico <- AICc.table.all(env.var.pico, matrix.char=distance(torun,"jaccard"),perm=999,method="jaccard", df=as(sample_data(torun),"data.frame"),comb.incl=c(4,10))
write.table(AICresults.pico, file="AIC_Table_pico_BD2.txt", append=T, sep="\t")

torun <- subset_samples(physeqb, size == "nano")
env.var.pico <- c( "NO2NO3", "water_temp", "conductivity" , "Bprod", "ice", "depth", "ZML_TS", "oxygen", "beam_trans", "bottom", "fluorescence")
AICresults.pico <- AICc.table.all(env.var.pico, matrix.char=distance(torun,"jaccard"),perm=999,method="jaccard", df=as(sample_data(torun),"data.frame"),comb.incl=c(4,10))
write.table(AICresults.pico, file="AIC_Table_nano_BD2.txt", append=T, sep="\t")

torun <- subset_samples(physeqb, size == "micro")
env.var.pico <- c( "NO2NO3", "water_temp", "conductivity" , "Bprod", "ice", "depth", "ZML_TS", "oxygen", "beam_trans", "bottom", "fluorescence")
AICresults.pico <- AICc.table.all(env.var.pico, matrix.char=distance(torun,"jaccard"),perm=999,method="jaccard", df=as(sample_data(torun),"data.frame"),comb.incl=c(4,10))
write.table(AICresults.pico, file="AIC_Table_micro_BD2.txt", append=T, sep="\t")
