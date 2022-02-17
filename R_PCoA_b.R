setwd('/home/tuk61790/Network/')
library('ape')
library('phytools')
library('phyloseqCompanion')
library("ggplot2")
library("plyr")
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

Adonis.results = adonis(formula = distance(physeq, "jaccard") ~ size * goup * latitude * layer, data = as(sample_data(physeq), "data.frame"), permutations = 9999) 
cat("Adonis_results", capture.output(Adonis.results), file="summary_of_Adonis_results021722.txt", sep="\n", append=TRUE)

#adonis2 for some parameters and their interaction (*)
Adonis2.results = adonis2(formula = distance(physeq, "jaccard") ~ size * group * latitude * layer, data = as(sample_data(physeq), "data.frame"), permutations = 9999) 
cat("Adonis_results", capture.output(Adonis2.results), file="summary_of_Adonis2_results021722.txt", sep="\n", append=TRUE)
# try term= NULL or "margin" or "term"
#tried mixing parameters
Adonis2.results = adonis2(formula = distance(physeq, "jaccard") ~ layer * latitude * size * group , data = as(sample_data(physeq), "data.frame"), permutations = 9999) 
cat("Adonis_results", capture.output(Adonis2.results), file="summary_of_Adonis2b_results021722.txt", sep="\n", append=TRUE)

Adonis2.results = adonis2(formula = distance(physeq, "jaccard") ~ group + layer + latitude + size  , data = as(sample_data(physeq), "data.frame"), permutations = 9999, by = "margin") 
cat("Adonis_results", capture.output(Adonis2.results), file="summary_of_Adonis2bmargin2_results021722.txt", sep="\n", append=TRUE)

Adonis3.results = adonis2(formula = distance(physeq, "jaccard") ~ longitude + latitude + time + bottom + depth + ZML_TS + Ze + ZCM + ice + air_temp + water_temp + conductivity + salinity + PAR + surf_PAR + oxygen + O_saturation + beam_trans + fluorescence , data = as(sample_data(physeq), "data.frame"), permutations = 9999, by = "margin") 
cat("Adonis_results", capture.output(Adonis3.results), file="summary_of_Adonis_results_env.txt", sep="\n", append=TRUE)


Adonis3.results = adonis2(formula = distance(physeq, "jaccard") ~ longitude + latitude + time + bottom + depth + ZML_TS + Ze + ZCM + ice + air_temp + water_temp + conductivity + salinity + PAR + surf_PAR + oxygen + O_saturation + beam_trans + NH4 + NO2NO3 + PO4 + fluorescence + Chla + Pprod_Sun + Pprod_PAR + Pprod_PAR_20 + Babun + Bprod, data = as(sample_data(physeq), "data.frame"), permutations = 9999, by = "margin") 
cat("Adonis_results", capture.output(Adonis3.results), file="summary_of_Adonis_results_env.txt", sep="\n", append=TRUE)

Adonis3.results = adonis2(formula = distance(physeq, "jaccard") ~ bottom * ice * airT * depth_m * size_num *  depSM * waterT1 * conductivity1 * oxygen1 *  fluorescence *  beamTrans * PAR1 *  latitude *  longitude * timeJ *  altM *  spar * salinity1 * oxygenSaturation * bprod * pprod_Sun * pprod_PAR, data = as(sample_data(physeq), "data.frame"), permutations = 9999) 
cat("Adonis_results", capture.output(Adonis.results), file="summary_of_Adonis_results.txt", sep="\n", append=TRUE)

#envfit for each parameters in my env_data table (+)
ordMDS <- metaMDS(t(otu_table(physeq)))
ordMDS.fit <- envfit(ordMDS ~ group + size + layer + longitude + latitude + time + bottom + depth + ZML_TS + Ze + ZCM + ice + air_temp + water_temp + conductivity + salinity + PAR + surf_PAR + oxygen + O_saturation + beam_trans + NH4 + NO2NO3 + PO4 + fluorescence + Chla + Pprod_Sun + Pprod_PAR + Pprod_PAR_20 + Babun + Bprod, data=sample.data.table(physeq), perm=999, na.rm = TRUE)
pdf("envfit_1.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()
cat("envfit_results", capture.output(ordMDS.fit), file="summary_of_envfit_results.txt", sep="\n", append=TRUE)

need to look by size or grouping the size