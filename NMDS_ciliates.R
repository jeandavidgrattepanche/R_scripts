setwd('/home/tuk61790/Network')
library('phytools')
library('phyloseqCompanion')
library("ggplot2")
library("plyr")

data.all <- read.table('OTUtable_ingroup_PE.txt',sep="\t",header=TRUE,row.names=1)
env.all <- read.table('NBP1910_envdata_v4.5_edited.txt',sep="\t",header=TRUE,row.names=1)
env.all.sc <- data.frame(cbind(env.all[c(0:7)],scale(env.all[c(8:41)]))) 
envdata = sample_data(env.all.sc)
mytable = otu_table(cbind(data.all[c(0)],data.all[c(10:105)]), taxa_are_rows=TRUE,errorIfNULL=TRUE)
envdata = sample_data(env.all)
testb <- as.matrix(data.all[c(0:10)])
TAX <- tax_table(testb)
tree <- phytools::read.newick('rename2.tree')
phy_tree(tree, errorIfNULL=TRUE)
physeq <- phyloseq(mytable, envdata, tree, TAX)
pico <- subset_samples(physeq, size == "pico")
nano <- subset_samples(physeq, size == "nano")
micro <- subset_samples(physeq, size == "micro")
dUNIFRAC = UniFrac(physeq, weighted=TRUE, normalized=FALSE, fast=TRUE)


cil <- subset_taxa(physeq, T4 == "Ciliophora")
cilb = prune_samples(sample_sums(cil) != 0, cil)
cilb = transform_sample_counts(cilb, function(x) 100 * x/sum(x)) 
pico <- subset_samples(cilb, size == "pico")
nano <- subset_samples(cilb, size == "nano")
micro <- subset_samples(cilb, size == "micro")

ordMDS <- metaMDS(t(otu_table(cilb)), distance = "jaccard")
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(cilb), perm=999, na.rm = TRUE)
pdf("NMDS_Ciliates.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()


ciliates:
> ordMDS.fit

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)    
Bprod        -0.58806  0.80882 0.2553  0.001 ***
Pprod_PAR    -0.59321  0.80505 0.1530  0.001 ***
ice           0.03859  0.99926 0.0098  0.640    
air_temp     -0.49316  0.86994 0.0520  0.072 .  
depth         0.60948 -0.79281 0.1711  0.001 ***
PAR          -0.62909  0.77733 0.1657  0.002 ** 
salinity      0.42026 -0.90740 0.0705  0.037 *  
oxygen       -0.85045  0.52605 0.0459  0.116    
O_saturation  0.52139 -0.85332 0.1283  0.002 ** 
fluorescence -0.84285  0.53815 0.0402  0.151    
beam_trans    0.60895 -0.79321 0.1998  0.001 ***
water_temp   -0.51473  0.85735 0.1636  0.001 ***
conductivity -0.54113  0.84094 0.0666  0.051 .  
NH4          -0.48842  0.87261 0.1425  0.002 ** 
NO2NO3        0.93910  0.34365 0.0116  0.619    
PO4           0.20501  0.97876 0.0062  0.778    
Chla         -0.59675  0.80242 0.1512  0.002 ** 
Chao1        -0.26678 -0.96376 0.2874  0.001 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

torun <- pico

ordMDS <- metaMDS(t(otu_table(torun)), distance = "jaccard")
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_Ciliates_pico.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()
> ordMDS.fit

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)    
Bprod        -0.91941  0.39331 0.2439  0.017 *  
Pprod_PAR    -0.96150  0.27481 0.1712  0.068 .  
ice          -0.98790  0.15512 0.0206  0.725    
air_temp     -0.38254  0.92394 0.1258  0.108    
depth         0.77014  0.63788 0.4875  0.001 ***
PAR          -0.95277 -0.30370 0.1931  0.047 *  
salinity      0.82033  0.57190 0.2003  0.038 *  
oxygen       -0.52100 -0.85356 0.1056  0.195    
O_saturation  0.19211 -0.98137 0.4241  0.003 ** 
fluorescence -0.18659  0.98244 0.0722  0.333    
beam_trans    0.46981 -0.88277 0.2867  0.016 *  
water_temp   -0.24685  0.96905 0.4539  0.003 ** 
conductivity -0.07326  0.99731 0.3701  0.006 ** 
NH4          -0.93556 -0.35316 0.2846  0.012 *  
NO2NO3        0.06688  0.99776 0.0388  0.554    
PO4           0.00262  1.00000 0.1029  0.216    
Chla         -0.53089  0.84744 0.1903  0.060 .  
Chao1         0.10001 -0.99499 0.1226  0.125    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

torun <- nano

ordMDS <- metaMDS(t(otu_table(torun)), distance = "jaccard")
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_Ciliates_nano.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()

> ordMDS.fit

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)   
Bprod        -0.79276  0.60953 0.3345  0.004 **
Pprod_PAR    -0.53946  0.84201 0.2138  0.030 * 
ice          -0.23080 -0.97300 0.3512  0.002 **
air_temp     -0.77663 -0.62996 0.1132  0.195   
depth         0.78050 -0.62515 0.2075  0.044 * 
PAR          -0.78889  0.61453 0.1828  0.063 . 
salinity      0.94681 -0.32180 0.0748  0.338   
oxygen       -0.13365  0.99103 0.1742  0.084 . 
O_saturation  0.84614 -0.53297 0.2750  0.019 * 
fluorescence -0.35875  0.93343 0.2275  0.029 * 
beam_trans    0.77418 -0.63297 0.4076  0.002 **
water_temp   -0.85352  0.52106 0.3265  0.011 * 
conductivity -0.83399  0.55178 0.1739  0.082 . 
NH4          -0.96083 -0.27712 0.1017  0.233   
NO2NO3       -0.10581 -0.99439 0.2819  0.008 **
PO4          -0.20496 -0.97877 0.2771  0.005 **
Chla         -0.59586  0.80309 0.3078  0.009 **
Chao1         0.58891  0.80820 0.4143  0.002 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

torun <- micro

ordMDS <- metaMDS(t(otu_table(torun)), distance = "jaccard")
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_Ciliates_micro.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()

 ordMDS.fit

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)    
Bprod        -0.71602 -0.69808 0.3102  0.003 ** 
Pprod_PAR    -0.49235 -0.87040 0.2620  0.010 ** 
ice          -0.32740  0.94489 0.0394  0.572    
air_temp     -0.33222  0.94320 0.1482  0.092 .  
depth         0.80838  0.58866 0.1738  0.067 .  
PAR          -0.89584 -0.44438 0.2382  0.017 *  
salinity      0.51026  0.86002 0.0871  0.271    
oxygen       -0.21709 -0.97615 0.1664  0.076 .  
O_saturation  0.37437 -0.92728 0.2945  0.004 ** 
fluorescence -0.16033 -0.98706 0.0700  0.358    
beam_trans    0.78150 -0.62390 0.1992  0.036 *  
water_temp   -0.41589  0.90942 0.3285  0.001 ***
conductivity -0.27899  0.96029 0.2264  0.026 *  
NH4          -0.99385 -0.11071 0.1763  0.068 .  
NO2NO3       -0.00067  1.00000 0.1643  0.082 .  
PO4          -0.27005  0.96285 0.0399  0.577    
Chla         -0.99927  0.03816 0.1243  0.148    
Chao1         0.78518 -0.61927 0.2971  0.006 ** 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

