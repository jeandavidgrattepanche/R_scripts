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

#### use all OTUs (PE) and not only the abundant OTUs (100)

torun <- pico
dUNIFRAC = UniFrac(torun, weighted=TRUE, normalized=TRUE, fast=TRUE)
ordMDS <- metaMDS(t(otu_table(torun)), comm = dUNIFRAC)
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_PE_UnifracW_Pico.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()

 ordMDS.fit
 ***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)   
Bprod         0.38482  0.92299 0.0243  0.705   
Pprod_PAR    -0.11183  0.99373 0.0895  0.264   
ice           0.79623 -0.60499 0.2590  0.017 * 
air_temp      0.98560  0.16907 0.1524  0.107   
depth        -0.85174 -0.52397 0.1139  0.199   
PAR           0.79170  0.61091 0.1797  0.062 . 
salinity      0.34679  0.93794 0.0059  0.923   
oxygen       -0.98755 -0.15732 0.0049  0.948   
O_saturation -0.57923 -0.81517 0.3596  0.004 **
fluorescence -0.85844  0.51291 0.1379  0.124   
beam_trans   -0.76571 -0.64319 0.2038  0.056 . 
water_temp    0.59192  0.80599 0.3855  0.004 **
conductivity  0.57929  0.81512 0.2875  0.010 **
NH4           0.84629  0.53273 0.0365  0.588   
NO2NO3        0.98652 -0.16365 0.1093  0.199   
PO4           0.91113 -0.41213 0.0253  0.672   
Chla          0.37815  0.92574 0.1082  0.193   
Chao1        -0.77306  0.63433 0.0004  0.997   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

dUNIFRAC = UniFrac(torun, weighted=FALSE, normalized=TRUE, fast=TRUE)

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)    
Bprod        -0.10957 -0.99398 0.2095  0.029 *  
Pprod_PAR    -0.01995 -0.99980 0.1299  0.127    
ice          -0.38347  0.92355 0.5204  0.001 ***
air_temp     -0.95488  0.29701 0.1727  0.051 .  
depth         0.68506  0.72848 0.1459  0.097 .  
PAR          -0.75913 -0.65094 0.2006  0.047 *  
salinity     -0.43113  0.90229 0.0229  0.711    
oxygen        0.30553 -0.95218 0.1619  0.071 .  
O_saturation  0.48515  0.87443 0.3304  0.002 ** 
fluorescence  0.48015 -0.87718 0.4144  0.003 ** 
beam_trans    0.26972  0.96294 0.2301  0.009 ** 
water_temp   -0.46527 -0.88517 0.3691  0.001 ***
conductivity -0.53742 -0.84331 0.2390  0.021 *  
NH4          -0.86929 -0.49430 0.1126  0.208    
NO2NO3       -0.56324  0.82629 0.3590  0.002 ** 
PO4          -0.62043  0.78426 0.0956  0.223    
Chla          0.00814 -0.99997 0.2632  0.006 ** 
Chao1        -0.81517  0.57922 0.0434  0.526    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

cil <- subset_taxa(pico, T4 == "Cil")
cil = prune_samples(sample_sums(cil) != 0, cil)
cil = transform_sample_counts(cil, function(x) 100 * x/sum(x)) 

torun <- cil
dUNIFRAC = UniFrac(torun, weighted=TRUE, normalized=TRUE, fast=TRUE)
ordMDS <- metaMDS(t(otu_table(torun)), comm = dUNIFRAC)
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_PE_UnifracW_PicoCil.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()
***VECTORS

                  NMDS1      NMDS2     r2 Pr(>r)  
Bprod        -0.0004174  1.0000000 0.0998  0.187  
Pprod_PAR    -0.0006040  1.0000000 0.0972  0.183  
ice           0.0002711 -1.0000000 0.0348  0.819  
air_temp      0.0000304  1.0000000 0.0454  0.342  
depth         0.0001869 -1.0000000 0.1844  0.043 *
PAR          -0.0000803  1.0000000 0.2089  0.075 .
salinity      0.0010819  1.0000000 0.0310  0.636  
oxygen       -0.0006733  1.0000000 0.1062  0.200  
O_saturation -0.0000237 -1.0000000 0.1400  0.134  
fluorescence -0.0012768 -1.0000000 0.0308  0.616  
beam_trans    0.0001270 -1.0000000 0.1290  0.124  
water_temp    0.0000014  1.0000000 0.1486  0.112  
conductivity  0.0000664  1.0000000 0.1152  0.173  
NH4          -0.0182003  0.9998300 0.0514  0.440  
NO2NO3        0.0066284  0.9999800 0.0472  0.592  
PO4           0.0007161 -1.0000000 0.0510  0.477  
Chla         -0.0002750  1.0000000 0.0896  0.178  
Chao1        -0.0002976 -1.0000000 0.0159  0.805  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999


dUNIFRAC = UniFrac(torun, weighted=FALSE, normalized=TRUE, fast=TRUE)
ordMDS <- metaMDS(t(otu_table(torun)), comm = dUNIFRAC)
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_PE_UnifracUnW_PicoCil.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()


***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)    
Bprod         0.37923 -0.92530 0.2351  0.029 *  
Pprod_PAR     0.67641 -0.73652 0.0475  0.493    
ice           0.99991 -0.01321 0.0704  0.361    
air_temp      0.97011  0.24267 0.0897  0.256    
depth        -0.26529  0.96417 0.4412  0.001 ***
PAR           0.81213 -0.58347 0.1983  0.046 *  
salinity     -0.17566  0.98445 0.1074  0.197    
oxygen       -0.23930 -0.97095 0.1215  0.148    
O_saturation -0.89803 -0.43994 0.2419  0.016 *  
fluorescence -0.70350 -0.71069 0.0240  0.686    
beam_trans   -0.91386  0.40602 0.1324  0.142    
water_temp    0.94131  0.33756 0.2671  0.009 ** 
conductivity  0.76531  0.64367 0.1950  0.037 *  
NH4           0.55326 -0.83301 0.2651  0.008 ** 
NO2NO3        0.59607  0.80293 0.1306  0.150    
PO4           0.52874  0.84878 0.1344  0.122    
Chla          0.77685 -0.62968 0.0576  0.410    
Chao1        -0.37262 -0.92798 0.0046  0.945    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

Mam <- subset_taxa(physeq, T4 == "Mam")
Mam = prune_samples(sample_sums(Mam) != 0, Mam)
Mam = transform_sample_counts(Mam, function(x) 100 * x/sum(x)) 

dUNIFRAC = UniFrac(torun, weighted=FALSE, normalized=TRUE, fast=TRUE)
ordMDS <- metaMDS(t(otu_table(torun)), comm = dUNIFRAC)
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_PE_UnifracUnW_Mam.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)   
Bprod        -0.76190  0.64769 0.0453  0.176   
Pprod_PAR    -0.85660  0.51599 0.0319  0.316   
ice           0.10446 -0.99453 0.1410  0.004 **
air_temp     -0.01556 -0.99988 0.0740  0.065 . 
depth         0.05330  0.99858 0.0881  0.035 * 
PAR          -0.05451 -0.99851 0.0641  0.088 . 
salinity      0.37684 -0.92628 0.0120  0.637   
oxygen       -0.06678  0.99777 0.1365  0.004 **
O_saturation  0.22184  0.97508 0.0507  0.157   
fluorescence -0.16442  0.98639 0.1342  0.007 **
beam_trans    0.95229  0.30518 0.0322  0.291   
water_temp   -0.24862 -0.96860 0.0576  0.118   
conductivity -0.17490 -0.98459 0.0348  0.295   
NH4          -0.18540 -0.98266 0.0167  0.535   
NO2NO3        0.09941 -0.99505 0.1615  0.002 **
PO4           0.13539 -0.99079 0.0144  0.586   
Chla         -0.59476  0.80390 0.0575  0.106   
Chao1         0.16413  0.98644 0.0216  0.440   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

Mam <- subset_taxa(pico, genus == "Micromonas")
Mam = prune_samples(sample_sums(Mam) != 0, Mam)
Mam = transform_sample_counts(Mam, function(x) 100 * x/sum(x)) 

dUNIFRAC = UniFrac(torun, weighted=FALSE, normalized=TRUE, fast=TRUE)
ordMDS <- metaMDS(t(otu_table(torun)), comm = dUNIFRAC)
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_PE_UnifracUnW_picoMicromonas.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)    
Bprod         0.84285  0.53815 0.2602  0.013 *  
Pprod_PAR     0.69800  0.71610 0.3403  0.003 ** 
ice          -0.41315 -0.91066 0.3663  0.002 ** 
air_temp     -0.00506 -0.99999 0.1010  0.214    
depth        -0.94029  0.34038 0.0306  0.652    
PAR           0.44942 -0.89332 0.0944  0.235    
salinity     -0.35540 -0.93471 0.1095  0.189    
oxygen        0.20478  0.97881 0.2182  0.025 *  
O_saturation -0.55919  0.82904 0.2699  0.007 ** 
fluorescence  0.36045  0.93278 0.4074  0.001 ***
beam_trans   -0.99726  0.07392 0.1944  0.030 *  
water_temp    0.62691 -0.77909 0.3039  0.001 ***
conductivity  0.43381 -0.90101 0.2009  0.030 *  
NH4           0.77620 -0.63049 0.0215  0.725    
NO2NO3       -0.24163 -0.97037 0.4279  0.002 ** 
PO4          -0.21431 -0.97677 0.1223  0.141    
Chla          0.81938  0.57325 0.3399  0.001 ***
Chao1        -0.72468 -0.68909 0.0856  0.296    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

dUNIFRAC = UniFrac(torun, weighted=TRUE, normalized=TRUE, fast=TRUE)
ordMDS <- metaMDS(t(otu_table(torun)), comm = dUNIFRAC)
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_PE_UnifracW_picoMicromonas.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)  
Bprod         0.83177  0.55512 0.0618  0.380  
Pprod_PAR    -0.00507  0.99999 0.0528  0.412  
ice           0.92003  0.39184 0.1706  0.041 *
air_temp      0.04409  0.99903 0.0449  0.309  
depth        -0.42917 -0.90322 0.2262  0.019 *
PAR          -0.15617  0.98773 0.3316  0.029 *
salinity     -0.99042 -0.13811 0.2143  0.035 *
oxygen        0.99966 -0.02607 0.0962  0.244  
O_saturation  0.70351 -0.71068 0.1375  0.113  
fluorescence  0.09356 -0.99561 0.0376  0.575  
beam_trans   -0.86215 -0.50666 0.0978  0.121  
water_temp   -0.60909  0.79310 0.1133  0.131  
conductivity -0.83529  0.54981 0.1604  0.105  
NH4           0.05096  0.99870 0.0854  0.276  
NO2NO3       -0.52444  0.85144 0.0653  0.428  
PO4          -0.86615 -0.49978 0.0151  0.811  
Chla          0.90957  0.41554 0.0442  0.397  
Chao1        -0.82259  0.56864 0.1449  0.115  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999
