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

#ordMDS <- metaMDS(t(otu_table(cilb)), distance = "jaccard")
ordMDS <- metaMDS(t(otu_table(cilb)), comm = dUNIFRAC)
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

mixo <- subset_taxa(physeq, Trophic =="mixo?" | Trophic == "Mixo" | Trophic == "Klepto?")
mixo = prune_samples(sample_sums(mixo) != 0, mixo)
mixo = transform_sample_counts(mixo, function(x) 100 * x/sum(x)) 
pico <- subset_samples(mixo, size == "pico")
nano <- subset_samples(mixo, size == "nano")
micro <- subset_samples(mixo, size == "micro")
torun <- micro
ordMDS <- metaMDS(t(otu_table(torun)), distance = "jaccard")
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_mixo_micro.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()
ordMDS.fit

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)    
Bprod        -0.55583  0.83129 0.4468  0.001 ***
Pprod_PAR    -0.45534  0.89032 0.2549  0.012 *  
ice          -0.28962 -0.95714 0.2708  0.007 ** 
air_temp     -0.92790  0.37283 0.0693  0.332    
depth         0.92944 -0.36898 0.1441  0.098 .  
PAR          -0.81332  0.58181 0.1805  0.048 *  
salinity      0.95637 -0.29216 0.0471  0.478    
oxygen       -0.09575  0.99541 0.2617  0.011 *  
O_saturation  0.48033 -0.87709 0.3512  0.004 ** 
fluorescence -0.12169  0.99257 0.4658  0.001 ***
beam_trans    0.33725 -0.94142 0.6316  0.001 ***
water_temp   -0.49349  0.86975 0.4111  0.002 ** 
conductivity -0.42467  0.90535 0.2407  0.020 *  
NH4          -0.76359 -0.64570 0.3493  0.004 ** 
NO2NO3       -0.15103 -0.98853 0.4513  0.001 ***
PO4          -0.20168 -0.97945 0.1211  0.137    
Chla         -0.24194  0.97029 0.6358  0.001 ***
Chao1         0.88045  0.47414 0.4601  0.001 ***
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

torun <- nano
ordMDS <- metaMDS(t(otu_table(torun)), distance = "jaccard")
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_nano_micro.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()
ordMDS.fit

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)    
Bprod        -0.49137  0.87095 0.1513  0.106    
Pprod_PAR    -0.20510  0.97874 0.1911  0.056 .  
ice          -0.53511 -0.84478 0.4054  0.001 ***
air_temp     -0.72611 -0.68758 0.1880  0.063 .  
depth         0.57559 -0.81774 0.4509  0.001 ***
PAR          -0.99890  0.04680 0.1030  0.224    
salinity      0.21111 -0.97746 0.2100  0.054 .  
oxygen        0.07500  0.99718 0.2499  0.019 *  
O_saturation  0.86667  0.49887 0.0770  0.326    
fluorescence  0.15438  0.98801 0.1288  0.181    
beam_trans    0.61152 -0.79123 0.1489  0.133    
water_temp   -0.97193 -0.23528 0.0864  0.281    
conductivity -0.58715 -0.80948 0.0717  0.362    
NH4          -0.99986  0.01666 0.1419  0.135    
NO2NO3       -0.31359 -0.94956 0.4579  0.002 ** 
PO4          -0.21459 -0.97670 0.4154  0.001 ***
Chla         -0.18348  0.98302 0.1786  0.097 .  
Chao1         0.15123  0.98850 0.2480  0.020 *  
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

torun <- pico
ordMDS <- metaMDS(t(otu_table(torun)), distance = "jaccard")
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_pico_micro.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()
ordMDS.fit
***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)    
Bprod         0.15260 -0.98829 0.3712  0.004 ** 
Pprod_PAR     0.08389 -0.99647 0.1985  0.048 *  
ice           0.84559 -0.53383 0.0830  0.268    
air_temp      0.83982  0.54286 0.2592  0.012 *  
depth        -0.26560  0.96408 0.5631  0.001 ***
PAR           0.99294  0.11860 0.1626  0.064 .  
salinity      0.06979  0.99756 0.6947  0.001 ***
oxygen       -0.11380 -0.99350 0.3350  0.006 ** 
O_saturation -0.55461 -0.83211 0.4071  0.001 ***
fluorescence -0.26800 -0.96342 0.0635  0.392    
beam_trans   -0.66571  0.74621 0.1754  0.048 *  
water_temp    0.65206  0.75817 0.3619  0.001 ***
conductivity  0.42414  0.90560 0.4845  0.001 ***
NH4           0.37415 -0.92737 0.1989  0.032 *  
NO2NO3        0.34320  0.93926 0.4308  0.002 ** 
PO4           0.30300  0.95299 0.1774  0.054 .  
Chla          0.21898 -0.97573 0.1450  0.095 .  
Chao1         0.09452  0.99552 0.0803  0.299    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999


phago <- subset_taxa(physeq, Trophic =="phagotroph" | Trophic == "Phagotroph" )
phago = prune_samples(sample_sums(phago) != 0, phago)
phago = transform_sample_counts(phago, function(x) 100 * x/sum(x)) 
pico <- subset_samples(phago, size == "pico")
nano <- subset_samples(phago, size == "nano")
micro <- subset_samples(phago, size == "micro")
torun <- micro
ordMDS <- metaMDS(t(otu_table(torun)), distance = "jaccard")
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_phago_micro.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()
ordMDS.fit

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)   
Bprod        -0.94667  0.32221 0.2599  0.010 **
Pprod_PAR    -0.63706  0.77082 0.1531  0.077 . 
ice          -0.26482 -0.96430 0.1451  0.110   
air_temp     -0.82853 -0.55994 0.0620  0.384   
depth         0.99861  0.05266 0.1872  0.041 * 
PAR          -0.88594  0.46380 0.1463  0.084 . 
salinity      0.98232 -0.18723 0.0640  0.331   
oxygen       -0.70849  0.70572 0.0810  0.262   
O_saturation  0.67922  0.73393 0.1170  0.162   
fluorescence -0.74724  0.66456 0.0378  0.544   
beam_trans    0.64723  0.76230 0.2249  0.014 * 
water_temp   -0.70418 -0.71002 0.1465  0.091 . 
conductivity -0.55408 -0.83246 0.0701  0.327   
NH4          -0.99564  0.09326 0.2325  0.016 * 
NO2NO3        0.06809 -0.99768 0.0955  0.203   
PO4          -0.41244 -0.91098 0.0202  0.714   
Chla         -0.98650 -0.16376 0.0880  0.256   
Chao1         0.98159 -0.19099 0.2176  0.020 * 
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

torun <- nano
ordMDS <- metaMDS(t(otu_table(torun)), distance = "jaccard")
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_phago_nano.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()
ordMDS.fit

***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)   
Bprod        -0.87167  0.49010 0.3129  0.004 **
Pprod_PAR    -0.64386  0.76515 0.1840  0.069 . 
ice          -0.24639 -0.96917 0.1372  0.158   
air_temp     -0.37410 -0.92739 0.1202  0.163   
depth         0.98629 -0.16504 0.2287  0.028 * 
PAR          -0.86279  0.50556 0.1759  0.077 . 
salinity      0.86134 -0.50803 0.1514  0.090 . 
oxygen       -0.25760  0.96625 0.1424  0.135   
O_saturation  0.90631  0.42262 0.1605  0.092 . 
fluorescence -0.50040  0.86579 0.1461  0.131   
beam_trans    0.99933  0.03649 0.2431  0.038 * 
water_temp   -0.93522 -0.35407 0.2013  0.058 . 
conductivity -0.78488 -0.61965 0.0899  0.250   
NH4          -0.91988  0.39219 0.2061  0.051 . 
NO2NO3        0.01034 -0.99995 0.2364  0.025 * 
PO4          -0.06183 -0.99809 0.1608  0.094 . 
Chla         -0.83048  0.55706 0.1773  0.075 . 
Chao1         0.31405  0.94941 0.3101  0.008 **
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999

torun <- pico
ordMDS <- metaMDS(t(otu_table(torun)), distance = "jaccard")
ordMDS.fit <- envfit(ordMDS ~ Bprod + Pprod_PAR + ice + air_temp + depth + PAR + salinity + oxygen + O_saturation + fluorescence + beam_trans + water_temp + conductivity + NH4 + NO2NO3 + PO4 + Chla + Chao1 , data=sample.data.table(torun), perm=999, na.rm = TRUE)
pdf("NMDS_phago_pico.pdf")
plot(ordMDS, dis="site")
plot(ordMDS.fit)
text(ordMDS, display="sites")
dev.off()
ordMDS.fit
***VECTORS

                NMDS1    NMDS2     r2 Pr(>r)    
Bprod        -0.79399 -0.60793 0.1728  0.070 .  
Pprod_PAR    -0.93283 -0.36032 0.1544  0.087 .  
ice          -0.31678 -0.94850 0.0584  0.410    
air_temp     -0.56233  0.82691 0.1141  0.140    
depth         0.63619  0.77153 0.2550  0.014 *  
PAR          -0.96429 -0.26485 0.1271  0.106    
salinity      0.35602  0.93448 0.1602  0.086 .  
oxygen       -0.36031 -0.93283 0.0496  0.474    
O_saturation  0.37576 -0.92672 0.3575  0.005 ** 
fluorescence -0.73487  0.67820 0.0163  0.784    
beam_trans    0.57360 -0.81913 0.2541  0.027 *  
water_temp   -0.42131  0.90692 0.3753  0.004 ** 
conductivity -0.27751  0.96072 0.3212  0.006 ** 
NH4          -0.42335 -0.90596 0.3900  0.001 ***
NO2NO3       -0.29684  0.95493 0.0260  0.689    
PO4          -0.54798  0.83649 0.0210  0.749    
Chla         -0.59554  0.80332 0.1339  0.144    
Chao1         0.21598 -0.97640 0.0645  0.375    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
Permutation: free
Number of permutations: 999



