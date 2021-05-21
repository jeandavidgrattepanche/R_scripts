setwd('/home/tuk61790/HPCr/')
library('ape')
library('phytools')
library('phyloseqCompanion')
library("ggplot2")
library("plyr")
library(stringr)
#library(factoextra)
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

#read distribution OTUs and Samples
readnumb = data.frame(readn = sort(taxa_sums(physeq), TRUE), sorted = 1:ntaxa(physeq), type = "OTUS")
readnumb = rbind(readnumb, data.frame(readn = sort(sample_sums(physeq), TRUE), sorted = 1:nsamples(physeq), type = "Samples")
titlenb = "read number distribution"
p = ggplot(readnumb, aes(x= sorted, y = readn)) + geom_bar(stat="identity")
p + ggtitle(titlenb) + scale_y_log10() + facet_wrap(~type, 1, scales = "free")


#check if empty sample - can be useful after next step
any(sample_sums(physeq) == 0)
physeq = prune_samples(sample_sums(physeq) >0, physeq)

#remove singleton ! be careful if samples are already rarefy
physeq.bc = physeq # backup in case 
physeq = prune_taxa(taxa_sums(physeq) > 1, physeq)

#subsample to reduce the number of OTU (crash R on the server)
#physeq.mil <- get_top_taxa(physeq, 1500000, relative = FALSE, discard_other= TRUE, other_label="Other")

# add a variable to reduce number of samples : merge station to group
sample_data(physeq)$GLS <- mapply(paste0, as.character(get_variable(physeq,"Group")), collaspe = '_', as.character(get_variable(physeq,"layer")),collaspe = '_', as.character(get_variable(physeq,"size")))
physeq.group <- merge_samples(physeq, "GLS")

#sample_data(physeq.group)$GLS <- levels(sample_data(physeq)$GLS)
renam <- str_split_fixed(sample_names(physeq.group), "_", 3)
sample_data(physeq.group)$G <- mapply(paste0, as.character(renam[,1]))
sample_data(physeq.group)$L <- mapply(paste0, as.character(renam[,2]))
sample_data(physeq.group)$S <- mapply(paste0, as.character(renam[,3]))


#convert abundance in percentage:
physeq.group = transform_sample_counts(physeq.group, function(x) 100 * x/sum(x))

#convert in presence/absence data
#mytable_pa = otu_table(cbind(data.all[c(0)],((data.all[c(16:111)]>0)*1L)), taxa_are_rows=TRUE,errorIfNULL=TRUE)
#physeq_pa <- phyloseq(mytable_pa, envdata, TAX)

# testing: 
# physeq.spec <- taxa_merge(physeq, taxa_names(physeq)[range of OTU or a list of OTU name]) not working for my purpose which is decreasing the data 
print("select taxa")
physeq.sp <- tax_glom(physeq.group, taxrank="Bsp", NArm = FALSE) 
# not sure about pruning OTU without Taxo
#pico <- subset_samples(physeq, size == "pico")
#nano <- subset_samples(physeq, size == "nano")
#micro <- subset_samples(physeq, size == "micro")
print('data ready to compute')

datatorun <- physeq.sp

#pb = plot_bar(datatorun,"Group", fill="Btaxo_rank4", facet_grid=layer~size)
pb = plot_bar(datatorun, "G",fill="Btaxo_rank4", facet_grid=L~S)

print('first plot running')
jpeg(filename='barplot_group_sp15m.jpeg', width=60, height=15, unit="cm", res=100)
plot <- pb + geom_bar(aes(color=Btaxo_rank4, fill=Btaxo_rank4), stat="identity", position="stack")
print(plot)
dev.off()
print('second plot legend done and saved')

print('done!')

#PCoA with unifrac !! need a tree !!
#dUNIFRAC = UniFrac(physeq, weighted=TRUE, normalized=FALSE, fast=TRUE)
#PCoA <- ordinate(physeq, "PCoA", distance = dUNIFRAC)
#ordplotP <- plot_ordination(physeq, PCoA, color="Lat", shape="size",label="sample", title="Allsize")
#ordplotPplus <- ordplotP + geom_point(size=3) + stat_ellipse(geom = "polygon", type="norm", alpha=0.4, aes(fill=size))
#quartz()
#plot(ordplotPplus)

#create PCoA plot for each size fraction, using 3 classical diversity indices. size can be replace by another varaible
sizes <- c("micro","nano","pico")
indices <- c("jaccard","bray","chao") #"bray","cao","raup","chisq"
list <- list()
i = 0
for(ind in indices){
	for(sizedata in sizes){
		i = i + 1
		named <- paste(ind,sizedata,sep = "-")
		print(named)
		PCoAm <- ordinate(subset_samples(datatorun, S == sizedata), "PCoA", distance = ind)
		ordplotPm <- plot_ordination(subset_samples(datatorun, S == sizedata), PCoAm, color="G", shape="L",label="G",title=named)
		ordplotPmplus <- ordplotPm + geom_point(size=2)
		assign(paste("p",i,sep=""), ordplotPm)
		list[[length(list)+1]] <- paste("p",i,sep="")
		plot(get(paste("p",i,sep="")))
	}
}
pdf('PCoAs_Bsp_test15m.pdf')
multiplot(p1 + theme(legend.position="none"),p2+ theme(legend.position="none"),p3+ theme(legend.position="none"),p4+ theme(legend.position="none"),p5+ theme(legend.position="none"),p6+ theme(legend.position="none"),p7+ theme(legend.position="none"),p8+ theme(legend.position="none"),p9+ theme(legend.position="none"),cols=3)
dev.off()

# same as above for presence absence data
#indices <- c("bray","jaccard","raup")#,"chisq"
#list <- list()
#i = 0
#for(ind in indices){
#	for(sizedata in sizes){
#		i = i + 1
#		named <- paste(ind,sizedata,sep = "-")
#		print(named)
#		PCoAm <- ordinate(subset_samples(physeq_pa, size == sizedata), "PCoA", distance = ind)
#		ordplotPm <- plot_ordination(subset_samples(physeq_pa, size == sizedata), PCoAm, color="Group", shape="layer",label="station",title=named)
#		ordplotPmplus <- ordplotPm + geom_point(size=3)
#		assign(paste("p",i,sep=""), ordplotPm)
#		list[[length(list)+1]] <- paste("p",i,sep="")
#		plot(get(paste("p",i,sep="")))
#	}
#}
#pdf('PCoAs_pa.pdf')
#multiplot(p1,p2,p3,p4,p5,p6,p7,p8,p9,cols=3)
#dev.off()

#NMDS plot ... some issue with creatign the merge plot
#list <- list()
#j = 0
#for(sizedata in sizes){
#	j = j + 1
#	print(sizedata)
#	ordMDS <- metaMDS(t(otu_table(subset_samples(physeq, size == sizedata))))
#	ordMDS.fit <- envfit(ordMDS ~ Chla + Babun + bottom + ice + airT + ZML_T + ZML_TS + ZML_TSP + Ze + ZCM + depth + water_temp + conductivity + salinity + oxygen + O_saturation + fluorescence + beam_trans + PAR + surf_PAR + Bprod, data=sample.data.table(subset_samples(physeq, size == sizedata)), perm=999, na.rm = TRUE)
#	part123 <- c(plot(ordMDS, dis="site") , plot(ordMDS.fit) , text(ordMDS, display="sites"))
#	qa <- plot(ordMDS, dis="site")
#	qb <- plot(ordMDS.fit)
#	qc <- text(ordMDS, display="sites")
#	assign(paste("q",j,"a",sep=""), qa)
#	assign(paste("q",j,"b",sep=""), qb)
#	assign(paste("q",j,"c",sep=""), qc)
#}
#pdf('NMDS_2.pdf')
#multiplot(q1a,q1b,q1c,q2a,q2b,q2c,q3a,q3b,q3c, cols=3) #,q4,q5,q6,q7,q8,q9,cols=3)
#dev.off()