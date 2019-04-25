# DAWG R Code for Spring 2019 #
# source("http://bioconductor.org/biocLite.R")
# install.packages("BiocManager")
# BiocManager::install("ShortRead")
# BiocManager::install("devtools")
# BiocManager::install("dada2", version = "3.8")

library(dada2)
#packageVersion("dada2")

##### 1. Set working directory ####

setwd("~/Downloads/Berts_ITS_Data/demultiplexed_fastq")

#Tell R to pull files from the path
path<- setwd("~/Downloads/Berts_ITS_Data/demultiplexed_fastq")

list.files(path)


#### 2. Sort the forward and reverse reads #####
# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern=".fastq", full.names = TRUE))
#You should see in your global env, fnFs chr[1:30] "and the path directory"

# Extract sample names
sample.names <- sapply(strsplit(basename(fnFs), ".fastq"), `[`, 1)


#### 3. Examine the quality profiles ####

## forward reads ##
quartz()
plotQualityProfile(fnFs[1:6])

# green line shows median quality score at each position
# red line shows how many reads extend to at least that position
# 

#### 4. Assign Filtering samples ####
# Assign the filenames for the filtered fastq.gz files

# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_filt.fastq.gz"))


# the commands included in this filter and trim are dependent on your specific data- 
# especially the truncLen- because you want to trim your seqs based on the quality profiles you generated- quality scores of 30 + are good
# so you should choose the trunclen to reflect the quality scores of your data 

out.2<- filterAndTrim(fnFs, filtFs, truncLen=400, 
                      minLen = 200, maxN=0, maxEE=5, truncQ = 2, 
                      multithread=TRUE, compress=TRUE)
# set multithread to FALSE if using windows
# this step takes awhile to perform
# setting trunQ=2 will remove sequences that have at least one read with a high probability of erroneous base assignment (>63%)
# minLen will remove sequences that are less than the specific value
# maxN = 0 will remove all sequences that have ambiguous base cales
# truncLen = 400 will cut off each sequence at base pair 400 and needs to be changed to account for every sequencing run
# look at quality profiles to determine what minLen and truncLen should be set to
# maxEE = 5 sets the maximum number of "expected errors" allowed in a read to be 5- in Illumina set to 2,2 as default but we are working with longer reads
# Compress = TRUE stores output as .gz files to save space
out.2
head(out.2)

#### 5. Learn the Error Rates ####

# learn error rates
# error models are run specific
errf<- learnErrors(filtFs, multithread=TRUE) #set multithread=FALSE if on windows


plotErrors(errf, nominalQ = TRUE)
# shows error rates of transitions from one nucleotide to another 
# black line shows the estimated error rates after convergence of the machine-learning algorithm
# red line shows the error rates expected under the nominal definition of the Q-score
# want the samples (dots) to track well with the black line


#### 6. Dereplicate the sequences ####

derepFs <- derepFastq(filtFs, verbose=TRUE)
# If there are 10 sequences that are the same, derepFs will contain 1 of the sequences
# speeds up downstream processing


names(derepFs) <- sample.names

##### 7. Call ESVs ####

# incorporates consensus quality profiles & abundances of unique sequences
# determines if the sequence is more likely to be of biological origin or spurious


dadaFs<-dada(derepFs, err=errf, HOMOPOLYMER_GAP_PENALTY=-1, BAND_SIZE=32, USE_QUALS=TRUE, multithread=TRUE)
# set multithread=FALSE if using windows
# homopolymer gap penalty of -1,  band size = 32, and setting USE_QUALS=TRUE is suggested by the dada2 team
# The cost of gaps in homopolymer regions (>=3 repeated bases). Default is NULL, which causes homopolymer gaps to be treated as normal gaps.
# The default value of BAND_SIZE is 16. If DADA is applied to sequencing technologies with high rates of indels, such as 454 sequencing, the BAND_SIZE parameter should be increased.
# USE_QUALS: If TRUE, the dada(...) error model takes into account the consensus quality score of the dereplicated unique sequences. If FALSE, quality scores are ignored


dadaFs[[2]]
# output shows X ESVs inferred from XXX input seqs from the first sample
# OMEGA-A is default parameter sets the level of “statistical evidence” (think p-value) required for inferences of a new ESV
# OMEGA-C The threshold at which unique sequences inferred to contain errors are corrected in the final output. The probability that each unique sequence is generated at its observed abundance from the center of its final partition is evaluated, and compared to OMEGA_C. If that probability is >= OMEGA_C, it is "corrected", i.e. replaced by the partition center sequence
# BandSize was set by the user in the filterAndTrim step


seqtab<-makeSequenceTable(dadaFs)
dim(seqtab)
#30 samples with 696 ESVs

#### 8. Remove Chimeric Sequences ####
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#147 seqs were removed because they were chimeras
dim(seqtab.nochim)
# 30 samples with 549 ESVs
sum(seqtab.nochim)/sum(seqtab)
# ~95% of the sequences in the dataset are real biological sequences


getN <- function(x) sum(getUniques(x))

track <- cbind(out.2, sapply(dadaFs, getN), rowSums(seqtab.nochim))

colnames(track) <- c("input", "filtered", "denoisedF", "nonchim")
rownames(track) <- sample.names
head(track)
# Check this to see how many of your samples were dropped and if there are any steps where a majority of your samples were dropped


##### 9. Assign Taxonomy with Unite ####
taxa <- assignTaxonomy(seqtab.nochim, "/Users/admin/Downloads/Unite/sh_general_release_dynamic_01.12.2017.fasta", multithread=TRUE)
taxa #check to make sure taxonomy looks appropriate





# Session 2 start here #
#### 10. Format data for analysis using Phyloseq object ####
source('http://bioconductor.org/biocLite.R')
biocLite('phyloseq')

# OR

install.packages("BiocManager")
BiocManager::install("phyloseq")
library(phyloseq)
install.packages("tibble")
library(tibble)
# three files for phyloseq, esv table (otu table), taxonomy table, metadata

#import the metadata file
Map<-import_qiime_sample_data("/Users/admin/Downloads/Berts_ITS_Data/BertsMetadata.txt")


# making phyloseq requires common linker to merge the files into one object
# take a look at each file individual
seqtab.nochim # you'll notice that the column names are sequences, not ESVs
taxa # you'll notice that the row names are sequences
Map # row names are sample names

# we need to change the sequences to ESVs  for both the seqtab.nochim (esv table) and taxa table

# A. Duplicate your final ESV table #
seqtab.nochim.mod<-seqtab.nochim


## B. Replace column names for ESV table with ESV numbers 
ncol(seqtab.nochim.mod)
#549 ESVs
fin.ESV.names<-paste("ESV",sep = "",rep(1:ncol(seqtab.nochim.mod)))

colnames(seqtab.nochim.mod)<-fin.ESV.names
# open up seqtab.nochim.mod, sequences are not replaced by ESVs

esv.table<-otu_table(seqtab.nochim.mod, taxa_are_rows=FALSE)
# this formats the esv table into an "otu_table", a requirement for phyloseq objects


## C. Add ESV names to Taxonomy table and then write to file
# We will add the ESVs from the ESV table to the taxonomy table 
# then move the column of ESVs into the rows

fin.ESV.names<-colnames(esv.table) # puts the ESVs into a character string

taxa.mod<-cbind(taxa, fin.ESV.names) # adds the character string of ESVs to the taxa table
# open up taxa.mod

taxa.df<-rownames_to_column(as.data.frame(taxa.mod), var="sequence") 
# puts sequences into a column, because we want to keep them 

taxa.evs.final<-column_to_rownames(taxa.df, var="fin.ESV.names")
# this moves the column "fin.ESV.names" to the rows

# now we have to format the taxa table into a tax_table
tax.evs.final.table<-tax_table(taxa.evs.final)
# we get an error saying ... coerce to matrix manually

tax.evs.final.mat<-as.matrix(taxa.evs.final)
tax.evs.final.table<-tax_table(tax.evs.final.mat)
head(taxa_names(tax.evs.final.table))

# D. Now we will build the phyloseq object
ps <- phyloseq(esv.table, tax.evs.final.table, Map)

#Now, let's make sure that the correct information is included in our phyloseq object


ps@otu_table #Should include ESV info
ps@tax_table #Should include taxonomic info, but we get a NULL output
ps@sam_data




#Now export the taxonomy table for records
write.table(tax.evs.final.table, file = "ESV_Taxonomy.txt")

# Write ESV table to file
write.table(seqtab.nochim.mod, file="ESV_Abund_Table.txt")


#### Remove ESVs not assigned to the Fungal Kingdom ####
badTaxa = c("ESV15", "ESV169", "ESV277", "ESV334", "ESV349", "ESV364", "ESV464")

allTaxa = taxa_names(ps)
allTaxa <- allTaxa[!(allTaxa %in% badTaxa)]
ex1 = prune_taxa(allTaxa, ps)
# new phyloseq object with just the taxa you kept.
ex1

#### 11. General processing and analysis using silva phyloseq object####

#####Now we can begin using the phyloseq object for analyses, but first, we need to rarefy the dataset using the following commands

sample_sums(ex1)
# the sample with the smallest sequence count is E28 with 6, so let's remove samples with fewer than 2,000 seqs
ex.pruned <- prune_samples(sample_sums(ex1)>=1400, ex1)
sample_sums(ex.pruned)

# 1 sample removed (E28)

# now we will rarefy the phyloseq object
set.seed(500)
ps.rare<-rarefy_even_depth(ex.pruned)
# 56 OTUs removed


sample_sums(ps.rare)
# Now 1,413 sequences are in each of the 29 samples

#Relative abundance calculation


ps.perc<-transform_sample_counts(ps.rare, function(x) x / sum(x)) #This calculates the percentage of each read in each sample


#### 12. Alpha-Diversity ####

# Make a data frame from the otu_data
# packages that you will need for this:
BiocManager::install("microbiome") 
#source('http://bioconductor.org/biocLite.R')
#biocLite('microbiome')
library(microbiome)
install.packages("lattice")
library(lattice)
install.packages("FSA")
library(FSA)

meta.df <- data.frame(sample_data(ps.perc))

alpha.div<-estimate_richness(ps.rare, measures=c("Shannon", "Simpson", "Observed"))
alpha.div


histogram(alpha.div$Shannon)
histogram(alpha.div$Simpson)

boxplot(alpha.div$Shannon ~ meta.df$Treatment)
# due to uneven sampling size and data being non-normal for shannon and simpson, use kruskal-wallis test
shannon.model<-kruskal.test(alpha.div$Shannon ~ meta.df$Treatment)
shannon.model
shan.pairwise<- dunnTest(alpha.div$Shannon ~ meta.df$Treatment, method = "bh")
shan.pairwise
capture.output(shan.pairwise, file = "ShannonPairwise_BHpAdjust.txt")

simpson.model<-kruskal.test(alpha.div$Simpson ~ meta.df$Treatment)
simpson.model

observed.model<-kruskal.test(alpha.div$Observed ~ meta.df$Treatment)
observed.model
observed.pairwise<- dunnTest(alpha.div$Observed ~ meta.df$Treatment, method = "bh")
observed.pairwise
capture.output(observed.pairwise, file = "ObservedPairwise_BHpAdjust.txt")



# now calculate evenness and perform test
alpha.div$evenness<-evenness(ps.rare, index = "all", zeroes = TRUE)

pileu<- kruskal.test(alpha.div$evenness$pielou ~ meta.df$Treatment)
pileu

jpeg(file="AlphaDivBoxPlot_ITS_Berts.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
plot_anova_diversity(ps.rare, method = c("richness", "shannon") , "Treatment",
                     pValueCutoff = 0.01)
dev.off()
#### 13. Beta-Diversity Estimates ####

# packages that you will need for this:
install.packages("ggplot2")
library(ggplot2)
install.packages("vegan")
library(vegan)
BiocManager::install("microbiomeSeq")
library(microbiomeSeq)


plot_heatmap(ps.rare, "PCoA", "bray", "Treatment")
# samples are grouped by their bray-curtis similarities, but using OTUs (ESVs) is somewhat useless


# use this to look at phylum level
ps.phylum.rare<- taxa_level(ps.rare, "Phylum") #this calculates the abunddances of each phylum across every sample into a new phyloseq object

jpeg(file="PhylaLevelHeatmap_ITS_Berts.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
plot_heatmap(ps.phylum.rare, "PCoA", "bray", "Treatment")
dev.off()

# Now look at Classes #

ps.class<- taxa_level(ps.rare, "Class") #this calculates the abunddances of each phylum across every sample into a new phyloseq object

jpeg(file="ClassLevelHeatmap_ITS_Berts.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
plot_heatmap(ps.class, "PCoA", "bray", "Treatment")
dev.off()

# Now look at Family #

ps.family <- taxa_level(ps.rare, "Family") #this calculates the abunddances of each phylum across every sample into a new phyloseq object

jpeg(file="FamilyLevelHeatmap_ITS_Berts.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
plot_heatmap(ps.family, "PCoA", "bray", "Treatment")
dev.off()

# Now look at Order #

ps.order <- taxa_level(ps.rare, "Order") #this calculates the abunddances of each phylum across every sample into a new phyloseq object

jpeg(file="OrderLevelHeatmap_ITS_Berts.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
plot_heatmap(ps.order, "PCoA", "bray", "Treatment")
dev.off()
# Use this to look at genera # 
ps.genera.perc<- taxa_level(ps.rare, "Genus")
jpeg(file="GenusLevelHeatmap_ITS_Berts.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
plot_heatmap(ps.genera.perc, "PCoA", "bray", "Treatment", na.value = "grey",  low= "#66CCFF", high = "#000033" )
dev.off()

plot_heatmap(ps.genera.perc, "PCoA", "jaccard", "Treatment", na.value = "grey",  low= "#66CCFF", high = "#000033" )

# Rhizopogon villosulus, Rhizopogon luteolus, Rhizopogon amylopogon, Rhizopogon fulvigleba, Pisolithus tinctorius, Scleroderma cepa, Scleroderma citrinum were the starting inoculom


# A. Start with a heatmap of bray

# B. Calculate Bray-Curtis
bc_rt_highcover_ns_abun <- phyloseq::distance(otu_table(ps.perc), "bray")

# C. Calculate Jaccard
jc_rt_highcover_ns_abun <- phyloseq::distance(otu_table(ps.perc), "jaccard")


# D. Set ggtheme
ggtheme <- theme(axis.title = element_text(colour="black",family = "Helvetica",
                                           size = rel(1.5)), 
                 axis.text = element_text(family = "Helvetica",colour = "black",
                                          size = rel(1)), 
                 axis.line = element_line(size = 0.5,colour = "black"), 
                 axis.ticks = element_blank(),
                 panel.grid.major = element_line(colour="grey",size = rel(0.25)), 
                 panel.grid.minor = element_blank(), 
                 panel.background = element_blank(),
                 plot.title = element_text(colour = "black", face = "bold",
                                           size = rel(2),family = "Helvetica",hjust = 0.5,
                                           geom_point(size = 15)))


## E. Perform PERMANOVA to check significance of categories: BC

# Use the vegan function adonis to conduct PERMANOVA. PERMANOVA is a nonparametric multivariate analysis, meant to test differences between groups the same way an ANOVA would do, but considering many variables, comparing centroids and dispersion rather than means, and using permutations to avoid possible biases. Rather than using our raw data, PERMANOVA takes distance matrices and tells you the the degree of effect on the distribution (beta diversity) for a grouping variable. Here, we use the Bray-Curtis matrix we generated in the previous step and select the variables we want to test.
bc.adon.trt <- adonis(bc_rt_highcover_ns_abun ~ Treatment, data = meta.df)
bc.adon.site <- adonis(bc_rt_highcover_ns_abun ~ Site, data = meta.df)

bc.adon.trt
bc.adon.site

capture.output(bc.adon.trt, bc.adon.site, file = "PERMANOVA_Results_Bray_Curtis_By_Site_Or_treatment.txt")
## All variables are found to have a significant effect on the distribution (i.e. grouping) of these microorganisms.


## F. Perform PERMANOVA to check significance of categories: JC

# PERMANOVA analyses using the Jaccard matrix.
adonis(jc_rt_highcover_ns_abun ~ Treatment, data = meta.df)
adonis(jc_rt_highcover_ns_abun ~ Site, data = meta.df)


### G. Posthoc: determining group differences in all data: (~ Treatment + Site): BC

#packages that you will need for this:
install.packages("RVAideMemoire")
library(RVAideMemoire)

set.seed(100)
# Pariwise Permutation MANOVAs
pair.perm.trt <- pairwise.perm.manova(bc_rt_highcover_ns_abun, meta.df$Treatment ,nperm=1000)
pair.perm.site <- pairwise.perm.manova(bc_rt_highcover_ns_abun, meta.df$Site ,nperm=1000)

capture.output(pair.perm.trt, pair.perm.site, file = "PairwisePERMANOVA_Results_Bray_Curtis_By_Site_Or_treatment.txt")

### H. Posthoc: determining group differences in all data: (~ Treatment + Site): JC
set.seed(1)
# Pariwise Permutation MANOVAs
pairwise.perm.manova(jc_rt_highcover_ns_abun, meta.df$Treatment,nperm=1000)
pairwise.perm.manova(jc_rt_highcover_ns_abun, meta.df$Site,nperm=1000)



### I. BRAY-CURTIS PCoA plot

#Calculate the PCoA on BC distance matrix: Treatment
rt.pcoa = ordinate(ps.perc, method="PCoA", distance=bc_rt_highcover_ns_abun)
plot_scree(rt.pcoa, "Screen plot")
# PCoA ordination
pcoa<-0
pcoa <- plot_ordination(ps.perc, rt.pcoa, "samples", color="Treatment")  + ggtheme 
quartz()
pcoa

# Create output folder
dir.create("~/Downloads/Berts_ITS_Data/plots/")

# Save figure as PDF
pdf("~/Downloads/Berts_ITS_Data/plots/PCoA_ps.perc_BC_Treatment.pdf",width=9,height=5)
pcoa
dev.off()

## J. JACCARD PCoA plot

#Calculate the PCoA on JC distance matrix: Treatment
rt.pcoa = ordinate(ps.perc, method="PCoA", distance=jc_rt_highcover_ns_abun)
plot_scree(rt.pcoa, "Screen plot")
# PCoA ordination
pcoa<-0
pcoa <- plot_ordination(ps.perc, rt.pcoa, "samples", color="Treatment") + ggtheme
pcoa


# Save figure as PDF
pdf("output/plots/PCoA_ps.perc_JC_Treatment.pdf",width=9,height=5)
pcoa
dev.off()

#### Session 3 Differential Abundance 

# Download phyloseq file
rt_raw<-load("Berts_ITS_RData4Analyses.RData")


## BiocAnyway(package name)##
## If your R version is older than 3.5, the function will use biocLite(), unless it will use BiocManager ##

BiocAnyway <- function(name_of_package){
  if (substr(version$version.strin, 11, 15) >= 3.5){
    install.packages("BiocManager")
    BiocManager::install(name_of_package)}
  else {source('http://bioconductor.org/biocLite.R')
    biocLite(name_of_package)}
}

#### Principle Component Analysis of Metals Data ####

library(devtools)

install_github("kassambara/factoextra")

library("factoextra")

library("FactoMineR")

install.packages("cowplot")
library(cowplot)

# pca requires numerical factors, so we will retrieve the metadata from the phyloseq object and then remove the categorical variables
metadata.t<- data.frame(sample_data(ps.rare)) 
metadata.t[,1:5]<-NULL # we can get rid of the categorical variables with this


res.pca <- PCA(metadata.t, graph = FALSE)
quartz()
plot(res.pca) #provides a pca plot of the metals data 


fviz_screeplot (res.pca, ncp=10) # % variance explained by each PC
fviz_pca_var(res.pca) # Variables factor map
var.PCA1.plot <-fviz_contrib(res.pca, choice = "var", axes = 1)  # set axes=1 to consider only PC1
var.PCA2.plot <- fviz_contrib(res.pca, choice = "var", axes = 2)  # set axes=1 to consider only PC1
pca.plot <- fviz_pca_biplot(res.pca)


log_plot <- plot_grid(var.PCA1.plot,var.PCA2.plot, nrow =3)

jpeg(file="PCAplots_ITS_Berts.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
log_plot
dev.off()


pca.biplot <- fviz_pca_biplot(res.pca, repel = TRUE, label = "var",
                              col.ind = meta.df$Treatment)
quartz()
pca.biplot

mtcars.pca <- prcomp(metadata.t, center = TRUE, scale = TRUE)
#this wont work because the "Hg" column has all 0's
# let's get rid of Hg first
metadata.t$Hg<-NULL

mtcars.pca <- prcomp(metadata.t, center = TRUE, scale = TRUE)
# now this works perfectly
pca.data <- as.data.frame(mtcars.pca$x)
meta.df$Site<- as.factor(meta.df$Site)

percentage <- round(mtcars.pca$sdev / sum(mtcars.pca$sdev) * 100, 2)
percentage <- paste( colnames(pca.data), "(", paste( as.character(percentage), "%", ")", sep="") )
pca.data$group <- sapply( strsplit(as.character(row.names(metadata.t)), "_"), "[[", 1 )
theme<-theme(panel.background = element_blank(),panel.border=element_rect(fill=NA),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),strip.background=element_blank(),axis.text.x=element_text(colour="black"),axis.text.y=element_text(colour="black"),axis.ticks=element_line(colour="black"),plot.margin=unit(c(1,1,1,1),"line"))
pca.plot<-ggplot(pca.data,aes(x=PC1,y=PC2, color = meta.df$Site)) + scale_fill_manual(breaks = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10")) + theme_classic() +geom_point()
pca.plot <-pca.plot + geom_point(aes(color = meta.df$Site, shape =meta.df$Treatment, size = 5)) + xlab(percentage[1]) + ylab(percentage[2]) + guides(colour = guide_legend(override.aes = list(size=6))) + guides(shape = guide_legend(override.aes = list(size = 6))) + theme(legend.text = element_text(size = 12))


quartz()
pca.plot

jpeg(file="PCAplot_Treatment.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
pca.plot
dev.off()

env.data.adon <- adonis(metadata.t ~ meta.df$Treatment, method = "euclidean")

cc.env <- vegdist(metadata.t, "manhattan")
cc.env

perm.env <- pairwise.perm.manova(cc.env, meta.df$Treatment, nperm=100)

capture.output(env.data.adon, perm.env, file = "MetalsData-Permanovas.txt")

mod.perm <- betadisper(cc.env, meta.df$Treatment)
mod.perm

## Perform test
anova(mod.perm)

## Permutation test for F
permutest(mod, pairwise = TRUE)


install.packages("PerformanceAnalytics")
library(PerformanceAnalytics)

quartz()
chart.Correlation(metadata.t, 
                  method="spearman",
                  histogram=TRUE,
                  pch=16)

install.packages("psych")
library(psych)

corr.test(metadata.t, 
          use = "pairwise",
          method="spearman",
          adjust="bonferroni",      # Can adjust p-values; see ?p.adjust for options
          alpha=.05)

##### look at correlations between taxa and env data #####

taxa_level <- function(physeq,which_level){
  #enforce orientation
  if(taxa_are_rows(physeq)){
    physeq <- t(physeq)
  }
  OTU <- otu_table(physeq)
  SAM <- sample_data(physeq)
  OTU_taxonomy <- tax_table(physeq)
  new_abund_table<-NULL
  if(which_level=="Otus"){
    OTU_tree <- phy_tree(physeq)
    new_abund_table<-OTU
  } else {
    list<-na.omit(unique(OTU_taxonomy[,which_level]))
    new_abund_table<-NULL
    for(i in list){
      rt <- na.omit(rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i])
      tmp<-data.frame(rowSums(OTU[,rt]))
      if(i==""){colnames(tmp)<-c("__Unknowns__")} else {colnames(tmp)<-paste("",i,sep="")}
      if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
    }
  }
  OTU<-as.data.frame(as(new_abund_table,"matrix"))
  #Convert the data to phyloseq format
  OTU = otu_table(as.matrix(OTU), taxa_are_rows = FALSE)
  TAX = tax_table(as.matrix(OTU_taxonomy))
  SAM = sample_data(SAM)
  #reconstruct the phyloseq object
  physeq<-NULL
  if(which_level=="Otus"){
    physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,midpoint(OTU_tree))
  } else {
    physeq<-merge_phyloseq(phyloseq(OTU),SAM)
  }
  return(physeq)
}

install_github("umerijaz/microbiomeSeq")  # Install the package
library(microbiomeSeq)
library(phyloseq)
library(ggplot2)

physeq <- taxa_level(ps.rare, "Phylum")
env.taxa.cor <- taxa.env.correlation(physeq, grouping_column = "Treatment", method = "pearson", 
                                     pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                     select.variables = NULL)

p <- plot_taxa_env(env.taxa.cor)
print(p)


# Correlations: Class
env.taxa.cor.class <- taxa.env.correlation(ps.class, grouping_column = "Treatment", method = "pearson", 
                                         pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                         select.variables = NULL)

p.5 <- plot_taxa_env(env.taxa.cor.class)

jpeg(file="ClassLevelMetalsCorrelations_ITS_Berts.jpeg", width=20, height=10, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
p.5
dev.off()
# Correlations: Order
physeq_ord<- taxa_level(ps.rare, "Order")
physeq_ord@sam_data$Site <- as.factor(physeq_ord@sam_data$Site)
env.taxa.cor.ord <- taxa.env.correlation(ps.order, grouping_column = "Treatment", method = "pearson", 
                                         pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                         select.variables = NULL)

p.4 <- plot_taxa_env(env.taxa.cor.ord)

jpeg(file="OrderLevelMetalsCorrelations_ITS_Berts.jpeg", width=20, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
p.4
dev.off()

# Correlations: Family
physeq_fam <- taxa_level(ps.rare, "Family")
physeq_fam@sam_data$Site <- as.factor(physeq_fam@sam_data$Site)

env.taxa.cor.fam <- taxa.env.correlation(ps.family, grouping_column = "Treatment", method = "pearson", 
                                         pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                         select.variables = NULL)

p.2 <- plot_taxa_env(env.taxa.cor.fam)
print(p.2)

jpeg(file="FamilyLevelMetalsCorrelations_ITS_Berts.jpeg", width=30, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
p.2
dev.off()

# Correlations: Genus
physeq_gen<- taxa_level(ps.rare, "Genus")
physeq_gen@sam_data$Site <- as.factor(physeq_gen@sam_data$Site)
env.taxa.cor.gen <- taxa.env.correlation(ps.genus, grouping_column = "Treatment", method = "pearson", 
                                         pvalue.threshold = 0.05, padjust.method = "BH", adjustment = 5, num.taxa = 50, 
                                         select.variables = NULL)

p.3 <- plot_taxa_env(env.taxa.cor.gen)

jpeg(file="GenusLevelMetalsCorrelations_ITS_Berts.jpeg", width=20, height=10, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
p.3
dev.off()



#### ANOVA for all metals, grouped by Treatment ####
plot_anova_env(ps.rare, grouping_column = "Treatment", pValueCutoff = 0.05, 
               select.variables = c("Al", "As", "B"))
plot_anova_env(ps.rare, grouping_column = "Treatment", pValueCutoff = 0.05, 
               select.variables = c("Ba", "Be", "Ca"))
plot_anova_env(ps.rare, grouping_column = "Treatment", pValueCutoff = 0.05, 
               select.variables = c("Cd", "Co", "Cr"))
plot_anova_env(ps.rare, grouping_column = "Treatment", pValueCutoff = 0.05, 
               select.variables = c("Cu", "Fe"))

# What happens when you try to plot them all together
jpeg(file="AnovaDiversity_A_MetalsData_ITS_Berts.jpeg", width=14, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
plot_anova_env(ps.rare, grouping_column = "Treatment", pValueCutoff = 0.05, 
               select.variables = c("Al", "As", "B", "Ba", "Be", "Ca", "Cd", "Co", "Cr", "Cu", "Fe", "K", "Li"))
dev.off()

jpeg(file="AnovaDiversity_B_MetalsData_ITS_Berts.jpeg", width=14, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
plot_anova_env(physeq, grouping_column = "Treatment", pValueCutoff = 0.05, 
               select.variables = c("Mg", "Mn", "Mo", "Na", "Ni", "P", "Pb", "S", "Se", "Sr", "Ti", "V", "Zn"))
dev.off()

#### 14.Differential Abundance Analysis ####
# packages you will need include
BiocAnyway(DESeq2)


# Deseq uses a generalized linear model and calculate the log2 fold changes, as well as adjusted p-values, for taxanomic groups based on an input variable
# Input variables can be both categorical or continuous
# We will start with continous variables


# Look for Genera that are differentially abundant across Fe concentrations
ps.genus<- taxa_level(ps.rare, which_level = "Genus")

gen.fe <- phyloseq_to_deseq2(ps.genus, ~ Fe) # converts the phyloseq object into a Deseq object
gen.fe.counts <- counts(gen.fe) # contains sequence counts for each genera
geoMeans = apply(gen.fe.counts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds.gen.fe <-estimateSizeFactors(gen.fe, geoMeans=geoMeans) 
dds.gen.fe <- DESeq(dds.gen.fe, test="Wald", fitType="parametric")

res.gen.fe = results(dds.gen.fe, cooksCutoff = FALSE) # extracts results table from DESeq analysis
alpha=0.05 # set alpha value 
sigtab.gen.fe = res.gen.fe[which(res.gen.fe$padj < alpha), ] # deseq output that excludes taxa that have p-values greater than alpha
sigtab.gen.fe #this file includes only the genera that were significantly abundant across the treatment specified--
# need to transform this into a vector that only includes the taxa differentially abundant to pull those out of the phyloseq object
capture.output(sigtab.gen.fe, file = "Genera_Level_DESeqResults_ByFe.txt")


esv.fe<- rownames_to_column(as.data.frame(sigtab.gen.fe), "Genera")
esv.fe[,2:17]<-NULL
esv.fe.vec<-as.vector(esv.fe$Genera)
esv.fe.ps<-prune_taxa(esv.fe.vec, ps.genus) # this removes the significant genera out of the main phyloseq object and makes a new phyloseq object with only those taxa

diff.gen.fe = data.frame(otu_table(esv.fe.ps)) # exports otu table from phyloseq object for downstream use
diff.gen.meta = data.frame(sample_data(esv.fe.ps)) # exports sample data from phyloseq object for downstream use
diff.gen.fe$Fe <- diff.gen.meta$Fe # places Fe from the sample data into the otu table file
diff.gen.fe$Treatment <- diff.gen.meta$Treatment

lm_eqn = function(m) {
  # Displays regression line equation and R^2 value on plot
  # Usage:
  # p + annotate("text", x=25, y=300, label=lm_eqn(lm(y ~ x, df)), parse=TRUE)
  
  l <- list(a = format(coef(m)[1], digits = 2),
            b = format(abs(coef(m)[2]), digits = 2),
            r2 = format(summary(m)$r.squared, digits = 3));
  
  if (coef(m)[2] >= 0)  {
    eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2,l)
  } else {
    eq <- substitute(italic(y) == a - b %.% italic(x)*","~~italic(r)^2~"="~r2,l)    
  }
  
  as.character(as.expression(eq));                 
} # function that calculates a regression and rsquared of the x-y values

jpeg(file="Schleroderma&Fe_regression.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
ggplot(diff.gen.fe, aes(x= Fe, y = g__Scleroderma)) + geom_point(aes(col=Treatment), size = 3) + geom_smooth(method = lm, se=TRUE, fullrange=TRUE) + labs(y = "Scleroderma Abundance", x = "Fe Concentration") + geom_text(aes(x = 11000, y = 2000, label = lm_eqn(lm(g__Scleroderma ~ Fe, diff.gen.fe))), size = 3, parse=TRUE) + theme_classic()
dev.off()


# now we can do the same for Order level
ps.order<- taxa_level(ps.rare, which_level = "Order")

ord.fe <- phyloseq_to_deseq2(ps.order, ~ Fe) # converts the phyloseq object into a Deseq object
counts.fe<- counts(ord.fe) 
geoMeans = apply(counts.fe, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
ddsLove.ord.fe <-estimateSizeFactors(ord.fe, geoMeans=geoMeans) 
dds.ord.fe <- DESeq(ddsLove.ord.fe, test="Wald", fitType="parametric") 

res.ord.fe = results(dds.ord.fe, cooksCutoff = FALSE)
alpha=0.05
sigtab.ord.fe = res.ord.fe[which(res.ord.fe$padj < alpha), ]
sigtab.ord.fe

ord.fe.df<- rownames_to_column(as.data.frame(sigtab.ord.fe), "Order")
ord.fe.df[,2:17]<-NULL
ord.fe.vec<-as.vector(ord.fe.df$Order)
ord.fe.ps<-prune_taxa(ord.fe.vec, ps.order) # this removes the significant genera out of the main phyloseq object and makes a new phyloseq object with only those taxa

diff.ord.fe = data.frame(otu_table(ord.fe.ps)) # exports otu table from phyloseq object for downstream use
diff.ord.meta = data.frame(sample_data(ord.fe.ps)) # exports sample data from phyloseq object for downstream use
diff.ord.fe$Fe <- diff.ord.meta$Fe # places Fe from the sample data into the otu table file
diff.ord.fe$Treatment <- diff.ord.meta$Treatment

ggplot(diff.ord.fe, aes(x= Fe, y = o__Boletales)) + geom_point(aes(col=Treatment), size = 3) + geom_smooth(method = lm, se=TRUE, fullrange=TRUE) + labs(y = "Boletales Abundance", x = "Fe Concentration") + geom_text(aes(x = 11000, y = 2000, label = lm_eqn(lm(o__Boletales ~ Fe, diff.ord.fe))), size = 3, parse=TRUE)

##### Canonical Correspondance Analysis ####

install.packages("vegan")
library(vegan)

esv.table.new<- data.frame(otu_table(ps.rare))

cca.all<- cca(esv.table.new ~ ., data = metadata.t) # build a CCA model with all variables

mod0<- cca(esv.table.new ~ 1, metadata.t) # create a null model
mod <- step(mod0, scope = formula(cca.all), test = "perm", perm.max = 100)

capture.output(step(mod0, scope = formula(cca.all), test = "perm", perm.max = 100)
, file = "ImportantEnvData.txt")
# compares the model with all the data to see which variables are most significant
# use the significant variables to build a new model
cca.var<- cca(formula=esv.table.new ~ Na + S + Mo + V + Be + Al, data = metadata.t)
# now we can compare how well each of them performs
anova(cca.all, by = "terms")
anova(cca.var) #significant p-values suggest model explains significantly more variation in data than a null model
# can also test whether there is significant variable explained by the axes
anova(cca.var, by = "axis", perm = 100) # indicates that CCA1 axis is significant 
anova(cca.var, by = "terms")
# this indicates that "V" and "Al" do not add significance to the model, so let's take it out

cca.var.2 <- cca(formula = esv.table.new ~ Na + S + Mo + Be, data = metadata.t)
anova(cca.var, cca.var.2) # this indicates that there is not a statistical difference between the two models

capture.output(anova(cca.var.2), file="CCA_Anova_Berts_ITS.txt")


plot(cca.all)

quartz()
#jpeg(file="CCA_WithAllImpMetals.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
plot(cca.all, type = "n", display = "sites")
points(cca.all, display = "sites", col = "black", pch=c(17,8,5)[meta.df$Treatment], cex = 0.9, font = 1)
points(cca.all, display = "bp", lwd = 2, col = "blue")
text(cca.all, display = "bp", col = "blue", font = 2)
#text(cca.all, display = "sites", col = "black")
with(meta.df, legend("topright", inset = c(0.1,0), legend=levels(meta.df$Treatment), pch=c(17,8,5), text.font=0.9, cex = 0.9, ncol = 1, bty = "n", title = "Treatments"))
#dev.off()

quartz()
tiff(file="CCA_wMostImpMetals.tiff", width=10, height=8, units ="in", res=300) #this makes a jpeg file named XXX in your output directory
plot(cca.var.2, type = "n", display = "sites")
points(cca.var.2, display = "sites", col = "black", pch=c(17,8,5)[meta.df$Treatment], cex = 0.9, font = 1)
points(cca.var.2, display = "bp", lwd = 2, col = "blue")
text(cca.var.2, display = "bp", col = "blue", font = 2)
with(meta.df, legend("topright", inset = c(0.1,0), legend=levels(meta.df$Treatment), pch=c(17,8,5), cex = 0.9, ncol = 1, bty = "n", title = "Treatments"))
dev.off()


jpeg(file="CCA_wMostImpMetals.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
plot(cca.var.2, type = "n", display = "sites")
points(cca.var.2, display = "sites", col = "black", pch=c(17,8,5)[meta.df$Treatment], cex = 0.9, font = 1)
points(cca.var.2, display = "bp", lwd = 2, col = "blue")
text(cca.var.2, display = "bp", col = "blue", font = 2)
with(meta.df, legend("topright", inset = c(0.1,0), legend=levels(meta.df$Treatment), pch=c(17,8,5), cex = 0.9, ncol = 1, bty = "n", title = "Treatments"))
dev.off()


##### Session 4 #####
#### Differential Abundance at Taxonomic Levels with All Variables
##Load phyloseq object and transfer into DESeq2 objectfor binomial test ##
## BiocAnyway(package name)##
## If your R version is older than 3.5, the function will use biocLite(), unless it will use BiocManager ##

BiocAnyway <- function(name_of_package){
  if (substr(version$version.strin, 11, 15) >= 3.5){
    install.packages("BiocManager")
    BiocManager::install(name_of_package)}
  else {source('http://bioconductor.org/biocLite.R')
    biocLite(name_of_package)}
}

BiocAnyway(DESeq2)
BiocAnyway(phyloseq)
install.packages("ggplot2")
library(phyloseq)
library(DESeq2)
library(ggplot2)

BiocManager::valid()
#if the response to the above is not "TRUE", you wil have problems running DESeq
# need to update the packages and then restart R



diagdds.flav <- phyloseq_to_deseq2(ps.rare, ~ Treatment)
counts.dds<- counts(diagdds.flav)
geoMeans = apply(counts.dds, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
ddsLove<-estimateSizeFactors(diagdds.flav, geoMeans=geoMeans)
estimatensizefactors<-as.matrix(counts(ddsLove, normalized=TRUE)) 
dds <- estimateDispersions(ddsLove)
estimatedispersions<-as.data.frame(counts(dds, normalized=TRUE))
dds.1 <- nbinomWaldTest(dds)
##testing log2 fold change ##
alpha=0.05
#DESeq does pairwise comparisons, so you need to test two levels at a time
res.nat_ino = results(dds.1, cooksCutoff = FALSE, contrast = c("Treatment", "Natural","Inoculated"))
sigtab.nat_ino = res.nat_ino[which(res.nat_ino$padj < alpha), ]
sigtab.nat_ino_df = data.frame(sigtab.nat_ino)
sigtab.nat_ino_df = cbind(ESV_name=row.names(sigtab.nat_ino_df), sigtab.nat_ino_df)
sigtab.nat_ino_df = sigtab.nat_ino_df[order(sigtab.nat_ino_df$log2FoldChange),]
sigtab.nat_ino_df # the output includes only the taxonomic groups that are differentially abundance (adjusted p<0.05)

# now, you'll want to create a graph of this information
d1 <-ggplot(sigtab.nat_ino_df, aes(x=log2FoldChange, y=reorder(ESV_name, log2FoldChange))) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) + coord_flip() +
  geom_point(size=3) +
  theme_set(theme_bw()) +
  labs(x ="log2 Fold Change", y ="ESV") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=10),axis.text.y= element_text(size=10), legend.position="none") +
  theme(axis.title= element_text(size=10))
d1 = d1 + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
d1

res.ref_ino = results(dds.1, cooksCutoff = FALSE, contrast = c("Treatment", "Ref","Inoculated"))
sigtab.ref_ino = res.ref_ino[which(res.ref_ino$padj < alpha), ]
sigtab.ref_ino_df = data.frame(sigtab.ref_ino)
sigtab.ref_ino_df = cbind(ESV_name=row.names(sigtab.ref_ino_df), sigtab.ref_ino_df)
sigtab.ref_ino_df = sigtab.ref_ino_df[order(sigtab.ref_ino_df$log2FoldChange),]
sigtab.ref_ino_df

d2 <-ggplot(sigtab.ref_ino_df, aes(x=log2FoldChange, y=reorder(ESV_name, log2FoldChange))) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) + coord_flip() +
  geom_point(size=3) +
  theme_set(theme_bw()) +
  labs(x ="log2 Fold Change", y ="ESV") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=10),axis.text.y= element_text(size=10), legend.position="none") +
  theme(axis.title= element_text(size=10))
d2 = d2 + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
d2

res.nat_ref = results(dds.1, cooksCutoff = FALSE, contrast = c("Treatment",  "Ref","Natural"))
sigtab.nat_ref = res.nat_ref[which(res.nat_ref$padj < alpha), ]
sigtab.nat_ref_df = data.frame(sigtab.nat_ref)
sigtab.nat_ref_df = cbind(ESV_name=row.names(sigtab.nat_ref_df), sigtab.nat_ref_df)
sigtab.nat_ref_df = sigtab.nat_ref_df[order(sigtab.nat_ref_df$log2FoldChange),]

sigtab.nat_ref_df

d3 <-ggplot(sigtab.nat_ref_df, aes(x=log2FoldChange, y=reorder(ESV_name, log2FoldChange, fill=ESV_name))) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) +
  geom_point(size=3) + coord_flip() +
  theme_set(theme_bw()) +
  labs(x ="log2 Fold Change", y ="ESV") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=10),axis.text.y= element_text(size=10), legend.position="none") +
  theme(axis.title= element_text(size=10))
d3 = d3 + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
d3

install.packages("cowplot")
library(cowplot)
log_plot <- plot_grid(d1,d2,d3, nrow =3)
log_plot


#jpeg(file="DESeq_Log2FoldChange_Plots.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
log_plot
#dev.off()


## Abundance Plot ##
# to graphically display these genera we need to convert the genera names to a vector
# From that vector we can copy the samples from the main phyloseq object to a new phyloseq object
install.packages("tibble")
library(tibble)
genera.esvs.gen1<-rownames_to_column(sigtab.nat_ino_df, "ESVs_gen")#this adds the ESVs from the rows to a column
genera.esvs.gen1[,2:17]<-NULL #this deletes the other columns from the file
genera.df1<-as.vector(genera.esvs.gen1$ESVs_gen) #this puts the ESVs into a vector
diff.genera.gen1<- prune_taxa(genera.df1, ps.rare)#this extracts the ESVS from the phyloseq relative abundance object
merged.gen1<- merge_samples(diff.genera.gen1, "Treatment")
merged.gen1<- subset_samples(merged.gen1, Treatment =="1"| Treatment=="2")
quartz()
p1<- plot_bar(merged.gen1, fill = "OTU", title = "Natural VS Inoculated")


genera.esvs.gen2<-rownames_to_column(sigtab.ref_ino_df, "ESVs_gen")#this adds the ESVs from the rows to a column
genera.esvs.gen2[,2:17]<-NULL #this deletes the other columns from the file
genera.df2<-as.vector(genera.esvs.gen2$ESVs_gen) #this puts the ESVs into a vector
diff.genera.gen2<- prune_taxa(genera.df2, ps.rare)#this extracts the ESVS from the phyloseq relative abundance object
merged.gen2<- merge_samples(diff.genera.gen1, "Treatment")
merged.gen2<- subset_samples(merged.gen2, Treatment =="1"| Treatment=="3")
p2 <- plot_bar(merged.gen2, fill = "OTU", title = "Inoculated VS Reference")
p2

genera.esvs.gen3<-rownames_to_column(sigtab.nat_ref_df, "ESVs_gen")#this adds the ESVs from the rows to a column
genera.esvs.gen3[,2:17]<-NULL #this deletes the other columns from the file
genera.df3<-as.vector(genera.esvs.gen3$ESVs_gen) #this puts the ESVs into a vector
diff.genera.gen3<- prune_taxa(genera.df3, ps.rare)#this extracts the ESVS from the phyloseq relative abundance object
merged.gen3<- merge_samples(diff.genera.gen3, "Treatment")
merged.gen3<- subset_samples(merged.gen3, Treatment =="2"| Treatment=="3")
p3 <- plot_bar(merged.gen3, fill = "OTU", title = "Natural VS Reference")
p3

ab_plot <- plot_grid(p1,p2,p3, nrow=3, ncol=1)

fin_plot <-plot_grid(log_plot, ab_plot,rel_widths = c(3,1), ncol=2)
fin_plot


#### Now we can do that at other taxonomic levels, too####
gen.diagdds <- phyloseq_to_deseq2(ps.genus, ~ Treatment)
counts.dds.gen<- counts(gen.diagdds)
geoMeans = apply(counts.dds.gen, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
ddsLove.gen <-estimateSizeFactors(gen.diagdds, geoMeans=geoMeans)
estimatensizefactors.gen<-as.matrix(counts(ddsLove.gen, normalized=TRUE)) 
dds.gen <- estimateDispersions(ddsLove.gen)
estimatedispersions.gen<-as.data.frame(counts(dds.gen, normalized=TRUE))
dds.gen <- nbinomWaldTest(dds.gen)

##testing log2 fold change ##
alpha=0.05
res.gen.nat_ino = results(dds.gen, cooksCutoff = FALSE, contrast = c("Treatment", "Natural","Inoculated"))
sigtab.gen.nat_ino = res.gen.nat_ino[which(res.gen.nat_ino$padj < alpha), ]
sigtab.gen.nat_ino_df = data.frame(sigtab.gen.nat_ino)
sigtab.gen.nat_ino_df = cbind(ESV_name=row.names(sigtab.gen.nat_ino_df), sigtab.gen.nat_ino_df)
sigtab.gen.nat_ino_df = sigtab.gen.nat_ino_df[order(sigtab.gen.nat_ino_df$log2FoldChange),]
sigtab.gen.nat_ino_df

d4 <-ggplot(sigtab.gen.nat_ino_df, aes(x=log2FoldChange, y=reorder(ESV_name, log2FoldChange))) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) + coord_flip() +
  geom_point(size=3) +
  theme_set(theme_bw()) +
  labs(x ="log2 Fold Change", y ="Genus") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=10),axis.text.y= element_text(size=10), legend.position="none") +
  theme(axis.title= element_text(size=10))
d4 = d4 + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
d4

res.gen.ref_ino = results(dds.gen, cooksCutoff = FALSE, contrast = c("Treatment", "Ref","Inoculated"))
sigtab.gen.ref_ino = res.gen.ref_ino[which(res.gen.ref_ino$padj < alpha), ]
sigtab.gen.ref_ino_df = data.frame(sigtab.gen.ref_ino)
sigtab.gen.ref_ino_df = cbind(ESV_name=row.names(sigtab.gen.ref_ino_df), sigtab.gen.ref_ino_df)
sigtab.gen.ref_ino_df = sigtab.gen.ref_ino_df[order(sigtab.gen.ref_ino_df$log2FoldChange),]
sigtab.gen.ref_ino_df

d5 <-ggplot(sigtab.gen.ref_ino_df, aes(x=log2FoldChange, y=reorder(ESV_name, log2FoldChange))) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) + coord_flip() +
  geom_point(size=3) +
  theme_set(theme_bw()) +
  labs(x ="log2 Fold Change", y ="Genus") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=10),axis.text.y= element_text(size=10), legend.position="none") +
  theme(axis.title= element_text(size=10))
d5 = d5 + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
d5

res.gen.nat_ref = results(dds.gen, cooksCutoff = FALSE, contrast = c("Treatment",  "Ref","Natural"))
sigtab.gen.nat_ref = res.gen.nat_ref[which(res.gen.nat_ref$padj < alpha), ]
sigtab.gen.nat_ref_df = data.frame(sigtab.gen.nat_ref)
sigtab.gen.nat_ref_df = cbind(ESV_name=row.names(sigtab.gen.nat_ref_df), sigtab.gen.nat_ref_df)
sigtab.gen.nat_ref_df = sigtab.gen.nat_ref_df[order(sigtab.gen.nat_ref_df$log2FoldChange),]

sigtab.gen.nat_ref_df

d6 <-ggplot(sigtab.gen.nat_ref_df, aes(x=log2FoldChange, y=reorder(ESV_name, log2FoldChange, fill=ESV_name))) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) +
  geom_point(size=3) + coord_flip() +
  theme_set(theme_bw()) +
  labs(x ="log2 Fold Change", y ="Genus") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=10),axis.text.y= element_text(size=10), legend.position="none") +
  theme(axis.title= element_text(size=10))
d6 = d6 + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
d6

library(cowplot)
log_plot.gen <- plot_grid(d4,d5,d6, nrow =3)
log_plot.gen


jpeg(file="DESeq_Log2FoldChange_Plots-Genus.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
log_plot.gen
dev.off()


## Abundance Plot ##
# to graphically display these genera we need to convert the genera names to a vector
# From that vector we can copy the samples from the main phyloseq object to a new phyloseq object
library(tibble)
genera.esvs.gen1<-rownames_to_column(sigtab.nat_ino_df, "ESVs_gen")#this adds the ESVs from the rows to a column
genera.esvs.gen1[,2:17]<-NULL #this deletes the other columns from the file
genera.df1<-as.vector(genera.esvs.gen1$ESVs_gen) #this puts the ESVs into a vector
diff.genera.gen1<- prune_taxa(genera.df1, ps.rare)#this extracts the ESVS from the phyloseq relative abundance object
merged.gen1<- merge_samples(diff.genera.gen1, "Treatment")
merged.gen1<- subset_samples(merged.gen1, Treatment =="1"| Treatment=="2")
p1 <- plot_bar(merged.gen1, fill = "OTU", title = "Natural VS Inoculated")
p1

genera.esvs.gen2<-rownames_to_column(sigtab.ref_ino_df, "ESVs_gen")#this adds the ESVs from the rows to a column
genera.esvs.gen2[,2:17]<-NULL #this deletes the other columns from the file
genera.df2<-as.vector(genera.esvs.gen2$ESVs_gen) #this puts the ESVs into a vector
diff.genera.gen2<- prune_taxa(genera.df2, ps.rare)#this extracts the ESVS from the phyloseq relative abundance object
merged.gen2<- merge_samples(diff.genera.gen1, "Treatment")
merged.gen2<- subset_samples(merged.gen2, Treatment =="1"| Treatment=="3")
p2 <- plot_bar(merged.gen2, fill = "OTU", title = "Inoculated VS Reference")
p2

genera.esvs.gen3<-rownames_to_column(sigtab.nat_ref_df, "ESVs_gen")#this adds the ESVs from the rows to a column
genera.esvs.gen3[,2:17]<-NULL #this deletes the other columns from the file
genera.df3<-as.vector(genera.esvs.gen3$ESVs_gen) #this puts the ESVs into a vector
diff.genera.gen3<- prune_taxa(genera.df3, ps.rare)#this extracts the ESVS from the phyloseq relative abundance object
merged.gen3<- merge_samples(diff.genera.gen3, "Treatment")
merged.gen3<- subset_samples(merged.gen3, Treatment =="2"| Treatment=="3")
p3 <- plot_bar(merged.gen3, fill = "OTU", title = "Natural VS Reference")
p3

quartz()
p3
ab_plot <- plot_grid(p1,p2,p3, nrow=3, ncol=1)

fin_plot <-plot_grid(log_plot, ab_plot,rel_widths = c(3,1), ncol=2)
fin_plot
