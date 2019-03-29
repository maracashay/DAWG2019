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


#### 11. General processing and analysis using silva phyloseq object####

#####Now we can begin using the phyloseq object for analyses, but first, we need to rarefy the dataset using the following commands

sample_sums(ps)
# the sample with the smallest sequence count is E28 with 6, so let's remove samples with fewer than 2,000 seqs
ps.pruned <- prune_samples(sample_sums(ps)>=1400, ps)
sample_sums(ps.pruned)

# 1 sample removed (E28)

# now we will rarefy the phyloseq object
set.seed(500)
ps.rare<-rarefy_even_depth(ps.pruned)
# 60 OTUs removed


sample_sums(ps.rare)
# Now 1,452 sequences are in each of the 29 samples

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

alpha.div<-estimate_richness(ps.rare, measures=c("Shannon", "Simpson"))
alpha.div


histogram(alpha.div$Shannon)
histogram(alpha.div$Simpson)

boxplot(alpha.div$Shannon ~ meta.df$Treatment)
# due to uneven sampling size and data being non-normal for shannon and simpson, use kruskal-wallis test
shannon.model<-kruskal.test(alpha.div$Shannon ~ meta.df$Treatment)
shannon.model
dunnTest(alpha.div$Shannon ~ meta.df$Treatment, method = "bonferroni")

simpson.model<-kruskal.test(alpha.div$Simpson ~ meta.df$Treatment)
simpson.model

# now calculate evenness and perform test
alpha.div$evenness<-evenness(ps.rare, index = "all", zeroes = TRUE)

pileu<- kruskal.test(alpha.div$evenness$pielou ~ meta.df$Treatment)
pileu

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
plot_heatmap(ps.phylum.rare, "PCoA", "bray", "Treatment")


# Use this to look at genera # 
ps.genera.perc<- taxa_level(ps.rare, "Genus")
plot_heatmap(ps.genera.perc, "PCoA", "bray", "Treatment", na.value = "grey",  low= "#66CCFF", high = "#000033" )
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
adonis(bc_rt_highcover_ns_abun ~ Treatment, data = meta.df)
adonis(bc_rt_highcover_ns_abun ~ Site, data = meta.df)

## All variables are found to have a significant effect on the distribution (i.e. grouping) of these microorganisms.


## F. Perform PERMANOVA to check significance of categories: JC

# PERMANOVA analyses using the Jaccard matrix.
adonis(jc_rt_highcover_ns_abun ~ Treatment, data = meta.df)
adonis(jc_rt_highcover_ns_abun ~ Site, data = meta.df)


### G. Posthoc: determining group differences in all data: (~ Treatment + Site): BC

#packages that you will need for this:
install.packages("RVAideMemoire")
library(RVAideMemoire)

set.seed(1)
# Pariwise Permutation MANOVAs
pairwise.perm.manova(bc_rt_highcover_ns_abun, meta.df$Treatment ,nperm=1000)
pairwise.perm.manova(bc_rt_highcover_ns_abun, meta.df$Site ,nperm=1000)


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
pcoa

# Create output folder
dir.create("~/Downloads/Berts_ITS_Data/demultiplexed_fastq/plots/")

# Save figure as PDF
pdf("~/Downloads/Berts_ITS_Data/demultiplexed_fastq/plots/PCoA_ps.perc_BC_Treatment.pdf",width=9,height=5)
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


#### 14.Differential Abundance Analysis ####
# packages you will need include
install.packages("cowplot")
library(cowplot)
BiocManager::install("DESeq2")
library(DESeq2)


# Look for ESVs that are differentially abundance across the treatments
diagdds.flav <- phyloseq_to_deseq2(ps.rare, ~ Treatment) # converts the phyloseq object into a Deseq object
counts.dds<- counts(diagdds.flav) 
geoMeans = apply(counts.dds, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
ddsLove<-estimateSizeFactors(diagdds.flav, geoMeans=geoMeans) 
#calculates the geometric mean of the rows and uses the median of the ratios as a size factor for each column
dds <- estimateDispersions(ddsLove)
dds.1 <- nbinomWaldTest(dds)

## testing log2 fold change ##
# only a single pairwise comparison can be made per analysis
alpha=0.05 # set alpha
res.nat_ino = results(dds.1, cooksCutoff = FALSE, contrast = c("Treatment", "Natural","Inoculated"))
# creates deseq result object with base mean, log2 fold change, lfcse, stat, p-value, and adjusted p-value scores
sigtab.nat_ino = res.nat_ino[which(res.nat_ino$padj < alpha), ]
# keeps only the ESVs that had an adjusted p-value below the value we set for alpha
sigtab.nat_ino_df = data.frame(sigtab.nat_ino) # convers the output to a data frame
sigtab.nat_ino_df = cbind(ESV_name=row.names(sigtab.nat_ino_df), sigtab.nat_ino_df)
sigtab.nat_ino_df = sigtab.nat_ino_df[order(sigtab.nat_ino_df$log2FoldChange),]
# orders the ESVs by the log 2 fold change
sigtab.nat_ino_df

# now that the ESVs that are differentiall abundance across the Natural and Inoculated samples, we can graph their log 2 fold changes
d1 <-ggplot(sigtab.nat_ino_df, aes(x=log2FoldChange, y=reorder(ESV_name, log2FoldChange))) + 
  geom_vline(xintercept = 0.0, color = "red", size = 1) +
  geom_point(size=3) + coord_flip() +
  theme_set(theme_bw()) +
  labs(x ="log2 Fold Change", y ="ESV") + 
  theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust=0.5, size=10),axis.text.y= element_text(size=10), legend.position="none") +
  theme(axis.title= element_text(size=10))
d1 = d1 + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
d1

# now we can calculate log 2 fold change between the Reference and the inoculated samples
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

# lastly, test the reference and the natural samples
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


log_plot <- plot_grid(d1,d2,d3, nrow =3)

#now, let's export the set of graphs
jpeg(file="Diff_Abundant_ESVs_LogChange.jpeg", width=10, height=8, units ="in", res=800) #this makes a jpeg file named XXX in your output directory
log_plot
dev.off()

## Abundance Plot ##
# to graphically display these genera we need to convert the genera names to a vector
# From that vector we can copy the samples from the main phyloseq object to a new phyloseq object
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

ab_plot <- plot_grid(p1,p2,p3, nrow=3, ncol=1)
ab_plot




