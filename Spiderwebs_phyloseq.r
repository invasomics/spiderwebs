# Following: https://benjjneb.github.io/dada2/tutorial_1_8.html

# Load in the required libraries
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(ggplot2)
library(plyr)
library(dplyr)
library(readr)
library(readxl)
library(vegan)
library(car)
library(tidyverse)
library(ape)
library(microbiome)
library(gridExtra)
library(ggpubr)
library(ggeasy)
library("MiscMetabar")

# Substitute in the input files based on the dataset/amplicon of choice;
# COI amplicon demonstrated here

########################################### 
#### PHYLOSEQ
# Ste the environment theme
theme_set(theme_bw())

# Read in the ASV counts file:
asv_mat<- read_tsv("ASVs_counts_bespokeCOI.tsv")

# Read in the taxonomy file
tax_mat<- read.csv("MergedTaxonomy_COI.csv", header=T)
samples_df <- read_excel("Spiderwebs_metadata_COI.xlsx")  

#Define the row names from the asv column
asv_mat <- asv_mat %>%
  tibble::column_to_rownames("ASV") 
#Sane for the two other matrices
tax_mat <- tax_mat %>% 
  tibble::column_to_rownames("ASV")
samples_df <- samples_df %>% 
  tibble::column_to_rownames("sample") 

# Transform into matrices asv and tax tables (sample table can be left as data frame)
asv_mat <- as.matrix(asv_mat)
tax_mat <- as.matrix(tax_mat)

# Construct a phyloseq object
ASV <- otu_table(asv_mat, taxa_are_rows = TRUE)
TAX <- tax_table(tax_mat)
samples <- sample_data(samples_df)

ps <- phyloseq(ASV, TAX, samples)
ps

########################################### 

# DATA EXPLORATION

sample_names(ps)
rank_names(ps)
sample_variables(ps)

# Check your phyla & filter for uncharacterised reads
table(tax_table(ps)[, "Phylum"], exclude = NULL)

ps0 <- subset_taxa(ps, !is.na(Phylum) & !Phylum %in% c("", "uncharacterized"))

ps0 #this new ps object should have less taxa if any of the Phyla results above returned NA

# Compute prevalence of each feature, store as data.frame 
prevdf = apply(X = otu_table(ps0), 
               MARGIN = ifelse(taxa_are_rows(ps0), yes = 1, no = 2), 
               FUN = function(x){sum(x > 0)})

# Add taxonomy and total read counts to this data.frame 
prevdf = data.frame(Prevalence = prevdf, 
                    TotalAbundance = taxa_sums(ps0),
                    tax_table(ps0))

#Visualise prevalence of each phyla
plyr::ddply(prevdf, "Phylum", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(prevdf, "Order", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})
plyr::ddply(prevdf, "Class", function(df1){cbind(mean(df1$Prevalence),sum(df1$Prevalence))})

#Check ps object & if everything looks OK, rename to ps
ps0 
ps <- ps0
ps

# Create a sample sum table to look at coverage metrics
sample_sums(ps)

sample_sum_df <- data.frame(sum=sample_sums(ps))
ggplot(sample_sum_df, aes(x = sum)) + 
  geom_histogram(color = "black", fill = "indianred", binwidth = 2500) +
  ggtitle("Distribution of sample sequencing depth") + 
  xlab("Read counts") +
  theme(axis.title.y = element_blank())

write.csv(sample_sums(ps),"COI_Reads_Per_Sample.csv") #outputs a coverage csv file

# Look for skew in coverage across sample types
set.seed(711)
DATA.2 <- ps  

df = as.data.frame(sample_data(DATA.2))
df$LibrarySize = sample_sums(DATA.2)
df = df[order(df$LibrarySize),]
df$Index = seq(nrow(df))

#Plot ordered library size & coloured by 'collection method'  
level_order2 <- c('Water', 'Web_artificial', 'Web_natural') 

ggplot(data=df, aes(x=Index, y=LibrarySize, colour= collection_method))+  
  geom_point(size =2.5)+
  facet_wrap(~ factor(collection_method, level = level_order2)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  scale_y_continuous(trans='sqrt') 

## Exploring the taxonomic diversity
table(phyloseq::tax_table(ps)[, "Phylum"])
table(phyloseq::tax_table(ps)[, "Class"])  
table(phyloseq::tax_table(ps)[, "Order"])
table(phyloseq::tax_table(ps)[, "Family"])
table(phyloseq::tax_table(ps)[, "Genus"])
table(phyloseq::tax_table(ps)[, "Species"])

#Basic bar graphs based on taxonomic group
plot_bar(ps, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

# Remove Hominidae reads
ps_nohomo <- (subset_taxa(ps, Family != "Hominidae"))
ps_nohomo
plot_bar(ps_nohomo, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

# IF SOME SAMPLES HAVE REALLY HIGH ABUNDANCE COMPARED TO OTHERS, MAY NEED TO RAREFY:
# this is just part of data exploration for now, to understand the data
ps_rare1 <- rarefy_even_depth(ps_nohomo, sample.size=100,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE) 
ps_rare1 # run this value to see how many taxa and samples are left when you reduced the dataset down to the abundance set by sample.size
# Replot and see what it looks like
plot_bar(ps_rare1, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
# Plot will now show abundance (at x reads) across samples, with x = the value you set for sample.size

## Top 20 hits in terms of relative abundance for unrarefied data:
ps_rel_abund <- phyloseq::transform_sample_counts(ps_nohomo, function(x){x / sum(x)})
ps_rel_abund
head(otu_table(ps_rel_abund))

top20 <- names(sort(taxa_sums(ps_rel_abund), decreasing=TRUE))[1:20]
ps.top20 <- transform_sample_counts(ps_rel_abund, function(OTU) OTU/sum(OTU))
ps.top20 <- prune_taxa(top20, ps.top20)
plot_bar(ps.top20, x="Sample", fill="Phylum") + facet_wrap(~collection_method, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Class") + facet_wrap(~collection_method, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Order") + facet_wrap(~collection_method, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Family") + facet_wrap(~collection_method, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Genus") + facet_wrap(~collection_method, scales="free_x") + labs(x="Sample names",y="Relative abundance")
plot_bar(ps.top20, x="Sample", fill="Species") + facet_wrap(~collection_method, scales="free_x") + labs(x="Sample names",y="Relative abundance")

########################################### 

#RAREFACTION CURVES # not normalised
rarecurve(t(asv_mat), step=100, col=mycol, lwd=2, ylab="ASVs", label=F, color=T, xlim=c(0,60000)) 

#############################################
### SACs    # not normalised

p <- accu_plot(ps_nohomo, "collection_method", add_nb_seq = TRUE, by.fact = TRUE, step = 10)
p +
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon'))

# see https://rdrr.io/github/adrientaudiere/MiscMetabar/man/accu_plot.html

########################################### 

#ALPHA DIVERSITY -- FULL DATASET, requires raw read counts/abundances for richness measures

########################################### 
#Visualise alpha-diversity:
alpha_div <- estimate_richness(ps_nohomo, measures = "Shannon")  
alpha_div

alpha_div$samplename <- rownames(alpha_div)
geo_data <- data.frame(sample_data(ps_nohomo))
geo_data$samplename = rownames(geo_data)
geo_data <- merge(geo_data, alpha_div, by = "samplename")

ggplot(geo_data, aes(x = collection_method, y = Shannon, color = collection_method)) +
  geom_point(size = 4, alpha = 0.5) +
  geom_smooth(method = "glm", se = FALSE, color = "lightgray", linetype = "dashed") + 
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  theme_minimal() +
  labs(title = "",
       x = "Collection method",
       y = "Alpha diversity (Shannon)") +
  theme(legend.position = "none")

################################################################

#ORDINATION  -- NORMALISED DATA 

################################################################
ps_rare1

#Generating and visualising the PCoA
vst_pcoa <- ordinate(ps_rare1, method="MDS", distance="bray")  

eigen_vals <- vst_pcoa$values$Eigenvalues

plot_ordination(ps_rare1, vst_pcoa, color="collection_method") + 
  geom_point(size = 5) +
  coord_fixed(sqrt(eigen_vals[2]/eigen_vals[1])) +
  theme_bw() + theme(legend.position="right") +
  scale_color_manual(values = c("Water" = "cornflowerblue","Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  stat_ellipse()

## Permanovas
dist.matrix = phyloseq::distance(ps_rare1, method = 'bray') 

# If no samples lost during rarefaction, run:
ps.perm = adonis2(dist.matrix ~ samples$collection_method, data = samples_df) 

# If samples lost during rarefaction, use:
samples_df2 <- subset(samples_df, sample_name !="NSW2" & sample_name != "NSW10")   
# need to update for your data by comparing sample_names(ps_rare1) to sample_names(ps_nohomo)
ps.perm = adonis2(dist.matrix ~ samples_df2$collection_method, data = samples_df2)

ps.perm

# Beta dispersion  
bd <- betadisper(dist.matrix, samples$collection_method)
# or
bd <- betadisper(dist.matrix, samples_df2$collection_method)
bd
mycol <- c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')
boxplot(bd, col=mycol)
permutest(bd, pairwise=T)

#############################################

#### CLAMTESTS -- normalised

#################################################
# First subset to just webs:
webs_only <- subset_samples(ps_nohomo, sample_data(ps_nohomo)$collection_method !="Water")  

plot_bar(webs_only, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
webs_rare <- rarefy_even_depth(webs_only, sample.size=1000,rngseed=7,replace=F,trimOTUs=TRUE,verbose=TRUE)
webs_rare   

webs_clam <- clamtest(t(otu_table(webs_rare)), (sample_data(webs_rare)$collection_method), alpha = 0.01, specialization = 2/3)

web_general <- webs_clam[webs_clam$Classes=="Generalist",]
web_rare <- webs_clam[webs_clam$Classes=="Too_rare",]
web_specialist <- webs_clam[webs_clam$Classes=="Specialist_presence",]
web_spec_abs <- webs_clam[webs_clam$Classes=="Specialist_absence",]
plot(webs_clam)

summary(webs_clam)

# Get species taxon names
web_specialist_nat <- webs_clam[webs_clam$Classes=="Specialist_Web_natural",]
ASV_specialist_natural <- web_specialist_nat$Species
tax_table_nat <- as.data.frame(tax_table(webs_only))
tax_subset_nat <- tax_table_nat[row.names(tax_table_nat) %in% ASV_specialist_natural, ]
write.csv(tax_subset_nat, file = "ITS_specialist_naturalwebs.csv")

web_specialist_art <- webs_clam[webs_clam$Classes=="Specialist_Web_artificial",]
ASV_specialist_artificial <- web_specialist_art$Species
tax_table_art <- as.data.frame(tax_table(webs_only))
tax_subset_art <- tax_table_art[row.names(tax_table_art) %in% ASV_specialist_artificial, ]
write.csv(tax_subset_art, file = "ITS_specialist_artificialwebs.csv")

#################################################################

# MORE COMPLEX ANALYSES ON DATA SUBSETS AND PLOTTING IN 
# MULTI-PANEL PLOTS; ONLY USED FOR COI DATASET

#################################################################

# 1. Subset data to Fungi, Chordates, Inverts

table(tax_table(ps)[, "Phylum"], exclude = NULL)

ps_fungi <- subset_taxa(ps, Phylum =="Ascomycota" | Phylum =="Basidiomycota"| Phylum =="Mucoromycota") 
ps_chordates <- subset_taxa(ps, Phylum =="Chordata") 
ps_inverts <- subset_taxa(ps, Phylum =="Annelida" | Phylum =="Arthropoda"| Phylum =="Cnidaria" | Phylum == "Gastrotricha" | Phylum =="Mollusca" | 
                            Phylum =="Nematoda" | Phylum =="Platyhelminthes" | Phylum =="Porifera" | Phylum == 'Rotifera' | Phylum == 'Tardigrada') 
ps_plantalg <- subset_taxa(ps, Phylum =="Bacillariophyta" | Phylum =="Chlorophyta"| Phylum =="Haptophyta" | Phylum =="Rhodophyta" | Phylum == 
                             "Streptophyta") 

# Look for skew in coverage across sample types
set.seed(2010)
DATA.1 <- ps_fungi 
DATA.2 <- ps_chordates 
DATA.3 <- ps_inverts  
DATA.4 <- ps_plantalg

df1 = as.data.frame(sample_data(DATA.1))
df1$LibrarySize = sample_sums(DATA.1)
df1 = df1[order(df1$LibrarySize),]
df1$Index = seq(nrow(df1))

df2 = as.data.frame(sample_data(DATA.2))
df2$LibrarySize = sample_sums(DATA.2)
df2 = df2[order(df2$LibrarySize),]
df2$Index = seq(nrow(df2))

df3 = as.data.frame(sample_data(DATA.3))
df3$LibrarySize = sample_sums(DATA.3)
df3 = df3[order(df3$LibrarySize),]
df3$Index = seq(nrow(df3))

df4 = as.data.frame(sample_data(DATA.4))
df4$LibrarySize = sample_sums(DATA.4)
df4 = df4[order(df4$LibrarySize),]
df4$Index = seq(nrow(df4))

#Plot ordered library size & coloured by 'collection method'

level_order2 <- c('Water', 'Web_artificial', 'Web_natural')

p1 <- ggplot(data=df1, aes(x=Index, y=LibrarySize, colour= collection_method))+
  geom_point(size =2.5)+
  facet_wrap(~ factor(collection_method, level = level_order2)) +
  ggtitle("Fungi") +
  ggeasy::easy_center_title() +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) 

p2 <- ggplot(data=df2, aes(x=Index, y=LibrarySize, colour= collection_method))+
  geom_point(size =2.5)+
  ggtitle("Chordates") +
  ggeasy::easy_center_title() +
  facet_wrap(~ factor(collection_method, level = level_order2)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) 

p3 <- ggplot(data=df3, aes(x=Index, y=LibrarySize, colour= collection_method))+
  geom_point(size =2.5)+
  ggtitle("Invertebrates") +
  ggeasy::easy_center_title() +
  facet_wrap(~ factor(collection_method, level = level_order2)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) 

p4 <- ggplot(data=df4, aes(x=Index, y=LibrarySize, colour= collection_method))+
  geom_point(size =2.5)+
  ggtitle("Plants/Algae") +
  ggeasy::easy_center_title() +
  facet_wrap(~ factor(collection_method, level = level_order2)) +
  theme(strip.background = element_blank(), strip.text.x = element_blank()) +
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) 

ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2, common.legend = TRUE, legend="right")

#Basic bar graphs based on taxonomic group
plot_bar(ps_fungi, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(ps_chordates, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(ps_inverts, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

plot_bar(ps_plantalg, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

## Top 20 relative abundance 
ps_rel_abund_F <- phyloseq::transform_sample_counts(ps_fungi, function(x){x / sum(x)})
ps_rel_abund_C <- phyloseq::transform_sample_counts(ps_chordates, function(x){x / sum(x)})
ps_rel_abund_I <- phyloseq::transform_sample_counts(ps_inverts, function(x){x / sum(x)})
ps_rel_abund_P <- phyloseq::transform_sample_counts(ps_plantalg, function(x){x / sum(x)})

top20_F <- names(sort(taxa_sums(ps_rel_abund_F), decreasing=TRUE))[1:20]
ps.top20_F <- transform_sample_counts(ps_rel_abund_F, function(OTU) OTU/sum(OTU))
ps.top20_F <- prune_taxa(top20_F, ps.top20_F)
plot_bar(ps.top20_F, x="Sample", fill="Phylum") + facet_wrap(~collection_method, scales="free_x") + labs(x="Sample names",y="Relative abundance")

top20_C <- names(sort(taxa_sums(ps_rel_abund_C), decreasing=TRUE))[1:20]
ps.top20_C <- transform_sample_counts(ps_rel_abund_C, function(OTU) OTU/sum(OTU))
ps.top20_C <- prune_taxa(top20_C, ps.top20_C)
plot_bar(ps.top20_C, x="Sample", fill="Phylum") + facet_wrap(~collection_method, scales="free_x") + labs(x="Sample names",y="Relative abundance")

top20_I <- names(sort(taxa_sums(ps_rel_abund_I), decreasing=TRUE))[1:20]
ps.top20_I <- transform_sample_counts(ps_rel_abund_I, function(OTU) OTU/sum(OTU))
ps.top20_I <- prune_taxa(top20_I, ps.top20_I)
plot_bar(ps.top20_I, x="Sample", fill="Phylum") + facet_wrap(~collection_method, scales="free_x") + labs(x="Sample names",y="Relative abundance")

top20_P <- names(sort(taxa_sums(ps_rel_abund_P), decreasing=TRUE))[1:20]
ps.top20_P <- transform_sample_counts(ps_rel_abund_P, function(OTU) OTU/sum(OTU))
ps.top20_P <- prune_taxa(top20_P, ps.top20_P)
plot_bar(ps.top20_P, x="Sample", fill="Phylum") + facet_wrap(~collection_method, scales="free_x") + labs(x="Sample names",y="Relative abundance")

# any errors are because some samples = 0; can ignore

########################################### 

#ALPHA DIVERSITY -- FULL DATASET, requires raw read counts/abundances for richness measures

########################################### 
#Visualise alpha-diversity:
alpha_div_F <- estimate_richness(ps_fungi, measures = "Shannon")  #ignore the warning message about singletons (they're removed deliberately!)
alpha_div_C <- estimate_richness(ps_chordates, measures = "Shannon")  #ignore the warning message about singletons (they're removed deliberately!)
alpha_div_I <- estimate_richness(ps_inverts, measures = "Shannon")  #ignore the warning message about singletons (they're removed deliberately!)
alpha_div_P <- estimate_richness(ps_plantalg, measures = "Shannon")  #ignore the warning message about singletons (they're removed deliberately!)

alpha_div_F$samplename <- rownames(alpha_div_F)
geo_data_F <- data.frame(sample_data(ps_fungi))
geo_data_F$samplename = rownames(geo_data_F)
geo_data_F <- merge(geo_data_F, alpha_div_F, by = "samplename")

alpha_div_C$samplename <- rownames(alpha_div_C)
geo_data_C <- data.frame(sample_data(ps_chordates))
geo_data_C$samplename = rownames(geo_data_C)
geo_data_C <- merge(geo_data_C, alpha_div_C, by = "samplename")

alpha_div_I$samplename <- rownames(alpha_div_I)
geo_data_I <- data.frame(sample_data(ps_inverts))
geo_data_I$samplename = rownames(geo_data_I)
geo_data_I <- merge(geo_data_I, alpha_div_I, by = "samplename")

alpha_div_P$samplename <- rownames(alpha_div_P)
geo_data_P <- data.frame(sample_data(ps_plantalg))
geo_data_P$samplename = rownames(geo_data_P)
geo_data_P <- merge(geo_data_P, alpha_div_P, by = "samplename")

p1_a <-ggplot(geo_data_F, aes(x = collection_method, y = Shannon, color = collection_method)) +
  geom_point(size = 4, alpha = 0.5) +
  ggtitle("Fungi") +
  ggeasy::easy_center_title() +
  geom_smooth(method = "glm", se = FALSE, color = "lightgray", linetype = "dashed") +  # Optional smoothing line
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  theme_minimal() +
  labs(x = "Collection method",
       y = "Alpha diversity (Shannon)") +
  theme(legend.position = "none")

p2_a <-ggplot(geo_data_C, aes(x = collection_method, y = Shannon, color = collection_method)) +
  geom_point(size = 4, alpha = 0.5) +
  ggtitle("Chordates") +
  ggeasy::easy_center_title() +
  geom_smooth(method = "glm", se = FALSE, color = "lightgray", linetype = "dashed") +  # Optional smoothing line
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  theme_minimal() +
  labs(x = "Collection method",
       y = "Alpha diversity (Shannon)") +
  theme(legend.position = "none")

p3_a <-ggplot(geo_data_I, aes(x = collection_method, y = Shannon, color = collection_method)) +
  geom_point(size = 4, alpha = 0.5) +
  ggtitle("Invertebrates") +
  ggeasy::easy_center_title() +
  geom_smooth(method = "glm", se = FALSE, color = "lightgray", linetype = "dashed") +  # Optional smoothing line
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  theme_minimal() +
  labs(x = "Collection method",
       y = "Alpha diversity (Shannon)") +
  theme(legend.position = "none")

p4_a <-ggplot(geo_data_P, aes(x = collection_method, y = Shannon, color = collection_method)) +
  geom_point(size = 4, alpha = 0.5) +
  ggtitle("Plants/Algae") +
  ggeasy::easy_center_title() +
  geom_smooth(method = "glm", se = FALSE, color = "lightgray", linetype = "dashed") +  # Optional smoothing line
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  theme_minimal() +
  labs(x = "Collection method",
       y = "Alpha diversity (Shannon)") +
  theme(legend.position = "none")

ggarrange(p1_a,p2_a,p3_a,p4_a, nrow = 2, ncol = 2)

################################################################ 

#ORDINATION  -- NORMALISE FIRST

################################################################
# First normalise using rarefy
ps_rare_F <- rarefy_even_depth(ps_fungi, sample.size=100,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE)
ps_rare_F   
ps_rare_C <- rarefy_even_depth(ps_chordates, sample.size=100,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE)
ps_rare_C  
ps_rare_I <- rarefy_even_depth(ps_inverts, sample.size=100,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE)
ps_rare_I  
ps_rare_P <- rarefy_even_depth(ps_plantalg, sample.size=100,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE)
ps_rare_P  

# Cannot do chordates or plant/algae as there is not enough data; just proceed with fungi and inverts
# Use geometric means:
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

diagdds <- phyloseq_to_deseq2(ps, ~ collection_method)  #on whole ps object
# calculate geometric means prior to estimate size factors

diagdds <- estimateSizeFactors(diagdds, type = "poscounts")

# Normalise counts
normalised_counts <- counts(diagdds, normalized = TRUE)
otu_table_norm <- otu_table(normalised_counts, taxa_are_rows = TRUE)

#Ordination_fungi
ps_fungi_norm <- phyloseq(otu_table_norm, tax_table(ps_fungi), sample_data(ps_fungi))

vst_pcoa_F_norm <- ordinate(ps_fungi_norm, method="MDS", distance="bray", na.rm = T) 
eigen_vals_F_norm <- vst_pcoa_F_norm$values$Eigenvalues

p_m1 <- plot_ordination(ps_fungi_norm, vst_pcoa_F_norm, color="collection_method") + 
  geom_point(size = 4) +
  coord_fixed(sqrt(eigen_vals_F_norm[2]/eigen_vals_F_norm[1])) +
  ylim(-1, 1) +
  xlim(-1, 1) +
  ggtitle("Fungi; Axis 1:2") +
  ggeasy::easy_center_title() +
  theme_bw() + theme(legend.position="right") +
  scale_color_manual(values = c("Water" = "cornflowerblue","Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  stat_ellipse()

#Plot axis 1 vs 3
p_m2 <- plot_ordination(ps_fungi_norm, vst_pcoa_F_norm, axes = c(1,3), color="collection_method") + 
  geom_point(size = 4) +
  coord_fixed(sqrt(eigen_vals_F_norm[2]/eigen_vals_F_norm[1])) +
  ylim(-1, 1) +
  xlim(-1, 1) +
  ggtitle("Fungi; Axis 1:3") +
  ggeasy::easy_center_title() +
  theme_bw() + theme(legend.position="right") +
  scale_color_manual(values = c("Water" = "cornflowerblue","Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  stat_ellipse()

# Inverts
ps_invert_norm <- phyloseq(otu_table_norm, tax_table(ps_inverts), sample_data(ps_inverts))

vst_pcoa_I_norm <- ordinate(ps_invert_norm, method="MDS", distance="manhattan", na.rm = T) 
eigen_vals_I_norm <- vst_pcoa_I_norm$values$Eigenvalues

p_m3 <- plot_ordination(ps_invert_norm, vst_pcoa_I_norm, color="collection_method") + 
  geom_point(size = 4) +
  coord_fixed(sqrt(eigen_vals_I_norm[2]/eigen_vals_I_norm[1])) +
  # ylim(-1, 1) +
  #  xlim(-1, 1) +
  ggtitle("Invertebrates; Axis 1:2") +
  ggeasy::easy_center_title() +
  theme_bw() + theme(legend.position="right") +
  scale_color_manual(values = c("Water" = "cornflowerblue","Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  stat_ellipse()

#Plot axis 1 vs 3
p_m4 <- plot_ordination(ps_invert_norm, vst_pcoa_I_norm, axes = c(1,3), color="collection_method") + 
  geom_point(size = 4) +
  coord_fixed(sqrt(eigen_vals_I_norm[2]/eigen_vals_I_norm[1])) +
  #ylim(-1, 1) +
  #  xlim(-1, 1) +
  ggtitle("Invertebrates; Axis 1:3") +
  ggeasy::easy_center_title() +
  theme_bw() + theme(legend.position="right") +
  scale_color_manual(values = c("Water" = "cornflowerblue","Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) +
  stat_ellipse()

ggarrange(p_m1,p_m2,nrow = 2, ncol =1, common.legend = TRUE, legend="bottom")
ggarrange(p_m3, p_m4, nrow = 2, ncol =1, common.legend = TRUE, legend="bottom")

## Permanovas
dist.matrix = phyloseq::distance(ps_fungi_norm, method = 'bray')
ps.perm = adonis2(dist.matrix~samples$collection_method, data = samples_df)
ps.perm
  
dist.matrix2 = phyloseq::distance(ps_invert_norm, method = 'bray')
ps.perm2 = adonis2(dist.matrix2~samples$collection_method, data = samples_df)
ps.perm2   
# Above doesn't work; missing values in results

# Beta dispersion
bd1 <- betadisper(dist.matrix, samples$collection_method)  
bd1

mycol <- c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')
boxplot(bd1, col=mycol)
permutest(bd1, pairwise=T)

bd2 <- betadisper(dist.matrix2, samples$collection_method)   # doesn't work
bd2
boxplot(bd2, col=mycol)
permutest(bd2, pairwise=T)

########################################### 

### SACs    # not normalised

########################################### 

p_fu <- accu_plot(ps_fungi, "collection_method", add_nb_seq = TRUE, by.fact = TRUE, step = 10)
p2_r <- p_fu +
  ylim(0, 1000) +
  ggtitle("Fungi") +
  ggeasy::easy_center_title() +
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) 

p_in <- accu_plot(ps_inverts, "collection_method", add_nb_seq = TRUE, by.fact = TRUE, step = 10)
p4_r <- p_in +
  ylim(0, 1000) +
  ggtitle("Invertebrates") +
  ggeasy::easy_center_title() +
  scale_color_manual(values = c("Water" = "cornflowerblue", "Web_artificial" = "brown1", "Web_natural" = 'darksalmon')) 

ggarrange(p2_r,p4_r,nrow = 1, ncol = 2, common.legend = TRUE, legend="right")

#############################################

#### CLAMTESTS 

#################################################
# First subset to just webs:
webs_only_F <- subset_samples(ps_fungi, sample_data(ps_fungi)$collection_method !="Water")  
webs_only_I <- subset_samples(ps_inverts, sample_data(ps_inverts)$collection_method !="Water")  

plot_bar(webs_only_F, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")
plot_bar(webs_only_I, fill = "Phylum") + 
  geom_bar(aes(color=Phylum, fill=Phylum), stat="identity", position="stack")

webs_rare_F <- rarefy_even_depth(webs_only_F, sample.size=2000,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE)
webs_rare_I <- rarefy_even_depth(webs_only_I, sample.size=100,rngseed=7,replace=FALSE,trimOTUs=TRUE,verbose=TRUE)

webs_rare_F  
webs_rare_I   

webs_clam_F <- clamtest(t(otu_table(webs_rare_F)), (sample_data(webs_rare_F)$collection_method), alpha = 0.01, specialization = 2/3)
webs_clam_I <- clamtest(t(otu_table(webs_rare_I)), (sample_data(webs_rare_I)$collection_method), alpha = 0.01, specialization = 2/3)

web_general <- webs_clam_F[webs_clam_F$Classes=="Generalist",]
web_rare <- webs_clam_F[webs_clam_F$Classes=="Too_rare",]
web_specialist_nat <- webs_clam_F[webs_clam_F$Classes=="Specialist_Web_natural",]
web_specialist_art <- webs_clam_F[webs_clam_F$Classes=="Specialist_Web_artificial",]
plot(webs_clam_F)
summary(webs_clam_F)

# Get species taxon names
ASV_specialist_natural <- web_specialist_nat$Species
tax_table_nat <- as.data.frame(tax_table(webs_only_F))
tax_subset_nat <- tax_table_nat[row.names(tax_table_nat) %in% ASV_specialist_natural, ]
write.csv(tax_subset_nat, file = "Fungi_specialist_naturalwebs.csv")

ASV_specialist_artificial <- web_specialist_art$Species
tax_subset_art <- tax_table_nat[row.names(tax_table_nat) %in% ASV_specialist_artificial, ]
write.csv(tax_subset_art, file = "Fungi_specialist_artificialwebs.csv")

# Repeat for inverts
web_general_I <- webs_clam_I[webs_clam_I$Classes=="Generalist",]
web_rare_I <- webs_clam_I[webs_clam_I$Classes=="Too_rare",]
web_specialist_nat_I <- webs_clam_I[webs_clam_I$Classes=="Specialist_Web_natural",]
web_specialist_art_I <- webs_clam_I[webs_clam_I$Classes=="Specialist_Web_artificial",]
plot(webs_clam_I)
summary(webs_clam_I)

# Get species taxon names
ASV_specialist_natural_I <- web_specialist_nat_I$Species
tax_table_nat_I <- as.data.frame(tax_table(webs_only_I))
tax_subset_nat_I <- tax_table_nat_I[row.names(tax_table_nat_I) %in% ASV_specialist_natural_I, ]
write.csv(tax_subset_nat_I, file = "Inverts_specialist_naturalwebs.csv")

ASV_specialist_artificial_I <- web_specialist_art_I$Species
tax_subset_art_I <- tax_table_nat_I[row.names(tax_table_nat_I) %in% ASV_specialist_artificial_I, ]
write.csv(tax_subset_art_I, file = "Inverts_specialist_artificialwebs.csv")
