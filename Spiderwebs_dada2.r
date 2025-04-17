# Running DADA2 pipeline on raw fastq files; script used for all three amplicons, COI amplicon demonstrated here
# See text for further details
# Raw fastq files (not demultiplexed to amplicon; follow Claident instructions in main text) avaiable on SRA:
# Bioproject number: PRJNA1251198

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")

#Load Libraries
library(dada2)

#Setting Path for input files
setwd('/insert library path here')
path='/insert library path here'
list.files(path)

#Specify forward & reverse read fastqs
fnFs <- sort(list.files(path, pattern="_mtcointf.forward.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_mtcointf.reverse.fastq", full.names = TRUE))

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 3) #split filenames with '_', pick first character string

sample.names

#Data quality metrics
plotQualityProfile(fnFs) 
plotQualityProfile(fnRs) 

#Filtering reads based on sequence quality scores and primers
filtFs <- file.path(path, "filtered_COI", paste0(sample.names, "_F_filt.fastq.gz")) #creating filtered folder
filtRs <- file.path(path, "filtered_COI", paste0(sample.names, "_R_filt.fastq.gz"))
#Pulling sample names from filtered fasta files
names(filtFs) <- sample.names 
names(filtRs) <- sample.names
#Trimming primers - set your primer sequence
FWD <- "GGWACWGGWTGAACWGTWTAYCCYCC"  #miCOIintF
REV <- "TAIACYTCIGGRTGICCRAARAAYCA"   #jgHCO2198 
trimLeft = c(FWD,REV)

#Use known primer sequences to trim from your amplicon sequences
out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(220,200),   
             maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,             # filter out all reads with > maxN=0 ambiguous nucleotides and >2 expected errors
              compress=TRUE, multithread=TRUE,trimLeft = c(25,26))         # remove first 26 nucleotides of F/R reads (length of primers?)
# rm.phix is default and removes reads that match the phiX genome
# truncQ=2 is deafult and truncates reads at first instance of a quality score less than or equal to 2

# Learning Error rates
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

# Plotting out the errors
png(filename="Error_COIbespoke_F.png")
plotErrors(errF, nominalQ=TRUE)
dev.off()
png(filename="Error_COIbespoke_R.png")
plotErrors(errR, nominalQ=TRUE)
dev.off()

# Dereplication
derepFs <- derepFastq(filtFs, verbose=TRUE)  # if throws an error, this will be because some reads didn't pass filter.
derepRs <- derepFastq(filtRs, verbose=TRUE)

# Name the derep-class objects by the sample names:
names(derepFs) <- sample.names
names(derepRs) <- sample.names

# Sample inference from forward reads
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)  
# Sample inference from reverse reads
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

# Merge paired reads
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose = TRUE)

# Making your ASV table. This is synonymous to OTU table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

# Length of ASV's and their frequencies
table(nchar(getSequences(seqtab)))
# The sequence table is a matrix with rows corresponding to (and named by) the samples, and columns corresponding to (and named by)

# Identifying and removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus",multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

# Track reads through the pipeline
# Good checkpoint to ensure you did not lose too many reads

getN <- function(x) sum(getUniques(x))
track<- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers,getN), rowSums(seqtab.nochim))
colnames(track) <- c("input","filtered","denoisedF","denoisedR","merged","nonchim")
rownames(track) <- sample.names
head(track)

write.csv(track, "SequencingStatistics_bespokeCOI.csv")

taxa <- assignTaxonomy(seqtab.nochim, "/path to/NCBI_DATABASE_taxon_17-03-25.fasta", 
                       tryRC=TRUE, multithread=TRUE)

# bespoke database created via CRABS; see https://github.com/gjeunen/reference_database_creator

# Check your taxonomy assignment
taxa.print <- taxa # Removing sequence rownames for display only
rownames(taxa.print) <- NULL
dim(taxa.print)
head(taxa.print)

# Giving our seq headers more manageable names (ASV_1, ASV_2...)
asv_seqs <- colnames(seqtab.nochim)
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

# making and writing out a fasta of our final ASV seqs:
    
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "bespokeCOI.fa")

# count table:
asv_tab <- t(seqtab.nochim)
row.names(asv_tab) <- sub(">", "", asv_headers)
write.table(asv_tab, "ASVs_counts_bespokeCOI.tsv", sep="\t", quote=F, col.names=NA)

##  Giving taxonomy table corresponding names as above (ASV_1, ASV_2...)

row.names(taxa.print) <- sub(">", "", asv_headers)
write.table(taxa.print, "ASVs_assigned_bespokeCOI.tsv", sep="\t", quote=F, col.names=NA)
