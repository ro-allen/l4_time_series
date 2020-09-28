# Pipeline based on only forward ITS reads, according to recommendations of Pauvert et al. (2019; Fungal Ecology). This technique recovers the a very high proportion of Basidiomycetes and Ascomycetes and does not lose taxa with an ITS1 region longer than read merging would allow (i.e. 600 bp)

# Install/Load Packages ---------------------------------------------------
install.packages("devtools")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("dada2", version = "3.9")
BiocManager::install("phyloseq", version = "3.9")
BiocManager::install("ShortRead", version = "3.9")
BiocManager::install("Biostrings", version = "3.9")
library(dada2)
library(phyloseq)
library(ShortRead)
library(Biostrings)

install.packages("ggplot2")
install.packages("cowplot")
library(ggplot2)
library(cowplot)


# Locate Files ------------------------------------------------------------
path <- ("")
list.files(path)

## retain only forward reads
fnFs <- sort(list.files(path, pattern = "Forward", full.names = TRUE))


# Remove Primers ----------------------------------------------------------
## filter sequences to remove ambigous bases
fnFs.filtN <- file.path(path, "filtN", basename(fnFs)) 
filterAndTrim(fnFs, fnFs.filtN, maxN = 0, multithread = TRUE)

## sequence of primers used
ITS1F <- "CTTGGTCATTTAGAGGAAGTAA"  
ITS2 <- "GCTGCGTTCTTCATCGATGC" 

## all possible primer orientations (biostrings)
allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer)  # convert primer sequence in DNAString for Biostrings package handling
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # convert back to character vector
}
ITS1F.orients <- allOrients(ITS1F)
ITS2.orients <- allOrients(ITS2)

## manually check sequences
ITS1F.orients
ITS2.orients

## remove primers using cutadapt (through terminal operating through RStudio)
path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))

ITS1F.RC <- dada2:::rc(ITS1F)
ITS2.RC <- dada2:::rc(ITS2)

R1.flags <- paste("-g", ITS1F, "-a", ITS2.RC) 

for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], # output files
                             fnFs.filtN[i])) # input files
} #for more information see arguments for cutadapt

# check the successful removal of primers (0 values should be returned on all counts)
rbind(FWD.ForwardReads = sapply(ITS1F.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ForwardReads = sapply(ITS2.orients, primerHits, fn = fnFs.cut[[1]]))


# Inspect Quality Profiles ------------------------------------------------
## forward and reverse fastq filenames already processd by cutadapt
cutFs <- sort(list.files(path.cut, pattern = ".fastq.gz", full.names = TRUE))

## extract sample names 
get.sample.name <- function(fname) strsplit(basename(fname), ".fastq")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

## plot forward read quality scores
plotQualityProfile(cutFs[1:1]) 


# Filter and Trim Sequences -----------------------------------------------
filtFs <- file.path(path.cut, "filtered", basename(cutFs))

## trimRight to eliminate low quality sequences at the end of each read
out <- filterAndTrim(cutFs, filtFs, maxN = 0, maxEE = 1, truncQ = 2, minLen = 50, trimRight = 2, rm.phix = TRUE, compress = TRUE, multithread = TRUE)  # trim right aims to remove low quality reads characteristic of the end of an illumina run, which may interfere with taxonomic assignment
head(out)


# Learn Error Rates -------------------------------------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)


# Infer ASVs --------------------------------------------------------------
## depreplicate reads to reduce computational load
derepFs <- derepFastq(filtFs, verbose = TRUE)

## name the derepilicated objects by the sample names
names(derepFs) <- sample.names

## infer amplicon sequence variants using the DADA2 algorithm
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)

## construct sequence table
seqtab.r1 <- makeSequenceTable(dadaFs) 
dim(seqtab.r1)


# Remove Chimeras ---------------------------------------------------------
## Pauvert et al. (2017) use this step as optional - however, this should not be optional. Chimera removal is essential for the accurate identification of funal taxa from high-throughput sequencing data.
seqtab.nochim.r1 <- removeBimeraDenovo(seqtab.r1, method = "consensus", multithread = TRUE, verbose = TRUE)


# Pipeline Evaluation -----------------------------------------------------
## evaluate sequence table
table(nchar(getSequences(seqtab.nochim.r1)))

## evaluate loss of sequences throughout pipeline
getN <- function(x) sum(getUniques(x))
track.r1 <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab.nochim.r1))
colnames(track.r1) <- c("input", "filtered", "denoisedF", 
                        "nonchim")
rownames(track.r1) <- sample.names
track.r1


# Taxonomic Assignment ----------------------------------------------------
## full path necessary if file is outside working directory
unite.ref <- "sh_general_release_dynamic_02.02.2019.fasta"
taxa.r1 <- assignTaxonomy(seqtab.nochim.r1, unite.ref, multithread = TRUE, tryRC = TRUE)

taxa.print.r1 <- taxa.r1  # Removing sequence rownames (to reduce clutter)
rownames(taxa.print.r1) <- NULL
head(taxa.print.r1)

taxa.r1.vaulot <- assignTaxonomy(seqtab.nochim.r1, unite.ref, minBoot = 0, outputBootstraps = TRUE, multithread = TRUE, tryRC = TRUE, verbose = TRUE)
head(taxa.r1.vaulot$boot)
head(taxa.r1.vaulot$tax) # pulls out your bootstrap information. Set minBoot to minimum desired threshold (default = 50, appropriate = 80)

taxa.print.r1.vaulot <- taxa.r1.vaulot  # Removing sequence rownames (to reduce clutter)
rownames(taxa.print.r1.vaulot$tax) <- NULL
head(taxa.print.r1.vaulot$tax)
rownames(taxa.print.r1.vaulot$boot) <- NULL
head(taxa.print.r1.vaulot$boot)

## plot taxonomic assignment
taxa.print.r1 = as.data.frame(taxa.print.r1)

## vector of assigned reads at each rank
reads.as.r1 = as.vector(c(sum(!is.na(taxa.print.r1$Kingdom)), sum(!is.na(taxa.print.r1$Phylum)),sum(!is.na(taxa.print.r1$Class)),sum(!is.na(taxa.print.r1$Order)),sum(!is.na(taxa.print.r1$Family)),sum(!is.na(taxa.print.r1$Genus)),sum(!is.na(taxa.print.r1$Species))))
## vector of unassgined reads at each rank
reads.un.r1 = as.vector(c(sum(is.na(taxa.print.r1$Kingdom)), sum(is.na(taxa.print.r1$Phylum)),sum(is.na(taxa.print.r1$Class)),sum(is.na(taxa.print.r1$Order)),sum(is.na(taxa.print.r1$Family)),sum(is.na(taxa.print.r1$Genus)),sum(is.na(taxa.print.r1$Species))))

## create plotting dataframe
rank = as.data.frame(c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
bind.as.r1 = cbind(reads.as.r1, as.vector(rank))
bind.as.r1$type = "Assigned"
colnames(bind.as.r1) = c("value", "rank", "type")
bind.un.r1 = cbind(reads.un.r1, as.vector(rank))
bind.un.r1$type = "Unassigned"
colnames(bind.un.r1) = c("value", "rank", "type")
df.pipe.r1 = rbind(bind.as.r1, bind.un.r1)

## level fators for plotting
df.pipe.r1$rank = factor(df.pipe.r1$rank, levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
df.pipe.r1$type = factor(df.pipe.r1$type, levels = c("Unassigned", "Assigned"))

## create and save plot
p.pipe.assign <- ggplot(df.pipe.r1, aes(x=rank, y=value, fill=type)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values = c("gray80","firebrick1")) +
  scale_y_continuous(expand = c(0,0), limits = c(0,200)) +
  xlab("") +
  ylab("ASVs") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

save_plot("pipeline.assignment.pauvert.png", p.pipe.assign, base_width = 8, base_aspect_ratio = 1.3)





# EXTRAS


#Upload metadata
sd.l4 = as.data.frame(rownames(seqtab.nochim.r1))
colnames(sd.l4) <- "sample"
dim(sd.l4)
head(sd.l4)

# create a vector for sample names
s.vec = as.vector(1:length(rownames(sd.l4)))  #number should reflect your total number of samples
s.nam = cbind("sample_", s.vec)
s.nam = as.data.frame(s.nam)
s.names = paste0(s.nam$V1, s.nam$s.vec)
s.names = as.data.frame(s.names)

# apply sample names to metadata 
row.names(sd.l4) = s.names$s.names
sd.l4 = as.data.frame(sd.l4)
head(sd.l4)

# apply sample names to sequence table
row.names(seqtab.nochim.r1) = s.names$s.names

# construct ASV names tables
dim(taxa.print.r1)
dim(seqtab.nochim.r1)
a.vec = as.vector(1:length(colnames(seqtab.nochim.r1)))  #number should reflect your total ASVs
a.nam = cbind("asv_", a.vec)
a.nam = as.data.frame(a.nam)
asv.names = paste0(a.nam$V1, a.nam$a.vec)
asv.names = as.data.frame(asv.names)

colnames(seqtab.nochim.r1) = asv.names$asv.names
taxa.l4 <- tax_table(taxa.print.r1)
rownames(taxa.l4) <- asv.names$asv.names
colnames(taxa.l4) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")

#Construct Phyloseq Object
l4.phy = phyloseq(tax_table(taxa.l4), otu_table(seqtab.nochim.r1, taxa_are_rows = FALSE), sample_data(sd.l4))

#Analysis ----

#reads remaining per sample
rowSums(otu_table(l4.phy))
mean(rowSums(otu_table(l4.phy))) 
min(rowSums(otu_table(l4.phy))) 
max(rowSums(otu_table(l4.phy))) 

# Restructuring tax_table to include 'best' classification
library(zoo)
library(tibble)
bc.t = t(as.data.frame(tax_table(l4.phy)))
bc.t[bc.t=="k__"] <- ""
bc.t[bc.t=="p__"] <- ""
bc.t[bc.t=="c__"] <- ""
bc.t[bc.t=="o__"] <- ""
bc.t[bc.t=="f__"] <- ""
bc.t[bc.t=="g__"] <- ""
bc.t[bc.t=="s__"] <- ""
bc.t[bc.t==""] <- NA
bc.fill = na.locf(bc.t, na.rm = TRUE)
t.bc.fill = as.data.frame(t(bc.fill))
head(t.bc.fill)
rnc.bc = rownames_to_column(t.bc.fill, "ASV")
rnc.bc = as.data.frame(rnc.bc)

## Creates a column with the best classification and the ASV
rnc.bc$taxa_ASV = paste(rnc.bc$Genus,rnc.bc$Species,rnc.bc$ASV)

## Bind this column back onto the original tax_table 
safe.bc = as.data.frame(tax_table(l4.phy))
safe.bc$taxa_ASV = paste(rnc.bc.new$taxa_ASV)
View(safe.bc)

# Setup object as tax_table
bc.tax = tax_table(safe.bc)
colnames(bc.tax) = colnames(safe.bc)
rownames(bc.tax) = rownames(safe.bc)
View(bc.tax)

## Update phyloseq object with new table
identical(bc.tax[1:length(colnames(seqtab.nochim.r1)),1:7], tax_table(l4.phy)) #should be true
tax_table(l4.phy) = bc.tax

###Rarefaction Curve
l4.curve <- ggrare(l4.phy, step = 1000, se = FALSE)
l4.curve <- l4.curve + theme_bw()
l4.curve
ggsave("l4_rarefaction_curve.jpeg", l4.curve, width = 15, height = 7.5, units = "cm", device = "jpeg")

###Rarefaction
set.seed(711)
l4.rare = rarefy_even_depth(l4.phy, sample.size = 0000, trimOTUs = TRUE) 
