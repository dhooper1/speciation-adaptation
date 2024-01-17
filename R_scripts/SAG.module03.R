# Load packages
library(tidyverse)
library(cowplot)
library(fastman) ##https://github.com/kaustubhad/fastman
library(normentR) ##https://danielroelfs.com/blog/how-i-make-qq-plots-using-ggplot/

##Note: You'll need to set your own working directory accordingly
setwd("/Users/danielhooper/Documents/Teaching/RGGS Workshop/R/Module03")
rm(list = ls())

# Read in data
gwas.results <- read_delim("gwas/beakcolor.H3.GT.wg.covariate.maf05.LD_prune.lmm.assoc.txt", delim = "\t", col_names = c("CHROM", "SNP", "POS", "MISS", "REF", "ALT", "MAF", "BETA", "BETA.SE", "LOGL_H1", "REMLE", "PWALD"), skip = 1)
# Add in a new variable for chromosome identity by separating the SNP column
gwas.results <- separate(gwas.results, SNP, c("chr", NA), remove = FALSE)
# How many SNPs are on each chromosome?
as.data.frame(table(gwas.results$chr))

#We might want to organize our results by chromosome in one way or another. For example,
#we may want to order our chromosomes from largest to smallest autosomes and then chrZ.
#One way that we can accomplish this is using a vector of chromosome names in our desired order.
order <- c("chr2", "chr1", "chr3", "chr1A", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr4A", "chr13", "chr14", "chr20", "chr15", "chr17", "chr18", "chr19", "chr21", "chr24", "chr26", "chr23", "chr28", "chr27", "chr22", "chr25", "chr29", "chrZ")
gwas.results$chr=match(gwas.results$chr, order)  ## 
gwas.results=gwas.results[order(gwas.results$chr),]

# Save the result as an RDS file so we don't need to do the above steps every time:
saveRDS(gwas.results, file = "gwas/beakcolor.H3.GT.wg.covariate.maf05.lmm.assoc.rds")
# Next time just read in the RDS file
gwas.results <- readRDS("gwas/beakcolor.H3.GT.wg.covariate.maf05.lmm.assoc.rds")

#If you want to subset your results in order to reduce the size of output files 
#you can filter out results that are the least significant 
gwas.reduced <- gwas.results |> filter(PWALD < 0.01)

###########################################
##MANHATTAN PLOT GENOME-WIDE GWAS RESULTS##
###########################################

#We are going to plot our GWAS results using the R package 'fastman'. This newer R package
#is faster and much easier to use than the prior go-to R package for GWAS called 'qqman'

# Create a vector of chromosome IDs in the same order as 'order' but without the "chr" prefix 
chrom.labels <- c("2", "1", "3", "1A", "4", "5", "6", "7", "8", "9", "10", "11", "12", "4A", "13", "14", "20", "15", "17", "18", "19", "21", "24", "26", "23", "28", "27", "22", "25", "29", "Z")

# Determine the genome-wide significance threshold based on the number of SNPs tested
sig.threshold <- (0.01 / length(gwas.results$SNP))

# Plot results as a Manhattan plot
fastman(gwas.results, chr="CHROM", bp="POS", p="PWALD", snp="SNP", chrlabs = chrom.labels, maxP=40, sortchr=FALSE, annotatePval = sig.threshold, annotationWinMb=1, colAbovePval = TRUE, cex = 0.3, cex.text = 0.6, annotationAngle= 60, suggestiveline = FALSE, genomewideline = FALSE)
fastman(gwas.reduced, chr="CHROM", bp="POS", p="PWALD", snp="SNP", chrlabs = chrom.labels, maxP=40, sortchr=FALSE, annotatePval = sig.threshold, annotationWinMb=1, colAbovePval = TRUE, cex = 0.3, cex.text = 0.6, annotationAngle= 60, suggestiveline = FALSE, genomewideline = FALSE)

# Examine results on a single specific chromosome in the data set
fastman(gwas.results, chr="CHROM", bp="POS", p="PWALD", snp="SNP", chrsubset = 8, chrlabs = chrom.labels, maxP=40, sortchr=FALSE, annotatePval = sig.threshold, annotationWinMb=1, colAbovePval = TRUE, cex = 0.3, cex.text = 0.01, annotationAngle= 60, suggestiveline = FALSE, genomewideline = FALSE)

###################################
##MAKE A Q-Q PLOT OF GWAS RESULTS##
###################################

#When performing GWAS it is common practice to evaluate the extent to which your p-values
#might be inflated or not. If a trait being modeled is itself correlated with genetic structure
#and this genetic structure is not accounted for when doing the GWAS, then one might recover many
#spurious associations (i.e., false positives). An essential tool for evaluating whether this might be
#occurring in your data is to make a quantile-quantile plot or 'Q-Q plot'.

#A 'Q-Q plot' plots the p-values for your set of SNPs (i.e., the OBSERVED data/distribution), ordered 
#from least significant to most-significant, against an idealized distribution going from 1 to 1/N 
#where N equals the number of SNPs tested (i.e., the EXPECTED data). 

#If your data contains no true associations the observed and expected distributions will fall along
#a 1-to-1 line but if you do have true associations you will observe a 'tail' of elevated observed
#values against the expected. If ALL or MUCH of your observed p-values are greater than your expected
#then this is evidence of demographic/structure producing erroenous associations (false positives).

#The go-to metric for evaluating Q-Q plots is 'lambda'. Lambda values ~1 indicate that there is little
#p-value inflation to worry about. Lambda values much greater than 1 indicate considerable inflation
#that needs to be accounted for.

# We can quickly make QQ-plots with fastman function 'fastqq'

# Make QQ plots to evaluate genomic inflation score
fastqq(gwas.results$PWALD, speedup=TRUE, lambda=TRUE, fix_zero=TRUE, cex=0.6, cex.axis=0.9)

##########################################
##MANHATTAN PLOT CHROMOSOME GWAS RESULTS##
##########################################

# Extract a specific chromosome to focus attention on
gwas.results.chr2 <- subset(gwas.results, CHROM == "2" & PWALD < 0.01)
gwas.results.chr8 <- subset(gwas.results, CHROM == "8" & PWALD < 0.01)
gwas.results.chrZ <- subset(gwas.results, CHROM == "chrZ" & PWALD < 0.01)

theme_set(theme_classic())

# Plot results for three chromosomes of interest

gwas.results.chr2 |> ggplot(aes(x=POS, y=-log10(PWALD))) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(sig.threshold), linetype = "dashed", colour = "red") +
  labs(x = "Chromosome 2", y = "-log10(p)")

gwas.results.chr8 |> ggplot(aes(x=POS, y=-log10(PWALD))) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(sig.threshold), linetype = "dashed", colour = "red") +
  labs(x = "Chromosome 8", y = "-log10(p)")

gwas.results.chrZ |> ggplot(aes(x=POS, y=-log10(PWALD))) +
  geom_point(size=1) +
  geom_hline(yintercept = -log10(sig.threshold), linetype = "dashed", colour = "red") +
  labs(x = "Chromosome Z", y = "-log10(p)")

###############################
##IDENTIFYING CANDIDATE GENES##
###############################

# Read in gff for chromosome of interest - note that file has already been filtered to only include genes
setwd("/Users/danielhooper/Documents/Teaching/RGGS Workshop/R/Module03/gwas/")
gff <- read.table("GCF_003957565.2_bTaeGut1.4.pri_genomic.chr8.gff", sep="\t", header=F)
# Subset and clear up the gff - add names
colnames(gff) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Arrange the gff
gff <- gff |> arrange(start, end)

# Make a gene mid point variable
gff <- gff |> mutate(mid = start + (end-start)/2)
# Make a variable that can used to plot genes
gff <- gff |> mutate(plot = 1.2*(start/start))

# Identify the 10 most significant loci on a chromosome
hits <- gwas.results.chr8 |> arrange(PWALD) |> head(n = 10)

# Find the nearest genes to our highest hit
x <- hits$POS[1]

# Find the set of genes within 150 kb of our top GWAS hit
gene_hits <- gff |> mutate(hit_dist = abs(mid - x)) |> arrange(hit_dist) |> filter(hit_dist < 150000)

# What are these genes?
gene_hits <- gene_hits |> dplyr::select(chr, start, end, attribute, hit_dist)

# Separate out the attribute column
gene_hits |> pull(attribute)
gene_hits |> pull(hit_dist)
