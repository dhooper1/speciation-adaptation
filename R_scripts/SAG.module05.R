# Load packages
library(tidyverse)
library(cowplot)

##Note: You'll need to set your own working directory accordingly
setwd("/Users/danielhooper/Documents/Teaching/RGGS Workshop/R/Module05")
rm(list = ls())

# Read in data
pca <- read_table("./module05/stitch.chrZ.repeatmask.filter.PL.mac1.sag_tutorial.pca.eigenvec", col_names = TRUE)
eigenval <- scan("./module05/stitch.chrZ.repeatmask.filter.PL.mac1.sag_tutorial.pca.eigenval")
sample.info <- read_delim("./module05/stitch.chrZ.repeatmask.filter.PL.mac1.sag_tutorial.fam", delim = "\t", col_names = c("FID", "IID", "PAT", "MAT", "SEX", "PHENO"))

#######################################
##ANNOTATING THE CHROMOSOME-WIDE DATA##
#######################################

#Note that our PCA results don't include sample sex but our .fam file DOES. 
#Let's join those two sets of information into a single df we can harness for plotting
#With just a single bit of information to add this is super easy to do as follows:

pca$SEX <- sample.info$SEX

#We will next convert our eigenval vector to percentage variance explained (i.e., pve)

pve <- data.frame(PC = 1:10, pve = eigenval/sum(eigenval)*100)

#Let's take a peak at the distribution of variance explained by PC
ggplot(pve, aes(PC, pve)) + 
  geom_bar(stat = "identity") + 
  ylab("Percentage variance explained")

#Wow - it looks like PC1 is explaining a ton of variance in our data

#Calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

#####################################
##PLOTTING THE CHROMOSOME-WIDE DATA##
#####################################

#We have information about source population, sex, and PC loadings within our pca df
#Let's utilize this information to create a visually informative plot

##Set ggplot theme
theme_set(theme_bw())

# Plot PC1 against PC2 and color-code populations and shape code sexes
pca |> ggplot(aes(PC1, PC2, col = FID, shape = as.character(SEX))) + 
  geom_point(size = 3) +
  scale_color_manual(values=c("#d53e4f", "#045a8d", "#8c6bb1")) +
  geom_text(aes(label = IID), vjust = -1.5, size = 2.5) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme(legend.position="none", panel.grid = element_blank())

# Plot PC1 against PC2 and color-code only on sex
pca |> ggplot(aes(PC1, PC2, col = as.character(SEX))) + 
  geom_point(shape = 1, size = 4) +
  scale_color_manual(values=c("#045a8d", "#f768a1")) +
  geom_text(aes(label = IID), vjust = -1.5, size = 2.5) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme(legend.position="none", panel.grid = element_blank())

#What can we infer from these results?
#Let's return back to Huxley now to generate sliding-window PCA data

######################################
##ANNOTATING THE SLIDING-WINDOW DATA##
######################################

##NOTE##
#Data for this section requires you having completed the sliding window PCA section
#of the lab notebook. 

# Read in data
sliding.pca <- read_delim("./module05/sliding_window_100kb_chrZ_full_results.txt", delim = "\t", col_names = c("FID", "IID", "PC1", "PC2", "PC3", "PC4", "WINDOW"))
# Remove a window with only 6 variants [uninformative]
sliding.pca <- subset(sliding.pca, WINDOW != "29.9")

#What are our expectations?

#We seek to identify regions of SIGNAL against a background of NOISE

#SIGNAL is produced by the presence of well differentiated haplotypes like those produced
#by an inversion polymorphism segregating as homozygous INV, homozygous WT, and heterozygous
#exactly in between.

#NOISE is produced by a lack of genetic structure between samples within the region investigated

# Plot the data by IID
sliding.pca |> ggplot(aes(x = WINDOW, y = PC1, colour = IID)) +
  geom_line(alpha = 0.4) +
  labs(x = "Chromosome Z [Mb]", y = "PC1 Loading") +
  theme(legend.position="none", panel.grid = element_blank())

# Plot the data by FID
sliding.pca |> ggplot(aes(x = WINDOW, y = PC1, colour = FID)) +
  geom_line(alpha = 0.4) +
  scale_color_manual(values=c("#d53e4f", "#045a8d", "#8c6bb1")) +
  labs(x = "Chromosome Z [Mb]", y = "PC1 Loading")

# Plot the data by FID - for individual populations
sliding.pca |> filter(FID == "EB") |> ggplot(aes(x = WINDOW, y = PC1, colour = FID)) +
  geom_line(alpha = 0.4) +
  scale_color_manual(values=c("#045a8d")) +
  labs(x = "Chromosome Z [Mb]", y = "PC1 Loading")
sliding.pca |> filter(FID == "CH") |> ggplot(aes(x = WINDOW, y = PC1, colour = FID)) +
  geom_line(alpha = 0.4) +
  scale_color_manual(values=c("#d53e4f")) +
  labs(x = "Chromosome Z [Mb]", y = "PC1 Loading")
sliding.pca |> filter(FID == "SP") |> ggplot(aes(x = WINDOW, y = PC1, colour = FID)) +
  geom_line(alpha = 0.4) +
  scale_color_manual(values=c("#8c6bb1")) +
  labs(x = "Chromosome Z [Mb]", y = "PC1 Loading")

# Plot the data by IID - for individual populations
sliding.pca |> filter(FID == "EB") |> ggplot(aes(x = WINDOW, y = PC1, colour = IID)) +
  geom_line(alpha = 0.4) +
  labs(x = "Chromosome Z [Mb]", y = "PC1 Loading")
sliding.pca |> filter(FID == "CH") |> ggplot(aes(x = WINDOW, y = PC1, colour = IID)) +
  geom_line(alpha = 0.4) +
  labs(x = "Chromosome Z [Mb]", y = "PC1 Loading")
sliding.pca |> filter(FID == "SP") |> ggplot(aes(x = WINDOW, y = PC1, colour = IID)) +
  geom_line(alpha = 0.4) +
  labs(x = "Chromosome Z [Mb]", y = "PC1 Loading")

#######################
##FIRST-PASS INSIGHTS##
#######################

#What can we infer from the behavior of PC1 loading as we slide along the Z chromosome?

#Observing PC1 polarity flipping in an individual as you slide along a chromosome indicates 
#two potential phenomena - one biologically interesting and one a bit of a nuisance:
#The first, is that a recombination event has occurred and this individual carries recombined chromosomes
#The second, is a byproduct of the sign of PC loading being more or less random window to window such that
#you might observe all individuals flip from being negatively [or positively] loaded in one window to
#being positively [or negarively] loaded in an adjacent window. If all samples flip polarity, this suggests
#that this is not a real recombination event. You could 'correct' these signatures by multiplying PC1 by -1
#for these windows

#Below I've produced a few lines of code that will correct some of these nuisance flips
#Try it out and then reproduce the plots above. Do they appear any 'cleaner' now?

#Step #1: Correct wrong-sign windows
pca.subset <- subset(sliding.pca, WINDOW == "19.1" | WINDOW == "27" | WINDOW == "27.8" | WINDOW == "52.8" | WINDOW == "56.9")
pca.subset$PC1 <- -1*(pca.subset$PC1)
#Step #2: Remove wrong-sign windows from sliding.pca df
remove <- subset(sliding.pca, WINDOW != "19.1" & WINDOW != "27" & WINDOW != "27.8" & WINDOW != "52.8" & WINDOW != "56.9")
#Step #3: Concatenate the corrected sign windows with the subset df
sliding.pca <- rbind(remove, pca.subset)

########################
##SECOND-PASS INSIGHTS##
########################

#Let's reconsider our results after correcting a few nuisance sign polarity flipping windows.
#What do we see?

#First, we should observe a NOISE-like signature indicating the absence of genetic structure
#in our samples between 0.0 and 11.3 Mb, between 58.2 and 64.0 Mb, and between ~71.0 and 75.4 Mb

# This plot neatly shows this phenomenon
sliding.pca |> ggplot(aes(x = WINDOW, y = PC1, colour = IID)) +
  geom_line(alpha = 0.4) +
  labs(x = "Chromosome Z [Mb]", y = "PC1 Loading") +
  theme(legend.position="none", panel.grid = element_blank())

#Second, we should observe an extended SIGNAL of consistently differentiated PC1 loading in a set of 
#samples from hybrid zone population "SP". Notably, there is NO evidence of recombination in these
#samples between 28.8 and 58.1 Mb. Indeed, we observe THREE well-structured groups in all windows
#between these two coordinates. This is exactly consistent with an inversion across this region

# This plot neatly shows this phenomenon
sliding.pca |> filter(FID == "SP") |> ggplot(aes(x = WINDOW, y = PC1, colour = IID)) +
  geom_line(alpha = 0.4) +
  geom_vline(xintercept = 28.8, linetype = "dashed", colour = "red") +
  geom_vline(xintercept = 58.1, linetype = "dashed", colour = "red") +
  labs(x = "Chromosome Z [Mb]", y = "PC1 Loading") +
  theme(legend.position="none", panel.grid = element_blank())

#We might conclude then that there is likely to be a ~30 Mb inversion on chrZ distinguishing
#our two long-tailed finch subspecies.

#Note, while there is also evidence of fairly well-differentiated trios in other parts of the 
#chromosome, for example between 11.3 and 28.8 Mb, you should notice occasional crossing-over
#of individuals samples; suggesting recombination is occurring in this region. If we plotted
#data for more samples from the hybrid zone this would be even more obvious.

#Let's return back to Huxley now to calculate heterozygosity in these samples within and outside
#of this inverted region. This is another population genetic metric that is helpful for identifying
#inversions

##################################################################
##COMPARING HETEROZYGOSITY WITHIN AND OUTSIDE OF INVERTED REGION##
##################################################################

##NOTE##
#Data for this section requires you having completed the heterozygosity section of the lab notebook. 

# load packages
library(gghybrid) #https://github.com/ribailey/gghybrid

##Read in dataset consisting of ancestry informative autosomal, chrZ, and chrM markers
inv_dat=read.data("module05/gghybrid/inverted.region.heterozygosity_check.recode.strct_in",precol.headers=0,nprecol=2,MISSINGVAL=NA,NUMINDS=18,ONEROW=1,PLOIDY=2,MAPDISTANCE=1,MARKERNAME=1)
uninv_dat=read.data("module05/gghybrid/uninverted.region.heterozygosity_check.recode.strct_in",precol.headers=0,nprecol=2,MISSINGVAL=NA,NUMINDS=18,ONEROW=1,PLOIDY=2,MAPDISTANCE=1,MARKERNAME=1)

##Prepare data 
inv_prepdata <- data.prep(data=inv_dat$data,
                      loci=inv_dat$loci,
                      alleles=inv_dat$alleles,
                      S0=c("2"), #POPID names for the first parental reference set#
                      S1=c("1"),     #POPID names for the second parental reference set#
                      precols=inv_dat$precols,
                      ###Filtering below###
                      max.S.MAF = 0.1,               #minor allele frq must be below this value in at least one of S0 or S1; loci with smallest MAF among parental reference sets above this will be removed#
                      min.diff = 0.1,                #loci with parental allele frequency difference < this value will be removed#
                      min.allele.copies.S0 = 5,     #loci with fewer non-missing allele copies in the S0 parental reference set will be removed#
                      min.allele.copies.S1 = 5,     #in this dataset the S1 sample size is much smaller than S0, so I'm being less strict with filtering for sample size#
                      AF.CIoverlap = TRUE,          #***RECOMMENDED*** filtering option IF parental reference samples are included - removes all loci for which there is overlap between S0 and S1 in the Bayesian posterior 95% credible intervals of allele frequency# 
                      ###Filtering above###
                      return.genotype.table=T,
                      return.locus.table=T)

uninv_prepdata <- data.prep(data=uninv_dat$data,
                          loci=uninv_dat$loci,
                          alleles=uninv_dat$alleles,
                          S0=c("2"), #POPID names for the first parental reference set#
                          S1=c("1"),     #POPID names for the second parental reference set#
                          precols=uninv_dat$precols,
                          ###Filtering below###
                          max.S.MAF = 0.1,               #minor allele frq must be below this value in at least one of S0 or S1; loci with smallest MAF among parental reference sets above this will be removed#
                          min.diff = 0.1,                #loci with parental allele frequency difference < this value will be removed#
                          min.allele.copies.S0 = 5,     #loci with fewer non-missing allele copies in the S0 parental reference set will be removed#
                          min.allele.copies.S1 = 5,     #in this dataset the S1 sample size is much smaller than S0, so I'm being less strict with filtering for sample size#
                          AF.CIoverlap = TRUE,          #***RECOMMENDED*** filtering option IF parental reference samples are included - removes all loci for which there is overlap between S0 and S1 in the Bayesian posterior 95% credible intervals of allele frequency# 
                          ###Filtering above###
                          return.genotype.table=T,
                          return.locus.table=T)

#How many SNPs?
inv_prepdata$Nloci.postfilter
uninv_prepdata$Nloci.postfilter

#Calculate heterozygosity
columns = c("INDLABEL1", "POPID", "Source", "HET")
inv.df = data.frame(matrix(nrow = length(inv_prepdata$geno.data$INDLABEL), ncol = 4))
colnames(inv.df) = columns
for(i in 1:length(inv_prepdata$geno.data$INDLABEL)) {
  inv.df$INDLABEL1[i] <- inv_prepdata$geno.data$INDLABEL[i]
  inv.df$POPID[i] <- inv_prepdata$geno.data$POPID[i]
  inv.df$Source[i] <- inv_prepdata$geno.data$Source[i]
  inv.df$HET[i] <- (sum(inv_prepdata$geno.data[i,] == "1") / inv_prepdata$Nloci.postfilter)
}

uninv.df = data.frame(matrix(nrow = length(uninv_prepdata$geno.data$INDLABEL), ncol = 4))
colnames(uninv.df) = columns
for(i in 1:length(uninv_prepdata$geno.data$INDLABEL)) {
  uninv.df$INDLABEL1[i] <- uninv_prepdata$geno.data$INDLABEL[i]
  uninv.df$POPID[i] <- uninv_prepdata$geno.data$POPID[i]
  uninv.df$Source[i] <- uninv_prepdata$geno.data$Source[i]
  uninv.df$HET[i] <- (sum(uninv_prepdata$geno.data[i,] == "1") / uninv_prepdata$Nloci.postfilter)
}

#Let's add in some additional information from our chromosome-wide PCA analysis 
#to help us plot our heterozygosity analysis results
tmp <- pca |> filter(SEX == 1) |> select(IID, FID, PC1)
inv.df$PC1 <- tmp$PC1
inv.df$FID <- tmp$FID
uninv.df$PC1 <- tmp$PC1
uninv.df$FID <- tmp$FID

#Plot results
a <- inv.df |> ggplot(aes(x = -PC1, y = HET, color = FID)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values=c("#d53e4f", "#045a8d", "#8c6bb1")) +
  coord_cartesian(ylim = c(0.0, 1.0)) +
  labs(x = "Ancestry [PC1]", y = "Heterozygosity") +
  ggtitle("ChrZ Inverted Region [28 - 58 Mb]")

b <- uninv.df |> ggplot(aes(x = -PC1, y = HET, color = FID)) +
  geom_point(size = 4, alpha = 0.7) +
  scale_color_manual(values=c("#d53e4f", "#045a8d", "#8c6bb1")) +
  coord_cartesian(ylim = c(0.0, 1.0)) +
  labs(x = "Ancestry [PC1]", y = "Heterozygosity") +
  ggtitle("ChrZ Uninverted Region [0 - 11 Mb]")

plot_grid(b, a, nrow = 1)

##What can be discern from the difference in heterozygosity between groups within and outside
##of our putative inverted region
