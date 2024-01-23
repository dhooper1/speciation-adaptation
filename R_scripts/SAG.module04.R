# install new packages
install.packages("rehh")

# load packages
library(rehh) #https://cran.r-project.org/web/packages/rehh/vignettes/rehh.html
library(tidyverse)
library(cowplot)

##Note: You'll need to set your own working directory accordingly
setwd("/Users/danielhooper/Documents/Teaching/RGGS Workshop/R/Module04")
rm(list = ls())

# Read in data for each Poephila subspecies/taxon
# acuticauda WT [i.e., yellow beaks]
paa_hh <- data2haplohh(hap_file = "vcfs/hapcut2.230815.chr2.20440774-22496205.paa.INFO_SCORE_0.6.vcf.gz", polarize_vcf = FALSE)
# hecki admixed [i.e., orange beaks]
pah_o_hh <- data2haplohh(hap_file = "vcfs/hapcut2.230815.chr2.20440774-22496205.pah_o.INFO_SCORE_0.6.vcf.gz", polarize_vcf = FALSE)
# hecki WT [i.e., red beaks]
pah_r_hh <- data2haplohh(hap_file = "vcfs/hapcut2.230815.chr2.20440774-22496205.pah_r.INFO_SCORE_0.6.vcf.gz", polarize_vcf = FALSE)

#Next, we will filter our data on a minor allele frequency or MAF. 
#This is really simple in rehh with the subset function:

# Filter on MAF - here 0.02 or 2%
paa_hh_f <- subset(paa_hh, min_maf = 0.02)
pah_o_hh_f <- subset(pah_o_hh, min_maf = 0.02)
pah_r_hh_f <- subset(pah_r_hh, min_maf = 0.02)

#This will remove some sites - and we’ll be ready to run our haplotype scans.

################################################
####Performing a haplotype genome scan - iHS####
################################################

#Before we can calculate the statistics we are interested in - iHS and xpEHH - 
#we need to calculate iES statistics. Luckily this is really easy using the rehh 
#function scan_hh.

# Perform scans with MAF filter
paa_scan <- scan_hh(paa_hh_f, polarized = FALSE)
pah_o_scan <- scan_hh(pah_o_hh_f, polarized = FALSE)
pah_r_scan <- scan_hh(pah_r_hh_f, polarized = FALSE)

# perform scans without MAF filter
#paa_scan2 <- scan_hh(paa_hh, polarized = FALSE)
#pah_r_scan2 <- scan_hh(pah_r_hh, polarized = FALSE)
#pah_o_scan2 <- scan_hh(pah_o_hh, polarized = FALSE)

#Note that we once again set the polarized argument to FALSE. Next we can use 
#the output of this scan to calculate iHS. We do this with the ihh2ihs function.

# Perform iHS analyses
paa_ihs <- ihh2ihs(paa_scan, freqbin = 1, min_maf = 0.02)
pah_o_ihs <- ihh2ihs(pah_o_scan, freqbin = 1, min_maf = 0.02)
pah_r_ihs <- ihh2ihs(pah_r_scan, freqbin = 1, min_maf = 0.02)

#Note that the freqbin = 1 argument is set again because we are not using 
#polarized data (i.e. we do not know which allele is ancestral or derived). If 
#we did, rehh can apply weights to different bins of allele frequencies in order 
#to test whether there is a significant deviation in the iHS statistic. However, 
#since this isn’t the case we bin everything into a single category.

#Check whether our iHS results show evidence of value inflation using QQ plot summaries
#Note that we expect to observe a small tail of deviations below and above the expected values
#when our data includes sites that exhibit evidence of partial sweeps
distribplot(paa_ihs$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(pah_o_ihs$ihs$IHS, xlab = "iHS", qqplot = TRUE)
distribplot(pah_r_ihs$ihs$IHS, xlab = "iHS", qqplot = TRUE)

#Having now calculated and checked the iHS statistic, let’s plot the results. 
#We can either plot the statistic itself like so:

theme_set(theme_bw())

# Plot iHS results
ggplot(paa_ihs$ihs, aes(POSITION/1e6, IHS)) + 
  geom_point(alpha = 0.4) +
  labs(x = "Chromosome 2 (Mb)", y = "iHS", title = "PAA [yellow]") +
  theme(legend.position="none")
ggplot(pah_o_ihs$ihs, aes(POSITION/1e6, IHS)) + 
  geom_point(alpha = 0.4) +
  labs(x = "Chromosome 2 (Mb)", y = "iHS", title = "PAH [orange]") +
  theme(legend.position="none")
ggplot(pah_r_ihs$ihs, aes(POSITION/1e6, IHS)) + 
  geom_point(alpha = 0.4) +
  labs(x = "Chromosome 2 (Mb)", y = "iHS", title = "PAH [red]") +
  theme(legend.position="none")

#Or we can plot the log P-value to test for outliers.
#Let's decide on a significance threshold based on the number of tests
#In this case the total number of tests are the sum of the number of SNPs in
#all three populations.

sig.threshold <- (0.01 / (length(paa_hh@positions) + length(pah_o_hh@positions) + length(pah_r_hh@positions)))

# Plot p-value associated with iHS
ggplot(paa_ihs$ihs, aes(POSITION/1e6, LOGPVALUE)) + 
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = -log10(sig.threshold), linetype = "dashed", colour = "red") +
  coord_cartesian(ylim = c(0.0, 10.0)) +
  labs(x = "Chromosome 2 (Mb)", y = "log10 p-value", title = "PAA [yellow]") +
  theme(legend.position="none")
ggplot(pah_o_ihs$ihs, aes(POSITION/1e6, LOGPVALUE)) + 
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = -log10(sig.threshold), linetype = "dashed", colour = "red") +
  coord_cartesian(ylim = c(0.0, 10.0)) +
  labs(x = "Chromosome 2 (Mb)", y = "log10 p-value", title = "PAH [orange]") +
  theme(legend.position="none")
ggplot(pah_r_ihs$ihs, aes(POSITION/1e6, LOGPVALUE)) + 
  geom_point(alpha = 0.4) +
  geom_hline(yintercept = -log10(sig.threshold), linetype = "dashed", colour = "red") +
  coord_cartesian(ylim = c(0.0, 10.0)) +
  labs(x = "Chromosome 2 (Mb)", y = "log10 p-value", title = "PAH [red]") +
  theme(legend.position="none")

#Here a log P-value of ~7.5 is equivalent to P = 10^-7.5 - which is quite a conservative 
#threshold for an outlier!

##################################################
####Performing a haplotype genome scan - xpEHH####
##################################################

#Next we will calculate xpEHH which is the cross-population EHH test. This is 
#essentially a test for the probability that if we randomly sampled haplotypes 
#from different populations, we would get different haplotypes. Again, rehh makes 
#this simple with the ies2xpehh function.

#Note that because this test requires a variant to be present in BOTH populations in
#order to be included we need to be cautious about filtering our data based on MAF if
#we want to include sites that are fixed differences between populations.

#We will consider the results using our MAF dataset below. On your own, try performing these
#analyses using a dataset that has NOT had a MAF filter applied. How do the results compare?

# Perform xp-ehh with MAF filtered dataset
paa_pah_r <- ies2xpehh(paa_scan, pah_r_scan,
                       popname1 = "acuticauda", popname2 = "hecki.r",
                       include_freq = T)
paa_pah_o <- ies2xpehh(paa_scan, pah_o_scan,
                       popname1 = "acuticauda", popname2 = "hecki.o",
                       include_freq = T)
pah_pah <- ies2xpehh(pah_o_scan, pah_r_scan,
                       popname1 = "hecki.o", popname2 = "hecki.r",
                       include_freq = T)

#Here we provide the names of our previous iES scans (paa_scan and pah_scan). 
#We can also provide the function with the names of our populations and finally, 
#if we set include_freq to TRUE, we get the frequencies of alleles in our output, 
#which might be useful if we want to see how selection is acting on a particular 
#position.

#Next, we can plot the xpEHH values, like so:

# Plot xpEHH results for each population pair
a1 <- ggplot(paa_pah_r, aes(POSITION/1e6, XPEHH_acuticauda_hecki.r)) + 
  geom_point(alpha = 0.6) + 
  coord_cartesian(ylim = c(-5.0, 11.0)) +
  labs(x = "Chromosome 2 (Mb)", y = "xpEHH", title = "PAA [yellow] versus PAH [red]")
a2 <- ggplot(paa_pah_o, aes(POSITION/1e6, XPEHH_acuticauda_hecki.o)) + 
  geom_point(alpha = 0.6) + 
  coord_cartesian(ylim = c(-5.0, 11.0)) +
  labs(x = "Chromosome 2 (Mb)", y = "xpEHH", title = "PAA [yellow] versus PAH [orange]")
a3 <- ggplot(pah_pah, aes(POSITION/1e6, XPEHH_hecki.o_hecki.r)) + 
  geom_point(alpha = 0.6) + 
  coord_cartesian(ylim = c(-5.0, 11.0)) +
  labs(x = "Chromosome 2 (Mb)", y = "xpEHH", title = "PAH [orange] versus PAH [red]")


#In the plots above,NEGATIVE values suggest selection in population 2 (hecki in 
#the first case) whereas POSITIVE values indicate selection in population 1 (acuticauda in this case). 
#Alternatively, like with iHS, we could plot the log P values.

b1 <- ggplot(paa_pah_r, aes(POSITION/1e6, LOGPVALUE)) + 
  geom_point(alpha = 0.6) + 
  coord_cartesian(ylim = c(0.0, 30.0)) +
  geom_hline(yintercept = -log10(sig.threshold), linetype = "dashed", colour = "red") +
  labs(x = "Chromosome 2 (Mb)", y = "log10 p-value", title = "PAA [yellow] versus PAH [red]")
b2 <- ggplot(paa_pah_o, aes(POSITION/1e6, LOGPVALUE)) + 
  geom_point(alpha = 0.6) + 
  coord_cartesian(ylim = c(0.0, 30.0)) +
  geom_hline(yintercept = -log10(sig.threshold), linetype = "dashed", colour = "red") +
  labs(x = "Chromosome 2 (Mb)", y = "log10 p-value", title = "PAA [yellow] versus PAH [orange]")
b3 <- ggplot(pah_pah, aes(POSITION/1e6, LOGPVALUE)) + 
  geom_point(alpha = 0.6) + 
  coord_cartesian(ylim = c(0.0, 30.0)) +
  geom_hline(yintercept = -log10(sig.threshold), linetype = "dashed", colour = "red") +
  labs(x = "Chromosome 2 (Mb)", y = "log10 p-value", title = "PAH [orange] versus PAH [red]")

# Let's organize our plots in a way that makes evaluation more straight forwards

plot_grid(a1, b1, nrow = 2) #acuticauda [yellow] versus allopatric hecki [red]
plot_grid(a2, b2, nrow = 2) #acuticauda [yellow] versus admixed hecki [orange]
plot_grid(a3, b3, nrow = 2) #Admixed hecki [orange] versus allopatric hecki [red]
plot_grid(a1, a2, a3, nrow = 3)
plot_grid(b1, b2, b3, nrow = 3)

##################################################################
####Examining haplotype structure around a target of selection####
##################################################################

#One other nice feature of rehh is that we can examine haplotype structure around 
#SNPs we think might be under selection. Before we do that, we need to identify 
#the SNP in our dataset with the strongest evidence of being an xpEHH outlier.

#In our examples above, one population shows the strongest evidence of a selective sweep
#That population is the acuticauda [yellow] population. What evidence do we have for this?

# Find the SNPs with highest -log10 p-values
hit <- paa_pah_o |> arrange(desc(LOGPVALUE)) |> top_n(10)
# Get SNP position
x <- hit$POSITION[1] #Select top hit

#Here we also set the position of our putative selection SNP as the object x. 
#This is because we need to identify where it occurs in our haplotype objects - 
#unfortunately we cannot use the position for this. In the code below, we find 
#the marker id for both our datasets.

# No MAF filter
marker_id_paa <- which(paa_hh@positions == x)
marker_id_pah.r <- which(pah_r_hh@positions == x)
marker_id_pah.o <- which(pah_o_hh@positions == x)

##Now we are ready to plot the bifurcation of haplotypes around our site of selection. 
##We do this like so:

# No MAF filter
paa_furcation <- calc_furcation(paa_hh, mrk = marker_id_paa)
pah_r_furcation <- calc_furcation(pah_r_hh, mrk = marker_id_pah.r)
pah_o_furcation <- calc_furcation(pah_o_hh, mrk = marker_id_pah.o)

#We can now plot each furcation to have a look at haplotype lengths:

par(mfrow=c(1,3)) ##[number of rows by number of columns]
plot(paa_furcation, xlim = c(21.42E+6, 21.49E+6))
plot(pah_o_furcation, xlim = c(21.42E+6, 21.49E+6))
plot(pah_r_furcation, xlim = c(21.42E+6, 21.49E+6))

#Calculating the furcation pattern also makes it possible to calculate the 
#haplotype length around our signature of selection.

paa_haplen <- calc_haplen(paa_furcation)
pah_r_haplen <- calc_haplen(pah_r_furcation)
pah_o_haplen <- calc_haplen(pah_o_furcation)

#With the haplotype length calculated, we can now plot this to see how haplotype 
#structure differs between each of our populations.

par(mfrow=c(1,3)) ##[number of rows by number of columns]
plot(paa_haplen, xlim = c(21.42E+6, 21.49E+6))
plot(pah_o_haplen, xlim = c(21.42E+6, 21.49E+6))
plot(pah_r_haplen, xlim = c(21.42E+6, 21.49E+6))

#Here we can see the red haplotype is much larger around this target and is 
#also far more frequent in the acuticauda than in hecki.

#computing EHH statistics for the focal SNP with name labeled above, i.e., "marker_id_paa"
#which displays a strong signal of selection
res1 <- calc_ehh(paa_hh, 
                 mrk = marker_id_paa, 
                 include_nhaplo = TRUE)
plot(res1, xlim = c(21.43E+6, 21.48E+6))
res1$ihh #Get values for the length of haplotypes carrying the ancestral (A) or derived (D) allele

res2 <- calc_ehh(pah_o_hh, 
                 mrk = marker_id_pah.o, 
                 include_nhaplo = TRUE)
plot(res2, xlim = c(21.43E+6, 21.48E+6))
res2$ihh #Get values for the length of haplotypes carrying the ancestral (A) or derived (D) allele

res3 <- calc_ehh(pah_r_hh, 
                 mrk = marker_id_pah.r, 
                 include_nhaplo = TRUE)
plot(res3, xlim = c(21.43E+6, 21.48E+6))
res3$ihh #Get values for the length of haplotypes carrying the ancestral (A) or derived (D) allele

#We can also visualize the haplotype structure directly
#These plots will look similar to the cartoon model examples of a selective sweep

hh_paa_subset <- subset(paa_hh, select.mrk = 25373:25573)
hh_pah_o_subset <- subset(pah_o_hh, select.mrk = 26397:26597)
hh_pah_r_subset <- subset(pah_r_hh, select.mrk = 25114:25314)

par(mfrow=c(3,1)) ##[number of rows by number of columns]
plot(
  hh_paa_subset,
  mrk = 101,
  group_by_allele = TRUE,
  ignore.distance = TRUE,
  col = c(NA, "red"),
  linecol = c("lightblue", "lightpink"),
  mrk.col = "black",
  cex = 0.1,
  pos.lab.hap = "none",
  pos.lab.mrk = "none"
)

plot(
  hh_pah_o_subset,
  mrk = 101,
  group_by_allele = TRUE,
  ignore.distance = TRUE,
  col = c(NA, "red"),
  linecol = c("lightblue", "lightpink"),
  mrk.col = "black",
  cex = 0.1,
  pos.lab.hap = "none",
  pos.lab.mrk = "none"
)

plot(
  hh_pah_r_subset,
  mrk = 101,
  group_by_allele = TRUE,
  ignore.distance = TRUE,
  col = c(NA, "red"),
  linecol = c("lightblue", "lightpink"),
  mrk.col = "black",
  cex = 0.1,
  pos.lab.hap = "none",
  pos.lab.mrk = "none"
)

#Plot summary figure for top hits
par(mfrow=c(3,3)) ##[number of rows by number of columns]

##Identifying candidate genes

# Read in gff for chromosome of interest - note that file has already been filtered to only include genes
setwd("/Users/danielhooper/Documents/Teaching/RGGS Workshop/R/Module03/gwas/")
gff <- read.table("GCF_003957565.2_bTaeGut1.4.pri_genomic.chr2.gff", sep="\t", header=F)

# Subset and clear up the gff - add names
colnames(gff) <- c("chr", "source", "feature", "start", "end", "score", "strand", "frame", "attribute")

# Remove lncRNAs
gff <- gff |> filter(!grepl("lncRNA", attribute))

# Arrange the gff
gff <- gff |> arrange(start, end)

# Make a gene mid point variable
gff <- gff |> mutate(mid = start + (end-start)/2)
# Make a variable that can used to plot genes
gff <- gff |> mutate(plot = 1.2*(start/start))

# Identify the highest peak of selection
hits <- paa_pah_o |> arrange(desc(LOGPVALUE)) |> top_n(10)
# Find the nearest genes to our highest hit
x <- hits$POSITION[1]

# Find hits within 150 Kb
gene_hits <- gff |> mutate(hit_dist = abs(mid - x)) |> arrange(hit_dist) |> filter(hit_dist < 150000)
# What are these genes?
gene_hits <- gene_hits |> select(chr, start, end, attribute, hit_dist)
# Separate out the attribute column
gene_hits |> pull(attribute) #Output attributes of each gene ordered from closest to furthest from selected SNP
gene_hits |> pull(hit_dist) #Output distance of each gene to the selected SNP

##Which gene is closest to the putative selective sweep? What does this gene do?
##This website is a great resource for looking up genes to examine what they do and what genes they interact with: https://string-db.org
##Doing a literature search for the gene name is another good approach

