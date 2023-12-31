# load packages
library(tidyverse)
library(cowplot)

##Note: You'll need to set your own working directory accordingly
setwd("/Users/danielhooper/Documents/Teaching/RGGS Workshop/R/Module02")
rm(list = ls())

##Set ggplot theme
theme_set(theme_bw())

# read in data
pca <- read_table2("./plink/cichlids.eigenvec", col_names = TRUE)
eigenval <- scan("./plink/cichlids.eigenval")

#######################
##ANNOTATING THE DATA##
#######################

#Note that plink output a family/pop ID named FID as well as an individual ID named IID
#In our dataset these values are identical but you would typically store population level
#information as FID. We can reconstitute this same demographic information from our sample
#names as follows:

#We would like to add a species, location and if required, a species x location vector. 
#We will do this using the R version of grep. We then use paste0 to combine the columns.

# sort out the individual species and pops
# spp
spp <- rep(NA, length(pca$IID))
spp[grep("PunPund", pca$IID)] <- "pundamilia"
spp[grep("PunNyer", pca$IID)] <- "nyererei"
# location
loc <- rep(NA, length(pca$IID))
loc[grep("Mak", pca$IID)] <- "makobe"
loc[grep("Pyt", pca$IID)] <- "python"
# combine - if you want to plot each in different colours
spp_loc <- paste0(spp, "_", loc)

#With these variables created, we can remake our data.frame like so. 
#Note the use of as.tibble to ensure that we make a tibble for easy summaries etc.

# remake data.frame
pca <- as_tibble(data.frame(pca, spp, loc, spp_loc))

#####################
##PLOTTING THE DATA##
#####################

#Now that we have done our housekeeping, we have everything in place to actually 
#visualise the data properly. First we will plot the eigenvalues. It is quite 
#straightforward to translate these into percentage variance explained (although note, 
#you could just plot these raw if you wished).

# first convert to percentage variance explained
pve <- data.frame(PC = 1:16, pve = eigenval/sum(eigenval)*100)

#With that done, it is very simple to create a bar plot showing the percentage of 
#variance each principal component explains.

# make plot
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + 
  ylab("Percentage variance explained")

#Cumulatively, they explain 100% of the variance but PC1, PC2 and possible PC3 
#together explain about 30% of the variance. We could calculate this with the 
#cumsum function, like so:

# calculate the cumulative sum of the percentage variance explained
cumsum(pve$pve)

#Next we move on to actually plotting our PCA. Given the work we did earlier to 
#get our data into shape, this doesn’t take much effort at all.

# plot pca
b <- ggplot(pca, aes(PC1, PC2, col = spp, shape = loc)) + geom_point(size = 3) +
  scale_colour_manual(values = c("red", "blue")) +
  coord_equal() +
  geom_text(aes(label = IID), vjust = -1.5, size = 2.5) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)"))

c <- ggplot(pca, aes(PC1, PC3, col = spp, shape = loc)) + geom_point(size = 3) +
  scale_colour_manual(values = c("red", "blue")) +
  coord_equal() +
  geom_text(aes(label = IID), vjust = -1.5, size = 2.5) +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC3 (", signif(pve$pve[3], 3), "%)"))

#Let's use program cowplot function plot_grid to create a single summary figure for our 
#principal component analysis

plot_grid(a, b, nrow = 1)
plot_grid(b, c, nrow = 1)

######################################
##GENETIC DIFFERENTITATION LANDSCAPE##
######################################

# read in the data
fst.younger <- read_table2("./vcftools/cichlid.python_pair.weir.fst", col_names = TRUE)
fst.older <- read_table2("./vcftools/cichlid.makobe_pair.weir.fst", col_names = TRUE)

#We need to clean up and merge our two datasets a little bit to ensure that they are comparable
#We will first create a new Fst dataframe and add in results for our younger and older comparisons.
fst <- fst.younger[1:2]
fst$FST.YOUNG <- fst.younger$WEIR_AND_COCKERHAM_FST
fst$FST.OLDER <- fst.older$WEIR_AND_COCKERHAM_FST
#We next want to remove any sites where the Fst estimate wasn't possible. Sites with result "-nan".
fst <- fst |> filter(FST.YOUNG != "-nan" & FST.OLDER != "-nan")
#We lastly need to change our Fst values to be stored as numeric data now that the character data is gone
fst$FST.YOUNG <- as.numeric(fst$FST.YOUNG)
fst$FST.OLDER <- as.numeric(fst$FST.OLDER)

#Let's first check some basis summary statistics of our Fst distributions
summary(fst$FST.YOUNG)
summary(fst$FST.OLDER)

#How many SNPs are located on each chromosome? Which chromosome has the most SNPs?
fst |> count(CHROM)
fst |> count(CHROM) |> arrange(desc(n))

# plot the results as a distribution
d <- fst |> ggplot(aes(x = FST.YOUNG)) + 
  geom_histogram(binwidth = 0.05, color="black", fill = "#bcbddc") +
  xlab("Fst between younger Python river populations")
e <- fst |> ggplot(aes(x = FST.OLDER)) + 
  geom_histogram(binwidth = 0.05, color="black", fill = "#756bb1") +
  xlab("Fst between older Makobe gulf populations")

plot_grid(d, e, nrow = 2)

#How many fixed differences exist between each population compared?
fst |> filter(FST.YOUNG == 1) |> count()
fst |> filter(FST.OLDER == 1) |> count()

#What if we wanted to store all fixed differences as a new dataframe?
fixed.diffs <- fst |> filter(FST.YOUNG == 1 | FST.OLDER == 1) |> select(CHROM, POS, FST.YOUNG, FST.OLDER)

#Is there any generalizable trend to the Fst distribution in the data?
fst |> ggplot(aes(x = FST.YOUNG, y = FST.OLDER)) +
  geom_point() +
  geom_abline(slope=1, intercept=0, color = "#feb24c") +
  xlab("Fst between younger Python river populations") +
  ylab("Fst between older Makobe gulf populations")

#There are many, many methods for performing outlier tests. Some use coalescent 
#simulations, others use Bayesian inference to estimate the posterior probability 
#that a SNP is an outlier and therefore putatively under selection.

#However, the simplest way to detect outliers is to use the empirical distribution 
#of Fst to look for SNPS that exceed an arbitrary threshold of differentiation. 
#This method has been used quite a lot in the literature but it is not without 
#caveats (and more on those later). Typically the threshold is set at either the 
#95th or 99th percentile of the empirical data.

#We can use the R function quantile to identify the 95th and 99th percentile of the Fst
#distribution for each comparison - note that below we use 0.975 and 0.995 because
#this is a one-tailed test (we're not interest in the lower tail of the Fst distribution.

quantile(fst$FST.YOUNG, c(0.975, 0.995), na.rm = T)
quantile(fst$FST.OLDER, c(0.975, 0.995), na.rm = T)

#We can use the 95th percentile threshold to select a set of shared outlier sites.
#Where are the shared highly-differentiated regions (i.e., HDRs) located?
shared.HDRs <- fst |> filter(FST.YOUNG >= 0.5 & FST.OLDER > 0.6)

#We can use count() to get a summary of how many shared HDRs are on each chromosome
shared.HDRs |> count(CHROM) |> arrange(desc(n))

##Let's next focus in on one chromosome of interest that contains some shared HDRs

# plot the results for one chromosome of interest as a Manhattan plot
fst |> filter(CHROM == "chr19") |> ggplot(aes(x = POS/1e6, y = FST.YOUNG)) +
  geom_point() +
  coord_cartesian(ylim = c(-0.2, 1.0)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#feb24c") +
  labs(x = "Chromosome 19 [Mb]", y = "Fst") +
  ggtitle("Python Pundamilia Pair [younger]")
fst |> filter(CHROM == "chr19") |> ggplot(aes(x = POS/1e6, y = FST.OLDER)) +
  geom_point() +
  coord_cartesian(ylim = c(-0.2, 1.0)) +
  geom_hline(yintercept = 0.6, linetype = "dashed", color = "#feb24c") +
  labs(x = "Chromosome 19 [Mb]", y = "Fst") +
  ggtitle("Makobe Pundamilia Pair [older]")

#Would we expect that some of these shared HDR regions might be associated with genes 
#releated to the phenotypic and ecological differences between pundamilia and nyererei?
