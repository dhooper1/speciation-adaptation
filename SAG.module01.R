# install new packages
install.packages("tidyverse") # if not already installed
install.packages("cowplot") # if not already installed

# load packages
library(tidyverse)
library(cowplot)

setwd("/Users/danielhooper/Documents/Teaching/RGGS Workshop/R/Module01")
rm(list = ls())

##Set ggplot theme
theme_set(theme_bw())

###################
##VARIANT QUALITY##
###################

var_qual <- read_delim("./summaries/cichlid_subset.lqual", delim = "\t", col_names = c("chr", "pos", "qual"), skip = 1)

a <- ggplot(var_qual, aes(qual)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)

#We can see that most sites have quite high quality scores [~1000]
#Let's set a minimum quality filter of 100 and filter more strongly on other features

#################
##VARIANT DEPTH##
#################

var_depth <- read_delim("./summaries/cichlid_subset.ldepth.mean", delim = "\t", col_names = c("chr", "pos", "mean_depth", "var_depth"), skip = 1)

b <- ggplot(var_depth, aes(mean_depth)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)

#Hmm - this plot is a bit misleading because clearly, there are very few variants with extremely 
#high coverage indeed. Let’s take a closer at the mean depth:

summary(var_depth$mean_depth)

#Since we all took different subsets, these values will likely differ slightly but clearly in 
#this case most variants have a depth of 17-20x whereas there are some extreme outliers. 
#We will redraw our plot to exclude these and get a better idea of the distribution of mean depth.

b <- b + theme_light() + xlim(0, 100)

#This gives a better idea of the distribution. We could set our minimum coverage at the 
#5 and 95% quantiles but we should keep in mind that the more reads that cover a site, 
#the higher confidence our basecall is. 10x is a good rule of thumb as a minimum cutoff 
#for read depth, although if we wanted to be conservative, we could go with 15x.

#What is more important here is that we set a good maximum depth cufoff. As the outliers show, 
#some regions clearly have extremely high coverage and this likely reflects mapping/assembly 
#errors and also paralogous or repetitive regions. We want to exclude these as they will bias 
#our analyses. Usually a good rule of thumb is something the mean depth x 2 - so in this case 
#we could set our maximum depth at 40x.

#So we will set our minimum depth to 10x and our maximum depth to 40x.

#######################
##VARIANT MISSINGNESS##
#######################

var_miss <- read_delim("./summaries/cichlid_subset.lmiss", delim = "\t", col_names = c("chr", "pos", "nchr", "nfiltered", "nmiss", "fmiss"), skip = 1)

c <- ggplot(var_miss, aes(fmiss)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)

#Our cichlid data has a very promising missingness profile - clearly most individuals have a 
#call at almost every site. Indeed if we look at the summary of the data we can see this even 
#more clearly.

summary(var_miss$fmiss)

#Most sites have almost no missing data. Although clearly, there are some (as the max value shows). 
#This means we can be quite conservative when we set our missing data threshold. We will remove 
#all sites where over 10% of individuals are missing a genotype. One thing to note here is that 
#vcftools inverts the direction of missigness, so our 10% threshold means we will tolerate 90% 
#missingness (yes this is confusing and counterintuitive… but that’s the way it is!). Typically 
#missingness of 75-95% is used.

##########################
##MINOR ALLELE FREQUENCY##
##########################

#Last of all for our per variant analyses, we will take a look at the distribution of allele 
#frequencies. This will help inform our minor-allele frequency (MAF) thresholds. As previously, 
#we read in the data:

var_freq <- read_delim("./summaries/cichlid_subset.frq", delim = "\t", col_names = c("chr", "pos", "nalleles", "nchr", "a1", "a2"), skip = 1)

#However, this is simply the allele frequencies. To find the minor allele frequency at each site, 
#we need to use a bit of dplyr based code.

# find minor allele frequency
var_freq$maf <- var_freq |> select(a1, a2) |> apply(1, function(z) min(z))

#Here we used apply on our allele frequencies to return the lowest allele frequency at each variant. 
#We then added these to our dataframe as the variable maf. Next we will plot the distribution.

d <- ggplot(var_freq, aes(maf)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3)

#The distribution might look a little odd - this is partly because of the low number of individuals
#we have in the dataset (i.e., N = 16), meaning there are only certain frequencies possible. 
#Nonetheless, it is clear that a large number of variants have low frequency alleles. We can also 
#look at the distribution in more detail:

summary(var_freq$maf)

#The upper bound of the distribution is 0.5, which makes sense because if MAF was more than this, 
#it wouldn’t be the MAF! How do we interpret MAF? It is an important measure because low MAF alleles 
#may only occur in one or two individuals. It is possible that some of these low frequency alleles 
#are in fact unreliable base calls - i.e. a source of error.

#With 16 individuals, there are 28 alleles for a given site. Therefore MAF = 0.04 is equivalent to 
#a variant occurring as one allele in a single individual (i.e. 28 * 0.04 = 1.12). Alternatively, 
#an MAF of 0.1 would mean that any allele would need to occur at least twice (i.e. 28 * 0.1 = 2.8).

#Setting MAF cutoffs is actually not that easy or straightforward. Hard MAF filtering (i.e. setting 
#a high threshold) can severely bias estimation of the site frequency spectrum and cause problems 
#with demographic analyses. Similarly, an excesss of low frequency, ‘singleton’ SNPs (i.e. only 
#occurring in one individual) can mean you keep many uniformative loci in your dataset that make 
#it hard to model things like population structure.

#Usually then, it is best practice to produce one dataset with a good MAF threshold and keep 
#another without any MAF filtering at all. For now however, we will set our MAF to 0.1

#Let's use program cowplot function plot_grid to create a single summary figure for these 
#five site based summaries

plot_grid(a, b, c, d, rows = 2)

###############################
##INDIVIDUAL BASED STATISTICS##
###############################

#As well as a our per variant statistics we generated earlier, we also calculated some individual 
#metrics too. We can look at the distribution of these to get an idea whether some of our individuals 
#have not sequenced or mapped as well as others. This is good practice to do with a new dataset. 
#A lot of these statistics can be compared to other measures generated from the data (i.e. principal 
#components as a measure of population structure) to see if they drive any apparent patterns in the data.

#############################
##MEAN DEPTH PER INDIVIDUAL##
#############################

ind_depth <- read_delim("./summaries/cichlid_subset.idepth", delim = "\t", col_names = c("ind", "nsites", "depth"), skip = 1)

e <- ggplot(ind_depth, aes(depth)) + geom_histogram(fill = "#fb6a4a", colour = "black", alpha = 0.3)

#Because we are only plotting data for 16 individuals, the plot looks a little disjointed. While 
#there is some evidence that some individuals were sequenced to a higher depth than others, 
#there are no extreme outliers. So this doesn’t suggest any issue with individual sequencing depth.

#############################################
##PROPORTION OF MISSING DATA PER INDIVIDUAL##
#############################################

ind_miss  <- read_delim("./summaries/cichlid_subset.imiss", delim = "\t", col_names = c("ind", "ndata", "nfiltered", "nmiss", "fmiss"), skip = 1)

f <- ggplot(ind_miss, aes(fmiss)) + geom_histogram(fill = "#fb6a4a", colour = "black", alpha = 0.3)

#Again this shows us, the proportion of missing data per individual is very small indeed. It ranges 
#from 0.01-0.16, so we can safely say our individuals sequenced well.

##########################################################
##HETEROZYGOSITY & INBREEDING COEFFICIENT PER INDIVIDUAL##
##########################################################

ind_het <- read_delim("./summaries/cichlid_subset.het", delim = "\t", col_names = c("ind","ho", "he", "nsites", "f"), skip = 1)

g <- ggplot(ind_het, aes(f)) + geom_histogram(fill = "#fb6a4a", colour = "black", alpha = 0.3)

#All individuals have a slightly negative inbreeding coefficient suggesting that we observed a bit 
#less heterozygote genotypes in these individuals than we would expect under Hardy-Weinberg 
#equilibrium. However, here we combined samples from four species and thus violate the assumption 
#of Hardy-Weinberg equilibrium. We would expect slightly negative inbreeding coefficients due to 
#the Wahlund-effect. Given that all individuals seem to show similar inbreeding coefficients, 
#we are happy to keep all of them. None of them shows high levels of allelic dropout (strongly 
#negative F) or DNA contamination (highly positive F).

#Let's use program cowplot function plot_grid to create a single summary figure for these 
#eight sites indivdiual based summaries

plot_grid(a, b, c, d, e, f, g, rows = 2)
