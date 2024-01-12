# install new packages
install.packages("devtools")
devtools::install_github("ribailey/gghybrid")

# load packages
library(tidyverse)
library(gghybrid) #https://github.com/ribailey/gghybrid
library(coda)
library(cowplot)

##Note: You'll need to set your own working directory accordingly
setwd("/Users/danielhooper/Documents/Teaching/RGGS Workshop/R/Module02")
rm(list = ls())

##Read in dataset consisting of ancestry informative autosomal, chrZ, and chrM markers
dat=read.data("gghybrid/all_chrom.repeatmask.filter.AIMs.wild.LD_prune.strct_in",precol.headers=0,nprecol=2,MISSINGVAL=NA,NUMINDS=983,ONEROW=1,PLOIDY=2,MAPDISTANCE=1,MARKERNAME=1)
dat

#Run data.prep including the different filtering options. Try hashing out or changing these 
#options to see how it affects the number of loci in the resulting analysis table (see 
#prepdata$Nloci.postfilter in the data.prep object).

prepdata <- data.prep(data=dat$data,
                      loci=dat$loci,
                      alleles=dat$alleles,
                      S0=c("BR","MH", "EL", "NL", "SV"), #POPID names for the first parental reference set#
                      S1=c("MY", "PR", "LS", "OC", "MA", "CH", "LH"),     #POPID names for the second parental reference set#
                      precols=dat$precols,
                      ###Filtering below###
                      max.S.MAF = 0.1,               #minor allele frq must be below this value in at least one of S0 or S1; loci with smallest MAF among parental reference sets above this will be removed#
                      min.diff = 0.8,                #loci with parental allele frequency difference < this value will be removed#
                      min.allele.copies.S0 = 50,     #loci with fewer non-missing allele copies in the S0 parental reference set will be removed#
                      min.allele.copies.S1 = 50,     #in this dataset the S1 sample size is much smaller than S0, so I'm being less strict with filtering for sample size#
                      AF.CIoverlap = TRUE,          #***RECOMMENDED*** filtering option IF parental reference samples are included - removes all loci for which there is overlap between S0 and S1 in the Bayesian posterior 95% credible intervals of allele frequency# 
                      ###Filtering above###
                      return.genotype.table=T,
                      return.locus.table=T)

#How many SNPs?
prepdata$Nloci.postfilter

#############################################################################
#############################################################################
#Running hybrid index estimation#############################################
#############################################################################
#############################################################################

?esth

###############################################
#In the presence of parental reference samples#
###############################################

hindlabel=esth(
  data.prep.object=prepdata$data.prep,
  read.data.precols=dat$precols,
  include.Source=TRUE,	#Leave at default TRUE if you want hybrid indices for the parental reference individuals, which is often useful#
  plot.ind = c("SP21","SP26","SP69","SP06","SP42","SP51"),  #Optionally plot some individuals in real time. Merely shows how well the adaptive burnin is working#
  plot.col = c("#FAB255", "#8CB271", "#38A68A", "#19859C", "#616A71", "#DD5129"),
  nitt=1000,    #Testing suggests nitt=3000,burnin=1000 are sufficient for accurate posterior estimation#
  burnin=500
)

##Save the result as an RDS file so we don't need to do the above step every time:
saveRDS(hindlabel, file = "gghybrid/all.repeatmask.filter.AIMs.wild.LD_prune.HI.rds")
##Next time just read in the RDS file
hindlabel <- readRDS("gghybrid/all.repeatmask.filter.AIMs.wild.LD_prune.HI.rds")

#There will be a warning if any NAs are produced for the best posterior estimate 
#(h_posterior_mode).

setkey(hindlabel$hi,beta_mean)#beta_mean is the mean of the posterior and is always >0 and <1, which is sometimes useful, but the mode is the best posterior estimate in the presence of a prior#
hindlabel#Run this line twice if hindlable$hi doesn't show up the first time#
hist(hindlabel$hi$h_posterior_mode, breaks = 50)

#If you want to see all the results#

View(hindlabel$hi)
write.csv(hindlabel$hi, file = "gghybrid/all.repeatmask.filter.AIMs.wild.LD_prune.HI.csv", row.names = FALSE)

#It is possible to get a result of NA or NaN for the posterior hybrid index estimate 
#(h_posterior_mode). If this happens it probably means the variance of the proposal 
#distribution during burnin was too high. The default is to use the variance 
#(1/(n_allele_copies_per_test_subject/2)/10 after the first 100 iterations. If this fails 
#for any individuals, add another zero and set the option 'init.var2' to this value. And/or 
#use the MCMC plots within esth to try and identify the problem. I often find that if I run an 
#individual a second time even without changing the settings, it works fine.

#nitt=3000, burnin=1000 should be sufficient for convergence of the posterior estimates across 
#multiple runs according to the Gelman-Rubin diagnostic, which should have a value < 1.2 
#to indicate convergence. Both nitt and burnin can be increased if the G-R diagnostic indicates 
#some individuals with poor convergence across runs. It's unlikely more than nitt=5000, 
#burnin=2000 would be necessary, and even nitt=2000 is very likely to be sufficient.

#################################################################
####################GEOGRAPHIC CLINE OF HYBRID INDEX#############
#################################################################

##Read in a csv of results previously computed by gghybrid function esth
hindlabel.complete <- read_delim("gghybrid/all.repeatmask.filter.AIMs.wild.LD_prune.HI.csv", delim = ",", col_names = c("SOURCE", "IID", "POPID", "H_POSTERIOR", "H_LOWER", "H_UPPER", "BETA_MEAN", "BETA_VAR", "SHAPE1", "SHAPE2"), skip = 1)
###OR###
##Create a new df of output from gghybrid function esth
hindlabel.complete <- hindlabel$hi

##Read in a csv containing sampling longitude and latitude for all examined samples
source.info <- read_delim("gghybrid/all_chrom.repeatmask.filter.AIMs.wild.LD_prune.sample_sites.csv", delim = ",", col_names = c("IID", "LONGITUDE", "LATITUDE"), skip = 1)

##Add latitude and longitude information to output from gghybrid function esth
hindlabel.complete$LAT <- source.info$LATITUDE
hindlabel.complete$LONG <- source.info$LONGITUDE

theme_set(theme_bw())

ggplot(hindlabel.complete, aes(x=LONG, y=h_posterior_mode, color=h_posterior_mode)) +
  scale_color_gradient2(midpoint = 0.5, mid = "#2ca25f", high = "#b2182b", low = "#2166ac") +
  theme(legend.position="none") + 
  geom_point(size=2, alpha = 0.6) + 
  coord_cartesian(xlim = c(122.5, 140), ylim = c(00, 1.0)) + 
  labs(x = "LONGITUDE", y = "Hybrid Index")

ggplot(hindlabel.complete, aes(x=LONG, y=LAT, color=h_posterior_mode)) +
  geom_point()

#Just how confident can we be in our inference of hybrid index?
#Run esth twice more on the same data and carry out the Gelman-Rubin diagnostic test
#This test allows us to evaluate how confident we can be in our hybrid index scores
#In practice, the G-R diagnostic evaluates convergence of independent MCMC runs to the same result

hindlabel2=esth(data.prep.object=prepdata$data.prep,read.data.precols=dat$precols,include.Source=TRUE,nitt=1000,burnin=500)
setkey(hindlabel2$hi,beta_mean)
hindlabel3=esth(data.prep.object=prepdata$data.prep,read.data.precols=dat$precols,include.Source=TRUE,nitt=1000,burnin=500)
setkey(hindlabel3$hi,beta_mean)

#Join the three results objects#

setkey(hindlabel$hi,INDLABEL);setkey(hindlabel2$hi,INDLABEL);setkey(hindlabel3$hi,INDLABEL)
hindall=hindlabel$hi[,.(Source,INDLABEL,POPID,beta_shape1,beta_shape2)][hindlabel2$hi[,.(INDLABEL,beta_shape1,beta_shape2)]][hindlabel3$hi[,.(INDLABEL,beta_shape1,beta_shape2)]]

hindall[,rn:=seq(1,.N)];

#Take a sample of N=10000 from the posterior of each run for each individual, logit transform 
#(the G-R diagnostic assumes normally distributed posteriors), and calculate the Gelman-Rubin 
#diagnostic per individual. This diagnostic evaluates convergence in hi estimation across runs.

library(coda)

hindall[,h_gelman:=coda::gelman.diag(
  mcmc.list(
    mcmc(qlogis(rbeta(10000,beta_shape1,beta_shape2))), 
    mcmc(qlogis(rbeta(10000,i.beta_shape1,i.beta_shape2))),
    mcmc(qlogis(rbeta(10000,i.beta_shape1.1,i.beta_shape2.1)))
  )
)$psrf[1],
by=rn];

#Plot. Ideally all values should be < 1.2#

hist(hindall$h_gelman, breaks = 100)

#Take a look at which are most extreme. With this dataset it's always those with hi estimates 
#very close to 0 or 1. This may be an issue with the logit-transform of these beta-distributed 
#posterior samples being asymmetrical and therefore breaking the assumptions of the G-R 
#diagnostic, rather than a problem with convergence.

hindall[order(h_gelman)]

#For example, plot a histogram of the logit-transformed posterior samples from the individual 
#with the most extreme G-R value. Adjust sample number below based on N of sample-set#

setkey(hindall,h_gelman);hindall
hist(hindall[983,qlogis(rbeta(10000,beta_shape1,beta_shape2))])

plot(hindlabel$hi$h_posterior_mode,hindlabel2$hi$h_posterior_mode)

#############################################################################
#############################################################################
#Quantify ancestry informative marker heterozygosity#########################
#############################################################################
#############################################################################

hi <- read.csv("gghybrid/all.repeatmask.filter.AIMs.wild.LD_prune.HI.csv", header = TRUE)
##OR##
hi <- hindlabel$hi
##NOTE - SAMPLE IDs MUST BE ORDERED BY INDLABEL 
hi <- hi[order(hi$INDLABEL), ]

#sum(prepdata$geno.data[1,] == "1") / prepdata$Nloci.postfilter

columns = c("INDLABEL1", "POPID", "Source", "HI", "HET")
df = data.frame(matrix(nrow = length(prepdata$geno.data$INDLABEL), ncol = 5))
colnames(df) = columns
for(i in 1:length(prepdata$geno.data$INDLABEL)) {
  df$INDLABEL1[i] <- prepdata$geno.data$INDLABEL[i]
  df$POPID[i] <- prepdata$geno.data$POPID[i]
  df$Source[i] <- prepdata$geno.data$Source[i]
  df$HI[i] <- hi$h_posterior_mode[i]
  df$HET[i] <- (sum(prepdata$geno.data[i,] == "1") / prepdata$Nloci.postfilter)
}

write.csv(df, file = "gghybrid/all.repeatmask.filter.AIMs.wild.LD_prune.HI_Ho.csv", row.names = FALSE)
df <- read.csv("gghybrid/all.repeatmask.filter.AIMs.wild.LD_prune.HI_Ho.csv", header = TRUE)

#Quick triangle plot
plot(df$HI,df$HET)
#Add a title and axis labels.
title(xlab="Hybrid index",ylab="Heterozygosity",cex.main=1.5,cex.lab=1.5)

#Make a somewhat nicer triangle plot of hybrid index against heterozygosity with ggplot
a <- ggplot(df, aes(x=HI, y=HET, color=Source)) + 
  theme(panel.grid = element_blank()) +
  geom_abline(intercept = 0, slope = 2, color="grey", size=0.5) + 
  geom_abline(intercept = 2, slope = -2, color="grey", size=0.5) + 
  geom_hline(yintercept=1.0) + geom_hline(yintercept=0.0) + 
  geom_vline(xintercept=0.0) + geom_vline(xintercept=1.0)
a <- a + geom_point(size=3, alpha=0.2) + 
  labs(x = "Hybrid Index", y = "Heterozygosity") + 
  theme(legend.position="none") + 
  coord_cartesian(ylim = c(0.0, 1.0))
a + scale_color_manual(values=c("#0571b0", "#ca0020", "#38A68A"))

#Hybrid index and heterozygosity histograms for non-parental populations
df2 <- subset(df, Source == "TEST")

b <- ggplot(df2, aes(x=HI, y=HET, color=Source)) + 
  theme(panel.grid = element_blank()) +
  geom_point(size=3, alpha=0.4) + 
  labs(x = "Hybrid Index", y = "Heterozygosity") +
  theme(legend.position="none") + 
  coord_cartesian(ylim = c(0.0, 1.0)) +
  scale_colour_manual(values=c("#38A68A"))

h2 <- ggplot(df2, aes(x=HI, color=Source, fill=Source)) + 
  geom_histogram(binwidth=0.05, position="dodge") + 
  theme(legend.position="none", panel.grid = element_blank()) + 
  scale_color_manual(values=c("#38A68A")) + scale_fill_manual(values=c("#38A68A"))

h3 <- ggplot(df2, aes(x=HET, color=Source, fill=Source)) + 
  geom_histogram(binwidth=0.05, position="dodge") + 
  theme(legend.position="none", panel.grid = element_blank()) + 
  scale_color_manual(values=c("#38A68A")) + scale_fill_manual(values=c("#38A68A"))


##########################
##########################
#Genomic cline estimation#
##########################
##########################

#Genomic cline analysis estimates 2 parameters, v and centre:

#centre = the genome-wide hybrid index at which the focal test.subject allele frequency is 
#halfway between those of the parental reference sets.

#v = cline steepness. It would be 1/cline width IF this were geographic cline analysis 
#(which has an unconstrained x axis of distance along a transect). With the x axis constrained 
#to [0,1], i.e. the set of genome-wide hybrid index estimates, v is not exactly 1/width but 
#remains a measure of cline steepness.

?ggcline

#See above for preparing all the objects needed for genomic cline analysis.

#The full set of arguments for ggcline including all defaults (comments indicate which entries 
#are not defaults).

tmp=ggcline(
  data.prep.object=prepdata$data.prep,    #Needs an entry#
  esth.object=hindlabel,                  #Needs an entry#
  esth.colname="h_posterior_mode",
  test.subject = "locus",
  poolv = FALSE,
  poolcentre = FALSE,
  include.Source = TRUE,                  #Default is FALSE#
  return.likmeans = TRUE,                 #Default is FALSE#
  read.data.precols=dat$precols,          #Needs an entry#
  fix.subject.v = FALSE,
  fix.value.v,
  fix.subject.centre = FALSE,
  fix.value.centre,
  plot.test.subject = c("chr8:3144828","chrZ:67498406"), #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.col = c("#fdae61","#91bfdb"),      #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.ylim = c(-3, 6),                   #Remove plots to speed up runtime, they are just to check the MCMC is working#
  plot.pch.v.centre = c(1, 3),            #Remove plots to speed up runtime, they are just to check the MCMC is working#
  prior.logitcentre = c(0, sqrt(50)),     #Default#
  prior.logv = c(0, sqrt(10)),            #Default#
  nitt=1000,                               #Needs an entry, 5000 should be sufficient for standard per-locus estimates#
  burnin=500,                             #Needs an entry, 2000 should be sufficient for standard per-locus estimates#
  start.v = NULL,
  start.centre = NULL,
  init.var.v = NULL,
  init.var.centre = NULL,
  init.cov.vcentre = NULL,
  print.k = 100
)

##Rate parameter \code{v} is plotted with OPEN CIRCLES, and logit(\code{centre}) with PLUS SIGNS, 
##both the same colour for the same test subject.

#The main output, gc$gc, includes rounded parameter estimates on the original scale 
#(exp_mean_log_v and invlogit_mean_logit_centre) and their 95% credible intervals, and the 
#(not rounded) posterior parameter means (mean_log_v and mean_logit_centre) and multivariate 
#normally distributed covariance matrix on the latent scale (var_log_v, var_logit_centre, 
#cov_log_v_logit_centre).

gc

##Save the result as an RDS file so we don't need to do the above step every time:
saveRDS(gc, file = "gghybrid/all.repeatmask.filter.AIMs.wild.LD_prune.GC.rds")
##Next time just read in the RDS file
gc <- readRDS("gghybrid/all.repeatmask.filter.AIMs.wild.LD_prune.GC.rds")

#ggcline includes Bayesian p values for the cline parameters but does not include false 
#discovery rate or other adjustments, which are left to the user. There are several R packages 
#for FDR calculations.

#Be aware that occasionally p is estimated as exactly 0 due to very narrow credible intervals, 
#and this may cause FDR calculations to fail.

#If you want to see all the results
setkey(gc$gc,centre_pvalue)
View(gc$gc)
write.csv(gc$gc, file = "gghybrid/all_chrom.repeatmask.filter.AIMs.wild.LD_prune.GC.csv", row.names = FALSE)
##Or read in results from a previously finished run
gc <- readRDS("gghybrid/all.repeatmask.filter.AIMs.wild.LD_prune.GC.rds")

##Let's plot the distribution of cline centre and steepness results
hist(gc$gc$exp_mean_log_v, breaks = 50, xlab = "v - cline steepness", main = "Histogram of cline steepness across loci")
hist(gc$gc$invlogit_mean_logit_centre, breaks = 50, xlab = "c - cline centre", main = "Histogram of cline centre across loci")
plot(gc$gc$exp_mean_log_v, gc$gc$invlogit_mean_logit_centre, ylab = "c - cline centre", xlab = "v - cline steepness")

#ggplot version
theme_set(theme_bw())
ggplot(gc$gc, aes(x=exp_mean_log_v, y=invlogit_mean_logit_centre)) + 
  geom_point(aes(size = 1/centre_pvalue), alpha = 0.85) + coord_cartesian(ylim=c(0, 1.0)) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  geom_vline(xintercept=1.0, linetype="dashed", color = "black") +
  labs(x = "cline steepness", y = "cline centre")

##What if we would like to evaluate our results based on chromosome identity?
##Let's first install/load a neat R color palette package named MetBrewer to help us out

install.packages("MetBrewer")
library(MetBrewer) ##https://github.com/BlakeRMills/MetBrewer
chrom_colors <- met.brewer(name="Hiroshige", n=9, type="continuous", direction=-1)

##Let's create a new dataframe for plotting with locus information for chromosome and SNP position
plot.gc <- gc$gc
plot.gc <- separate(plot.gc, locus, into = c("chrom", "pos"), sep = ":", remove = FALSE)

##We can combine this information to plot in a way that allows us to color-code by chromosome
plot.gc |> ggplot(aes(x=exp_mean_log_v, y=invlogit_mean_logit_centre, colour=chrom)) + 
  coord_cartesian(ylim=c(0, 1.0)) +
  geom_point(aes(fill=chrom), colour="white", pch=21, size=4) +
  scale_fill_manual(values=chrom_colors) +
  geom_hline(yintercept=0.5, linetype="dashed", color = "black") +
  geom_vline(xintercept=1.0, linetype="dashed", color = "black") +
  labs(x = "cline steepness", y = "cline centre")

#########################
#########################
#Plotting genomic clines#
#########################
#########################

?plot_clinecurve

##Plot genomic clines for representative SNPs from two BILL COLOUR associated regions on the same plot
##Figure 4B for bill colour manuscript
plot_clinecurve(
  ggcline.object=gc$gc,
  cline.locus=c("chr8:3144828", "chrZ:67498406"),
  locus.column="locus",
  cline.col=c("#fdae61", "#91bfdb"),
  cline.centre.line=c("chr8:3144828", "chrZ:67498406"),
  cline.centre.col=c("#fdae61", "#91bfdb")
)

##Plot some null curve results
plot_clinecurve(
  ggcline.object=gc$gc,
  cline.locus=random,
  locus.column="locus",
  cline.col="#bdbdbd",
  cline.centre.line=c("chr8:3144828", "chrZ:67498406"),
  cline.centre.col=c("#fdae61", "#91bfdb")
)

#Add a title and axis labels.
title(xlab="Hybrid index",ylab="Locus allele frequency",cex.main=1,cex.lab=1)


###################################################
#Add sampled cline curves to represent uncertainty#
###################################################

#There is no option for shading credible intervals within plot_clinecurve, so here is a suggestion for how to add some: 
#sample parameter values from the joint posterior, and add a grey curve to the plot for each posterior sample. 
#This will overlay the original cline curve plotted above, so it needs to be plotted again.

#First generate random samples from the posterior using gghybrid's rtmvnormDT3 function.

gcsamp=gc$gc[locus%in%c("chr8:3144828"),
             rtmvnormDT3(
               nsamp=1000,                                #default is 1 sample per locus#
               meanlogv=mean_log_v,                       #no default, column must be specified#
               meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
               varlogv=var_log_v,                         #no default, column must be specified#
               varlogitcentre=var_logit_centre,           #no default, column must be specified#
               covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
               lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
               upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
             ),                                           #end of rtmvnormDT3 function#
             by="locus"];                                 #indicate a grouping variable for the sampling#

#Add the best posterior estimates to gcsamp as the final row, so it's not obscured by the mass of grey lines.
locipost=gc$gc[locus=="chr8:3144828",.(locus,mean_log_v,mean_logit_centre)];setnames(locipost,c("locus","log_v","logit_centre"));
gcsamp=rbind(gcsamp,locipost)

#Add the individual curves to the plot in a 'for' loop.
for(i in 1:nrow(gcsamp)){
  v = gcsamp[i,exp(log_v)]
  u = gcsamp[i,logit_centre*exp(log_v)]
  S1.prop_1 = gc$gc[locus=="chr8:3144828",S1.prop_1];
  S0.prop_1 = gc$gc[locus=="chr8:3144828",S0.prop_1];
  
  par(new=T);
  
  if(i < nrow(gcsamp)){
    curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
          from=0,to=1,axes=F,xlab="",ylab="", col="#fed976", alpha=0.1, lwd=0.5,ylim=c(0,1));
  }else{                                                                            #Add the best-fitting curve to the plot again, so it's not obscured#
    curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
          from=0,to=1,axes=F,xlab="",ylab="", col="#fdae61",lwd=3,ylim=c(0,1));
  };
};#end of for loop#

####
###NEXT ADD THE TTC39B RESULTS
###

#First generate random samples from the posterior using gghybrid's rtmvnormDT3 function.

gcsamp=gc$gc[locus%in%c("chrZ:67498406"),
             rtmvnormDT3(
               nsamp=1000,                                #default is 1 sample per locus#
               meanlogv=mean_log_v,                       #no default, column must be specified#
               meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
               varlogv=var_log_v,                         #no default, column must be specified#
               varlogitcentre=var_logit_centre,           #no default, column must be specified#
               covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
               lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
               upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
             ),                                           #end of rtmvnormDT3 function#
             by="locus"];                                 #indicate a grouping variable for the sampling#

#Add the best posterior estimates to gcsamp as the final row, so it's not obscured by the mass of grey lines.
locipost=gc$gc[locus=="chrZ:67498406",.(locus,mean_log_v,mean_logit_centre)];setnames(locipost,c("locus","log_v","logit_centre"));
gcsamp=rbind(gcsamp,locipost)

#Add the individual curves to the plot in a 'for' loop.
for(i in 1:nrow(gcsamp)){
  v = gcsamp[i,exp(log_v)]
  u = gcsamp[i,logit_centre*exp(log_v)]
  S1.prop_1 = gc$gc[locus=="chrZ:67498406",S1.prop_1];
  S0.prop_1 = gc$gc[locus=="chrZ:67498406",S0.prop_1];
  
  par(new=T);
  
  if(i < nrow(gcsamp)){
    curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
          from=0,to=1,axes=F,xlab="",ylab="", col="#c6dbef",alpha=0.1,lwd=0.5,ylim=c(0,1));
  }else{                                                                            #Add the best-fitting curve to the plot again, so it's not obscured#
    curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
          from=0,to=1,axes=F,xlab="",ylab="", col="#91bfdb",lwd=3,ylim=c(0,1));
  };
};#end of for loop#

###Create SNP specific supplemental figures with sample GTs

####
####CYP2J19
####

plot_clinecurve(
  ggcline.object=gc$gc,
  cline.locus="chr8:3144828",
  locus.column="locus",
  cline.col="#fdae61",
  cline.centre.line="chr8:3144828",
  cline.centre.col="#fdae61"
)

#Add a title and axis labels.
title(main = "chr8:3144828 [CYP2J19]",xlab="Hybrid index",ylab="Locus allele frequency",cex.main=1,cex.lab=1)

gcsamp=gc$gc[locus%in%c("chr8:3144828"),
             rtmvnormDT3(
               nsamp=1000,                                #default is 1 sample per locus#
               meanlogv=mean_log_v,                       #no default, column must be specified#
               meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
               varlogv=var_log_v,                         #no default, column must be specified#
               varlogitcentre=var_logit_centre,           #no default, column must be specified#
               covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
               lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
               upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
             ),                                           #end of rtmvnormDT3 function#
             by="locus"];                                 #indicate a grouping variable for the sampling#

#Add the best posterior estimates to gcsamp as the final row, so it's not obscured by the mass of grey lines.
locipost=gc$gc[locus=="chr8:3144828",.(locus,mean_log_v,mean_logit_centre)];setnames(locipost,c("locus","log_v","logit_centre"));
gcsamp=rbind(gcsamp,locipost)

#Add the individual curves to the plot in a 'for' loop.
for(i in 1:nrow(gcsamp)){
  v = gcsamp[i,exp(log_v)]
  u = gcsamp[i,logit_centre*exp(log_v)]
  S1.prop_1 = gc$gc[locus=="chr8:3144828",S1.prop_1];
  S0.prop_1 = gc$gc[locus=="chr8:3144828",S0.prop_1];
  
  par(new=T);
  
  if(i < nrow(gcsamp)){
    curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
          from=0,to=1,axes=F,xlab="",ylab="", col="#969696",lwd=0.5,ylim=c(0,1));
  }else{                                                                            #Add the best-fitting curve to the plot again, so it's not obscured#
    curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
          from=0,to=1,axes=F,xlab="",ylab="", col="#fdae61",lwd=3,ylim=c(0,1));
  };
};#end of for loop#

######################################################
#Add data to the plot as genotypes scaled from 0 to 1#
######################################################

ploidy=2

#Create a data.table containing individual reference, hybrid index, locus name and genotype. See the ggcline examples for creating prepdata and hindlabel.
setkey(prepdata$data.prep,INDLABEL);setkey(hindlabel$hi,INDLABEL);
genodat=hindlabel$hi[,.(INDLABEL,h_posterior_mode)][prepdata$data.prep[locus=="chr8:3144828",sum(Source_allele)/ploidy,by=c("INDLABEL","locus")]]

#Add the resulting points to the cline plot.
points(genodat[,V1]~genodat[,h_posterior_mode],pch=16,cex=0.7,col="grey")
###
dev.off()#END OF PLOT CODE FOR LOCUS CYP2J19************************************************


####
####TTC39B
####

plot_clinecurve(
  ggcline.object=gc$gc,
  cline.locus="chrZ:67498406",
  locus.column="locus",
  cline.col="#91bfdb",
  cline.centre.line="chrZ:67498406",
  cline.centre.col="#91bfdb"
)

#Add a title and axis labels.
title(main = "chrZ:67498406 [TTC39B]",xlab="Hybrid index",ylab="Locus allele frequency",cex.main=1,cex.lab=1)

gcsamp=gc$gc[locus%in%c("chrZ:67498406"),
             rtmvnormDT3(
               nsamp=1000,                                #default is 1 sample per locus#
               meanlogv=mean_log_v,                       #no default, column must be specified#
               meanlogitcentre=mean_logit_centre,         #no default, column must be specified#
               varlogv=var_log_v,                         #no default, column must be specified#
               varlogitcentre=var_logit_centre,           #no default, column must be specified#
               covlogvlogitcentre=cov_log_v_logit_centre, #no default, column must be specified#
               lower=c(-Inf,-Inf),                        #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
               upper=c(Inf,Inf)                           #default, not needed for sampling from the ggcline posterior but included in case it's useful in other contexts#
             ),                                           #end of rtmvnormDT3 function#
             by="locus"];                                 #indicate a grouping variable for the sampling#

#Add the best posterior estimates to gcsamp as the final row, so it's not obscured by the mass of grey lines.
locipost=gc$gc[locus=="chrZ:67498406",.(locus,mean_log_v,mean_logit_centre)];setnames(locipost,c("locus","log_v","logit_centre"));
gcsamp=rbind(gcsamp,locipost)

#Add the individual curves to the plot in a 'for' loop.
for(i in 1:nrow(gcsamp)){
  v = gcsamp[i,exp(log_v)]
  u = gcsamp[i,logit_centre*exp(log_v)]
  S1.prop_1 = gc$gc[locus=="chrZ:67498406",S1.prop_1];
  S0.prop_1 = gc$gc[locus=="chrZ:67498406",S0.prop_1];
  
  par(new=T);
  
  if(i < nrow(gcsamp)){
    curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
          from=0,to=1,axes=F,xlab="",ylab="", col="#bdbdbd",lwd=0.5,ylim=c(0,1));
  }else{                                                                            #Add the best-fitting curve to the plot again, so it's not obscured#
    curve(S0.prop_1 + (x^v/(x^v + (1 - x)^v*exp(u)))*(S1.prop_1 - S0.prop_1), #The logit-logistic cline function#
          from=0,to=1,axes=F,xlab="",ylab="", col="#91bfdb",lwd=3,ylim=c(0,1));
  };
};#end of for loop#

######################################################
#Add data to the plot as genotypes scaled from 0 to 1#
######################################################

ploidy=2

#Create a data.table containing individual reference, hybrid index, locus name and genotype. See the ggcline examples for creating prepdata and hindlabel.
setkey(prepdata$data.prep,INDLABEL);setkey(hindlabel$hi,INDLABEL);
genodat=hindlabel$hi[,.(INDLABEL,h_posterior_mode)][prepdata$data.prep[locus=="chrZ:67498406",sum(Source_allele)/ploidy,by=c("INDLABEL","locus")]]

#Add the resulting points to the cline plot.
points(genodat[,V1]~genodat[,h_posterior_mode],pch=16,cex=0.7,col="grey")
###
dev.off()#END OF PLOT CODE FOR LOCUS TTC39B************************************************
