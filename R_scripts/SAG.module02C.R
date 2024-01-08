# install new packages
install.packages("devtools")

# in order to install geographic cline package BAHZ, you'll need to install RStan on your computer
# see https://github.com/tjthurman/BAHZ for details on how to install RStan for Mac, Windows, or Linux users.

# once Rstan is installed you can run:
devtools::install_github("tjthurman/BAHZ")

# load packages
library(tidyverse)
library(cowplot)
library(bahz) #https://github.com/tjthurman/BAHZ

##Note: You'll need to set your own working directory accordingly
setwd("/Users/danielhooper/Documents/Teaching/RGGS Workshop/R/Module02")
rm(list = ls())

##Load dataframe of cline data
data <- read.csv(file = "bahz/billcolor.snps.csv", header = TRUE)
data <- data |> filter(auto.2N >= 10)
data$cyp.freq <- (data$chr8.3144828 / data$auto.2N)
data$ttc.freq <- (data$chrZ.67547840 / data$chrZ.2N)

##Create SNP specific dfs 
cyp.data <- data |> select(distance,chr8.3144828,auto.2N)
cyp.data <- cyp.data |> rename(transectDist = distance, nFocalAllele = chr8.3144828, nTotalAlleles = auto.2N)
ttc.data <- data |> select(distance,chrZ.67547840,chrZ.2N)
ttc.data <- ttc.data |> rename(transectDist = distance, nFocalAllele = chrZ.67547840, nTotalAlleles = chrZ.2N)
inv.data <- data |> select(distance,chrZ.INV1,chrZ.2N)
inv.data <- inv.data |> rename(transectDist = distance, nFocalAllele = chrZ.INV1, nTotalAlleles = chrZ.2N)


##Prepare your genetic cline data for loading into Stan
cyp.cline <- prep_geno_data(cyp.data, type = "bi")
ttc.cline <- prep_geno_data(ttc.data, type = "bi")
inv.cline <- prep_geno_data(inv.data, type = "bi")

##Specify priors
make_prior_config(name = "prior_config_template.yaml", overwrite = F) ##Might need to not have some tidyverse packages loaded

##Fit the model
set.seed(11238)

cyp.cline.fit <- fit_geno_cline(data = cyp.data, prior_file = "bahz/prior_config_CYP2J19.yaml",
                            type = "bi", tails = "none")

ttc.cline.fit <- fit_geno_cline(data = ttc.data, prior_file = "bahz/prior_config_TTC39B.yaml",
                                type = "bi", tails = "none")

inv.cline.fit <- fit_geno_cline(data = inv.data, prior_file = "bahz/prior_config_INV1.yaml",
                                type = "bi", tails = "none")

##Analyze the results
cline_summary(cyp.cline.fit)
cline_summary(ttc.cline.fit)
cline_summary(inv.cline.fit)

##Basic cline plotting
plot_cline(stanfit = cyp.cline.fit, data = cyp.data,
           add.obs = T, confidence = T, cline.col = "orange", point.col = "orange",
           xlab = "distance along transect", ylab = "chr8:3144828 allele frequency")
plot_cline(stanfit = ttc.cline.fit, data = ttc.data,
           add.obs = T, confidence = T, cline.col = "lightblue", point.col = "lightblue",
           xlab = "distance along transect", ylab = "chrZ.67547840 allele frequency")
plot_cline(stanfit = inv.cline.fit, data = inv.data,
           add.obs = T, confidence = T, cline.col = "pink", point.col = "pink",
           xlab = "distance along transect", ylab = "chrZ.INV1 allele frequency")

##Higher quality plotting
cyp_pred_cline <- predict_cline(stanfit = cyp.cline.fit, 
                            distance = 0:1825, 
                            confidence = T)
ttc_pred_cline <- predict_cline(stanfit = ttc.cline.fit, 
                            distance = 0:1825, 
                            confidence = T)
inv_pred_cline <- predict_cline(stanfit = inv.cline.fit, 
                                distance = 0:1825, 
                                confidence = T)

theme_set(theme_classic())
ggplot(data = cyp_pred_cline, aes(x = transectDist, y = p_mean,
                              ymin = low_0.95_HPDI,
                              ymax = up_0.95_HPDI)) +
  geom_ribbon(fill = alpha("#fdae61", 0.5)) +
  geom_line(linewidth = 1.0, color = "#fdae61") +
  xlab("distance along transect") +
  ylab("allele frequency [CYP2J19]") +
  theme_classic()

ggplot(data = ttc_pred_cline, aes(x = transectDist, y = p_mean,
                                  ymin = low_0.95_HPDI,
                                  ymax = up_0.95_HPDI)) +
  geom_ribbon(fill = alpha("#91bfdb", 0.5)) +
  geom_line(linewidth = 1.0, color = "#91bfdb") +
  xlab("distance along transect") +
  ylab("allele frequency [TTC39B]") +
  theme_classic()

ggplot(data = inv_pred_cline, aes(x = transectDist, y = p_mean,
                                  ymin = low_0.95_HPDI,
                                  ymax = up_0.95_HPDI)) +
  geom_ribbon(fill = alpha("#E76254", 0.5)) +
  geom_line(linewidth = 1.0, color = "#E76254") +
  xlab("distance along transect") +
  ylab("allele frequency [INV1]") +
  theme_classic()

##Plot both clines in one figure

tmp <- ttc_pred_cline[2:5]
tmp <- tmp |> rename(ttc.p_mean = p_mean, ttc.p_median = p_median, ttc.low_0.95_HPDI = low_0.95_HPDI, ttc.up_0.95_HPDI = up_0.95_HPDI)
tmp2 <- inv_pred_cline[2:5]
tmp2 <- tmp2 |> rename(inv.p_mean = p_mean, inv.p_median = p_median, inv.low_0.95_HPDI = low_0.95_HPDI, inv.up_0.95_HPDI = up_0.95_HPDI)
merge <- cbind(cyp_pred_cline, tmp, tmp2)

##Plot best-fit clines as a single panel
theme_set(theme_classic())
plot <- ggplot(data = merge, aes(x = transectDist)) +
  geom_ribbon(aes(y = p_mean, ymin = low_0.95_HPDI, ymax = up_0.95_HPDI), fill = alpha("#fdae61", 0.5)) +
  geom_line(aes(y = p_mean), linewidth = 1.0, color = "#fdae61") +
  geom_ribbon(aes(y = ttc.p_mean, ymin = ttc.low_0.95_HPDI, ymax = ttc.up_0.95_HPDI), fill = alpha("#91bfdb", 0.5)) +
  geom_line(aes(y = ttc.p_mean), linewidth = 1.0, color = "#91bfdb") +
  xlab("distance along transect") +
  ylab("allele frequency") +
  theme_classic()

##Add population level data to plot
plot + geom_point(aes(x=0, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=304.3, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=470.49, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=508.94, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=541.08, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=563.3, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=576.86, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=583.06, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=603.88, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=618.3, y=0.988), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=628.23, y=0.991), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=629.62, y=0.971), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=633.79, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=639.87, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=641.69, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=646.92, y=0.986), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=677.89, y=0.983), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=672.45, y=1.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=803.05, y=0.650), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=858.05, y=0.533), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=869.8, y=0.533), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=881.97, y=0.160), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=927.57, y=0.50), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=946.15, y=0.40), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=966.23, y=0.423), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=1138.57, y=0.067), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=1294.23, y=0.035), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=1401.95, y=0.077), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=1783.76, y=0.025), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=1814.7, y=0.00), colour = "#fdae61", shape = 3, alpha = 0.9) +
  geom_point(aes(x=0, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=304.3, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=470.49, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=508.94, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=541.08, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=563.3, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=576.86, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=583.06, y=0.962), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=603.88, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=618.3, y=0.972), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=628.23, y=0.989), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=629.62, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=633.79, y=0.990), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=639.87, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=641.69, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=646.92, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=677.89, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=672.45, y=0.987), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=803.05, y=0.952), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=858.05, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=869.8, y=0.864), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=881.97, y=0.900), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=927.57, y=1.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=946.15, y=0.935), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=966.23, y=0.889), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=1138.57, y=0.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=1294.23, y=0.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=1401.95, y=0.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=1783.76, y=0.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_point(aes(x=1814.7, y=0.00), colour = "#91bfdb", shape = 3, alpha = 0.9) +
  geom_vline(xintercept = 866.28, linewidth = 0.8, linetype = "dashed", colour = "#fdae61") +
  geom_vline(xintercept = 1022.02, linewidth = 0.8, linetype = "dashed", colour = "#91bfdb")
  
##What about if you want to compare against the geographic cline for INV1?
plot2 <- ggplot(data = merge, aes(x = transectDist)) +
  geom_ribbon(aes(y = p_mean, ymin = low_0.95_HPDI, ymax = up_0.95_HPDI), fill = alpha("#fdae61", 0.5)) +
  geom_line(aes(y = p_mean), linewidth = 1.0, color = "#fdae61") +
  geom_ribbon(aes(y = ttc.p_mean, ymin = ttc.low_0.95_HPDI, ymax = ttc.up_0.95_HPDI), fill = alpha("#91bfdb", 0.5)) +
  geom_line(aes(y = ttc.p_mean), linewidth = 1.0, color = "#91bfdb") +
  geom_ribbon(aes(y = inv.p_mean, ymin = inv.low_0.95_HPDI, ymax = inv.up_0.95_HPDI), fill = alpha("#E76254", 0.5)) +
  geom_line(aes(y = inv.p_mean), linewidth = 1.0, color = "#E76254") +
  xlab("distance along transect") +
  ylab("allele frequency") +
  geom_vline(xintercept = 647.54, linewidth = 0.8, linetype = "dashed", colour = "#E76254") +
  geom_vline(xintercept = 866.28, linewidth = 0.8, linetype = "dashed", colour = "#fdae61") +
  geom_vline(xintercept = 1022.02, linewidth = 0.8, linetype = "dashed", colour = "#91bfdb") +
  theme_classic()

##Just plot the population level data
coordinates <- data |> select(distance, cyp.freq, ttc.freq)
coordinates |> ggplot(aes(x = distance)) +
  geom_point(aes(y = cyp.freq), colour = "#fdae61", shape = 3, size = 2, alpha = 0.9) +
  geom_point(aes(y = ttc.freq), colour = "#91bfdb", shape = 3, size = 2, alpha = 0.9) +
  coord_cartesian(xlim = c(0, 1826)) +
  xlab("distance along transect") +
  ylab("allele frequency")
  
