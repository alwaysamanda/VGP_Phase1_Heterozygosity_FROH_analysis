#### Script to plot MSMC analysis results of VGP species ####
## Date: 20250508 (May 8th, 2025)
## Author: Amanda Gardiner
## Version 1
## GOAL: To create a script which will take the output file of MSMC analysis and plot the results
## NOTES:

####

## Load in necessary packages
library(dplyr)
args <- commandArgs()

## Load in variables and data
clade <- args[6]
spec_name <- args[7]
mu <- as.numeric(args[8]) ## Mutation rate
gen_time <- as.numeric(args[9]) ## Generation time
output_file_name <- args[10]

msmc_dat <- args[11]
dat <- read.table(msmc_dat, header=TRUE)

## Create png for plot results
png(file = output_file_name)

### Plot population change ###
## Intial plot
plot(dat$left_time_boundary/mu*gen_time, (1/dat$lambda)/(2*mu), 
ylim=c(0,1250000), xlim=c(0,1200000), type="n", xlab="Years before present", ylab="Effective Population Size")

lines(dat$left_time_boundary/mu*gen_time, (1/dat$lambda)/(2*mu), type="s", col="#2734A3", lwd=1.5)

## Plot bootstrapped data
for (i in 0:29) {
    file_name <- paste0(clade, "/", spec_name, "/MSMC/Bootstrap_results/", spec_name, "_Bootstrapping_", i, ".msmc2.final.txt")
    bootstrap_dat <- read.table(file_name, header=TRUE)
    lines(bootstrap_dat$left_time_boundary/mu*gen_time, 
    ylim=c(0,1250000), xlim=c(0,1200000), 
    (1/bootstrap_dat$lambda)/(2*mu), type="s", col="#6A8BDC", lty=2)
}


dev.off()
