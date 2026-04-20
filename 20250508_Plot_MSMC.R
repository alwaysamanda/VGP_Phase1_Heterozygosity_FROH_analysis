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
## Intial plot by years
# plot(dat$left_time_boundary/mu*gen_time, (1/dat$lambda)/(2*mu), 
# ylim=c(0,200000), xlim=c(0,4000000), type="n", xlab="Years before present", ylab="Effective Population Size")

## Intiial plot by number of generations
plot(dat$left_time_boundary/mu, (1/dat$lambda)/(2*mu), 
ylim=c(0,300000), xlim=c(0,800000), type="n", xlab="Generations", ylab="Effective Population Size")

lines(dat$left_time_boundary/mu, (1/dat$lambda)/(2*mu), type="s", col="#2734A3", lwd=1.5)

dev.off()

