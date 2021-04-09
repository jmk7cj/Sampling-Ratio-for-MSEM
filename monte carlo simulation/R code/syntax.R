#----------------------------------------------------------------------------------------------#
# Author: Joseph Kush (jkush1@jhu.edu) (jmk7cj@virginia.edu)
# 
# Title: The sampling ratio in multilevel structural equation models: Considerations to inform
#        study design R code for Monte Carlo simulations
#
# Date: 4/20/2021
#
# Purpose: Master .R file to set up and run a 2-step Monte Carlo simulation study in Mplus
#          Step 1: Generate facets for Monte Carlo study
#          Step 2: Generate input files for population data generation 
#          Step 3: Use MplusAutomation to run all input files
#          Step 4: Sample from each population datafile according to a sampling ratio
#          Step 5: Create a replist for each sampled datafile
#          Step 6: Generate input files for analytic sample 
#          Step 7: Use MplusAutomation to run all input files on analytic sample
#          Step 8: Read in results, compute bias, RMSE, etc. 
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Step 1: Generate facets for Monte Carlo study 
#----------------------------------------------------------------------------------------------#
# Load necessary packages
library("MplusAutomation"); library("splitstackshape"); library("data.table"); 
library("doParallel"); library("weights"); library("foreach"); library("future")

# Remove working environment, close any open connections
rm(list = ls()); closeAllConnections()

# Set desired path to create new folders
setwd("~/Desktop/my_folder"); x <- getwd()

# Facets to vary: 3 x 3 x 2 x 2 x 4 = 144 conditions
l2_groups <- c(50, 100, 500) # number of clusters
nT <- c(20, 100, 1000) # number of units per cluster
ICC <- c(.05, .25) # ICCs
lambda <- c(.5, .8) # standardized factor loadings
SR <- c(.05, .2, .5, .8) # sampling ratios

# Total number of replications 
number_of_reps <- 1 # 1 as example

# Create a dataframe for each combination of facets
facets <- as.data.frame(expand.grid(l2_groups, nT, ICC, lambda, SR))
colnames(facets) <- c("j", "tot_ij", "ICC", "std_loading", "SR")
facets[,ncol(facets)+1] <- (facets$`tot_ij` * facets$`SR`)
colnames(facets)[ncol(facets)] <- c("samp_ij")
facets$l2_resid_var = round(facets$ICC*(1-facets$std_loading^2)/(facets$std_loading^2), digits=2)
facets$l1_resid_var = round((1-facets$ICC)*(1-facets$std_loading^2)/(facets$std_loading^2), digits=2)

# Create a series of folders and subfolders to store data for specific conditions
for(a in l2_groups) { 
  for(b in nT) { 
    for(c in lambda) { 
      for(d in ICC) { 
        for(e in SR){
dir.create(paste0(x,"/j_",a,"/nt_",b,"/load_",c,"/icc_",d,"/sr_",e), recursive=T)
}
}
}
}
}
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Simulation Step 2: Generate input files for population data generation 
#----------------------------------------------------------------------------------------------#
# Determine number of processors available -1
my_processors <- detectCores() - 1

# Create .inp file within unique folder for each condition
for(i in 1:(nrow(facets))){
j <- facets[i,1]
tot_ij <- facets[i,2]
icc <- facets[i,3]
load <- facets[i,4]    
sr <- facets[i,5]
samp_ij <- facets[i,6]
l2_resid_var <- facets[i,7]
l1_resid_var <- facets[i,8]

input <- paste(
"title: 
Step 1 - generate population data using the doubly latent model 

montecarlo: 
seed = 123; 
names = x1-x4; 
ncsizes = 1; 
csizes = ", j,"(", tot_ij,"); 
nobservations = ", j*tot_ij, "; 
nreps = ",number_of_reps,"; 
repsave = all; 
save = rep*.dat;

model population:
%within%
x1-x4*",l1_resid_var,";
Factor_w BY x1-x4*1;
Factor_w*",(1-icc),";

%between%
x1-x4*",l2_resid_var,";
Factor_b BY x1-x4*1;
Factor_b*",icc,";

analysis:
type = twolevel basic;
processors = " ,my_processors,";
", sep="") 

setwd(paste0(x,"/j_",j,"/nt_",tot_ij,"/load_",load,"/icc_",icc))
write.table(input, "step1.inp", quote=F, row.names=F, col.names=F)
}
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Simulation Step 3: Use MplusAutomation to run all input files
#----------------------------------------------------------------------------------------------#
processors <- parallel::makeCluster(my_processors, setup_strategy = "sequential")
registerDoParallel(processors)

foreach(j = l2_groups, .combine=rbind) %:%
foreach(nt = nT, .combine=rbind) %:%
foreach(load = lambda, .combine=rbind) %:%
foreach(icc = ICC, .combine=rbind) %dopar% {
library("MplusAutomation")
runModels(target=paste0(x,"/j_",j,"/nt_",nt,"/load_",load,"/icc_",icc), logFile=NULL)
}
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Simulation Step 4: Sample from each population datafile according to a sampling ratio
#----------------------------------------------------------------------------------------------#
foreach(j = l2_groups, .combine=rbind) %:%
foreach(nt = nT, .combine=rbind) %:%
foreach(load = lambda, .combine=rbind) %:%
foreach(icc = ICC, .combine=rbind) %:% 
foreach(sr = SR, .combine=rbind) %:%
foreach(rep = 1:number_of_reps, .combine=rbind) %dopar% {
    
library("MplusAutomation"); library("data.table"); library('splitstackshape')

setwd(paste0(x,"/j_",j,"/nt_",nt,"/load_",load,"/icc_",icc))
data <-fread(paste0("rep",rep,".dat"))
length <- nrow(data)/j
colnames(data) <- c("x1","x2","x3","x4","schid")
randomdraw = stratified(as.data.frame(data), group="schid", size=sr*length)
rm(data)
setwd(paste0(x,"/j_",j,"/nt_",nt,"/load_",load,"/icc_",icc,"/sr_",sr))
write.table(randomdraw, paste0("part2_data",rep,".csv"), col.names=F)
rm(randomdraw)
}
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Simulation Step 5: Create a replist for each sampled datafile
#----------------------------------------------------------------------------------------------#
for(j in l2_groups) {
for(nt in nT) {
for(load in lambda) {
for(icc in ICC) {
for(sr in SR) {

replist <- paste("part2_data",seq(from=1, to=number_of_reps),".csv", sep="") 
setwd(paste0(x,"/j_",j,"/nt_",nt,"/load_",load,"/icc_",icc,"/sr_",sr))
write.table(replist, "replist.txt", quote=F, row.names=F, col.names=F)  
}
}
}
}
}
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Simulation Step 6: Generate input files for analytic sample 
#----------------------------------------------------------------------------------------------#
for(i in 1:(nrow(facets))) {
j <- facets[i,1]
tot_ij <- facets[i,2]
icc <- facets[i,3]
load <- facets[i,4]    
sr <- facets[i,5]
samp_ij <- facets[i,6]
l2_resid_var <- facets[i,7]
l1_resid_var <- facets[i,8]

input <- paste(
"title: 
Step 2 - analyze sampled data using the doubly latent model   

data:
file = replist.txt;
type = montecarlo; 

variable:
names = i x1-x4 schid;
usevariables = x1-x4 schid;
cluster = schid;

model: 
%within%
x1-x4*",l1_resid_var,";
Factor_w BY x1*1(1)
x2*1(2)
x3*1(3)
x4*1(4);
Factor_w@",(1-icc),";

%between%
x1-x4*",l2_resid_var,";
Factor_b BY x1*1(1)
x2*1(2)
x3*1(3)
x4*1(4);
Factor_b*",icc,";

analysis:
type = twolevel;
processors = " ,my_processors,";
estimator = mlr;
", sep="") 

setwd(paste0(x,"/j_",j,"/nt_",tot_ij,"/load_",load,"/icc_",icc,"/sr_",sr))
write.table(input, "step2.inp", quote=F, row.names=F, col.names=F)
}
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Simulation Step 7: Use MplusAutomation to run all input files on analytic sample
#----------------------------------------------------------------------------------------------#
foreach(j = l2_groups, .combine=rbind) %:%
foreach(nt = nT, .combine=rbind) %:%
foreach(load = lambda, .combine=rbind) %:%
foreach(icc = ICC, .combine=rbind) %:% 
foreach(sr = SR, .combine=rbind) %dopar% {

library("MplusAutomation")
runModels(target=paste0(x,"/j_",j,"/nt_",nt,"/load_",load,"/icc_",icc,"/sr_",sr,"/step2.inp"), logFile=NULL)
}
#----------------------------------------------------------------------------------------------#


#----------------------------------------------------------------------------------------------#
# Simulation Step 8: Read in results, compute bias, RMSE, etc. 
#----------------------------------------------------------------------------------------------#
results <- as.data.frame(expand.grid(l2_groups, nT, ICC, lambda, SR))
colnames(results) <- c("j", "tot_ij", "icc", "std_loading", "SR")

for(i in 1:(nrow(results))) {
j <- results[i,1]
tot_ij <- results[i,2]
sr <- results[i,5]
icc <- results[i,3]
load <- results[i,4]

requested <- as.numeric(strsplit(summary(as.data.frame(
strsplit(as.vector(read.delim(paste0(x,"/j_",j,"/nt_",tot_ij,"/load_",load,"/icc_",icc,"/sr_",sr,"/step2.out"))
[38,]), " ")))[2], " ")[[1]][1])
results[i,6] <- requested

completed <- as.numeric(strsplit(summary(as.data.frame(
strsplit(as.vector(read.delim(paste0(x,"/j_",j,"/nt_",tot_ij,"/load_",load,"/icc_",icc,"/sr_",sr,"/step2.out"))
[39,]), " ")))[2], " ")[[1]][1])
results[i,7] <- completed

convergence <- completed/requested

r <- readModels(paste0(x,"/j_",j,"/nt_",tot_ij,"/load_",load,"/icc_",icc,"/sr_",sr,"/step2.out"))$parameters

x1_avg <- ifelse(is.character(r$unstandardized[10,4]) == T, NA, r$unstandardized[10,4]) 
x2_avg <- ifelse(is.character(r$unstandardized[11,4]) == T, NA, r$unstandardized[11,4]) 
x3_avg <- ifelse(is.character(r$unstandardized[12,4]) == T, NA, r$unstandardized[12,4]) 
x4_avg <- ifelse(is.character(r$unstandardized[13,4]) == T, NA, r$unstandardized[13,4]) 

x1_95cov <- ifelse(is.character(r$unstandardized[10,8]) == T, NA, r$unstandardized[10,8]) 
x2_95cov <- ifelse(is.character(r$unstandardized[11,8]) == T, NA, r$unstandardized[11,8]) 
x3_95cov <- ifelse(is.character(r$unstandardized[12,8]) == T, NA, r$unstandardized[12,8]) 
x4_95cov <- ifelse(is.character(r$unstandardized[13,8]) == T, NA, r$unstandardized[13,8]) 

x1_rmse <- ifelse(is.character(r$unstandardized[10,7]) == T, NA, r$unstandardized[10,7]) 
x2_rmse <- ifelse(is.character(r$unstandardized[11,7]) == T, NA, r$unstandardized[11,7]) 
x3_rmse <- ifelse(is.character(r$unstandardized[12,7]) == T, NA, r$unstandardized[12,7]) 
x4_rmse <- ifelse(is.character(r$unstandardized[13,7]) == T, NA, r$unstandardized[13,7]) 

x1_rel_bias <- ifelse(is.character(r$unstandardized[10,4]) == T, NA, 
                      (r$unstandardized[10,4] - r$unstandardized[10,3])/r$unstandardized[10,3])
x2_rel_bias <- ifelse(is.character(r$unstandardized[11,4]) == T, NA, 
                      (r$unstandardized[11,4] - r$unstandardized[11,3])/r$unstandardized[11,3])
x3_rel_bias <- ifelse(is.character(r$unstandardized[12,4]) == T, NA, 
                      (r$unstandardized[12,4] - r$unstandardized[12,3])/r$unstandardized[12,3])
x4_rel_bias <- ifelse(is.character(r$unstandardized[13,4]) == T, NA, 
                      (r$unstandardized[13,4] - r$unstandardized[13,3])/r$unstandardized[13,3])

results[i,8] <- convergence
results[i,9] <- x1_avg
results[i,10] <- x2_avg
results[i,11] <- x3_avg   
results[i,12] <- x4_avg 
results[i,13] <- x1_95cov  
results[i,14] <- x2_95cov
results[i,15] <- x3_95cov  
results[i,16] <- x4_95cov
results[i,17] <- sqrt(x1_rmse)   
results[i,18] <- sqrt(x2_rmse)
results[i,19] <- sqrt(x3_rmse) 
results[i,20] <- sqrt(x4_rmse) 
results[i,21] <- x1_rel_bias*100   
results[i,22] <- x2_rel_bias*100 
results[i,23] <- x3_rel_bias*100
results[i,24] <- x4_rel_bias*100
}

results <- results[with(results, order(j, tot_ij, icc, std_loading, SR)),]

colnames(results)[6:24] <- c("requested","completed","convergence", 
                             "x1_avg", "x2_avg", "x3_avg", "x4_avg",
                             "x1_95cov", "x2_95cov", "x3_95cov", "x4_95cov",
                             #"x1_power", "x2_power", "x3_power", "x4_power",
                             "x1_rmse", "x2_rmse", "x3_rmse", "x4_rmse",
                             "x1_rel_bias", "x2_rel_bias", "x3_rel_bias", "x4_rel_bias")

write.csv(results, paste0(x,"/results.csv"))
summary(results)
#----------------------------------------------------------------------------------------------#
