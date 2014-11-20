#########################################################
### A) Installing and loading required packages
#########################################################

if (!require("gplots")) {
   install.packages("gplots", dependencies = TRUE)
   library(gplots)
   }
if (!require("RColorBrewer")) {
   install.packages("RColorBrewer", dependencies = TRUE)
   library(RColorBrewer)
   }
libs <- c("qgraph", "Hmisc", "igraph", "colorspace")
lapply(libs, require, character.only=T)

#########################################################
### B) Installing and loading required custom functions
#########################################################
source("mcl.R");
source("heatmap.R");

#########################################################
### C) Check overrriden values (results and figures)
#########################################################
if(!exists("dir.data.name") || dir.data.name == ""){
	dir.data.name <- "data"
}
if(!exists("dir.results.name") || dir.results.name == ""){
	dir.results.name <- "results"
}
if(!exists("batch_name")){
	batch_name <- ""
}

if(!exists("inf_lev1")){
	inf_lev1 <- 1.5
}
if(!exists("inf_lev2")){
	inf_lev2 <- 1.5
}

if(!exists("minClusSize_L1")){
	minClusSize_L1 <- 5
}
if(!exists("minClusSize_L2")){
	minClusSize_L2 <- 5
}
if(!exists("p_value")){
	p_value <- 0.05
}


#if(!exists("thresholds")){
#	thresholds <- c(-0.4,0.4);
#}

#Clustering level 1
#ex_from <- c(-1.0, thresholds[1], thresholds[1]+0.01, 0, thresholds[2]-0.01, thresholds[2], 1.0);
#ex_to <- c(-1.0, thresholds[1], -0.01, 0, 0.01, thresholds[2], 1.0);
#extrap_values_1 <- array(c(ex_from,ex_to), dim=c(length(ex_from),2));

#Clustering level 2
#ex_from <- c(-1.0, thresholds[2]-0.01, thresholds[2], 1.0);
#ex_to <-   c( 0.0, 0.01, thresholds[2], 1.0);
#extrap_values_2 <- array(c(ex_from,ex_to), dim=c(length(ex_from),2));

#########################################################
### D) Local files and folders (results and figures)
#########################################################
dir.course.local <- getwd();
dir.data <- file.path(dir.course.local,dir.data.name); # directory containing the datasets
dir.create(dir.data, showWarnings=F,recursive=T);
dir.results <- file.path(dir.course.local,dir.results.name); # directory containing the results
dir.create(dir.results, showWarnings=F,recursive=T);
dir.figures <- file.path(dir.results,paste('figures',batch_name, sep="")); # directory containing the final results
dir.create(dir.figures, showWarnings=F,recursive=T);