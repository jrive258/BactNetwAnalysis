# Clear R workspace
rm(list=ls(all=TRUE))

# Set a working directory
setwd("C:/Users/Tatiana Fetecua/Documents/01. Dropbox/Dropbox/02. Summer 2014/Research/Bacteria/NAV2");
source("config.R")
source("auxFunctions.R")

# Pick the minimum number of reads needed for an OTU to be considered
minreads <- 100

# Read in data frames
abund_raw <- read.csv(file.path(dir.data, "All Reads Raw.csv"))
#abund_norm <- read.csv(file.path(dir.data, "All Reads Normalized.csv"))
temp <- abund_raw[,3:ncol(abund_raw)] #* abund_raw[,3:ncol(abund_raw)]
temp_sum <- rowSums(temp)
temp <- temp / temp_sum
#abund_norm <- sqrt(temp)
abund_norm <- temp
abund_norm <- cbind(Group=abund_raw[,"Group"], Smoking = abund_raw[,"Smoking.Status"], abund_norm)

write.table(abund_norm, file=file.path(dir.data, "All Reads Normalized.csv"), sep=",", append=FALSE, col.names=TRUE, na="", row.names=FALSE)
# Discard low abundance OTUs and subjects
keep.cols <- abund_norm[3:length(abund_norm)][which(colSums(abund_raw[3:length(abund_raw)])>minreads)]

# Pick groups to keep

# Select subjects (Based on column 2)
group <- c("Active")
input <- keep.cols[which(abund_raw[2]==group),]
input <- cbind(Status = 'Active', input[sort(colnames(input))])
write.table(input, file=file.path(dir.data, "normFileActive.csv"), sep=",", append=FALSE, col.names=TRUE, na="", row.names=FALSE)
group <- c("Former")
input <- keep.cols[which(abund_raw[2]==group),]
input <- cbind(Status = 'Former', input[sort(colnames(input))])
write.table(input, file=file.path(dir.data, "normFileFormer.csv"), sep=",", append=FALSE, col.names=TRUE, na="", row.names=FALSE)
group <- c("Never")
input <- keep.cols[which(abund_raw[2]==group),]
input <- cbind(Status = 'Never', input[sort(colnames(input))])
write.table(input, file=file.path(dir.data, "normFileNever.csv"), sep=",", append=FALSE, col.names=TRUE, na="", row.names=FALSE)

group <- c("Smoker")
input <- keep.cols[which(abund_raw[1]==group),]
input <- cbind(Status = 'Non COPD', input[sort(colnames(input))])
write.table(input, file=file.path(dir.data, "normFileNonCOPD.csv"), sep=",", append=FALSE, col.names=TRUE, na="", row.names=FALSE)
group <- c("COPD")
input <- keep.cols[which(abund_raw[1]==group),]
input <- cbind(Status = 'COPD', input[sort(colnames(input))])
write.table(input, file=file.path(dir.data, "normFileCOPD.csv"), sep=",", append=FALSE, col.names=TRUE, na="", row.names=FALSE)

#Run the comparisons
#batch_name = "_Active_Former"
#load_inputFile1 = "rawFileActive.csv"
#load_inputFile2 = "rawFileFormer.csv"
#source("NetworkCompare.R")

#batch_name = "_Active_Never"
#load_inputFile1 = "rawFileActive.csv"
#load_inputFile2 = "rawFileNever.csv"
#source("NetworkCompare.R")

#batch_name = "_Former_Never"
#load_inputFile1 = "rawFileFormer.csv"
#load_inputFile2 = "rawFileNever.csv"
#source("NetworkCompare.R")

#COPD
#batch_name = "_COPD_NonCOPD"
#load_inputFile1 = "rawFileCOPD.csv"
#load_inputFile2 = "rawFileNonCOPD.csv"
#source("NetworkCompare.R")

#batch_name = "_COPD_Never"
#load_inputFile1 = "rawFileCOPD.csv"
#load_inputFile2 = "rawFileNever.csv"
#source("NetworkCompare.R")

#batch_name = "_NonCOPD_Never"
#load_inputFile1 = "rawFileNonCOPD.csv"
#load_inputFile2 = "rawFileNever.csv"
#source("NetworkCompare.R")


