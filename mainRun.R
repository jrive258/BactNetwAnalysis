setwd("C:/Users/Tatiana Fetecua/Documents/01. Dropbox/Dropbox/02. Summer 2014/Research/Bacteria/NAV2");

configName="config.R";
inf_lev1 <- 2.0
inf_lev2 <- 2.0
minClusSize_L1 = 2;
minClusSize_L2 = 1;
p_value = 0.01;


#for(il in c(1.2, 1.5, 2.5, 3.5, 4.5, 5.5)){
#for(il in c(1.5, 3.5)){
#inf_lev1 <- il
#inf_lev2 <- il
batch_name=paste("_Former_Size_",minClusSize_L1,"_",minClusSize_L2,"_Inf_",inf_lev1,"_",inf_lev2,"_P_",p_value, sep="");
inputFile="normFileFormer.csv";
source("runBacteria.R");

batch_name=paste("_Active_Size_",minClusSize_L1,"_",minClusSize_L2,"_Inf_",inf_lev1,"_",inf_lev2,"_P_",p_value, sep="");
inputFile="normFileActive.csv";
source("runBacteria.R");

batch_name=paste("_Never_Size_",minClusSize_L1,"_",minClusSize_L2,"_Inf_",inf_lev1,"_",inf_lev2,"_P_",p_value, sep="");
inputFile="normFileNever.csv";
source("runBacteria.R");

batch_name=paste("_COPD_Size_",minClusSize_L1,"_",minClusSize_L2,"_Inf_",inf_lev1,"_",inf_lev2,"_P_",p_value, sep="");
inputFile="normFileCOPD.csv";
source("runBacteria.R");

#batch_name=paste("_TestRun","_", sep="");
#inputFile="normFileNever.csv";
#source("runBacteria.R");

#}

