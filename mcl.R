################################################################
## This scripts shows how run the MCL clustering algorithm
## developed by Stijn van Dongen with R
## source('~/Documents/enseignement/bioinformatics_courses/statistics_bioinformatics/R-files/network_topology/sylvain/mcl.R')
##
## Authors: Sylvain Brohee
## Date:    October 2008

########################################
## Functions
# Add the identity matrix to the matrix
add.one <- function (M) {
  for (i in 1:dim(M)[1]) {
    if (M[i,i] == 0) {
      M[i,i] <- M[i,i] + 1;
    }
  }
  return (M);
}

# Inflation step of MCL
inflate <- function (M,
		     inf) {
  M <- M^(inf);
  return (M);
}

# Normalize the matrix by column
norm <- function (M) {
  colum.sum <- apply(M,2,sum)
  M <- t(M) / colum.sum
  return (t(M))
}

# MCL procedure
mcl <- function (M, 	# Matrix
		 inf, 	# Inflation value
		 iter, 	# Number of iterations
		 verbose = F,
		 heatmaps = F
		 ) { 
  if(heatmaps){
	plot_st.heatmap(M, paste("mcl_", "_Iter.png", sep=""));
  }
  for (i in 1:iter) {
    old.M <- M;
    M.norm <- norm(M);
	if(heatmaps && i==1){
	  plot_st.heatmap(M.norm, paste("mcl_", "_Iter0.png", sep=""));
    }
    M <- M.norm%*%M.norm;
    M <- inflate(M, inf);
    M <- norm(M);
	if(heatmaps){
		plot_st.heatmap(M, paste("mcl_", "_Iter", i, ".png", sep=""));
	}
    if (sum(old.M == M) == dim(M)[1]*dim(M)[2]) {
      break;
    }
    if (verbose) {
      print (paste ("iteration", i));
    } 
  }
  return (M);
}


# collect.mcl.clusters
# having a MCL output, makes a data frame having 
# in the first column the node id
# in the second column the cluster to which this node belongs
# run through the columns of the mcl result matrix and looks for all nodes that compose a cluster
collect.mcl.clusters <- function (M 	# Matrix (mcl result)
		                 ) {
  cluster.list <- list();
  cluster.number <- 1;
  M.names <- row.names(M);
  clustered.nodes <- vector(mode = "logical", length = dim(M)[1])
  for (i in 1:dim(M)[1]) {
    nodes <- M.names[which(M[i,] != 0)];
    if (length(nodes) > 0 && !clustered.nodes[which(M[i,] != 0)]) {
		#print (nodes);
		#Add the nodes to the cluster list
		cluster.list[[cluster.number]] <- nodes;
		
		clustered.nodes[which(M[i,] != 0)] = T;
		cluster.number <- cluster.number + 1;
    }
  }
  return(cluster.list);
}

# collect.mcl.clusters
# having a MCL output, makes a data frame having 
# in the first column the node id
# in the second column the cluster to which this node belongs
# run through the columns of the mcl result matrix and looks for all nodes that compose a cluster
# Only selects clusters with more than the indicated number of elements
collect.mcl.clusters2 <- function (M, size 	# Matrix (mcl result)
		                 ) {
  cluster.list <- list();
  cluster.number <- 1;
  M.names <- row.names(M);
  clustered.nodes <- vector(mode = "logical", length = dim(M)[1])
  for (i in 1:dim(M)[1]) {
    nodes <- M.names[which(M[i,] != 0)];
    if (length(nodes) >= size && !clustered.nodes[which(M[i,] != 0)]) {
		#print (nodes);
		#Add the nodes to the cluster list
		cluster.list[[cluster.number]] <- nodes;
		
		clustered.nodes[which(M[i,] != 0)] = T;
		cluster.number <- cluster.number + 1;
    }
  }
  return(cluster.list);
}

collect.mcl.leaders <- function (M 	# Matrix (mcl result)
		                 ) {
  
  M.names <- row.names(M);
  leaders.list <- vector();
  clustered.nodes <- vector(mode = "logical", length = dim(M)[1])
  for (i in 1:dim(M)[1]) {
    nodes <- M.names[which(M[i,] != 0)];
    if (length(nodes) >= 1 && !clustered.nodes[which(M[i,] != 0)]) {
		#Add the nodes to the cluster list
		leaders.list <- c(leaders.list, M.names[i]);
		clustered.nodes[which(M[i,] != 0)] = T;
	}
  }
  return(leaders.list);
}


###################################
## Example

#file <- file.path(dir.data,"data2.txt");
#mcl.data <- read.delim(file, header = T, row.names = 1);
#mcl.data <- as.matrix(mcl.data);
#mcl.data <- add.one(mcl.data);
#inf <- 2

# Launch MCL
#mcl.clusters <- mcl(mcl.data,inf,200, verbose = T);
#collect.mcl.clusters(mcl.clusters);

#fileOut <- file.path(dir.data,"clusters.out.csv");
#write.table(mcl.clusters, file=fileOut, sep=",", append=TRUE, col.names=NA, na="");
