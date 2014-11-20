###############################################################################
#Aux_ Functions
#This R file only makes sense in the context of NetworkCompare_R
###############################################################################


###############################################################################
#Cosine Difference
###############################################################################
cosineDif <- function	(V1,	#First Vector
						V2	#Second Vector
						) {
  return( sum(V1 * V2) / ( sqrt(sum(V1 * V1)) * sqrt(sum(V2*V2)) ) )
}

###############################################################################
#Get the path of a result file, takes into account the batch_name
###############################################################################
getResFile <- function	(in_file_name, in_file_extension, in_file_folder){
	return (file.path(in_file_folder, paste(in_file_name, batch_name, ".", in_file_extension, sep="")))
}

###############################################################################
#Store the cluster information in a csv file
###############################################################################
storeClusters <- function	(in_file_name, in_cluster_list, in_levels = 2) {
	if(in_levels == 2){
		for(i in 1:length(in_cluster_list)){
			cat(paste("Cluster ", i), "\n", file=in_file_name, append=TRUE)
			for(j in 1:length(in_cluster_list[[i]])){
				cat(paste("Cluster ", i, "-", j), "\n", file=in_file_name, append=TRUE)
				write.table(in_cluster_list[[i]][[j]], file=in_file_name, sep=",", append=TRUE, col.names=FALSE, na="")
			}
		}
	}else{
		for(i in 1:length(in_cluster_list)){
			cat(paste("Cluster ", i), "\n", file=in_file_name, append=TRUE)
			write.table(in_cluster_list[[i]], file=in_file_name, sep=",", append=TRUE, col.names=FALSE, na="")
		}
	}
}

storeClustersMatrix <- function	(in_file_name, in_cluster_list, in_data) {
	#Get the split sections
	final_matrix <- in_data[rle(unlist(in_cluster_list))$values, rle(unlist(in_cluster_list))$values]
	sections <- list()
	br_point <- 0
	for(i in 1:(length(in_cluster_list)-1)){
		br_point <- br_point + length(in_cluster_list[[i]])
		final_matrix <- cbind(final_matrix[,1:br_point],c(rep.int(NA,nrow(final_matrix))),final_matrix[,(br_point+1):ncol(final_matrix)])
		final_matrix <- rbind(final_matrix[1:br_point,],c(rep.int(NA,ncol(final_matrix))),final_matrix[(br_point+1):nrow(final_matrix),])
		br_point <- br_point + 1
	}
	write.table(final_matrix, file=in_file_name, sep=",", append=TRUE, col.names=NA, na="")
	return(final_matrix)
}


###############################################################################
# extrapolate
# Changes the values of the matrix extrapolating 
# based on the values of E.
# E contains the extrapolation values (from, to) in pairs
###############################################################################
extrapolate <- function (M, E 	
		                 ) {
	M2 <- M;
	first_idx <- E[1,];
	for(i in 2:nrow(E)){
		second_idx = E[i,];
		a = first_idx[1];
		b = first_idx[2];
		c = second_idx[1];
		d = second_idx[2];
		indexes = which(M>=a & M<=c, arr.ind=TRUE);
		M2[indexes] <- ( ( M[indexes] - a ) * (d-b) / (c-a) ) + b;
		first_idx <- second_idx;
	}
	return(M2);
}

###############################################################################
#MCL clustering for the differences (two level clustering)
###############################################################################
clusterMatrix <- function	(in_matrix){
	#Map (values in the config file)
	in_matrix_abs <- abs(extrapolate(in_matrix, extrap_values_1))
	mcl_clusters <- mcl(in_matrix_abs,mclc_inf_level1,2000, verbose = F, heatmaps=F);
	mcl_list <- collect.mcl.clusters2(mcl_clusters,mclc_minClusSize_level1);
	
	mcl_listFull <- list();
	for(i in 1:length(mcl_list)){
		if(length(mcl_list[[i]]) >= mclc_minClusSize_level1 && length(mcl_list[[i]]) > 1){
			new_data <- in_matrix[mcl_list[[i]], mcl_list[[i]]];
			new_data <- extrapolate(new_data, extrap_values_2);
			mcl_clusters2 <- mcl(new_data,mclc_inf_level2,2000, verbose = F);
			mcl_leaders2 <- collect.mcl.leaders(mcl_clusters2);
			mcl_list2 <- collect.mcl.clusters2(mcl_clusters2, mclc_minClusSize_level2);
			#mcl_listFull <- c(mcl_listFull, mcl_list2);
			mcl_listFull[[length(mcl_listFull)+1]] <- mcl_list2
		}
	}	
	return (mcl_listFull)
}

###############################################################################
#MCL clustering for the differences (single level clustering)
###############################################################################
clusterMatrixSingleLevel <- function	(in_matrix){
	#Map (values in the config file)
	in_matrix_abs <- abs(extrapolate(in_matrix, extrap_values_1))
	mcl_clusters <- mcl(in_matrix_abs,mclc_inf_level1,2000, verbose = F, heatmaps=F);
	mcl_list <- collect.mcl.clusters2(mcl_clusters,mclc_minClusSize_level1);
	return (mcl_list)
}