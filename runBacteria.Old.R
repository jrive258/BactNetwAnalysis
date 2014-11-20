if(!exists("configName") || configName==""){
	configName <- "config.R";
}

#setwd("C:/Users/Tatiana Fetecua/Documents/01. Dropbox/Dropbox/02. Summer 2014/Research/Bacteria/NAV2");
source(configName);

####################
#Log the parameters
####################
parameters.list <- list();
#parameters.list[[length(parameters.list)+1]] <- "Threshold 1";
#parameters.list[[length(parameters.list)+1]] <- ex_values_1;
#parameters.list[[length(parameters.list)+1]] <- "Threshold 2";
#parameters.list[[length(parameters.list)+1]] <- ex_values_2;
parameters.list[[length(parameters.list)+1]] <- "P value";
parameters.list[[length(parameters.list)+1]] <- p_value;
parameters.list[[length(parameters.list)+1]] <- "Min Cluster Size Level 1";
parameters.list[[length(parameters.list)+1]] <- minClusSize_L1;
parameters.list[[length(parameters.list)+1]] <- "Min Cluster Size Level 2";
parameters.list[[length(parameters.list)+1]] <- minClusSize_L2;
parameters.list[[length(parameters.list)+1]] <- "Inflation 1";
parameters.list[[length(parameters.list)+1]] <- inf_lev1;
parameters.list[[length(parameters.list)+1]] <- "Inflation 2";
parameters.list[[length(parameters.list)+1]] <- inf_lev2;
fileOut <- file.path(dir.figures,paste("parameters",batch_name,".txt", sep=""));
lapply(parameters.list, cat, "\n", file=fileOut, append=TRUE);

#########################################################
#First step, calculate the correlations from Input File
#########################################################
# Read in source files
if(!exists("inputFile") || inputFile == ""){
	inputFile <- "inputFile.csv";
}
pc <- read.csv(file.path(dir.data, inputFile));

# Create correlation matrices
#mcl.data <- as.matrix(cor(pc[,-1], use="pairwise.complete.obs"))
correlations <- rcorr(as.matrix(pc[,-1]))
mcl.data <- as.matrix(correlations$r)

# Replace NA values with 0
mcl.data[is.na(mcl.data)] <- 0

#Write correlations before cut-off
fileOut <- file.path(dir.figures,"correlation.output.csv");
write.table(mcl.data, file=fileOut, sep=",", append=TRUE, col.names=NA, na="") #col.names=NA is used to align column names correctly

#Eliminate non-significant correlations (based on p-value given by the pearson correlation
mcl.data.original <- mcl.data
mcl.data[which(correlations$P>p_value)] <- 0

# Write to output file
# Using write.table instead of write.csv in order to append all
fileOut <- file.path(dir.figures,"correlation.output.pvalued.csv");
write.table(mcl.data, file=fileOut, sep=",", append=TRUE, col.names=NA, na="") #col.names=NA is used to align column names correctly


#########################################################
#Second step, calculate the first clustering
#########################################################
#Do mapping
#Map (values in the config file)
#mcl.data.ext <- extrapolate(mcl.data, ex_values_1);
mcl.data.ext <- mcl.data;

#Absolute Values
mcl.data.abs <- abs(mcl.data.ext);

#Remove columns and rows with only 0s, normalization cannot be done in mcl if a column value is 0
mcl.data.abs <- mcl.data.abs[,which(colSums(mcl.data.abs[,])!=0)]
mcl.data.abs <- mcl.data.abs[which(rowSums(mcl.data.abs[,])!=0),]

# Launch MCL
mcl.clusters <- mcl(mcl.data.abs,inf_lev1,2000, verbose = F, heatmaps=F);
mcl.list <- collect.mcl.clusters2(mcl.clusters,minClusSize_L1);
fileOut <- file.path(dir.figures,paste("Level1Clusters",batch_name,".csv", sep=""));
for(j in 1:length(mcl.list)){
	write.table(mcl.list[[j]], file=fileOut, sep=",", append=TRUE, col.names=NA, na="");
}

#First heat map of the clusters
plot_st.heatmap(mcl.data[rle(unlist(mcl.list))$values, rle(unlist(mcl.list))$values], paste("Level1Clusters",batch_name,".png", sep=""));

#########################################################
#Third step, create heatmap for each cluster found
#########################################################
#for(i in 1:length(mcl.list)){
#	if(length(mcl.list[[i]]) > 1){
#		fname <- paste("Level1",batch_name,"Cluster", i, ".png", sep="");
#		#print(fname);
#		plot_st.heatmap(mcl.data[mcl.list[[i]], mcl.list[[i]]], fname);
#	}
#}

#########################################################
#Fourth step, create level 2 clusters
#########################################################
mcl.listFull <- list();
for(i in 1:length(mcl.list)){
	if(length(mcl.list[[i]]) >= minClusSize_L1 && length(mcl.list[[i]]) > 1){
		fname <- paste("Level1",batch_name,"Cluster", i, "_Sorted.png", sep="");
		fnameEP <- paste("Level1",batch_name,"Cluster", i, "_NewEP.png", sep="");
		fnameCsv <- paste("Level1",batch_name,"Cluster", i, ".csv", sep="");
		fnameCsvMean <- paste("Level2",batch_name,"Clusters", "_means.csv", sep="");
		fnameCsvStd <- paste("Level2",batch_name,"Clusters", "_std.csv", sep="");
		fnameCsvLeaders <- paste("Level2",batch_name,"Clusters", "_leaders.csv", sep="");
		
		new.data <- mcl.data[mcl.list[[i]], mcl.list[[i]]];
		
		#new.data <- extrapolate(new.data, ex_values_2);
		#Make all negative values = 0
		new.data[which(new.data<0)] <- 0
		
		mcl.clusters2 <- mcl(new.data,inf_lev2,2000, verbose = F);
		mcl.leaders2 <- collect.mcl.leaders(mcl.clusters2);
		mcl.list2 <- collect.mcl.clusters2(mcl.clusters2, minClusSize_L2);
		mcl.listFull <- c(mcl.listFull, mcl.list2);
		
		mean_values <- matrix(data=NA,nrow=length(mcl.list2),ncol=length(mcl.list2));
		std_values <- matrix(data=NA,nrow=length(mcl.list2),ncol=length(mcl.list2));
		
		fileOut <- file.path(dir.figures,fnameCsv);
		if(length(mcl.list2) > 0){
			for(j in 1:length(mcl.list2)){
				#write.table(mcl.list2[[j]], file=fileOut, sep=",", append=TRUE, col.names=NA, na="");
			}
			#plot_st.heatmap(mcl.data[rle(unlist(mcl.list2))$values, rle(unlist(mcl.list2))$values], fname);
			
			
			#Calculate the mean and standard deviation of intercluster correlations
			for(j in 1:length(mcl.list2)){
				for(k in j:length(mcl.list2)){
					corr_info <- mcl.data[mcl.list2[[j]], mcl.list2[[k]]];
					if(j==k){
						rel_data <-	corr_info[lower.tri(corr_info, diag=FALSE)];
					}else{
						rel_data <- corr_info;
					}
					mean_values[j,k] = mean(rel_data);
					std_values[j,k] = sd(rel_data);
				}
			}
			fileOut <- file.path(dir.figures,fnameCsvMean);
			write.table(mean_values, file=fileOut, sep=",", append=TRUE, col.names=NA, na="");
			fileOut <- file.path(dir.figures,fnameCsvStd);
			write.table(std_values, file=fileOut, sep=",", append=TRUE, col.names=NA, na="");
			fileOut <- file.path(dir.figures,fnameCsvLeaders);
			write.table(mcl.leaders2, file=fileOut, sep=",", append=TRUE, col.names=NA, na="");
		}
	}
}
fileNameOutput <- paste("Level1Clusters",batch_name,"_InnerSortedCorrelations.csv", sep="");
fileNameOutput <- file.path(dir.figures,fileNameOutput);

write.table(mcl.data[rle(unlist(mcl.listFull))$values, rle(unlist(mcl.listFull))$values], file=fileNameOutput, sep=",", append=TRUE, col.names=NA, na="");

fileOut <- file.path(dir.figures,paste("Level1Clusters",batch_name,"_InnerSorted.csv", sep=""));
for(j in 1:length(mcl.listFull)){
	write.table(mcl.listFull[[j]], file=fileOut, sep=",", append=TRUE, col.names=NA, na="");
}

# Assign heat colors and shapes
heat.data <- mcl.data[rle(unlist(mcl.listFull))$values, rle(unlist(mcl.listFull))$values]
heat.colors <- c()
node_names_heat <- colnames(heat.data)
#color_list <- rainbow(length(mcl.listFull))
color_list <- rainbow(31, alpha=1)
for(j in 1:ncol(heat.data)){
	heat.colors[j] <- "#BBBBBB99"
	#node.shapes[j] <- shape_list[10]
	for(k in 1:length(mcl.listFull)){
		if(length(which(mcl.listFull[[k]]==node_names_heat[j]))>0){
			col_index <- ( ( k * 6 ) %% 31 ) + 1
			heat.colors[j] <- color_list[col_index]
		}
	}
	#for(k in 1:length(mcl.list)){
		#if(length(which(mcl.list[[k]]==node_names_graph[j]))>0){
			#node.shapes[j] <- shape_list[k]
			#node.names[j] <- paste(k, "_", node.names[j])
		#} 
	#}
}

plot_st.heatmap_w_colors(mcl.data[rle(unlist(mcl.listFull))$values, rle(unlist(mcl.listFull))$values], paste("Level1Clusters",batch_name,"_InnerSorted.png", sep=""), noFolder=FALSE, colorList=heat.colors);


#Now we print the graph
#mcl.data has all the correlations
xn <- as.matrix(pc[,-1])
xm <- mcl.data
xc <- colSums(xn[,])

### Node Sizes
xcfix <- xc
minxcfix <- min(xc[which(xc>0)])
xcfix[which(xc==0)] <- minxcfix

scale <- ceiling(max(xcfix)/min(xcfix))
if (max(xcfix) %% min(xcfix) == 0) {
	node.size <- log10(10^(nchar(scale)-1)*xcfix)
} else {
	node.size <- log10(10^(nchar(scale))*xcfix)             
}
node.size <- 1.7*node.size

### Get short names for node names
xnames <- colnames(xn)
xnames <- gsub("Phylum", "P", xnames)
xnames <- gsub("Class", "C", xnames)
xnames <- gsub("Order", "O", xnames)
xnames <- gsub("Family", "F", xnames)
names.front <- substr(xnames, 1, 5)
substrRight <- function(x, n){
substr(x, nchar(x)-n+1, nchar(x))
}
names.back <- substrRight(xnames, 2)
node.names <- paste0(names.front, ".", names.back)

logt <- function(x) x^3
xm <- sapply(xm, logt)
xm <- matrix(xm, ncol=length(node.names))
colnames(xm) <- node.names
rownames(xm) <- node.names

# Assign node shapes
#shape_list <- c("circle","square","triangle","diamond","heart","star","circle","square","triangle","diamond","heart","star","circle","square","triangle","diamond","heart","star","circle","square","triangle","diamond","heart","star","circle","square","triangle","diamond","heart","star","circle","square","triangle","diamond","heart","star")
node.shapes <- "circle"

# Assign node colors and shapes
node.shapes <- c()
node.colors <- c()
node_names_graph <- colnames(xn)
#color_list <- rainbow(length(mcl.listFull))
color_list <- rainbow(31)
for(j in 1:length(node_names_graph)){
	node.colors[j] <- "#BBBBBBFF"
	#node.shapes[j] <- shape_list[10]
	for(k in 1:length(mcl.listFull)){
		if(length(which(mcl.listFull[[k]]==node_names_graph[j]))>0){
			col_index <- ( ( k * 6 ) %% 31 ) + 1
			node.colors[j] <- color_list[col_index]
		}
	}
	#for(k in 1:length(mcl.list)){
		#if(length(which(mcl.list[[k]]==node_names_graph[j]))>0){
			#node.shapes[j] <- shape_list[k]
			#node.names[j] <- paste(k, "_", node.names[j])
		#} 
	#}
}
#node.colors <- "#999999FF"

fileOut <- file.path(dir.figures,paste("Level1ClustersGraph",batch_name,".png", sep=""));
png(file=fileOut, height=2160, width=3840, res = 300, pointsize = 7)
Q <- qgraph(xm)
qgraph(Q, layout="spring",
       border.width=0.65, 
       vsize=node.size,
       labels=node.names,
       shape=node.shapes,
       color=node.colors,
       label.scale=FALSE, 
       label.cex=1.2, 
       label.position=5, 
       legend=FALSE, 
       legend.cex=0.4,
       diag=TRUE,
       ###Modify maximum, minimum, cut, and colFactor to tweak image display
       maximum=1,
       minimum=0.02,
       #minimum=0.1,
       cut=.5,
       colFactor=.5,
       #Green edges only
       #negCol="#FFFFFF00",
       #Red edges only
       #posCol="#FFFFFF00",
       #edge.width=0.2,
	   usePCH = T,
       mar=c(1,1,1,1))
dev.off()

#Analysis of types of triads
z <- nrow(mcl.data)
count <- 0
total_triads <- c(0, 0, 0, 0)
triads <- list(c(), c(), c(), c())
print("Calculating p-valued triads....")
for(i in (1:(z-2))){
	for(j in ((i+1):(z-1))){
		for(k in ((j+1):z)){
			count <- 0
			# i vs j
			if(mcl.data[i,j]>0){
				count <- count + 1
			}
			if(mcl.data[i,j]==0){
				count <- -10
			}
			# i vs k
			if(mcl.data[i,k]>0){
				count <- count + 1
			}
			if(mcl.data[i,k]==0){
				count <- -10
			}
			
			# j vs k
			if(mcl.data[j,k]>0){
				count <- count + 1
			}
			if(mcl.data[j,k]==0){
				count <- -10
			}
			
			if(count >= 0){
				total_triads[count + 1] <- total_triads[count + 1] + 1
				triads[[count + 1]] <- rbind(triads[[count + 1]], c(i,j,k))
			}
		}
	}
}

#Analysis of types of triads
z <- nrow(mcl.data.original)
count.original <- 0
total_triads.original.original <- c(0, 0, 0, 0)
triads.original <- list(c(), c(), c(), c())
print("Calculating original triads....")
for(i in (1:(z-2))){
	#print(paste("Calculating originals for ", i))
	for(j in ((i+1):(z-1))){
		for(k in ((j+1):z)){
			count.original <- 0
			# i vs j
			if(mcl.data.original[i,j]>0){
				count.original <- count.original + 1
			}
			if(mcl.data.original[i,j]==0){
				count.original <- -10
			}
			# i vs k
			if(mcl.data.original[i,k]>0){
				count.original <- count.original + 1
			}
			if(mcl.data.original[i,k]==0){
				count.original <- -10
			}
			
			# j vs k
			if(mcl.data.original[j,k]>0){
				count.original <- count.original + 1
			}
			if(mcl.data.original[j,k]==0){
				count.original <- -10
			}
			
			if(count.original >= 0){
				total_triads.original.original[count.original + 1] <- total_triads.original.original[count.original + 1] + 1
				#triads.original[[count.original + 1]] <- rbind(triads.original[[count.original + 1]], c(i,j,k))
			}
		}
	}
}
print("Saving values....")
results.list <- list();
results.list[[length(results.list)+1]] <- "P-valued Values";
results.list[[length(results.list)+1]] <- paste(total_triads[1],total_triads[2],total_triads[3],total_triads[4], sep=",");
results.list[[length(results.list)+1]] <- "Original values";
results.list[[length(results.list)+1]] <- paste(total_triads.original.original[1],total_triads.original.original[2],total_triads.original.original[3],total_triads.original.original[4], sep=",");


#Analysis of types of triads
print("Calculating cliques")
z <- nrow(mcl.data)
clique_count <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0)
clique_members <- list(c(), c(), c(), c(), c(), c(), c(), c(), c(), c(), c(), c(),c(),c(),c(),c(),c())
c_idx <- c()
for(c_idx_1 in (1:(z-2))){
	#print(paste("Cliques calc", c_idx_1))
	for(c_idx_2 in ((c_idx_1+1):(z-1))){
		if(mcl.data[c_idx_1, c_idx_2] > 0){
		for(c_idx_3 in ((c_idx_2+1):z)) {
			if(mcl.data[c_idx_1, c_idx_3] > 0 && mcl.data[c_idx_2, c_idx_3] > 0){
				clique_count[3] <- clique_count[3] + 1;
			if(c_idx_3 < z){
			for(c_idx_4 in ((c_idx_3+1):z)) {
				if(mcl.data[c_idx_1, c_idx_4] > 0 && mcl.data[c_idx_2, c_idx_4] > 0 && mcl.data[c_idx_3, c_idx_4] > 0){
					clique_count[4] <- clique_count[4] + 1;
					clique_members[[4]] <- rbind(clique_members[[4]], c(c_idx_1, c_idx_2, c_idx_3, c_idx_4))
				if(c_idx_4 < z){
				for(c_idx_5 in ((c_idx_4+1):z)) {
					if(mcl.data[c_idx_1, c_idx_5] > 0 && mcl.data[c_idx_2, c_idx_5] > 0 && mcl.data[c_idx_3, c_idx_5] > 0
					   && mcl.data[c_idx_4, c_idx_5] > 0){
						clique_count[5] <- clique_count[5] + 1;
						clique_members[[5]] <- rbind(clique_members[[5]], c(c_idx_1, c_idx_2, c_idx_3, c_idx_4, c_idx_5))
					if(c_idx_5 < z){
					for(c_idx_6 in ((c_idx_5+1):z)) {
						if(mcl.data[c_idx_1, c_idx_6] > 0 && mcl.data[c_idx_2, c_idx_6] > 0 && mcl.data[c_idx_3, c_idx_6] > 0
						   && mcl.data[c_idx_4, c_idx_6] > 0 && mcl.data[c_idx_5, c_idx_6] > 0){
							clique_count[6] <- clique_count[6] + 1;
							clique_members[[6]] <- rbind(clique_members[[6]], c(c_idx_1, c_idx_2, c_idx_3, c_idx_4, c_idx_5, c_idx_6))
						if(c_idx_6 < z){
						for(c_idx_7 in ((c_idx_6+1):z)) {
							if(mcl.data[c_idx_1, c_idx_7] > 0 && mcl.data[c_idx_2, c_idx_7] > 0 && mcl.data[c_idx_3, c_idx_7] > 0
							   && mcl.data[c_idx_4, c_idx_7] > 0 && mcl.data[c_idx_5, c_idx_7] > 0 && mcl.data[c_idx_6, c_idx_7] > 0){
								clique_count[7] <- clique_count[7] + 1;
								clique_members[[7]] <- rbind(clique_members[[7]], c(c_idx_1, c_idx_2, c_idx_3, c_idx_4, c_idx_5, c_idx_6, c_idx_7))
							if(c_idx_7 < z){
							for(c_idx_8 in ((c_idx_7+1):z)) {
								if(mcl.data[c_idx_1, c_idx_8] > 0 && mcl.data[c_idx_2, c_idx_8] > 0 && mcl.data[c_idx_3, c_idx_8] > 0
								   && mcl.data[c_idx_4, c_idx_8] > 0 && mcl.data[c_idx_5, c_idx_8] > 0 && mcl.data[c_idx_6, c_idx_8] > 0
								   && mcl.data[c_idx_7, c_idx_8] > 0){
									clique_count[8] <- clique_count[8] + 1;
									clique_members[[8]] <- rbind(clique_members[[8]], c(c_idx_1, c_idx_2, c_idx_3, c_idx_4, c_idx_5, c_idx_6, c_idx_7, c_idx_8))
								if(c_idx_8 < z){
								for(c_idx_9 in ((c_idx_8+1):z)) {
									if(mcl.data[c_idx_1, c_idx_9] > 0 && mcl.data[c_idx_2, c_idx_9] > 0 && mcl.data[c_idx_3, c_idx_9] > 0
									   && mcl.data[c_idx_4, c_idx_9] > 0 && mcl.data[c_idx_5, c_idx_9] > 0 && mcl.data[c_idx_6, c_idx_9] > 0
									   && mcl.data[c_idx_7, c_idx_9] > 0 && mcl.data[c_idx_8, c_idx_9] > 0){
										clique_count[9] <- clique_count[9] + 1;
										clique_members[[9]] <- rbind(clique_members[[9]], c(c_idx_1, c_idx_2, c_idx_3, c_idx_4, c_idx_5, c_idx_6, c_idx_7, c_idx_8, c_idx_9))
									if(c_idx_9 < z){
									for(c_idx_10 in ((c_idx_9+1):z)) {
										if(mcl.data[c_idx_1, c_idx_10] > 0 && mcl.data[c_idx_2, c_idx_10] > 0 && mcl.data[c_idx_3, c_idx_10] > 0
										   && mcl.data[c_idx_4, c_idx_10] > 0 && mcl.data[c_idx_5, c_idx_10] > 0 && mcl.data[c_idx_6, c_idx_10] > 0
										   && mcl.data[c_idx_7, c_idx_10] > 0 && mcl.data[c_idx_8, c_idx_10] > 0 && mcl.data[c_idx_9, c_idx_10] > 0){
											clique_count[10] <- clique_count[10] + 1;
											clique_members[[10]] <- rbind(clique_members[[10]], c(c_idx_1, c_idx_2, c_idx_3, c_idx_4, c_idx_5, c_idx_6, c_idx_7, c_idx_8, c_idx_9, c_idx_10))
										if(c_idx_10 < z){
										for(c_idx_11 in ((c_idx_10+1):z)) {
											if(mcl.data[c_idx_1, c_idx_11] > 0 && mcl.data[c_idx_2, c_idx_11] > 0 && mcl.data[c_idx_3, c_idx_11] > 0
											   && mcl.data[c_idx_4, c_idx_11] > 0 && mcl.data[c_idx_5, c_idx_11] > 0 && mcl.data[c_idx_6, c_idx_11] > 0
											   && mcl.data[c_idx_7, c_idx_11] > 0 && mcl.data[c_idx_8, c_idx_11] > 0 && mcl.data[c_idx_9, c_idx_11] > 0
											   && mcl.data[c_idx_10, c_idx_11] > 0){
												clique_count[11] <- clique_count[11] + 1;
												clique_members[[11]] <- rbind(clique_members[[1]], c(c_idx_1, c_idx_2, c_idx_3, c_idx_4, c_idx_5, c_idx_6, c_idx_7, c_idx_8, c_idx_9, c_idx_10, c_idx_11))
											if(c_idx_11 < z){
											
											}
											}
										}
										}
										}
									}
									}
									}
								}
								}
								}
							}
							}
							}
						}
						}
						}
					}
					}
					}
				}
				}
				}
			}
			}
			}
		}
		}
	}
}
results.list[[length(results.list)+1]] <- "Cliques of,3,4,5,6,7,8,9,10,11"
results.list[[length(results.list)+1]] <- paste("",clique_count[3], clique_count[4], clique_count[5], clique_count[6], clique_count[7], clique_count[8], clique_count[9], clique_count[10], clique_count[11], sep=",")

print("Calculating cliques 2")

next_level_clique <- function (corr_data, level, clique_count, clique_members, c_idx){
	z <- nrow(corr_data)
	for(c_idx_lcl in ((c_idx[level-1]+1):z)) {
		c_idx[level]<-c_idx_lcl
		is_clique <- TRUE
		cluster_list <- c()
		for(i in (1:(level-1))){
			is_clique <- is_clique && (corr_data[c_idx[i], c_idx[level]] > 0)
			cluster_list <- c(cluster_list, c_idx[i])
		}
		cluster_list <- c(cluster_list, c_idx[level])
		if(is_clique){
			if(is.na(clique_count[level])){clique_count[level]<-0}
			clique_count[level] <- clique_count[level] + 1;
			#if(is.na(clique_members[[level]])){clique_members[[level]]<-c()}
			if(length(clique_members)>=level)
				clique_members[[level]] <- rbind(clique_members[[level]], cluster_list)
			else
				clique_members[[level]] <- cluster_list
			if(c_idx[level] < z){
				r <- next_level_clique(corr_data, level+1, clique_count, clique_members, c_idx)
				clique_count <- r$clique_count 
				clique_members <- r$clique_members 
				c_idx <- r$c_idx
			}
		}
	}
	return(list(clique_count=clique_count, clique_members=clique_members, c_idx=c_idx))
}

z <- nrow(mcl.data)
clique_count <- c(0)
clique_members <- list(c())
c_idx <- c(1)
for(j in (1:(z-1))){
	c_idx[1]<-j
	r <- next_level_clique(mcl.data, 2, clique_count, clique_members, c_idx)
	clique_count <- r$clique_count 
	clique_members <- r$clique_members 
	c_idx <- r$c_idx
}

str_clique <- "Cliques new of"
str_count <- " "
for(j in 1:length(clique_count)){
	str_clique <- paste(str_clique, j, sep=",")
	str_count <- paste(str_count, clique_count[j], sep=",") 
}
results.list[[length(results.list)+1]] <- str_clique
results.list[[length(results.list)+1]] <- str_count
fileOut <- file.path(dir.figures,paste("TriadsData",batch_name,".csv", sep=""));
lapply(results.list, cat, "\n", file=fileOut, append=TRUE);


