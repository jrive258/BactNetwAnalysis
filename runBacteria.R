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


if(1 == 2){ #Eliminate if you want to count the triads
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


results.list <- list();
results.list[[length(results.list)+1]] <- "P-valued Triad Counts";
results.list[[length(results.list)+1]] <- paste(total_triads[1],total_triads[2],total_triads[3],total_triads[4], sep=",");

#Calculating Cliques of any size
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
			if(level > 2){
				if(length(clique_members)>=level)
					clique_members[[level]] <- rbind(clique_members[[level]], cluster_list)
				else
					clique_members[[level]] <- cluster_list
			}
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

#Add the number of cliques to the result list
str_clique <- "Cliques new of"
str_count <- " "
for(j in 1:length(clique_count)){
	str_clique <- paste(str_clique, j, sep=",")
	str_count <- paste(str_count, clique_count[j], sep=",") 
}
results.list[[length(results.list)+1]] <- str_clique
results.list[[length(results.list)+1]] <- str_count

}




#Finding the bicliques
analyze_biclique2_2 <- function(corr_data, total_bicliques, bicliques, idx1, idx2, idx3, idx4){
	#2_2 biclique
	if(corr_data[idx1,idx3] > 0 && corr_data[idx1,idx4] > 0 && corr_data[idx2,idx3] > 0 && corr_data[idx2,idx4] > 0){
		if(corr_data[idx1,idx2] <= 0 && corr_data[idx3, idx4] <= 0){
			#Pure 2_2
			total_bicliques[1] <- total_bicliques[1] + 1
			bicliques[[1]] <- rbind(bicliques[[1]], c(idx1, idx2, idx3, idx4))
		}else if(corr_data[idx1,idx2] > 0 && corr_data[idx3, idx4] > 0){
			#4 clique do nothing
		}else{
			#Dirty 2_2
			total_bicliques[2] <- total_bicliques[1] + 1
			bicliques[[2]] <- rbind(bicliques[[2]], c(idx1, idx2, idx3, idx4))
		}
	}
	#2_2 Rival
	if(corr_data[idx1,idx3] < 0 && corr_data[idx1,idx4] < 0 && corr_data[idx2,idx3] < 0 && corr_data[idx2,idx4] < 0){
		total_bicliques[5] <- total_bicliques[5] + 1
		bicliques[[5]] <- rbind(bicliques[[5]], c(idx1, idx2, idx3, idx4))
	}
	return (list(total=total_bicliques, lists=bicliques));
}

analyze_biclique3_1 <- function(corr_data, total_bicliques, bicliques, idx1, idx2, idx3, idx4){
	#3_1 biclique
	if(corr_data[idx1,idx4] > 0 && corr_data[idx2,idx4] > 0 && corr_data[idx3,idx4] > 0){
		if(corr_data[idx1,idx2] <= 0 && corr_data[idx1, idx3] <= 0 && corr_data[idx2, idx3] <= 0){
			#Pure 3_1
			total_bicliques[3] <- total_bicliques[3] + 1
			bicliques[[3]] <- rbind(bicliques[[3]], c(idx1, idx2, idx3, idx4))
		}else if(corr_data[idx1,idx2] > 0 && corr_data[idx1, idx3] > 0 && corr_data[idx2, idx3] > 0){
			#4 clique do nothing
		}else{
			#Dirty 3_1
			total_bicliques[4] <- total_bicliques[4] + 1
			bicliques[[4]] <- rbind(bicliques[[4]], c(idx1, idx2, idx3, idx4))
		}
	}
	#3_1 Rival
	if(corr_data[idx1,idx4] < 0 && corr_data[idx2,idx4] < 0 && corr_data[idx3,idx4] < 0){
		total_bicliques[6] <- total_bicliques[6] + 1
		bicliques[[6]] <- rbind(bicliques[[6]], c(idx1, idx2, idx3, idx4))
	}
	return (list(total=total_bicliques, lists=bicliques));
}

if(1 == 2){
z <- nrow(mcl.data)
count <- 0
total_bicliques <- c(0, 0, 0, 0, 0, 0)  #2_2 pure, 2_2 dirty, 3_1 pure, 3_1 dirty, 2_2 rival, 3_1 rival
bicliques <- list(c(), c(), c(), c(), c(), c())
print("Calculating bicliques....")
for(i in (1:(z-3))){
	print(paste("Biclique ", i))
	for(j in ((i+1):(z-2))){
		for(k in ((j+1):(z-1))){
			for(l in ((k+1):z)){
				r <- analyze_biclique2_2(mcl.data, total_bicliques, bicliques, i,j,k,l)
				total_bicliques <- r$total
				bicliques <- r$lists
				
				r <- analyze_biclique2_2(mcl.data, total_bicliques, bicliques, i,k,j,l)
				total_bicliques <- r$total
				bicliques <- r$lists
				
				r <- analyze_biclique2_2(mcl.data, total_bicliques, bicliques, i,l,j,k)
				total_bicliques <- r$total
				bicliques <- r$lists
				
				r <- analyze_biclique3_1(mcl.data, total_bicliques, bicliques, i,j,k,l)
				total_bicliques <- r$total
				bicliques <- r$lists
				
				r <- analyze_biclique3_1(mcl.data, total_bicliques, bicliques, i,j,l,k)
				total_bicliques <- r$total
				bicliques <- r$lists
				
				r <- analyze_biclique3_1(mcl.data, total_bicliques, bicliques, i,k,l,j)
				total_bicliques <- r$total
				bicliques <- r$lists
				
				r <- analyze_biclique3_1(mcl.data, total_bicliques, bicliques, j,k,l,i)
				total_bicliques <- r$total
				bicliques <- r$lists
			}
		}
	}
}

#Add the number of bicliques to the result list
str_clique <- "Bicliques,2 2 Pure,2 2 Dirty,3 1 Pure,3 1 Dirty,2 2 Rival,3 1 Rival"
str_count <- " "
for(j in 1:length(total_bicliques)){
	str_count <- paste(str_count, total_bicliques[j], sep=",") 
}
results.list[[length(results.list)+1]] <- str_clique
results.list[[length(results.list)+1]] <- str_count


#Add the list of cliques
for(j in 3:length(clique_members)){
	results.list[[length(results.list)+1]] <- paste("Cliques of size ", j)
	for(k in 1:clique_count[j]){
		str_clique <- "";
		str_clique_name <- "";
		if(clique_count[j] == 1){
			for(l in 1:length(clique_members[[j]])){
				str_clique <- paste(str_clique, clique_members[[j]][l], sep=",")
				str_clique_name <- paste(str_clique_name, colnames(mcl.data)[clique_members[[j]][l]], sep=",")
			}
		}
		else{
			for(l in 1:ncol(clique_members[[j]])){
				str_clique <- paste(str_clique, clique_members[[j]][k,l], sep=",")
				str_clique_name <- paste(str_clique_name, colnames(mcl.data)[clique_members[[j]][k,l]], sep=",")
			}
		}
		results.list[[length(results.list)+1]] <- str_clique
		results.list[[length(results.list)+1]] <- str_clique_name
	}
}

#Add the list of bicliques
types <- c("2 2 Pure", "2 2 Dirty", "3 1 Pure", "3 1 Dirty", "2 2 Rival", "3 1 Rival")
for(j in 1:6){
	
	results.list[[length(results.list)+1]] <- paste("BiCliques of type ", types[j])
	for(k in 1:total_bicliques[j]){
		str_clique <- "";
		str_clique_name <- "";
		if(total_bicliques[j] == 0){
		}else if(total_bicliques[j] == 1){
			for(l in 1:length(bicliques[[j]])){
				str_clique <- paste(str_clique, bicliques[[j]][l], sep=",")
				str_clique_name <- paste(str_clique_name, colnames(mcl.data)[bicliques[[j]][l]], sep=",")
			}
		}
		else{
			for(l in 1:ncol(bicliques[[j]])){
				str_clique <- paste(str_clique, bicliques[[j]][k,l], sep=",")
				str_clique_name <- paste(str_clique_name, colnames(mcl.data)[bicliques[[j]][k,l]], sep=",")
			}
		}
		results.list[[length(results.list)+1]] <- str_clique
		results.list[[length(results.list)+1]] <- str_clique_name
	}
}

#Storing the results
fileOut <- file.path(dir.figures,paste("TriadsData",batch_name,".csv", sep=""));
lapply(results.list, cat, "\n", file=fileOut, append=TRUE);
}

#Check for structural balance
#1. Initialization
z <- nrow(mcl.data)
group <- c(1:z) - c(1:z)
n_stack <- group
n_stack[1] <- 1
c_stack <- 2
group[1] <- 1
subgroup <- c(1:z) - c(1:z)  #Indicates if the subgroup is balanced or not
subgroup[1] <- 1
c_subg <- 1

for( i_stack in 1:z ){
	if(n_stack[i_stack] == 0){ #No nodes in the stack (yet not every node has been visited)
		#I have to find the next node that has not been visited
		c_subg <- c_subg + 1
		subgroup[c_subg] <- 1
		for(k in 1:z){
			if(group[k] == 0){
				n_stack[i_stack] <- k
				c_stack <- i_stack + 1
				group[k] <- (c_subg * 2) - 1
				break
			}
		}
	}
	
	#If the stack was empty we made sure we filled it before we got here!!!
	if( n_stack[i_stack] != 0){
		#Look at the children of node n_stack[i_stack]
		for( j in 1:z){
			if(j != n_stack[i_stack] && mcl.data[n_stack[i_stack],j]>0){
				#The relation is positive, nodes should be on the same group
				if(group[j] == 0){
					group[j] <- group[n_stack[i_stack]]
					n_stack[c_stack] <- j
					c_stack <- c_stack + 1
				}else if(group[j] !=0){
					if(group[j] != group[n_stack[i_stack]]){
						subgroup[c_subg] <- -1
						print(paste("Failed at ", j, n_stack[i_stack]))
					}
				}
		 	}else if(j != n_stack[i_stack] && mcl.data[n_stack[i_stack],j]<0){
				#The relation is positive, nodes should be on the same group
				if(group[j] == 0){
					group[j] <- ( ((group[n_stack[i_stack]] - (c_subg * 2))+1) %% 2 ) * (-1) + (c_subg * 2)
					n_stack[c_stack] <- j
					c_stack <- c_stack + 1
				}else if(group[j] !=0){
					calc_group <- ( ((group[n_stack[i_stack]] - (c_subg * 2))+1) %% 2 ) * (-1) + (c_subg * 2)
					if(group[j] != calc_group){
						subgroup[c_subg] <- -1
						print(paste("Failed at ", j, n_stack[i_stack]))
					}
				}
			}
		}
	}
}

plot.graph.split( paste("SplitGraph",batch_name,".png", sep=""), mcl.data, group, subgroup )

left <- c()
right <- c()
for(i in 1:c_subg){
	aa <- (i*2) - 1
	bb <- (i*2)
	left <- c(left, colnames(mcl.data)[which(group==aa)])
	right <- c(right, colnames(mcl.data)[which(group==bb)])
}
split_file <- matrix(nrow=max(length(left),length(right)),ncol=2)
split_file[1:length(left),1] <- left
split_file[1:length(right),1] <- right

fileOut <- file.path(dir.figures,paste("SplitFile_", batch_name,".csv", sep=""))
write.table(split_file, file=fileOut, sep=",", append=TRUE, col.names=NA, na="") #col.names=NA is used to align column names correctly


for(i in 1:c_subg){
	aa <- (i*2) - 1
	bb <- (i*2)
	ww <- c(which(group==aa), which(group==bb))
	if(length(ww) > 1)
		plot_st.heatmap(mcl.data[ww, ww], paste("StructureBalance_",batch_name,"_", i, ".png", sep=""), noFolder=FALSE);
}

