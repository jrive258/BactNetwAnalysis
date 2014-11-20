plot_st.heatmap <- function(M, name, noFolder=FALSE) {
#########################################################
### Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(-1,0,length=100),  # for red
#  seq(0,0.5,length=100),              # for yellow
#  seq(0.5,1,length=100))              # for green
#col_breaks = c(seq(min(0,min(M)), max(M), length=300));
col_breaks = c(seq(-1, 1, length=300));

# creates a 5 x 5 inch image
if(noFolder){
	file <- name;
}else{
file <- file.path(dir.figures, name);
}

colores <- rainbow(ncol(M))

png(file,                # create PNG for the heat map        
  width = 10*300,        # 5 x 300 pixels
  height = 10*300,
  res = 300,            # 300 pixels per inch
  pointsize = 6)        # smaller font size

heatmap.2(M, 
  #cellnote = mat_data,  # same data set for cell labels
  main = "Correlation", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns on trace lines inside the heat map
  margins =c(15,15),    # widens margins around plot
  col=my_palette,       # use on color palette defined earlier 
  breaks=col_breaks,    # enable color transition at specified limits
  keysize=0.5,
  #ColSideColors = colores,
  #RowSideColors = colores,
  dendrogram="none",    # no dendograms
  Rowv="NA",            # turn off row clustering
  Colv="NA")            # turn off column clustering

dev.off()               # close the PNG device
}

plot_st.heatmap_w_colors<- function(M, name, noFolder=FALSE, colorList) {
#########################################################
### Customizing and plotting the heat map
#########################################################

# creates a own color palette from red to green
my_palette <- colorRampPalette(c("red", "black", "green"))(n = 299)

# (optional) defines the color breaks manually for a "skewed" color transition
#col_breaks = c(seq(-1,0,length=100),  # for red
#  seq(0,0.5,length=100),              # for yellow
#  seq(0.5,1,length=100))              # for green
#col_breaks = c(seq(min(0,min(M)), max(M), length=300));
col_breaks = c(seq(-1, 1, length=300));

# creates a 5 x 5 inch image
if(noFolder){
	file <- name;
}else{
file <- file.path(dir.figures, name);
}

colores <- colorList

png(file,                # create PNG for the heat map        
  width = 10*300,        # 5 x 300 pixels
  height = 10*300,
  res = 300,            # 300 pixels per inch
  pointsize = 6)        # smaller font size

heatmap.2(M, 
  #cellnote = mat_data,  # same data set for cell labels
  main = "Correlation", # heat map title
  notecol="black",      # change font color of cell labels to black
  density.info="none",  # turns off density plot inside color legend
  trace="none",         # turns on trace lines inside the heat map
  margins =c(15,15),    # widens margins around plot
  col=my_palette,       # use on color palette defined earlier 
  breaks=col_breaks,    # enable color transition at specified limits
  keysize=0.5,
  #ColSideColors = colores,
  RowSideColors = colores,
  dendrogram="none",    # no dendograms
  Rowv="NA",            # turn off row clustering
  Colv="NA")            # turn off column clustering

  #lines(c(0,0.1,0.2,0.3,0.4,0.5), c(0,0.1,0.2,0.3,0.4,0.5), type="b", col="#0000FFFF")
  
dev.off()               # close the PNG device
}

plot.network <- function( name, mcl.data, n_color, n_size, n_names, node_pos ){
	require(plotrix)
	file <- file.path(dir.figures, name);
	png(file,                # create PNG for the heat map        
	  width = 10*300,        # 5 x 300 pixels
	  height = 10*300,
	  res = 300,            # 300 pixels per inch
	  pointsize = 6)        # smaller font size
	 
	par(mar=c(1,1,1,1)) #Better margins
	plot(c(0,1),c(0,1),col="white",axes=FALSE,ann=FALSE)
	 
	for(i in (1:(ncol(mcl.data)-1))){
		for(j in ((i+1):ncol(mcl.data))){
			n_alpha = abs(mcl.data[i,j])
			l_color = "white"
			if(mcl.data[i,j] > 0)
				l_color <- rgb(0,1,0,alpha=n_alpha)
			else
				l_color <- rgb(1,0,0,alpha=n_alpha)
			if(mcl.data[i,j] != 0)
				lines(c(node_pos[i,1],node_pos[j,1]), c(node_pos[i,2],node_pos[j,2]), type="l", col=l_color)
		}
	}
	
	for(i in 1:nrow(node_pos)){
		draw.circle(node_pos[i,1],node_pos[i,2],n_size[i], col=n_color[i], border="black")
		text(node_pos[i,1],node_pos[i,2], labels = n_names[i])
	}
	
	
	dev.off()

}

plot.graph.split <- function( name, mcl.data, group, subgroup ){
	#First calculate the position of every node
	#Calculate the number of rows needed
	MAX_COLS <- 8
	GAP <- 3 / ((MAX_COLS-1)*10)
	tot_left <- 0
	tot_right <- 0
	for(i in 1:length(subgroup)){
		if(subgroup[i]!=0){
			tot_left <- tot_left + ceiling(length(which(group==((i*2)-1))) / MAX_COLS)
			tot_right <- tot_right + ceiling(length(which(group==(i*2))) / MAX_COLS)
			#print(paste(tot_left, tot_right, sep=" - "))
		}
	}
	node_pos <- matrix(nrow=nrow(mcl.data), ncol=2)
	#Now calculate the x,y positions for each left item
	curr_top_left <- 0
	curr_top_right <- 0
	for(i in 1:length(subgroup)){
		if(subgroup[i]!=0){
			items_left <- which(group==((i*2)-1))
			if(length(items_left)>0){
				need_rows <- ceiling(length(which(group==((i*2)-1))) / MAX_COLS)
				for(j in 1:length(items_left)){
					posx <- ceiling(j/need_rows)
					posy <- ( (j-1) %% need_rows ) + 1
					node_pos[items_left[j], 1] <- (0.35 - ((posx-1)*GAP))
					node_pos[items_left[j], 2] <- (1 - ((curr_top_left + posy)*(1/tot_left)))
				}
			}
			items_right <- which(group==(i*2))
			if(length(items_right)>0){
				need_rows <- ceiling(length(which(group==(i*2))) / MAX_COLS)
				for(j in 1:length(items_right)){
					posx <- ceiling(j/need_rows)
					posy <- ( (j-1) %% need_rows ) + 1
					node_pos[items_right[j], 1] <- (0.65 + ((posx-1)*GAP))
					node_pos[items_right[j], 2] <- (1 - ((curr_top_right + posy)*(1/tot_right)))
				}
			}
			curr_top_left <- curr_top_left + ceiling(length(which(group==((i*2)-1))) / MAX_COLS)
			curr_top_right <- curr_top_right + ceiling(length(which(group==(i*2))) / MAX_COLS)
		}
	}
	
	color_list <- rainbow(31)
	n_color <- c(1:nrow(node_pos)) - c(1:nrow(node_pos))
	for(j in 1:length(group)){
		n_color[j] <- color_list[( ( group[j] * 6 ) %% 31 ) + 1]
	}
	n_size <- c(1:nrow(node_pos)) - c(1:nrow(node_pos)) + 0.015
	
	xnames <- colnames(mcl.data)
	names.front <- substr(xnames, 1, 5)
	substrRight2 <- function(x, n){
	substr(x, nchar(x)-n+1, nchar(x))
	}
	names.back <- substrRight2(xnames, 2)
	n_names <- paste0(names.front, ".", names.back)
	
	plot.network(name, mcl.data, n_color, n_size, n_names, node_pos)
}

