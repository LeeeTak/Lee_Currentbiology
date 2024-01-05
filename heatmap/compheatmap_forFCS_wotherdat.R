#!/usr/bin/Rscript
#script written by Tak Lee, tlee@mpipz.mpg.de
#used for "hierarchical" clustering a set of data and plotting a heatmap
#usage: Rscript compheatmap_withclustering.R input_matrix output_filename heatmap_title(optional)

library(pheatmap)
library(RColorBrewer)
library(NbClust)
library(dplyr)
library(stringr)
library(ComplexHeatmap)
library(circlize)

args <- commandArgs(TRUE)
#read data
dat <- read.table("wtspot_DEGs_FCS.txt",sep='\t',header=T,row.names=1,check.names=F)
#finding the best number of clusters using silhouette index with NbClust, clustering distance and method could be modified
#reference:https://cran.r-project.org/web/packages/NbClust/NbClust.pdf
spdist <- as.dist(1-cor(t(dat),method="spearman"))
cl <- NbClust(dat,distance=NULL,diss=spdist,method="kmeans",min.nc=1,max.nc=30,index="silhouette")
clnum <- cl$Best.nc[1]

#determine the number of k with gap statistics:
#gap_stat <- clusGap(filtdat, FUN = kmeans, nstart = 4,K.max = 40, B = 500)
#k <- maxSE(f = gap_stat$Tab[, "gap"], SE.f = gap_stat$Tab[, "SE.sim"])

clusters <- as.data.frame(cl$Best.partition)
colnames(clusters) <- c("cluster")
ordcl <- clusters %>% arrange(cluster)

#write out genelists of each clusters
dir.create("clusters")
for(x in seq(1,clnum)){
    clgs <- as.data.frame(rownames(filter(ordcl, cluster == x)))
    write.table(clgs,file=paste("clusters/","wtspotinoc_cluster",x,".txt",sep=""),sep="\t",row.names=F,quote=F,col.names=F)
}
#print(clnum)

#determine color range based on min,max values
print(max(dat))
print(min(dat))
lgdcols <- colorRamp2(c(min(dat),-1.5,0,1.5,max(dat)),c("blue","#0080ff","grey45","yellow","red"))
lgd <- Legend(col_fun = lgdcols, title = "log2FC",direction = "horizontal")
#colsplit <- rep(c("WT","lsh1","lsh1/\nlsh2","pLjUBI:\nLSH1/2"),c(2,2,2,1))
#colsplit_factor <- factor(colsplit,levels=unique(colsplit))
#rowsplit <- rep(c("auxin signalling","cytokinin signalling","cell cycle"),c(8,7,3))
#rowsplit_factor <- factor(rowsplit,levels=unique(rowsplit))
#k <- 17
#output pdf file info
#adjust the size and margin of the pdf file
pdf(args[2],height=24,width=10)
#plot heatmap by adjusting the parameters
ht <- Heatmap(
	as.matrix(dat),
	col=lgdcols,
	#rect_gp = gpar(col = "white"),
	cluster_columns = F,
	cluster_rows = T,
	clustering_distance_rows = "spearman",
	clustering_method_rows="complete",
	#column_split = colsplit_factor,
	width = ncol(dat)*unit(10, "mm"), 
	height = nrow(dat)*unit(0.02, "mm"),
	#row_km = k,
	#row_split = rowsplit_factor,
	row_split=clusters$cluster,
	show_row_names = F,
	#column_order = colorder,
	#border=F,
	border_gp = gpar(col = "black"),
	#column_gap = unit(2.2,"mm"),
	heatmap_legend_param = list(title="log2FC"),
	show_heatmap_legend = F,
	use_raster = F,
	row_dend_reorder = T
	)

dat2 <- read.table("wtspotDEGs_lsh1lsh2_FCS.txt",sep='\t',header=T,row.names=1,check.names=F)
ht2 <- Heatmap(
	       as.matrix(dat2),
	       col=lgdcols,
	       column_labels = c("lsh1lsh2_24hrs","lsh1lsh2_72hrs"),
	       show_row_names=F,
	       cluster_columns = F,
	       border_gp = gpar(col = "black"),
	       show_heatmap_legend = F,
	       use_raster = F
	       )
dat3 <- read.table("wtspotDEGs_LSHchiptargets.txt",sep='\t',header=T,row.names=1,check.names=F)
#chipcolors <- structure(c("grey45","#b3b300","#ffff66" ),names = c("<2","2","3"))
ht3 <- Heatmap(
	       as.matrix(dat3),
	       show_row_names=F,
	       border_gp = gpar(col = "black"),
	       col = c("grey45","#66ff66","#ff66ff"),
	       heatmap_legend_param = list(title="ChIP\ntargets\nby confidence",labels = c("low/\nnon-targets","moderate","high")),
	       column_labels = c("LSH1_ChIP"),
	       use_raster = F
	       )

dat4 <- read.table("wtspotDEGs_LSHOE_FCS.txt",sep="\t",header=T,row.names=1,check.names=F)

ht4 <- Heatmap(
	       as.matrix(dat4),
	       col=lgdcols,
	       show_row_names=F,
	       cluster_columns = F,
	       border_gp = gpar(col = "black"),
	       show_heatmap_legend = F,
	       use_raster = F
	       )


draw(ht+ht3+ht2+ht4,padding = unit(c(10, 10, 10, 10), "mm"))
draw(lgd,x = unit(0.8,"npc"), y = unit(0.12,"npc"))
dev.off()
#klusters <- row_order(hm)
#dir.create("clusters")
#for (c in seq(1,k)) {
#    genes <- rownames(dat[klusters[[c]],])
#    fn <- paste("clusters/wtspotinoc_cluster",c,"_genes.txt",sep="")
#    write.table(genes,file=fn,col.names=F,row.names=F,quote=F)
#}




