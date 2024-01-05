#!/usr/bin/Rscript
#Author: Tak Lee, tlee@mpipz.mpg.de
#script used to identify Differentially expressed genes with DESeq2

library(sva)
library(DESeq2)
library(ggplot2)
library(tidyverse)
library(limma)
library(foreach)
library(ggforce)

cnts <- read.table("counts.txt", sep='\t',header=T,row.names=1,check.names=F)
cnts$Length <- NULL
coldata <- read.table("coldata.txt", sep='\t',header=T,row.names=1,check.names=F)

ddsMat <- DESeqDataSetFromMatrix(countData = cnts,
                                 colData = coldata,
                                 design = ~ condition)

dds <- estimateSizeFactors(ddsMat)
dat  <- counts(dds, normalized = TRUE)

###############filtering the lowly expressed genes################
geneidx <- rowSums(dat >= 10) >= 3
#geneidx  <- rowMeans(dat) > 1
dat <- dat[geneidx,]

###############finding the surrogate variables################
mod  <- model.matrix(~ condition, colData(dds))
mod0 <- model.matrix(~ 1,colData(dds))
#manual detection of surrogate variable
#svseq <- svaseq(dat, mod, mod0, n.sv = 6)

dmicnt <- rep(0,50)
if(FALSE){
#auto detection of surrogate variables
for (i in 1:5){
	svobj <- sva(dat, mod, mod0)
	dmicnt[svobj$n.sv] <- dmicnt[svobj$n.sv] +1
}
svnum <- which.max(dmicnt)
}

svnum <- 4
ddssva <- ""
for (k in 4:svnum){
	ddssva <- dds 
	fmla <- "~" 
	if (k > 0){
	for (j in 1:k){
	    svseq <- svaseq(dat, mod, mod0, n.sv = k)
	    newcol <- svseq$sv[,j]
	    colData(ddssva) <- cbind(colData(ddssva), newcol)
	    names(colData(ddssva))[length(names(colData(ddssva)))] <- paste0("SV",j)
	    if (j == 1){fmla=paste0(fmla,paste0("SV",j))}
	    else {fmla=paste(fmla,paste0("SV",j),sep="+")}
		}
	fmla <- paste(fmla,"condition",sep="+")
	} else {fmla <- paste(fmla,"condition",sep="")}
	print(fmla)
	design(ddssva) <- as.formula(fmla)
	samples <- coldata$condition
	design <- model.matrix(~samples)

	################plotting PCA#################
	#vst normalization of data
	vsd <- vst(dds)[geneidx,]
	samples <- coldata$condition
	if (k>0){
	design <- model.matrix(~samples)
	rldData <- assay(vsd) %>%
		  removeBatchEffect(covariates = svseq$sv, design = design)
	}
	else {rldData <- assay(vsd) }
	sampleIdx <- 1:length(samples)
	group <- c("R108_spotinoc")#"all","subset1","subset2")
	top <- 500
	allgn <- length(rownames(rldData))
	allidx <- list(sampleIdx)#allsampleIdx,group1idx,group2idx)
	for(i in 1:length(allidx)){
		idx <- allidx[[i]]
		pca <- prcomp(t(rldData[, idx]))
		length(rownames(rldData))
		for (n in c(top,allgn)){
		#for top 500 genes that contribute the most to the variance of PC1 and PC2.
			p1 <- order(abs(pca$rotation[,1]),decreasing=TRUE)[1:n]
			p2 <- order(abs(pca$rotation[,2]),decreasing=TRUE)[1:n]
			g1 <- names(pca$rotation[,1][p1])
			g2 <- names(pca$rotation[,2][p2])
			p1p2 <- union(g1,g2)
			pca2 <- prcomp(t(rldData[p1p2, idx]))
			pc1 <- pca2$x[,1]
			pc2 <- pca2$x[,2]
			percentVar <- pca2$sdev^2/sum(pca2$sdev^2)
			percentVar <- round(100 * percentVar)

			pcaData <- data.frame(PC1 = pc1, PC2 = pc2, Group = colData(dds)[idx, 1], ID = rownames(colData(dds))[idx])
			ggplot(pcaData, aes(x = PC1, y = PC2, colour = Group, label = ID)) +
				geom_point(size = 3) +
				xlab(paste0("PC1: ",percentVar[1],"% variance")) +
				ylab(paste0("PC2: ",percentVar[2],"% variance")) +
				geom_mark_ellipse(aes(color = Group, label=Group)) +
				#scale_colour_manual(name = 'SynCom',values = cols[colorIdx]) +
				#stat_ellipse(aes(x = PC1, y = PC2, group = Group), type = 't', linetype = 2, level = 0.8) +
				coord_fixed(1) +
				theme_classic() +
				geom_text(aes(label=ID),vjust=2) +
				theme(plot.title = element_text(hjust = 0.5, size = 12, face = 'bold'),
				      legend.text.align = 0,
				      axis.text = element_text(size = 13),
				      axis.title = element_text(size = 14),
				      legend.text=element_text(size= 13),
				      legend.title = element_text(size = 14))
			if(n == allgn){
				outpdf <- paste(group[i],"_PCA_sva",k,".pdf",sep="")
				ggsave(outpdf, width = 30,height=30)
			} else {
				outpdf <- paste(group[i],"_PCA_top",n,"_sva",k,".pdf",sep="")
				ggsave(outpdf, width = 30,height=30)
			}
			#outjpg <- paste(group[i],"_PCA_top500.jpg",sep="")
			#ggsave(outjpg, width = 30,height=30)
		}

	}
}
#if(FALSE){
##################Differential Expression of genes###################
ddssva <- DESeq(ddssva[geneidx, ])
contrasts <- list(
		  c("mock_12hrs","spotinoc_12hrs"),
		  c("mock_24hrs","spotinoc_24hrs"),
		  c("mock_72hrs","spotinoc_72hrs")
		  )

for(i in 1:length(contrasts)){
	ctrl <- contrasts[[i]][1]
	trt <- contrasts[[i]][2]
	res <- results(ddssva, contrast=c("condition",trt,ctrl))
	outf <- paste("Mtr5.0_R108_",trt,"_diff.txt",sep="")
	write.table(res,file=outf,quote=F, sep="\t")
	#print(outf)
	}
#}
