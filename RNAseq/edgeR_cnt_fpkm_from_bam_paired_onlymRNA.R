#!/usr/bin/Rscript
#Author: Tak Lee, tlee@mpipz.mpg.de
#using edgeR to produce the raw count and rpkm/fpkm files
#usage: Rscript edgeR_cnt_rpkm_from_bam.R listfile threadN

options(echo=TRUE) 

#required R packages. They can be downloaded from source("http://bioconductor.org/biocLite.R")
library(Rsubread); library(limma); library(edgeR); library(Rsamtools)
options(digits=2)
args <-commandArgs(TRUE)

#read in all input bam files
#inputs <- dir(".", "bam$")
inputs <- scan('bamlist.txt',what='character')
#inputs <- as.vector(inputs)
#number of threads
#tn <- args[1]

#do featureCounts
fc <- featureCounts(
		    files=inputs,
		    annot.ext="/netscratch/dep_psl/grp_rgo/taklee/genomes/Mtruncatula/MtrunA17r5.0-ANR/20220708_MtrunA17r5.0-ANR-EGN-r1.9_FunctionalAnnotation/annotation/MtrunA17r5.0-ANR-EGN-r1.9_onlymRNA.gtf",
		    isGTFAnnotationFile = TRUE,
		    GTF.featureType = "exon",
		    GTF.attrType = "gene_id",
		    nthreads=args[1],
		    isPairedEnd=TRUE
		    )
#don't add pseudocount
#fc$counts <- fc$counts + 1 
x <- DGEList(counts=fc$counts, genes=fc$annotation[,c("GeneID","Length")])

#print out results
outd <- "EdgeR_RESULT_mRNA"
dir.create(outd)
cntout <- paste(outd, "/counts.txt", sep="")
rpkmout <- paste(outd, "/rpkm.txt", sep="")

write.table(x=data.frame(fc$annotation[,c("GeneID","Length")],fc$counts,stringsAsFactors=FALSE),cntout,quote=FALSE,sep="\t",row.names=FALSE)

x_rpkm <- rpkm(x,x$genes$Length)
write.table(x_rpkm, rpkmout, col.names=NA, quote=FALSE, sep="\t")





