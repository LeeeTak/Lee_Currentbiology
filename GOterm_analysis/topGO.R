#!/usr/bin/Rscript
#Author: Tak Lee, tlee@mpipz.mpg.de
#script to run topGO for GO term enrichment analysis

library(topGO)
library(ggplot2)

g2go <- readMappings(file = "Mtr5.0_r1.9_combined_gene2go.map")
geneNames <- names(g2go)
go2g <- inverseList(g2go)

files <- scan(inputfiles.list,character())
for (f in files){
	degs <- rownames(read.csv(f,header=T,sep="\t",row.names=1))
	geneList <- factor(as.integer(geneNames %in% degs))
	names(geneList) <- geneNames
	if (length(levels(geneList)) > 1){
	    for (categ in c("BP")){
		    GOdata <- new("topGOdata",ontology=categ,allGenes=geneList,annot = annFUN.gene2GO, gene2GO = g2go)
		    testout <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
		    allRes <- GenTable(GOdata, Fis = testout,topNodes=300)
		    if (dim(allRes[allRes$Fis == "< 1e-30",])[1] != 0){
		    allRes[allRes$Fis == "< 1e-30",]$Fis <- 1e-30
		    }
		    allRes$Fis <- as.numeric(allRes$Fis)
		    pfilt <- allRes[allRes$Fis < 0.01,]
		    names(pfilt) <- c("GOid","GOterm","Annotated","NgenesInTerm","Expected","SignificanceP")
		    pfilt$GenesInTermRatio <- pfilt$NgenesInTerm/pfilt$Annotated
		    outfilen <- gsub(".tsv",paste("_GOenrichment_",categ,".tsv",sep=""),f)
		    allGO <- genesInTerm(GOdata)
		    genesinTerm <- lapply(allGO,function(x) x[x %in% degs] )
		    annotdegs <- c()
		    for (gg in pfilt$GOid){
			    annotg <- unlist(genesinTerm[gg])
			    ll <- length(annotg)
			    annotgenes <- paste(annotg,collapse=",")
			    annotdegs <- append(annotdegs,annotgenes)
		    }
		    pfilt$GOannotdegs <- annotdegs
		    write.table(pfilt,file=outfilen,sep="\t",quote=F,row.names=F)
		    title = ''
		    if (categ == "CC") { title = "Cellular_Component"} else if (categ == "BP") { title = "Biological_Process"} else { title = "Molecular_Function" }
		    m = (min(-log10(pfilt$SignificanceP)) + max(-log10(pfilt$SignificanceP)))/2
		    pfilt$GOs <- paste(pfilt$GOid,pfilt$GOterm)
		    pfilt$GOs <- factor(pfilt$GOs, levels = pfilt$GOs[order(pfilt$GenesInTermRatio)])
		    outpdf = gsub(".tsv",paste("_GOenrichment_",categ,".pdf",sep=""),f)
		    pdf(file=outpdf,width=16,height=dim(pfilt)[1]*0.2)
		    plot <- ggplot(pfilt, aes(x = GenesInTermRatio, y = GOs, color = -log10(SignificanceP), size = log10(NgenesInTerm))) +
		      geom_point() +
		      scale_color_gradient2(low="#ffcc00",mid="#ff6600",high="#e60000",midpoint=m) +
		      xlab("GenesInTermRatio") +
		      ylab("GOterms") +
		      ggtitle(paste(title,"enrichments")) +
		      scale_size_continuous(range = c(1, 10)) +
		      theme_classic()
		    print(plot)
		    dev.off()

	    }
	}
}
