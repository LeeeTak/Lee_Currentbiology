#!/usr/bin/bash
#Author: Tak Lee
#bash script used to read the paired reads and aligning to the genome by STAR. Aligned BAM files will be produced by this script

N=$(wc -l samps_72hrs | cut -d' ' -f1)
FILES=$(cat samps_72hrs)
INPUTS=(${FILES//' '/ })
#echo $FILES
for ((i=0;i<$N;i++));
do
	ulimit -n 10000
	#modify here to change the name of the input files
        f1="${INPUTS[i]}_R1_trimmed.fastq.gz"
	f2="${INPUTS[i]}_R2_trimmed.fastq.gz"
	STAR \
	--outReadsUnmapped Fastx \
	--genomeDir MtrunA17r5.0-ANR/gindex_gtf_149 \
	--readFilesCommand zcat \
	--readFilesIn $f1 $f2 \
	--outSAMtype BAM SortedByCoordinate \
	--runThreadN 40 \
	--outFileNamePrefix bam/${INPUTS[i]}_ 

done

