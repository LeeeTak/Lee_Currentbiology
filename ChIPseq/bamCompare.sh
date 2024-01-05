#!/bin/bash

#Author: Tak Lee, tlee@mpipz.mpg.de
#running deeptools bamcompare to identify chipseq regions that are high in abundance compared to the control
mkdir -p bamcomp

readarray files1 < $1
threads=$2
for ((i=0;i<${#files1[@]};i++));
do
    IFS=$'\t' read -ra ff <<< "${files1[i]}"
    ip=${ff[0]}
    ctrl=${ff[1]}
    #samtools index -@ $threads $input
    #samtools index -@ $threads $ip
    outname="bamcomp/${ff[2]}_bamCompare.bw"
    bamCompare -p $threads \
	-b1 $ip -b2 $ctrl \
	--binSize 20 \
	--scaleFactorsMethod None \
	--normalizeUsing BPM \
	--centerReads \
	--smoothLength 60 \
	-o $outname

done
