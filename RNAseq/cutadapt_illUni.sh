#!/bin/bash
#Author: Tak Lee
#script used to trim illumina universal adapters from raw fastq files

readarray files < samps_72hrs
cores=$1
for ((i=0;i<${#files[@]};i++));
do
        name=${files[i]//[$'\t\r\n']}
        infile1="${name}_R1_all.fastq.gz"
        infile2="${name}_R2_all.fastq.gz"
        #forout=(${infile//'_'/ })
        outfile1="${name}_R1_trimmed.fastq.gz"
        outfile2="${name}_R2_trimmed.fastq.gz"
        echo $infile1
        cutadapt --cores=$cores -a "AGATCGGAAGAG" -A "AGATCGGAAGAG" -m 10 -o $outfile1 -p $outfile2 $infile1 $infile2
done
