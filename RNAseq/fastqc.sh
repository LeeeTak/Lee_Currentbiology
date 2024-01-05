#!/bin/bash
#Author: Tak Lee
#bash script used to run fastqc to check the quality of fastq files
threadn=$1

fastqc *.gz --threads $threadn --outdir fastqc/
