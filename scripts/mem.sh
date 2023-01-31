#!/bin/bash


for SRR in $(cat SRR_Acc_List3.txt);\
do bwa mem chr7.fa $SRR.fastq | samtools sort > bams/$SRR.bam;
done
