#!/bin/bash

for SRR in $(cat SRR_Acc_List1.txt);\
do fastq-dump $SRR;\
done