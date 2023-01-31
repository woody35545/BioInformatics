#!/bin/bash

for SRR in $(cat SRR_Acc_List0.txt);\
do samtools index bams/$SRR.bam;\
bcftools mpileup -Ou -f chr7.fa -r chr7:55000000-55200000 bams/$SRR.bam | bcftools call -mv -Ob -o var/calls_$SRR.bcf;\
bcftools view -I var/calls_$SRR.bcf | bgzip -c > var/calls_$SRR.vcf.gz;\
bcftools index var/calls_$SRR.vcf.gz;\
bcftools norm -f chr7.fa var/calls_$SRR.vcf.gz -Ob -o var/calls_$SRR.norm.bcf;\
bcftools filter --IndelGap 5 var/calls_$SRR.norm.bcf -Ob -o var/calls_$SRR.norm.flt-indels.bcf;\
cat chr7.fa | bcftools consensus --mark-del D var/calls_$SRR.vcf.gz > consensus/consensus_$SRR.fa;\
done

for SRR in $(cat SRR_Acc_List0.txt);\
do cat consensus/consensus_$SRR.fa | tr -d '\n' > consensus/consensus_$SRR.txt;\
done
