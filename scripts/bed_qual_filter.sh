#!/bin/bash

vcf=$1
bed="./bed/PGRN-seq_ALL_target_Final_cononical.bed"


## filter on min QUAL and on target
vcftools --vcf $vcf --bed $bed --recode --minQ 10 --out ${vcf%%.vcf}.ontarget.minQ10


## filter non-variant sites
vcffixup ${vcf%%.vcf}.ontarget.minQ10.recode.vcf > tmp
vcffilter -f "AC > 0" tmp > tmp.acFilter.vcf
mv tmp.acFilter.vcf ${vcf%%.vcf}.ontarget.minQ10.vcf
rm tmp
rm ${vcf%%.vcf}.ontarget.minQ10.recode.vcf
