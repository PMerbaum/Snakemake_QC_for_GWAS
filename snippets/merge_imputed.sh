#!/bin/bash

#SBATCH --time 4:00:00
#SBATCH --mem 15G

sample = {'file_name'}

for chr in {1..22}; do 
    bcftools merge chr${chr}.dose.vcf.gz -O z -o ${sample}.vcf.gz ;
done