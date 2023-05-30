# Snakemake_QC_for_GWAS-
A Snakemake pipline to perform a pre-QC and QC with visualization for binary genome files (bim, bed, fam).
The pipeline performs harmonization to HRC, lifting to h38 (if needed) and basic quality control steps via plink. 
The output of the code are "clean" binary genome files which could be easily converted into a vcf file.
