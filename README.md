# Snakemake_QC_for_GWAS

This Snakemake pipeline performs all necessary pre-QC and QC steps with visualization for binary genome files (bim, bed, fam) that should be stored in data/raw.
The pipeline harmonizes SNP IDs to HRC, lifts to the required genome build, and runs basic quality control steps via PLINK All parameters might be tailored with a config file.
The code output is "clean" binary genome files which could be easily converted into a .vcf file.
