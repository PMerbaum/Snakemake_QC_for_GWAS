configfile: "config/config.yaml"
bimfile = config["BIMFILE"]

def matching(bimfile):
    output = [expand("data/processed/{{sample}}_HRCmatched.bim", 
        sample=config["SAMPLE"])]
    
    return output

rule all:
    input:
        matching(bimfile)
        

#check which HRC genotype reference to use - 37 or 38! 
rule hrc_harmonization:
    input:
        "data/raw/{sample}.bim"
    output:
        "data/processed/{sample}_snps_exclude.txt",
        "data/processed/{sample}_snps_flip.txt",
        "data/processed/{sample}_snps_update.txt",
        "data/processed/{sample}_snps_swap.txt"
    threads: 1
    resources: mem_mb=15000, time=600
    script:
        "scripts/match_bim2hrc_sm.R" #37 by default

rule exclude_non_hrc:
    input:
        expand("data/raw/{{sample}}.{ext}", ext=["bim","bed","fam"]),
        expand("data/processed/{{sample}}_{ext}", ext=["snps_exclude.txt", "snps_flip.txt", "snps_update.txt"])
    output:
        expand("data/processed/{{sample}}_hrc_ex.{ext}", ext=["bim","bed","fam"])
    params:
        in_prefix = lambda wildcards, input: input[0][:-4],
        out_prefix = lambda wildcards, output: output[0][:-4],
    shell:
        "plink --bfile {params.in_prefix} --list-duplicate-vars suppress-first --out data/processed/{wildcards.sample}  ; " 
        "plink --bfile {params.in_prefix} --exclude data/processed/{wildcards.sample}.dupvar  --make-bed --out {params.out_prefix}"


rule flip_hrc:
    input:
        expand("data/processed/{{sample}}_hrc_ex.{ext}", ext=["bim","bed","fam"]),
        expand("data/processed/{{sample}}_{ext}", ext=["snps_exclude.txt", "snps_flip.txt", "snps_update.txt"])
    output:
        expand("data/processed/{{sample}}_flipped.{ext}", ext=["bim","bed","fam"])
    params:
        in_prefix = lambda wildcards, input: input[0][:-4],
        out_prefix = lambda wildcards, output: output[0][:-4]
    shell:
        "plink --bfile {params.in_prefix} --exclude {input[3]} --flip {input[4]} --make-bed --out {params.out_prefix}"

rule update_hrc:
    input:
        expand("data/processed/{{sample}}_flipped.{ext}", ext=["bim","bed","fam"]),
        expand("data/processed/{{sample}}_{ext}", ext=["snps_exclude.txt", "snps_flip.txt", "snps_update.txt", 'snps_swap.txt'])
    output:
        expand("data/processed/{{sample}}_HRCmatched.{ext}", ext=["bim","bed","fam"])
    params:
        in_prefix = lambda wildcards, input: input[0][:-4],
        out_prefix = lambda wildcards, output: output[0][:-4]
    shell:
        "plink --bfile {params.in_prefix} --update-name {input[5]} --update-alleles {input[6]} --make-bed --out {params.out_prefix}"

# rule lift:
#     input:
#         expand("{{sample}}_upd.{ext}", ext=["bim"])
#     output:
#         "{sample}_lifted.txt",
#         "{sample}_not_lifted.txt",
#         "{sample}_to_upd.txt"
#     log:
#         'log/{sample}_lift.log'
#     params:
#         GChr_version = "hg19", #change to hg18 if needed 
#         in_prefix = lambda wildcards, input: input[0][:-4],
#         out_prefix = lambda wildcards, output: output[0][:-4]
#     script:
#         "lift_sm.py" 

# rule update_after_lift:
#     input:
#         expand("{{sample}}_upd.{ext}", ext=["bim","bed","fam"]), 
#          "{sample}_lifted.txt",
#         "{sample}_not_lifted.txt"
#     output:
#         expand("{{sample}}_h38.{ext}", ext=["bim","bed","fam"])
#     params:
#         in_prefix = lambda wildcards, input: input[0][:-4],
#         out_prefix = lambda wildcards, output: output[0][:-4]
#     shell:
#         "awk -F '\t' '{{print $2}}' {input[0]} | sort | uniq -d > multi.txt ; "
#         "plink --bfile {params.in_prefix} --exclude multi.txt --make-bed --out {params.in_prefix}_dup_out ; "
#         "plink --bfile {params.in_prefix}_dup_out --exclude {input[4]} --update-chr {input[3]} 1 2 --update-map {input[3]} 4 2 --make-bed --out {params.out_prefix}"


#here your actuall QC starts
rule exclude_missingness:
    input:
        expand("data/raw/{{sample}}.{ext}", ext=["bim","bed","fam"])
    output:
        expand("data/processed/{{sample}}_nomiss.{ext}", ext=["bim","bed","fam"])
    params:
        mind = config["mind"],
        in_prefix = lambda wildcards, input: input[0][:-4],
        out_prefix = lambda wildcards, output: output[0][:-4]
    shell:
        "plink --bfile {params.in_prefix} --missing --out {params.in_prefix} ; "
        "plink --bfile {params.in_prefix} --mind {params.mind} --make-bed --out {params.out_prefix}"

#Would not work if no sex chr are present 
# rule check_sex:
#     input:
#         expand("data/processed/{{sample}}_nomiss.{ext}", ext=["bim","bed","fam"])
#     output:
#         expand("data/processed/{{sample}}_sexchecked.{ext}", ext=["bim","bed","fam"])
        
#     params:
#         in_prefix = lambda wildcards, input: input[0][:-4],
#         out_prefix = lambda wildcards, output: output[0][:-4]
#     shell:
#         "plink --bfile {params.in_prefix} --check-sex ;" 
#         "grep PROBLEM plink.sexcheck > sexcheck_errors.txt ; "
#         "plink --bfile {params.in_prefix} --remove sexcheck_errors.txt --make-bed --out {params.out_prefix}"

rule maf_hwe_geno: 
    input:
        expand("data/processed/{{sample}}_nomiss.{ext}", ext=["bim","bed","fam"])
    output:
        expand("data/processed/{{sample}}_maf.{ext}", ext=["bim","bed","fam"])
    params:
        maf = config["maf"],
        hwe = config['hwe'],
        geno = config['geno'],
        in_prefix = lambda wildcards, input: input[0][:-4],
        out_prefix = lambda wildcards, output: output[0][:-4]
    shell:
        "plink --bfile {params.in_prefix} --maf {params.maf} --hwe {params.hwe} --geno {params.geno} --make-bed --out {params.out_prefix}"

rule relatedness:
    input:
        expand("data/processed/{{sample}}_maf.{ext}", ext=["bim","bed","fam"])
    output:
        expand("data/processed/{{sample}}_related_filter.{ext}", ext=["bim","bed","fam"]),
        "data/processed/{sample}.het",
        "data/processed/{sample}_maf.genome"
    params:
        genome = config['genome'],
        PI_HAT = config['PI_HAT'],
        in_prefix = lambda wildcards, input: input[0][:-4],
        out_prefix = lambda wildcards, output: output[0][:-4]
    shell:
        "plink --bfile {params.in_prefix} --genome --max {params.genome} --out {params.in_prefix} ; "
        "awk '$10 > {params.PI_HAT}' {params.in_prefix}.genome > data/processed/{wildcards.sample}_relatives.txt ; "
        "awk '$10 < {params.genome}' {params.in_prefix}.genome > data/processed/{wildcards.sample}_keep_relatedness.txt ; "
        "plink --bfile {params.in_prefix} --keep data/processed/{wildcards.sample}_keep_relatedness.txt --make-bed --out {params.out_prefix} ; "
        "plink --bfile {params.out_prefix} --het --out {params.out_prefix}"

rule heterozygosity_script:
    input:
        expand("data/processed/{{sample}}_related_filter.{ext}", ext=["bim","bed","fam"])
    output:
        expand("data/processed/{{sample}}_hetfail.{ext}", ext=["txt"])
    script:
        "R-heterozygosity.R"

rule exclude_hetfail:
    input:
        expand("data/processed/{{sample}}_hetfail.{ext}", ext=["txt"]),
        expand("data/processed/{{sample}}_related_filter.{ext}", ext=["bim","bed","fam"])
    output:
        expand("data/processed/{{sample}}_nohetfail.{ext}", ext=["bim","bed","fam"])
    params:
        in_prefix = lambda wildcards, input: input[1][:-4],
        out_prefix = lambda wildcards, output: output[0][:-4]
    shell:
        "plink --bfile {params.in_prefix} --exclude {input[0]} --make-bed --out {params.out_prefix}"

# #visualization and extracting clean data. Be sure R with all packages is available in the environment  
rule visualize_QC:
    input:
        "data/processed/{sample}_h38.bim",
        "data/processed/{sample}.het",
        "data/processed/{sample}_maf.genome"
    output:
        "data/processed/{sample}_fail_het_mis_visualization.png",
        "data/processed/{sample}_h38_individuals_failedQC.txt"
    params:
        imissTh = config["imissTh"],
        highIBDTh = config['highIBDTh'],
        hetTh= config['hetTh'],
        maleTh=config['maleTh'],
        femaleTh=config['femaleTh']
    script: 
        "visualization_plinkQC_sm.R"

rule marker_visualization:
    input:
        "data/processed/{sample}_h38.bim",
        "data/processed/{sample}_h38_individuals_failedQC.txt"
    output:
        "data/processed/{sample}_h38_snp_failedQC.txt"
    params:
        lmissTh=config["lmissTh"],
        hweTh=config["hweTh"],
        macTh=config["macTh"],
        mafTh=config["mafTh"]
    script: 
        "marker_vis.R"

rule ancestry_check:
    input:
        expand("data/processed/{{sample}}_nohetfail.{ext}", ext=["bim"]) #file after all other QC steps 
    output:
        expand("data/processed/{{sample}}_nohetfail.hapmap3_r3_b38_dbsnp150_illumina_fwd.consensus.qc.poly.{ext}", ext=["eigenval","eigenvec"])
    params:
        in_prefix = lambda wildcards, input: input[0][:-4],
        out_prefix = lambda wildcards, output: output[0][:-9]
    shell:
        "bash scripts/hapmap_anc_merge_sm.sh {params.in_prefix} {params.out_prefix}; " #creates a merged sample + hapmap file 
        "sed -i 's/hm3_//' {output[1]}" 


rule anestry_graph:
    input:
        expand("data/processed/{{sample}}_nohetfail.hapmap3_r3_b38_dbsnp150_illumina_fwd.consensus.qc.poly.{ext}", ext=["eigenval","eigenvec"]),
        expand("data/processed/{{sample}}_nohetfail.{ext}", ext=["bim"])
    output:
        "data/processed/{sample}_nohetfail_ancestry.png",
        "data/processed/{sample}_nohetfail_to_remove.txt"
    params:
        europeanTh=config['europeanTh']
    script:
        "create_ancestry_data_sm.R"


rule clean_data:
    input:
        expand("data/processed/{{sample}}_nohetfail.{ext}", ext=["bim","bed","fam"]),
        #"data/processed/{sample}_nohetfail_to_remove.txt"
    output:
         expand("data/result/{{sample}}_clean.{ext}", ext=["bim","bed","fam"])
    params:
        in_prefix = lambda wildcards, input: input[0][:-4],
        out_prefix = lambda wildcards, output: output[0][:-4]
    shell:
        "plink --bfile {params.in_prefix} --make-bed --out {params.out_prefix}" #--remove {input[3]} "
