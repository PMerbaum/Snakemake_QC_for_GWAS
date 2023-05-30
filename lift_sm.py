#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 22 10:43:01 2018
@author: wrheene2
"""
import pandas as pd
from pyliftover import LiftOver

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f
    
    lo = LiftOver(snakemake.params.GChr_version, "hg38")
    #lo = LiftOver("hg19", "hg38")

    def lift(pos):
        if str(pos[0]) == '23':
            pos[0]= "X"
        elif str(pos[0]) == '24':
            pos[0]= "Y"
        else:
            pass

        pos_new = lo.convert_coordinate("chr" + str(pos[0]), int(pos[1]))
        if pos_new == None or len(pos_new) == 0:
            bp_new = "NA"
            chrom_new = "NA"
        else:
            chrom_new = pos_new[0][0].replace("chr", "")
            bp_new = pos_new[0][1]
            
        return [chrom_new, bp_new]
    
    data = pd.read_csv(snakemake.input[0], sep="\s+",
                    names=["chr", "id", "cm", "bp", "a1", "a2"])

    old_pos = list(zip(data["chr"], data["bp"]))
    new_pos = [lift(list(pos)) for pos in old_pos]

    data["chr"], data["bp"] = zip(*new_pos)
    data["new_id"] = "chr" + data["chr"].astype(str) + ":" + data["bp"].astype(str) + ":" + data["a1"] + ":" + data["a2"]

    chromosomes = [str(x) for x in range(1,23)]
    chromosomes.append("X")
    chromosomes.append("Y")
    
    out_lifted = data[data['chr'].isin(chromosomes)]
    out_notlifted = data[~data['chr'].isin(chromosomes)]["id"]
    out_updateID = data[["id","new_id"]]

    out_lifted.to_csv(snakemake.output[0], sep="\t", index=False, header=False)
    out_notlifted.to_csv(snakemake.output[1], sep="\t", index=False, header=False)
    out_updateID.to_csv(snakemake.output[2], sep="\t", index=False, header=False)