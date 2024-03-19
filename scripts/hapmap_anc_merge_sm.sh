qcdir=$(pwd)
name=$1 
refdir='/hpc/hers_en/references/hapmap/'
refname='hapmap3_r3_b38_dbsnp150_illumina_fwd.consensus.qc.poly'
highld='/hpc/hers_en/pmerbaum/tools/R/x86_64-pc-linux-gnu-library/3.6/plinkQC/extdata/high-LD-regions-hg38-GRCh38.txt'

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $qcdir/$name.bim  > \
    $qcdir/$name.ac_gt_snps

awk 'BEGIN {OFS="\t"}  ($5$6 == "GC" || $5$6 == "CG" \
                        || $5$6 == "AT" || $5$6 == "TA")  {print $2}' \
    $refdir/$refname.bim  > \
    $qcdir/$refname.ac_gt_snps
   
plink --bfile  $refdir/$refname \
      --exclude $qcdir/$refname.ac_gt_snps \
      --make-bed \
      --out $qcdir/$refname.no_ac_gt_snps

plink --bfile  $qcdir/$name \
      --exclude $qcdir/$name.ac_gt_snps \
      --make-bed \
      --out $qcdir/$name.no_ac_gt_snps

plink --bfile  $qcdir/$name.no_ac_gt_snps \
      --exclude range  $highld \
      --indep-pairwise 50 5 0.2 \
      --out $qcdir/$name.no_ac_gt_snps

plink --bfile  $qcdir/$name.no_ac_gt_snps \
      --extract $qcdir/$name.no_ac_gt_snps.prune.in \
      --make-bed \
      --out $qcdir/$name.pruned

plink --bfile  $refdir/$refname \
      --extract $qcdir/$name.no_ac_gt_snps.prune.in \
      --make-bed \
      --out $qcdir/$refname.pruned

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$1; next} \
    ($2 in a && a[$2] != $1)  {print a[$2],$2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim | \
    sed -n '/^[XY]/!p' > $qcdir/$refname.toUpdateChr

plink --bfile $qcdir/$refname.pruned \
      --update-chr $qcdir/$refname.toUpdateChr 1 2 \
      --make-bed \
      --out $qcdir/$refname.updateChr

awk 'BEGIN {OFS="\t"} FNR==NR {a[$2]=$4; next} \
    ($2 in a && a[$2] != $4)  {print a[$2],$2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
    $qcdir/${refname}.toUpdatePos

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5)  {print $2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.pruned.bim > \
    $qcdir/$refname.toFlip

plink --bfile $qcdir/$refname.updateChr \
      --update-map $qcdir/$refname.toUpdatePos 1 2 \
      --flip $qcdir/$refname.toFlip \
      --make-bed \
      --out $qcdir/$refname.flipped

awk 'BEGIN {OFS="\t"} FNR==NR {a[$1$2$4]=$5$6; next} \
    ($1$2$4 in a && a[$1$2$4] != $5$6 && a[$1$2$4] != $6$5) {print $2}' \
    $qcdir/$name.pruned.bim $qcdir/$refname.flipped.bim > \
    $qcdir/$refname.mismatch

plink --bfile $qcdir/$refname.flipped \
      --exclude $qcdir/$refname.mismatch \
      --make-bed \
      --out $qcdir/$refname.clean

plink --bfile $qcdir/$name.pruned  \
      --bmerge $qcdir/$refname.clean.bed $qcdir/$refname.clean.bim \
         $qcdir/$refname.clean.fam --allow-no-sex \
      --make-bed \
      --out $qcdir/$name.merge.$refname

plink --bfile $qcdir/$name.merge.$refname \
      --pca --geno 0.01 --allow-no-sex\
      --out $qcdir/$name.$refname