library(data.table)
library(dplyr)
bimfile = snakemake@input[[1]]
fix_alleles = function(a1_ref, a2_ref, a1_match, a2_match){
  m_alleles = as.matrix(cbind(a1_ref, a2_ref, a1_match, a2_match))
  m_num = matrix(NA, nrow=nrow(m_alleles), ncol=ncol(m_alleles))
  m_num[m_alleles == "A"] = -1
  m_num[m_alleles == "C"] = -2
  m_num[m_alleles == "T"] = 1
  m_num[m_alleles == "G"] = 2
  m_flip = m_num * -1
  palin = ifelse(m_num[,3] == m_flip[,4], 1, 0)
  asis  = ifelse(m_num[,1] == m_num[,3] &
                 m_num[,2] == m_num[,4], 1, 0)
  swap  = ifelse((m_num[,1] == m_num[,4] &
                     m_num[,2] == m_num[,3]) |
                    (m_flip[,1] == m_num[,4] &
                     m_flip[,2] == m_num[,3]), 1, 0)
  flip  = ifelse((m_flip[,1] == m_num[,3] &
                     m_flip[,2] == m_num[,4]) |
                    (m_flip[,1] == m_num[,4] &
                     m_flip[,2] == m_num[,3]), 1, 0)
  action = ifelse(palin == 1 | is.na(m_num[,3]) | is.na(m_num[,4]), "excl",
           ifelse(asis == 1, "asis",
           ifelse(swap == 1, "swap",
           ifelse(swap == 1 & flip == 1, "flipswap",
           ifelse(flip == 1, "flip", "error")))))
   return(action)
}
# read in data
ref = fread("zcat /hpc/hers_en/references/HRC/HRC.r1-1.GRCh37.wgs.mac5.sites.tab.gz", header=T)
colnames(ref)[1:2] = c("chr","bp")
bmatch = fread(file=bimfile, header=F)
colnames(bmatch) = c("chr","id_match","cm","bp","a1_match","a2_match")
bmatch$chr = as.character(bmatch$chr)
# snps to exclude (not in ref)
exclude1 = anti_join(bmatch, ref, by=c("chr","bp")) %>%
            select(id_match)
# snp action:
shared = left_join(bmatch, ref, by=c("chr","bp"), relationship="many-to-many") %>%
          mutate(action = fix_alleles(REF, ALT, a1_match, a2_match)) %>%
          mutate(exclude = ifelse(ID=="." | action %in% c("error", "excl"), 1, 0))
# snps to exclude (palindromic, incompatible alleles, non-dbSNP)
exclude2 = filter(shared, exclude==1) %>%
            select(id_match)
exclude = rbind(exclude1, exclude2)
# snps to flip strand
flip = filter(shared, exclude==0 & action %in% c("flip", "flipswap")) %>%
         select(id_match)
# snp IDs to update
update = filter(shared, exclude==0) %>%
         select(id_match, ID)
# write output
fwrite(exclude, snakemake@output[[1]], col.names=F, row.names=F, quote=F, sep="\t")
fwrite(flip, snakemake@output[[2]], col.names=F, row.names=F, quote=F, sep="\t")
fwrite(update, snakemake@output[[3]], col.names=F, row.names=F, quote=F, sep="\t")