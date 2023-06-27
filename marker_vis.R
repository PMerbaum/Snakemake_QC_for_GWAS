.libPaths('/hpc/hers_en/pmerbaum/tools/miniconda3/envs/snakemake-tutorial/lib/R/library') 

library(ggplot2)
library(magrittr)
library(tidyverse)

fixMixup=FALSE
interactive=FALSE
verbose=FALSE
label_fail=TRUE
highlight_samples = NULL
highlight_type = c("text", "label", "color", "shape")
highlight_text_size = 3
highlight_color = "#c51b8a"
highlight_shape = 17
highlight_legend = FALSE
path2plink='/hpc/local/CentOS7/hers_en/software/plink-1.90/plink'
keep_individuals=NULL
remove_individuals=NULL
exclude_markers=NULL
extract_markers=NULL

#--SNP's overview--
legend_text_size = 5
legend_text_size = 5
legend_title_size = 7
axis_text_size = 5
axis_title_size = 7
title_size = 9
subplot_label_size = 9
lmissTh=snakemake@params[[1]]
hweTh=snakemake@params[[2]]
hweTh <- as.numeric(hweTh)
macTh=snakemake@params[[3]]
mafTh=snakemake@params[[4]]

name <- snakemake@input[1]
prefix <- gsub("\\.bim$", "", name)

#--SNP missingness--

    failids <- paste0(prefix, "_individuals_failedQC.txt")
    failids_df <- read.table(paste0(prefix, "_individuals_failedQC.txt"), header = TRUE, sep = "\t")

    if (!file.exists(failids) | file.size(failids) == 0){
        if (!file.exists(failids)) {
            message("File with individuals that failed perIndividualQC: ",
                    failids, " does not exist. Continue ",
                    "check_SNP_missingness for all samples in ", prefix,
                    ".fam")
        } else {
            message("No individuals failed perIndividualQC (",
                    failids, " is empty). Continue ",
                    "check_SNP_missingness for all samples in ", prefix,
                    ".fam")
        }
        suffix <- ""
        sys::exec_wait(path2plink,
                       args=c("--bfile", prefix,
                              "--missing",
                              "--freq",
                              "--out", paste(prefix, suffix, sep=""),
                              args_filter))

    } else {
        allsamples <- data.table::fread(paste(prefix, ".fam", sep=""),
                                        data.table=FALSE,
                                        stringsAsFactors=FALSE,
                                        header=FALSE)

        failsamples <-  failids_df$ID
        if(all(allsamples[,2] %in% failsamples)) {
            stop("All samples are contained in the", failids, "file ",
                 "from perIndividualQC, no samples remaining for check_maf")
        }
        
        suffix <- ".no_failIDs"
        sys::exec_wait(path2plink,
                       args=c("--bfile", prefix,
                              "--remove", failids,
                              "--missing",
                              "--freq",
                              "--out", paste(prefix, suffix, sep="")))
    
    }

    lmiss <- read.table(paste(prefix, suffix, ".lmiss",sep=""), as.is=TRUE,
                        header=TRUE)

    frq <- read.table(paste(prefix, suffix, ".frq", sep=""), header=TRUE,
                      as.is=TRUE)
    lmiss_frq <- merge(lmiss, frq)
    lmiss_frq$MAF_bin <- ifelse(lmiss_frq$MAF < mafTh, 1, 0)
    p_highMAF <- ggplot(dplyr::filter(lmiss_frq, .data$MAF_bin == 0),
                        aes_string('F_MISS'))
    p_highMAF <- p_highMAF + geom_histogram(binwidth = 0.005,
                                            fill="#66a61e") +
        ylab("Number of SNPs") +
        xlab("Proportion of missing data") +
        ggtitle("SNPs with MAF > ", mafTh) +
        geom_vline(xintercept=lmissTh, lty=2, col="red") +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = legend_title_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))
    p_lowMAF <- ggplot(dplyr::filter(lmiss_frq, .data$MAF_bin == 1),
                       aes_string('F_MISS'))
    p_lowMAF <- p_lowMAF + geom_histogram(binwidth = 0.005,
                                          fill="#e6ab02") +
        ylab("Number of SNPs") +
        xlab("Proportion of missing data") +
        ggtitle("SNPs with MAF < ", mafTh) +
        geom_vline(xintercept=lmissTh, lty=2, col="red") +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = legend_title_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))
    p_histo <- cowplot::plot_grid(p_lowMAF, p_highMAF)
    title <- cowplot::ggdraw() +
        cowplot::draw_label("Marker missingness rate", size=title_size)
    p_lmiss <- cowplot::plot_grid(title, p_histo, ncol = 1,
                                rel_heights = c(0.1, 1))

    ggsave(p_lmiss,file=paste0(prefix,"_lmiss.png"), device='png')
    
    fail_missingness <- lmiss[lmiss$F_MISS > lmissTh,]
    fail_miss <- fail_missingness$SNP
    fail_miss <- data.frame(fail_miss)
    colnames(fail_miss) <- c("SNP")
    if (!is.null(fail_miss) && nrow(fail_miss) > 0) {
        fail_miss$flag <- paste0("high missingness")   
    }  

#---Identification of SNPs showing a significant deviation from Hardy-Weinberg-equilibrium (HWE)---
    sys::exec_wait(path2plink,
                       args=c("--bfile", prefix, "--remove", failids,
                              "--hardy", "--out", paste(prefix, suffix, sep="")))
                            

    hwe <- read.table(paste(prefix, suffix, ".hwe", sep=""), header=TRUE,
                      as.is=TRUE)
    hwe <- hwe[grepl("ALL", hwe$TEST),]
    hwe$P_bin <- ifelse(hwe$P < hweTh, 1, 0)
    hwe$minus_log10P <- -log10(hwe$P)
    p_allP <- ggplot(hwe, aes_string('minus_log10P'))
    p_allP <- p_allP + geom_histogram(binwidth = 0.5,
                                      fill="#66a61e") +
        ylab("Number of SNPs") +
        xlab(expression(-log[10](HWE~exact~test~p-value))) +
        ggtitle(expression(All~p-value[HWE])) +
        geom_vline(xintercept=-log10(hweTh), lty=2, col="red") +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = legend_title_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))
    p_lowP <- ggplot(dplyr::filter(hwe, .data$P_bin == 1),
                     aes_string('minus_log10P'))
    p_lowP <- p_lowP + geom_histogram(binwidth = 0.5,
                                            fill="#e6ab02") +
        ylab("Number of SNPs") +
        xlab(expression(-log[10](HWE~exact~test~p-value))) +
        ggtitle(expression(p-value[HWE]<0.01)) +
        geom_vline(xintercept=-log10(hweTh), lty=2, col="red") +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = legend_title_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))
    p_histo <- cowplot::plot_grid(p_allP, p_lowP)
    title <- cowplot::ggdraw() +
        cowplot::draw_label(expression(Distribution~of~-log[10](p-value[HWE])),
                            size=title_size)
    p_hwe <- cowplot::plot_grid(title, p_histo, ncol = 1,
                                  rel_heights = c(0.1, 1))
    
    ggsave(p_hwe,file=paste0(prefix,"_hwe.png"), device='png')

    fail_hwe <- hwe$SNP[hwe$P < hweTh]
    fail_hwe <- data.frame(fail_hwe)
    colnames(fail_hwe) <- c("SNP")
    if (!is.null(fail_hwe) && nrow(fail_hwe) > 0) {
        fail_hwe$flag <- paste0("deviation from HWE")   
    }  

#---Identification of SNPs with low minor allele frequency---
    # sys::exec_wait(path2plink,
    #                    args=c("--bfile", prefix, "--freq", "--out",
    #                           paste(prefix, suffix, sep="")))

    maf <- read.table(paste(prefix, suffix, ".frq",sep=""),
                       header=TRUE, as.is=TRUE)

    

    if (is.null(mafTh) && is.null(macTh)) {
        stop("Either mafTh or macTh need to be provided")
    }


    p_maf <- ggplot(maf, aes_string('MAF'))
    p_maf <- p_maf + geom_histogram(binwidth = 0.01,
                                            fill="#999999") +
        ylab("Number of SNPs") +
        xlab("Minor allele frequency") +
        ggtitle("Minor allele frequency distribution") +
        geom_vline(xintercept=mafTh, lty=2, col="red") +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = title_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))

    ggsave(p_maf, file=paste0(prefix,"_maf.png"), device='png')
   
    fail_maf <- maf$SNP[maf$MAF < mafTh]
    fail_maf<- data.frame(fail_maf)
    colnames(fail_maf) <- c("SNP")
    if (!is.null(fail_maf) && nrow(fail_maf) > 0) {
        fail_maf$flag <- paste0("low maf")   
    } 

    comb_snp <- NULL

    if ("flag" %in% colnames(fail_miss) && !is.null(fail_miss$flag)) {
        comb_snp <- rbind(try(fail_miss %>% select(SNP, flag), silent = TRUE), comb_snp)
    }

    if ("flag" %in% colnames(fail_hwe)) {
        fail_hwe_select <- fail_hwe %>% select(SNP, flag)
        for (i in 1:nrow(fail_hwe_select)) {
            snp <- fail_hwe_select$SNP[i]
            flag <- fail_hwe_select$SNP[i]
            if (snp %in% comb_snp$SNP) {
                comb_snp$flag[comb_snp$SNP == snp] <- paste(comb_snp$flag[comb_snp$SNP == snp], flag, sep = ",")
            } else {
                comb_snp <- rbind(comb_snp, fail_hwe_select[i, ])
            }
        }
    }

    if ("flag" %in% colnames(fail_maf)) {
        fail_maf_select <- fail_maf %>% select(SNP, flag)
        for (i in 1:nrow(fail_maf_select)) {
            snp <- fail_maf_select$SNP[i]
            flag <- fail_maf_select$SNP[i]
            if (snp %in% comb_snp$SNP) {
                comb_snp$flag[comb_snp$SNP == snp] <- paste(comb_snp$flag[comb_snp$SNP == snp], flag, sep = ",")
            } else {
                comb_snp <- rbind(comb_snp, fail_maf_select[i, ])
            }
        }
    }



    # if ("flag" %in% colnames(fail_miss) && !is.null(fail_miss$flag)) {
    #     comb_snp <- rbind(try(fail_miss %>% select(SNP, flag), silent = TRUE), comb_snp)
    # }

    # if ("flag" %in% colnames(fail_hwe) && !is.null(fail_hwe$flag)) {
    #     comb_snp <- rbind(try(fail_hwe %>% select(SNP, flag), silent = TRUE), comb_snp)
    # }

    # if ("flag" %in% colnames(fail_maf) && !is.null(fail_maf$flag)) {
    #     comb_snp <- rbind(try(fail_maf %>% select(SNP, flag), silent = TRUE), comb_snp)
    # }

    if (!is.null(comb_snp)) {
        saveRDS(comb_snp, file = paste0(prefix, "_snp_failedQC.rds"))
        write.table(comb_snp, file = paste0(prefix, "_snp_failedQC.txt"), col.names = FALSE, row.names = FALSE, quote = FALSE)
    } else {
        print("No SNP's were excluded")
    }
