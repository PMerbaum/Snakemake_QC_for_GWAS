#create seperate png files for individual plinkQC 
#output 3 plots, het/imiss, gender mismatch , relatedness 
.libPaths('/hpc/hers_en/pmerbaum/tools/miniconda3/envs/snakemake-tutorial/lib/R/library') 

library(ggplot2)
library(magrittr)
library(tidyverse)

#HET and imiss ------

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
path2plink=NULL
keep_individuals=NULL
remove_individuals=NULL
exclude_markers=NULL
extract_markers=NULL
legend_text_size = 5
legend_title_size = 7
axis_text_size = 5
axis_title_size = 7
title_size = 9

imissTh=0.05 #(NOT DEFAULT)
hetTh=3 #(NOT DEFAULT)
name <- snakemake@input[1]
prefix <- gsub("\\.bim$", "", name) #{sample}_h38
sample <- gsub("_h38\\.bim$", "", name)

names_imiss <- c("FID", "IID", "MISS_PHENO", "N_MISS", "N_GENO", "F_MISS")
imiss <- read.table(paste(prefix, ".imiss", sep=""), header=TRUE,
                        as.is=TRUE)

if (!all(names_imiss == names(imiss))) {
        stop("Header of ", prefix, ".imiss is not correct. Was your
             file generated with plink --imiss?")
    }
fail_imiss <- imiss[imiss$F_MISS > imissTh,]

names_het <- c("FID", "IID", "O.HOM.", "E.HOM.", "N.NM.", "F")
het <- read.table(paste0(sample,'.het'), header=TRUE, as.is=TRUE)
if (!all(names_het == names(het))) {
    stop("Header of ",snakemake@input[2], " is not correct. Was your
            file generated with plink --het?")
    }
fail_het <- het[het$F < (mean(het$F)  - hetTh*sd(het$F)) |
                        het$F > (mean(het$F) + hetTh*sd(het$F)),]


nr_samples <- nrow(imiss)
imiss$logF_MISS <- log10(imiss$F_MISS)
het_imiss <- merge(imiss, het, by="IID")
fail_het_imiss <- het_imiss[which(het_imiss$IID %in%
                                          union(fail_het$IID, fail_imiss$IID)),]
if (nrow(fail_het_imiss) == 0) {
        fail_het_imiss <- NULL
    }
##
    het_imiss$type <- "pass"
    het_imiss$type[het_imiss$IID %in% fail_het$IID] <- "fail het"
    het_imiss$type[het_imiss$IID %in% fail_imiss$IID] <- "fail miss"
    het_imiss$type[het_imiss$IID %in%
                       intersect(fail_het$IID, fail_imiss$IID)] <- "fail het + miss"

    minus_sd <- mean(het_imiss$F) - 1:5*(sd(het_imiss$F))
    plus_sd <- mean(het_imiss$F) + 1:5*(sd(het_imiss$F))

    colors <- c("#666666", "#1b9e77", "#d95f02", "#7570b3")
    names(colors) <- c("pass", "fail het", "fail miss", "fail het + miss" )

    het_imiss$shape <- "general"
    shape_guide <- FALSE

    if(!is.null(highlight_samples)) {
        if (!all(highlight_samples %in% het_imiss$IID)) {
            stop("Not all samples to be highlighted are present in the",
                 "prefixMergedDataset")
        }
        highlight_type <- match.arg(highlight_type, several.ok = TRUE)
        if (all(c("text", "label") %in% highlight_type)) {
            stop("Only one of text or label highlighting possible; either ",
                 "can be combined with shape and color highlighting")
        }

        if ("shape" %in% highlight_type) {
            het_imiss$shape[het_imiss$IID %in% highlight_samples] <- "highlight"
            shape_guide <- highlight_legend
        }
        if ("color"  %in% highlight_type && highlight_legend) {
            het_imiss$type[het_imiss$IID %in% highlight_samples] <- "highlight"
            colors <- c(colors,  highlight_color)
            names(colors)[length(colors)] <- "highlight"
        }
    }
##
    het_imiss$type <- factor(het_imiss$type, levels=names(colors))
    het_imiss$shape <- as.factor(het_imiss$shape)

    p_het_imiss <- ggplot()
    p_het_imiss <- p_het_imiss + geom_point(data=het_imiss,
                                            aes_string(x='logF_MISS', y='F',
                                                       color='type',
                                                       shape="shape")) +
        scale_shape_manual(values=c(16, highlight_shape), guide="none") +
        scale_color_manual(values=colors) +
        labs(x = "Proportion of missing SNPs",
             y = "heterozygosity rate (and sd)",
             color = "Marker",
             title = "heterozygosity by Missingness across samples") +
        geom_hline(yintercept=c(minus_sd[1:3], plus_sd[1:3]), lty=2,
                   col="azure4") +
        scale_y_continuous(labels=c("-5", "-4", "-3" ,"+3", "+4", "+5"),
                           breaks=c(minus_sd[3:5], plus_sd[3:5])) +
        scale_x_continuous(labels=c(0.0001, 0.001, 0.01, 0.03, 0.05, 0.01, 1),
                           breaks=c(-4,-3,-2, log10(0.03), log10(0.05),-1,0)) +
        geom_hline(yintercept=mean(het_imiss$F) - (hetTh*sd(het_imiss$F)),
                   col="#e7298a", lty=2) +
        geom_hline(yintercept=mean(het_imiss$F) + (hetTh*sd(het_imiss$F)),
                   col="#e7298a", lty=2) +
        geom_vline(xintercept=log10(imissTh), col="#e7298a", lty=2)

    if (!is.null(fail_het_imiss) && label_fail) {
        p_het_imiss <-
            p_het_imiss + ggrepel::geom_label_repel(
                data=data.frame(x=fail_het_imiss$logF_MISS,
                                y=fail_het_imiss$F,
                                label=fail_het_imiss$IID),
                aes_string(x='x', y='y', label='label'),
                size=highlight_text_size)
    }
    ##
    highlight_data <- dplyr::filter(het_imiss, .data$IID %in% highlight_samples)
    if (!is.null(highlight_samples)) {
        if ("text"  %in% highlight_type) {
            p_het_imiss <- p_het_imiss +
                ggrepel::geom_text_repel(data=highlight_data,
                                         aes_string(x='logF_MISS', y='F',
                                                    label="IID"),
                                         size=highlight_text_size)
        }
        if ("label"  %in% highlight_type) {
            p_het_imiss <- p_het_imiss +
                ggrepel::geom_label_repel(data=highlight_data,
                                          aes_string(x='logF_MISS', y='F',
                                                     label="IID"),
                                          size=highlight_text_size)
        }

        if ("color"  %in% highlight_type && !highlight_legend) {
            p_het_imiss <- p_het_imiss +
                geom_point(data=highlight_data,
                           aes_string(x='logF_MISS', y='F', shape='shape'),
                           color=highlight_color,
                           show.legend=highlight_legend)
        }
        if ("shape"  %in% highlight_type && highlight_legend) {
            p_het_imiss <- p_het_imiss +
                guides(shape = "legend")
        }
    }
    ##
    p_het_imiss <- p_het_imiss +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              axis.text = element_text(size = axis_text_size),
              title = element_text(size = title_size),
              axis.title = element_text(size = axis_title_size))

ggsave(p_het_imiss,  file=paste0(sample,"_fail_het_mis_visualization.png"), device='png')


# #sex ------
maleTh=0.8
femaleTh=0.2
externalSex=NULL

fixMixup=FALSE
verbose=FALSE
label_fail=TRUE
highlight_samples = NULL
highlight_type = c("text", "label", "color", "shape")
highlight_text_size = 3
highlight_color = "#c51b8a"
highlight_shape = 17
highlight_legend = FALSE
path2plink=NULL
keep_individuals=NULL
remove_individuals=NULL
exclude_markers=NULL
extract_markers=NULL
legend_text_size = 5
legend_title_size = 7
axis_text_size = 5
axis_title_size = 7
title_size = 9

names_sexcheck <- c("FID", "IID", "PEDSEX", "SNPSEX", "STATUS", "F")

#  sexcheck <- read.table(paste(prefix, ".sexcheck",sep=""),
   sexcheck <- read.table('plink.sexcheck',
                           header=TRUE, stringsAsFactors=FALSE)
sexcheck %<>% filter(SNPSEX != 0)
fail_sex <- sexcheck %>% filter(STATUS != "PROBLEM")

sexcheck$LABELSEX <- "Unassigned"
sexcheck$LABELSEX[sexcheck$PEDSEX == 1] <- "Male"
sexcheck$LABELSEX[sexcheck$PEDSEX == 2] <- "Female"

sexcheck$shape <- "general"
shape_guide <- FALSE

 colors <- c("#999999", "#377eb8", "#e41a1c")
    names(colors) <- c("Unassigned", "Male", "Female")

sexcheck$LABELSEX <- factor(sexcheck$LABELSEX, levels=names(colors))
sexcheck$PEDSEX <- as.factor(sexcheck$PEDSEX)
sexcheck$shape <- as.factor(sexcheck$shape)

p_sexcheck <- ggplot()
    p_sexcheck <- p_sexcheck + geom_point(data=sexcheck,
                                          aes_string(x='PEDSEX', y='F',
                                                     color='LABELSEX',
                                                     shape='shape')) +
        scale_shape_manual(values=c(16, highlight_shape), guide="none") +
        scale_color_manual(values=colors, name="Sex") +
        labs(title="Check assigned sex versus SNP sex",
             x="Reported Sex (PEDSEX)",
             y="ChrX heterozygosity") +
        geom_segment(data=data.frame(x=0.8, xend=1.2, y=maleTh,
                                     yend=maleTh),
                     aes_string(x='x', xend='xend', y='y', yend='yend'), lty=2,
                     color="#e7298a") +
        geom_segment(data=data.frame(x=1.8, xend=2.2, y=femaleTh,
                                     yend=femaleTh), lty=2,
                     aes_string(x='x', xend='xend', y='y', yend='yend'),
                     color="#e7298a")
    if (!is.null(fail_sex) && label_fail) {
        p_sexcheck <- p_sexcheck +
            ggrepel::geom_label_repel(
                data=dplyr::filter(sexcheck, .data$IID %in% fail_sex$IID),
                aes_string(x='PEDSEX',
                           y='F',
                           label='IID'),
                size=highlight_text_size)
    }

    if (!is.null(highlight_samples)) {
        highlight_data <- dplyr::filter(sexcheck, .data$IID %in% highlight_samples)
        if ("text"  %in% highlight_type) {
            p_sexcheck <- p_sexcheck +
                ggrepel::geom_text_repel(data=highlight_data,
                                         aes_string(x='PEDSEX', y='F',
                                                    label="IID"),
                                         size=highlight_text_size)
        }
        if ("label"  %in% highlight_type) {
            p_sexcheck <- p_sexcheck +
                ggrepel::geom_label_repel(data=highlight_data,
                                          aes_string(x='PEDSEX', y='F',
                                                     label='IID'),
                                          size=highlight_text_size)
        }

        if ("color"  %in% highlight_type && !highlight_legend) {
            p_sexcheck <- p_sexcheck +
                geom_point(data=highlight_data,
                           aes_string(x='PEDSEX', y='F', shape='shape'),
                           color=highlight_color)

        }
        if ("shape"  %in% highlight_type && highlight_legend) {
            p_sexcheck <- p_sexcheck +
                guides(shape = "legend") +
                labs(shape = "Individual")
        }
    }
    p_sexcheck <- p_sexcheck +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = title_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))

ggsave(p_sexcheck,  file=paste0(prefix,"_fail_sex_visualization.png"), device='png')


#---- PI-HAT:
legend_text_size = 5
legend_title_size = 7
axis_text_size = 5
axis_title_size = 7
title_size = 9

highIBDTh = 0.9

genome <- read.table(paste0(sample, '_maf.genome'), header=TRUE,
                         as.is=TRUE, stringsAsFactors=FALSE)


    genome$PI_HAT_bin <- ifelse(genome$PI_HAT > 0.05, 0, 1)
    p_allPI_HAT <- ggplot(genome, aes_string('PI_HAT'))
    p_allPI_HAT <- p_allPI_HAT + geom_histogram(binwidth = 0.005,
                                                fill="#66a61e") +
        ylab("Number of pairs") +
        xlab("Estimated pairwise IBD (PI_HAT)") +
        ggtitle("IBD for all sample pairs") +
        geom_vline(xintercept=highIBDTh, lty=2, col="#e7298a") +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = legend_text_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))
    p_highPI_HAT <- ggplot(dplyr::filter(genome, .data$PI_HAT_bin == 0),
                           aes_string('PI_HAT'))
    p_highPI_HAT <- p_highPI_HAT + geom_histogram(binwidth = 0.005,
                                                  fill="#e6ab02") +
        ylab("Number of pairs") +
        xlab("Estimated pairwise IBD (PI_HAT)") +
        ggtitle("IBD for sample pairs with PI_HAT >0.1") +
        geom_vline(xintercept=highIBDTh, lty=2, col="#e7298a") +
        theme_bw() +
        theme(legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = legend_title_size),
              title = element_text(size = legend_text_size),
              axis.text = element_text(size = axis_text_size),
              axis.title = element_text(size = axis_title_size))
 

 ggsave(p_highPI_HAT,  file=paste0(prefix,"relatedness.png"), device='png')
 