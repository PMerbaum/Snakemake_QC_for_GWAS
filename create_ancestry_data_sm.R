library(dplyr)
library(magrittr)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)


datadir=getwd()
refname<-'hapmap3_r3_b38_dbsnp150_illumina_fwd.consensus.qc.poly'
name<- snakemake@input[[3]]
name <- gsub("\\.bim$", "", name)
prefix <- snakemake@input[[2]]
prefix <- gsub("\\.eigenvec$", "", prefix)
refSamples <- data.table::fread('/hpc/hers_en/pmerbaum/FTD_GWAS/HapMap.consensus.poly_Pop', header=TRUE, stringsAsFactors=FALSE,
                                 data.table=FALSE)
refColors <- data.table::fread('/hpc/hers_en/pmerbaum/FTD_GWAS/HapMap_PopColors.txt', header=TRUE, stringsAsFactors=FALSE,
                                 data.table=FALSE)
refSamplesFile <- '/hpc/hers_en/pmerbaum/FTD_GWAS/HapMap.consensus.poly_Pop'
refColorsFile <- '/hpc/hers_en/pmerbaum/FTD_GWAS/HapMap_PopColors.txt'
defaultRefSamples = c("HapMap", "1000Genomes")
refPopulation = c("CEU", "TSI")
refSamplesIID = "IID"
refSamplesPop = "Pop"
refColorsColor = "Color"
refColorsPop = "Pop"
studyColor = "#2c7bb6"
legend_labels_per_row = 6
europeanTh = 10

samples <- data.table::fread(paste(name, ".fam", sep=""),
                                 header=FALSE, stringsAsFactors=FALSE,
                                 data.table=FALSE)[,1:2]
    colnames(samples) <- c("FID", "IID")

    if (!file.exists(paste(prefix, ".eigenvec", sep=""))){
        stop("plink --pca output file: ", prefix, ".eigenvec does not exist.")
    }
    # testNumerics(numbers=c(europeanTh, legend_labels_per_row),
    #              positives=c(europeanTh, legend_labels_per_row))
    pca_data <- data.table::fread(paste(prefix, ".eigenvec", sep=""),
                                  stringsAsFactors=FALSE, data.table=FALSE)
    colnames(pca_data) <- c("FID", "IID", paste("PC",1:(ncol(pca_data)-2),
                                                sep=""))
    if (!any(samples$IID %in% pca_data$IID)) {
        stop("There are no ", prefix, ".fam samples in the .eigenvec file")
    }
    if (!all(samples$IID %in% pca_data$IID)) {
        stop("Not all ", prefix, ".fam samples are present in the",
             ".eigenvec file")
    }

    if (is.null(refSamples) && is.null(refSamplesFile)) {
        if (any(!defaultRefSamples %in% c("1000Genomes", "HapMap"))){
            stop("defaultRefSamples should be one of 'HapMap' or '1000Genomes'",
                 " but ", defaultRefSamples," provided")
        }
        defaultRefSamples <- match.arg(defaultRefSamples)
        if (defaultRefSamples == "HapMap") {
            refSamplesFile <- system.file("extdata", "HapMap_ID2Pop.txt",
                                          package="plinkQC")
            if(is.null(refColorsFile) && is.null(refColors)) {
                refColorsFile <-  system.file("extdata", "HapMap_PopColors.txt",
                                              package="plinkQC")
            }
        } else  {
            refSamplesFile <- system.file("extdata", "Genomes1000_ID2Pop.txt",
                                          package="plinkQC")
            if(is.null(refColorsFile) && is.null(refColors)) {
                refColorsFile <-  system.file("extdata",
                                              "Genomes1000_PopColors.txt",
                                              package="plinkQC")
            }
        }

        if (verbose) {
            message("Using ", defaultRefSamples, " as reference samples.")
        }
    }
    if (!is.null(refSamplesFile) && !file.exists(refSamplesFile)) {
        stop("refSamplesFile file", refSamplesFile, "does not exist.")
    }
    if (!is.null(refSamplesFile)) {
        refSamples <- read.table(refSamplesFile, header=TRUE,
                                 stringsAsFactors=FALSE)
    }
    if (!(refSamplesIID  %in% names(refSamples))) {
        stop(paste("Column", refSamplesIID, "not found in refSamples."))
    }
    if (!(refSamplesPop %in% names(refSamples))) {
        stop(paste("Column", refSamplesPop, "not found in refSamples."))
    }
    names(refSamples)[names(refSamples) == refSamplesIID] <- "IID"
    names(refSamples)[names(refSamples) == refSamplesPop] <- "Pop"
    refSamples <- dplyr::select(refSamples, .data$IID, .data$Pop)
    refSamples$IID <- as.character(refSamples$IID)
    refSamples$Pop <- as.character(refSamples$Pop)

    if (!is.null(refColorsFile) && !file.exists(refColorsFile)) {
        stop("refColorsFile file", refColorsFile, "does not exist.")
    }
    if (!is.null(refColorsFile)) {
        refColors <- read.table(refColorsFile, header=TRUE,
                                stringsAsFactors=FALSE)
    }
    if (!is.null(refColors)) {
        if (!(refColorsColor  %in% names(refColors))) {
            stop(paste("Column", refColorsColor, "not found in refColors."))
        }
        if (!(refColorsPop %in% names(refColors))) {
            stop(paste("Column", refColorsPop, "not found in refColors."))
        }
        names(refColors)[names(refColors) == refColorsColor] <- "Color"
        names(refColors)[names(refColors) == refColorsPop] <- "Pop"
        refColors <- dplyr::select(refColors, .data$Pop, .data$Color)
        refColors$Color <- as.character(refColors$Color)
        refColors$Pop <- as.character(refColors$Pop)
    } else {
        refColors <- data.frame(Pop=unique(as.character(refSamples$Pop)),
                                stringsAsFactors=FALSE)
        refColors$Color <- 1:nrow(refColors)
    }
    if (!all(refSamples$Pop %in% refColors$Pop)) {
        missing <- refSamples$Pop[!refSamples$Pop %in% refColors$Pop]
        stop("Not all refSamples populations found in population code of
             refColors; missing population codes: ", paste(missing,
                                                           collapse=","))
    }
    if (!all(refPopulation %in% refColors$Pop)) {
        missing <- refPopulation[!refPopulation %in% refColors$Pop]
        stop("Not all refPopulation populations found in population code of
             refColors; missing population codes: ", paste(missing,
                                                           collapse=","))
    }
    refSamples <- merge(refSamples, refColors, by="Pop", all.X=TRUE)

    ## Combine pca data and population information ####
    data_all <- merge(pca_data, refSamples, by="IID", all.x=TRUE)
    data_all$Pop[data_all$IID %in% samples$IID] <- name
    data_all$Color[data_all$IID %in% samples$IID] <- studyColor
    data_all <- na.omit(data_all)
    if (any(is.na(data_all))) {
        stop("There are samples in the prefixMergedDataset that cannot be found
             in refSamples or ", prefix, ".fam")
    }

    colors <-  dplyr::select(data_all, .data$Pop, .data$Color)
    colors <- colors[!duplicated(colors$Pop),]
    colors <- colors[order(colors$Color),]

    ## Find mean coordinates and distances of reference Europeans ####
    all_european <- dplyr::filter(data_all, .data$Pop %in% refPopulation)
    euro_pc1_mean <- mean(all_european$PC1)
    euro_pc2_mean <- mean(all_european$PC2)

    all_european$euclid_dist <- sqrt((all_european$PC1 - euro_pc1_mean)^2 +
                                         (all_european$PC2 - euro_pc2_mean)^2)

    max_euclid_dist <- max(all_european$euclid_dist)

    ## Find samples' distances to reference Europeans ####
    data_name <- dplyr::filter(data_all, .data$Pop == name)
    data_name$euclid_dist <- sqrt((data_name$PC1 - euro_pc1_mean)^2 +
                                      (data_name$PC2 - euro_pc2_mean)^2)
    non_europeans <- dplyr::filter(data_name, .data$euclid_dist >
                                        (max_euclid_dist * europeanTh))
    fail_ancestry <- dplyr::select(non_europeans, .data$FID, .data$IID)
    legend_rows <- round(nrow(colors)/legend_labels_per_row)

    data_all$shape <- "general"
    shape_guide <- FALSE
    
    colors$Pop <- factor(colors$Pop, levels=unique(colors$Pop))
    data_all$Pop <- factor(data_all$Pop, levels=levels(colors$Pop))

    saveRDS(data_all, paste0(name,'_ancestry.rds'))

ancestry <- readRDS(file = paste0(name,'_ancestry.rds'))
summary <- ancestry %>% filter(Pop == "CEU") %>% summarise(PC1_mean =mean(PC1,na.rm = T), PC1_sd=sd(PC1,na.rm = T), PC2_mean =mean(PC2,na.rm = T), PC2_sd =sd(PC2,na.rm = T))

#outliers:
ancestry %<>% mutate(flag=case_when(
                       ancestry$PC1 > (summary$PC1_mean + (25*summary$PC1_sd)) ~ TRUE,
                       ancestry$PC1 < (summary$PC1_mean - (25*summary$PC1_sd)) ~ TRUE,
                       ancestry$PC2 > (summary$PC2_mean + (25*summary$PC2_sd)) ~ TRUE,
                       ancestry$PC2 < (summary$PC2_mean - (25*summary$PC2_sd)) ~ TRUE,
                       TRUE ~ FALSE))

outliers <-ggplot(ancestry, aes(x=PC1,y=PC2, col=Color)) + geom_point() + theme_bw()


display.brewer.all(colorblindFriendly = TRUE)

rm <- ancestry %>% filter(flag == TRUE) %>% filter(Pop == paste0(name, "samples"))
P1 <- ggplot(ancestry) + geom_point(aes(x=PC1,y=PC2, col=Pop)) +
  theme_bw()+ geom_point(data=rm, aes(x=PC1,y=PC2, size= 2,shape=8, fill=Pop),show.legend = F) + 
  scale_color_manual(values = c(ancestry %>% distinct(Color,Pop) %>% arrange(Pop) %>% pull(Color), "red1")) +scale_shape_identity()+ggtitle("25SD from CEU mean")

ggsave(P1, file=paste0(name,"_ancestry.png"),width = 10, height=10)
write.table(rm %>% select(IID,FID),file=paste0(name,"_to_remove.txt"), col.names =F, quote=F, row.names = F,sep="\t")
