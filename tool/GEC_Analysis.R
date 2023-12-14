#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 3) {
  stop("Three arguments must be supplied \n", call. = FALSE)
}

print("=====================================")
print("=====         GEC Module        =====")
print("=====================================")

if (!exists("aoi_seurat")) {
    source("../fraction_utils.R")
    user_id <- args[3]
    setwd(file.path(user_id))
    aoi_seurat <- qread("data/aoi_seurat.qs")
    param <- args[1]
    goi <- suppressWarnings(read.table(args[2], header = FALSE, sep = ";"))
    goi <- unlist(goi, use.names = FALSE)
    output_dir <- "GEC/"
    suppressWarnings(dir.create(output_dir, recursive = TRUE))
}

#' Gene Expression Comparison Analysis
#' An integrated function to convert .csv/.h5ad/.rds objects to BPCell object
#' and finally to Seurat V5 object or Scanpy object.
#' 
#' @aoi_seurat an integrated seurat object with UMAP
#' @goi gene of interest, can be one or more genes
#' @param metadata column to compare gene expression
#' @output_dir directory to store output results, should always be "GEC/"
#' @threshold gene expression threshold in data matrix 
#' @floor number of cell threshold
#' 

GEC_Analysis <- function(aoi_seurat, goi,
                    param = "DISEASE", output_dir = "GEC/",
                    threshold = 0.5,
                    floor = 50) {
    aoi_seurat$orig.ident <- aoi_seurat$accession
    suppressWarnings(dir.create(output_dir, recursive = TRUE))
  # align format of each gene
    goi <- toupper(goi)
    goi <- gsub("\\-|\\_", "", goi)
    # check whether all goi in seurat features
    goi <- intersect(goi, rownames(aoi_seurat))
    # fraction table and fraction plot
    all_folders <- list.files(output_dir)
    frac_table_path <- paste0(output_dir, paste(unlist(goi), collapse = "_"),
            "_fractionTable.csv")
    if (file.exists(frac_table_path)) {
        cat("Previous Fraction Table Calculation has already been stored \n")
        frac_table <- read.csv(frac_table_path)
    } else {
        frac_table <- fraction_param_table(aoi_seurat,
                                    threshold = threshold,
                                    param = "celltype", #will not change
                                    goi = goi)
        write.csv(frac_table, frac_table_path)
    }
    ps <- fraction_plot(frac_table = frac_table, goi = goi, param = "groups",
                        threshold = threshold, floor = floor)
    feature_violin_folders <- list.files(output_dir,
            pattern = "featureDimPlot.png|ViolinComparePlot.png")
    # featureDim and ViolinCompare
    for (i in 1:length(goi)) {
        if (length(grep(goi[i], feature_violin_folders)) < 2) {
            cat("FeatureDim Plot \n")
            p <- featureDim_plot(aoi_seurat = aoi_seurat, goi[i],
                            reduction = "umap", param = param)
            p_title <- paste0(output_dir, goi[i],
                "_featureDimPlot.png")
            ggsave(filename = p_title,
                    plot = p, width = 12, height = 8, units = "in")
            cat("ViolinCompare Plot \n")
            p <- ViolinComparePlot(aoi_seurat = aoi_seurat, goi[i],
                                reduction = "umap", param = param)
            p_title <- paste0(output_dir, goi[i],
                "_ViolinComparePlot.png")
            ggsave(filename = p_title,
                    plot = p, width = 12, height = 8, units = "in")
        }
        # cat("Fraction Plot \n")
        p_title <- paste0(output_dir,
            goi[i], "_fractionPlot.png")
        ggsave(filename = p_title,
                plot = ps[[i]], width = 12, height = 8, units = "in")
    }
    cat("Successfully Run Gene Expression Comparison Analysis \n")
    # return(NULL)
}

suppressWarnings(GEC_Analysis(aoi_seurat,
                    param = param,
                    goi = goi,
                    output_dir = output_dir))