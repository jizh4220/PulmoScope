#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 3) {
  stop("Three arguments must be supplied \n", call. = FALSE)
}

print("=====================================")
print("=====         Global_DEG_Display Explorer         =====")
print("=====================================")


if (!exists("aoi_seurat")) {
        source("../fraction_utils.R")
        source("../BPCells_utils.R")
        library(qs)
        user_id <- args[3]
        setwd(file.path(user_id))
        # should be named gene_info.txt
        gene_info <- suppressWarnings(read.csv(args[2],
                header = TRUE))
        # Celltype  Disease  Gene
        all_disease <- "/mnt/root/lungDB_backend/global/disease/"
        cur_disease <- args[1]
        cur_celltype <- unique(gene_info$Celltype)
        #     cat("Current seurat path:", file.path(disease_header,
        #             cur_disease, "data/aoi_seurat.rds"))
        all_diseaseDEG <- list.files(file.path("~/lungDB_backend/tool/display/disease", cur_disease, "Step1/"), "*_groupwise_DEG_table.csv", full.names = TRUE)
        all_diseaseDEG <- grep("sig_", all_diseaseDEG, invert = TRUE, value = TRUE)
        ct_disease_path <- grep(gsub("\\+.*", "", cur_celltype), all_diseaseDEG, value = TRUE)
        aoi_seurat <- qread(file.path(all_disease,
                cur_disease, "data/aoi_seurat.qs"))
        if (length(ct_disease_path) == 0) {
                bp_counts <- open_matrix_dir(file.path(all_disease,
                                cur_disease, "data/BPCounts"))
                ct_disease_deg <- Groupwise_BPCells_DEG(
                        aoi_counts = bp_counts, output_dir = "Step1/",
                        celltype_options = cur_celltype,
                        log2fc = 0.5, pval = 0.05,
                        aoi_seurat = aoi_seurat, top_num = 200)
        } else {
                ct_disease_deg <- read.csv(ct_disease_path)
        }
        output_dir <- "DEG_Display/"
        suppressWarnings(dir.create(output_dir, recursive = TRUE))
}

ViolinComparePlot <- function(aoi_seurat, goi,
                                 reduction = "umap", param = "disease") {
  p <- VlnPlot(aoi_seurat,
                features = goi,
                group.by = "celltype",
                split.by = param,
                # filp = TRUE,
                add.noise = TRUE,
                pt.size = 0.05,
                split.plot = TRUE,
                same.y.lims = TRUE,
                ) + theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title = element_blank(),
                          axis.text = element_text(size = 20),
                          axis.text.x = element_text(vjust = 0.9, angle = 90))
  return(p)
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

DEG_Display <- function(aoi_seurat, gene_info,
                ct_disease_deg,
                param = "disease", output_dir = "GEC/",
                threshold = 0.5,
                floor = 50) {
        goi <- gene_info$Gene
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
        ps <- fraction_plot(frac_table = frac_table, goi = goi,
                                param = "disease",
                                threshold = threshold, floor = floor)
        feature_violin_folders <- list.files(output_dir,
                pattern = "featureDimPlot.png|ViolinComparePlot.png")
        # featureDim and ViolinCompare
        for (i in 1:length(goi)) {
                if (length(grep(goi[i], feature_violin_folders)) == 0) {
                cat("FeatureDim Plot \n")
                p <- featureDim_plot(aoi_seurat = aoi_seurat, goi[i],
                                reduction = "umap", param = param)
                p_title <- paste0(output_dir, goi[i],
                        "_featureDimPlot.png")
                ggsave(filename = p_title,
                        plot = p, width = 16, height = 12, units = "in")
                cat("ViolinCompare Plot \n")
                p <- ViolinComparePlot(aoi_seurat = aoi_seurat, goi[i],
                                        reduction = "umap", param = param)
                p_title <- paste0(output_dir, goi[i],
                        "_ViolinComparePlot.png")
                ggsave(filename = p_title,
                        plot = p, width = 16, height = 12, units = "in")
                }
                cat("Fraction Plot \n")
                p_title <- paste0(output_dir,
                goi[i], "_fractionPlot.png")
                ggsave(filename = p_title,
                        plot = ps[[i]], width = 16, height = 12, units = "in")
        }
        print(getwd())
        cat("Successfully Run DEG Display Module\n")
        return(NULL)
}

suppressWarnings(DEG_Display(aoi_seurat,
                ct_disease_deg = ct_disease_deg,
                gene_info = gene_info,
                output_dir = output_dir))