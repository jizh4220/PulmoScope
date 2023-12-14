#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 5) {
  stop("Five arguments must be supplied \n", call. = FALSE)
}

print("=====================================")
print("=====         GroupDEG Module        =====")
print("=====================================")


if (!exists("aoi_seurat")) {
    user_id <- args[5]
    setwd(file.path(user_id))
    # sample_info.txt
    log2fc <- args[1]
    pval <- args[2]
    top_num <- args[3]
    celltype_options <- args[4]
    library(Seurat)
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(future)
    library(BPCells)
    plan("multicore", workers = 8)
    print(getwd())
    source("../../BPCells_utils.R")
    aoi_seurat <- qread("data/aoi_seurat.qs")
}


if (!file.exists(file.path("data", "BPCounts"))) {
  bp_counts <- write_matrix_dir(aoi_seurat[["RNA"]]$counts,
            file.path("data", "BPCounts"), overwrite = TRUE)
} else {
  bp_counts <- open_matrix_dir(file.path("data", "BPCounts"))
}
print(table(aoi_seurat$groups))
deg_table <- Groupwise_BPCells_DEG(
            aoi_counts = bp_counts, output_dir = "Step1/",
            celltype_options = celltype_options,
            log2fc = as.numeric(log2fc), pval = as.numeric(pval),
            aoi_seurat = aoi_seurat, top_num = as.numeric(top_num))

library(reticulate)
source_python("../../csv2json.py")
path <- paste0("Step1/", celltype_options, "_groupwise_DEG_table.csv")
bulk_csv2json(path)
cat("Successfully finish Groupwise DEG \n")