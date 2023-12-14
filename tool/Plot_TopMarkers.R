#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 4) {
  stop("Four arguments must be supplied \n", call. = FALSE)
}

print("=====================================")
print("=====         TopMarkers Explorer         =====")
print("=====================================")

cat("Overall Top Markers Time Cost is:", system.time({
    source("../BPCells_utils.R")
    user_id <- args[4]
    setwd(file.path(user_id))
    # sample_info.txt
    log2fc <- args[1]
    pval <- args[2]
    topnum <- args[3]
    library(Seurat)
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(future)
    library(BPCells)
    plan("multicore", workers = 8)
    print(getwd())
    aoi_seurat <- qs::qread("data/aoi_seurat.qs")

  if (!file.exists(file.path("data", "BPCounts"))) {
    bp_counts <- write_matrix_dir(aoi_seurat[["RNA"]]$counts,
              file.path("data", "BPCounts"), overwrite = TRUE)
  } else {
    bp_counts <- open_matrix_dir(file.path("data", "BPCounts"))
  }

  deg_table <- BPCells_Topmarker(output_dir = "Step1/",
              aoi_counts = bp_counts,
              pval = as.numeric(pval),
              log2fc = as.numeric(log2fc),
              aoi_seurat = aoi_seurat,
              top_num = as.numeric(topnum))
  library(reticulate)
  source_python("../../csv2json.py")
  path <- 'Step1/TopMarkers_*.csv'
  bulk_csv2json(path)
}), "\n")