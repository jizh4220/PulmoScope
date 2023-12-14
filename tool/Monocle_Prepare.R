#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 3) {
  stop("Two arguments must be supplied \n", call. = FALSE)
}

print("=====================================")
print("=====       Monocle_Preparation Module        =====")
print("=====================================")

library(dplyr)
user_id <- args[3]
setwd(file.path(user_id))
celltype_options <- args[1]
param <- args[2]
source("../../Monocle_utils.R")
library(monocle)
dir.create("Monocle/")
output_dir <- "Monocle/"
if (file.exists(paste0(output_dir, celltype_options, "_cds.RData"))) {
        load(paste0(output_dir, celltype_options, "_cds.RData"))
        cat("Load Current CDS object and Plot Cell Trajectory\n")
        cell_trajectory_paramPlot(cds_ct_aoi = cds_ct_aoi,
                                param = param,
                                output_dir = output_dir)
} else {

        file.remove(list.files(output_dir,
            pattern = "_cds.RData", full.names = TRUE))
        library(qs)
        aoi_seurat <- qread("data/aoi_seurat.qs")
        # celltype_options <- "Macrophages"
        aoi_seurat <- subset(aoi_seurat, celltype %in% celltype_options)
        df <- table(aoi_seurat$disease, aoi_seurat$celltype)
        if (ncol(aoi_seurat) >= 2000) {
          print("Downsample to only 2000 cells")
          Seurat::Idents(aoi_seurat) <- aoi_seurat$groups
          aoi_seurat <- subset(aoi_seurat, downsample = 1000)
        }

        print(table(aoi_seurat$celltype))
        cds_ct_aoi <- monocle_preparation(aoi_seurat = aoi_seurat, celltype_options = celltype_options, output_dir = output_dir)
        cat("Save Current CDS object and Plot Cell Trajectory\n")
        cell_trajectory_paramPlot(cds_ct_aoi = cds_ct_aoi,
                                param = param,
                                output_dir = output_dir)
}
