#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 3) {
  stop("Three arguments must be supplied \n", call. = FALSE)
}

print("=====================================")
print("=====         GSEA Module        =====")
print("=====================================")

if (!exists("aoi_seurat")) {
    user_id <- args[3]
    # print(getwd())
    source("../Enrichment_utils.R")
    source("../BPCells_utils.R")
    # suppressWarnings(dir.create(file.path( user_id), recursive = TRUE))
    setwd(file.path(user_id))
    print(getwd())
    celltype_options <- args[1]
    # aoi_groups <- args[2]
    output_dir <- "GSEA/"
    kegg_go <- args[2]
    suppressWarnings(dir.create(output_dir, recursive = TRUE))
}




cat("Load in GroupDEG.\n")
groupwiseDEG_path <- paste0( "Step1/", celltype_options, "_groupwise_DEG_table.csv")
if (file.exists(groupwiseDEG_path)) {
        deg_table <- read.csv(groupwiseDEG_path)
        deg_table <- subset(deg_table, celltype %in% celltype_options)
} else {
        library(qs)
        library(BPCells)
        aoi_seurat <- qread("data/aoi_seurat.qs")
        if (!file.exists(file.path("data", "BPCounts"))) {
                bp_counts <- write_matrix_dir(aoi_seurat[["RNA"]]$counts,
                file.path("data", "BPCounts"), overwrite = TRUE)
        } else {
                bp_counts <- open_matrix_dir(file.path("data", "BPCounts"))
        }
        aoi_seurat <- subset(aoi_seurat, celltype == celltype_options)
        Idents(aoi_seurat) <- aoi_seurat$groups

        deg_table <- GroupBPDEG_Calculation(
                aoi_counts = bp_counts, output_dir = "Step1/",
                celltype_options = celltype_options,
                log2fc = 0.5, pval = 0.05,
                aoi_seurat = aoi_seurat, top_num = 1000)

}


res <- GSEA_allplots(deg_table, celltype_options, output_dir = "GSEA/", kegg_go)

library(reticulate)
source_python("../../csv2json.py")
path <- 'GSEA/*.csv'
bulk_csv2json(path)

cat("Finish GSEA Analysis \n")
