#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 4) {
  stop("Four arguments must be supplied \n", call. = FALSE)
}

if (!exists("aoi_seurat")) {
    library(monocle)
    library(dplyr)
    user_id <- args[4]
    branch_options <- as.numeric(args[1])
    param <- args[2]
    setwd(file.path("temp", user_id))
    goi <- read.table(args[3])
    print(getwd())
    output_dir <- "Monocle/"
}

library(Seurat)

goi_pseudotime_plot <- function(cds_ct_aoi,
                                param = param,
                                goi,
                                branch_options) {
    p <- plot_genes_branched_pseudotime(cds_ct_aoi[goi, ],
                            branch_point = 1,
                            color_by = param,
                            ncol = 1)
    p_title <- paste0(output_dir,
        "pseudotime_branchpoint",
        branch_options, "_branchgenes_goi",
        ".pdf")
    ggsave(filename = p_title,
        plot = p, width = 6, height = 6, units = "in")
}
beam_header <- "beam_branch_point"
beam_branch_path <- paste0(output_dir,
                param, beam_header,
                branch_options, "_cds.RData")
if (file.exists(beam_branch_path)) {
    load(beam_branch_path)
    cat("Pseudotimeplot of goi \n")
    goi <- goi[, 1]
    goi <- intersect(rownames(cds_ct_aoi), goi)
    goi_pseudotime_plot(goi = goi, cds_ct_aoi = cds_ct_aoi,
            param = param, branch_options = branch_options)
} else {
    cat("None CDS object found, please generate CDS object first \n")
}
