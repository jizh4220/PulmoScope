#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 5) {
  stop("Five arguments must be supplied \n", call. = FALSE)
}


library(monocle)
library(dplyr)
library(ggsci)
user_id <- args[5]
celltype_options <- args[1]
param <- args[2]
branch_options <- as.numeric(args[3])
goi <- suppressWarnings(read.table(args[4], header = FALSE, sep = ";"))
goi <- unlist(goi, use.names = FALSE)
setwd(file.path(user_id))
print(getwd())
output_dir <- "Monocle/"

branching_GEC <- function(cds_ct_aoi, branch_options = 1,
                    goi, param,
                    output_dir = "Monocle/") {
    # align format of each gene
    goi <- toupper(goi)
    goi <- gsub("\\-|\\_", "", goi)
    # check whether all goi in seurat features
    goi <- intersect(goi, rownames(cds_ct_aoi))
    cat("Presenting PseudotimePlot of Genes of Interest\n")
    for (i in 1:length(goi)) {
        pData(cds_ct_aoi)[goi[i]] <- log2(exprs(cds_ct_aoi)[goi[i], ] + 1)
        p <- plot_cell_trajectory(cds_ct_aoi,
                color_by = goi[i],
                cell_name_size = 0.5,
                legend.key.size = unit(1, "cm"),
                legend.text = element_text(size = 3),
                cell_size = 0.5) +
                facet_wrap(~param, nrow = 1)+ scale_color_gsea()
        p_title <- paste0(output_dir, 
            goi[i], "_pseudotimeplot.pdf")
        ggsave(filename = p_title, 
                plot = p,
                width = 8, height = 4, units = "in")
    }

    p <- plot_genes_branched_pseudotime(cds_ct_aoi[top_branch_genes[1:5], ],
                            branch_point = branch_options,
                            color_by = param,
                            ncol = 1)
    p_title <- paste0(output_dir,
        "pseudotime_branch_genes_branchpoint",
        branch_options,
        ".pdf")
    ggsave(filename = p_title,
        plot = p, width = 10, height = 14, units = "in")
}



if (file.exists(paste0(output_dir, celltype_options, "_cds.RData"))) {
    load(paste0(output_dir, celltype_options, "_cds.RData"))
    cat("Calculating branching genes based on BEAM \n")
    branched_gene <- branching_aoi(cds_ct_aoi = cds_ct_aoi,
                    branch_options = as.numeric(branch_options),
                    output_dir = output_dir)
    library(reticulate)
    source_python("../../csv2json.py")
    path <- 'Monocle/*.csv'
    bulk_csv2json(path)
} else {
    cat("None CDS object found, please generate CDS object first \n")
}
