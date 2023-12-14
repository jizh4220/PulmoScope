#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 3) {
  stop("Four arguments must be supplied \n", call. = FALSE)
}

print("=====================================")
print("=====       Monocle_BEAM Module        =====")
print("=====================================")

library(monocle)
library(dplyr)
user_id <- args[3]
param <- args[2]
branch_options <- as.numeric(args[1])
setwd(file.path(user_id))
print(getwd())
output_dir <- "Monocle/"

valid_gene <- function(deg_table) {
        ribo_mt <- "^MT-|^RPL|^RPS|^ATP|^RP.*\\-|LINC"
        deg_table <- deg_table[grep(ribo_mt,
                        rownames(deg_table), invert = TRUE), ]
        return(deg_table)
}


branching_aoi <- function(cds_ct_aoi, branch_options = 1,
                    topnum = 50,
                    output_dir = "Monocle/") {
    cat("Calculate BEAM on Current User's Selected Branch Point\n")
    
    branch_beam <- BEAM(cds = cds_ct_aoi,
                progenitor_method = "sequential_split",
                branch_point = branch_options, cores = 8)

    sig_branch_beam <- subset(branch_beam, qval < 1e-4 & num_cells_expressed >= 0.01 * ncol(cds_ct_aoi))

    write.csv(sig_branch_beam, paste0(output_dir,
            "beam_branch_point", branch_options,
            "_branchingDEG.csv"))
    # branch_beam <- valid_gene(branch_beam)
    top_branch_genes <- top_n(sig_branch_beam, n = topnum,
                    desc(qval)) %>% pull(gene_short_name) %>% as.character()
    p <- plot_genes_branched_heatmap(
                                cds_ct_aoi[
                                    top_branch_genes,
                                ],
                                branch_point = branch_options,
                                show_rownames = TRUE,
                                return_heatmap = TRUE
                            )

    p_title <- paste0(output_dir, 
            "beam_branch_point", branch_options, "_heatmap.png")
    ggsave(p_title, p$ph_res,
            width = 8, height = 6, units = "in")

    cat("Plot Branching Gene Heatmap \n")
    p <- plot_genes_branched_pseudotime(cds_ct_aoi[top_branch_genes[1:5], ],
                            branch_point = branch_options,
                            color_by = param,
                            ncol = 1)
    p_title <- paste0(output_dir,
        "pseudotime_branch_genes_branchpoint",
        branch_options,
        ".png")
    ggsave(filename = p_title,
        plot = p, width = 8, height = 6, units = "in")
}

cds_path <- grep("_cds.RData", list.files(output_dir, full.names = TRUE), value = TRUE)
load(cds_path)

cat("Downsample CDS data for better calculation performance\n")


cat("Calculating branching genes based on BEAM \n")
branched_gene <- branching_aoi(cds_ct_aoi = cds_ct_aoi,
                branch_options = as.numeric(branch_options),
                output_dir = output_dir)
library(reticulate)
source_python("../../csv2json.py")
path <- 'Monocle/*.csv'
bulk_csv2json(path)
