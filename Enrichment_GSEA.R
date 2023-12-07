#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 4) {
  stop("Four arguments must be supplied \n", call. = FALSE)
}

if (!exists("aoi_seurat")) {
    user_id <- args[4]
    # print(getwd())
    source("../Enrichment_utils.R")
    # suppressWarnings(dir.create(file.path( user_id), recursive = TRUE))
    setwd(file.path(user_id))
    print(getwd())
    celltype_options <- args[1]
    aoi_groups <- args[2]
    output_dir <- "GSEA/"
    kegg_go <- args[3]
    suppressWarnings(dir.create(output_dir, recursive = TRUE))
}


deg_table <- read.csv("Step1/sig_groupwise_DEG_table.csv")
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(DOSE)


df_res <- data.frame(res)
write.table(df_res, file = paste0(output_dir,
                  groupname,
                "GSEA.csv"), row.names = FALSE)

p <- gseaplot2(df_res, geneSetID = 1:5, pvalue_table = T)
# go_enrichment_goi(deg_table = deg_table, output_dir = output_dir,
#                 ct = celltype_options, kegg_go = kegg_go,
#                 aoi_groups = aoi_groups)
p_title <- paste0(output_dir,
                groupname,
                "Enrich_GSEAplot.pdf")
ggsave(filename = p_title,
        plot = p, width = 8, height = 8, units = "in")
pdf(p_title)
p
dev.off()

cat("Finish GSEA Analysis \n")
