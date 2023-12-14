#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
print(args)
# test if there is at least one argument: if not, return an error
if (length(args) != 4) {
  stop("Four arguments must be supplied \n", call. = FALSE)
}

print("=====================================")
print("=====        Correlate_GEC Module        =====")
print("=====================================")

if (!exists("aoi_seurat")) {
    library(Seurat)
    library(tidyr)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    user_id <- args[4]
    source("../fraction_utils.R")
    print(getwd())
    setwd(file.path(user_id))
    aoi_seurat <- qread("data/aoi_seurat.qs")
    if (length(unique(aoi_seurat$groups)) == 1 || length(unique(aoi_seurat$groups)) > 2) {
        aoi_seurat$groups <- aoi_seurat$disease
    }
    gene1 <- as.character(args[1])
    gene2 <- as.character(args[2])
    aoi_celltype <- as.character(args[3])
    output_dir <- "Correlated_GEC/"
    suppressWarnings(dir.create(output_dir, recursive = TRUE))
}


Correlated_GEC <- function(aoi_seurat, gene1,
                    gene2, aoi_celltype = "Macrophages",
                    output_dir = "corr_GEC/") {
    dir.create(output_dir, recursive = TRUE)
    meta <- aoi_seurat@meta.data
  # align format of each gene
    gene1 <- toupper(gene1)
    gene1 <- gsub("\\-|\\_", "", gene1)
    gene2 <- toupper(gene2)
    gene2 <- gsub("\\-|\\_", "", gene2)
    # check whether all gene1 in seurat features
    # gene2 <- intersect(gene2, rownames(aoi_seurat))
    CorrFeature_folders <- list.files(output_dir,
            pattern = "CorrFeatureplot.png")
    goi_dat <- FetchData(aoi_seurat, vars = c(gene1, gene2))
    aoi_seurat$gene1_pos <- paste0(gene1, "-")
    aoi_seurat$gene1_pos[goi_dat[, 1] > 0] <- paste0(gene1, "+")
    aoi_seurat$gene2_pos <- paste0(gene2, "-")
    aoi_seurat$gene2_pos[goi_dat[, 2] > 0] <- paste0(gene2, "+")
    aoi_seurat$goi_group <- paste0(aoi_seurat$gene1_pos, aoi_seurat$groups)
    
    cur_gene_exp <- FetchData(aoi_seurat, vars = c(gene1, gene2))
    cur_gene_exp$gene1_pos <- 0
    cur_gene_exp$gene1_pos[which(cur_gene_exp[, 1] > 0)] <- 1
    cur_gene_exp$gene2_pos <- 0
    cur_gene_exp$gene2_pos[which(cur_gene_exp[, 2] > 0)] <- 1
    cur_gene_exp$celltype <- aoi_seurat$celltype

    gene1_pct <- sum(cur_gene_exp$gene1_pos) / nrow(cur_gene_exp)
    gene2_pct <- sum(cur_gene_exp$gene2_pos) / nrow(cur_gene_exp)
    groupname <- paste0(gene1, "_", gene2)
    cat("Current gene expression Gene1: ", gene1_pct, " Gene2: ", gene2_pct, "\n")
    # print(cur_gene_pct)
    if (gene1_pct < 0.05 && gene2_pct < 0.05) {
        png(filename = paste0(output_dir, groupname, "_CorrFeatureplot.png"))
        plot(1, 1, type="n", ann=FALSE, axes=FALSE)
        # Add your text
        text(x = 1, y = 1, paste0("Current gene expression Gene1: ", gene1_pct, " Gene2: ", gene2_pct, " which are too low for visualization\n"))
        # Close the PNG device
        dev.off()
        png(filename = paste0(output_dir, groupname, "_CorrViolinplot.png"))
        plot(1, 1, type="n", ann=FALSE, axes=FALSE)
        # Add your text
        text(x = 1, y = 1, paste0("Current gene expression Gene1: ", gene1_pct, " Gene2: ", gene2_pct, " which are too low for visualization\n"))
        # Close the PNG device
        dev.off()
        png(filename = paste0(output_dir, groupname, "_CorrFraction.png"))
        plot(1, 1, type="n", ann=FALSE, axes=FALSE)
        # Add your text
        text(x = 1, y = 1, paste0("Current gene expression Gene1: ", gene1_pct, " Gene2: ", gene2_pct, " which are too low for visualization\n"))
        # Close the PNG device
        dev.off()
        png(filename = paste0(output_dir, groupname, "_FeatureScatterPlot.png"))
        plot(1, 1, type="n", ann=FALSE, axes=FALSE)
        # Add your text
        text(x = 1, y = 1, paste0("Current gene expression Gene1: ", gene1_pct, " Gene2: ", gene2_pct, " which are too low for visualization\n"))
        # Close the PNG device
        dev.off()
        next
    }


    cat("Plotting CorrFeaturePlot of gene:", groupname, "\n")
    # Corr features
    p <- FeaturePlot(aoi_seurat,
        features = c(gene1, gene2),
        blend = TRUE,
        split.by = "groups",
        pt.size = 0.1) +
                theme(axis.title = element_text(size = 12),
                axis.text = element_text(size = 11, angle = 45),
                legend.key.size = unit(4, "cm"),
                legend.text = element_text(size = 26),
                plot.title = element_text(size = 20, hjust = 0.5))
    p_title <- paste0(output_dir, groupname, "_CorrFeatureplot.png")
    ggsave(filename = p_title,
        plot = p, width = 12, height = 10, units = "in")

    ct_gene1_exp <- table(aoi_seurat$gene1_pos, aoi_seurat$celltype)
    ct_gene2_exp <- table(aoi_seurat$gene2_pos, aoi_seurat$celltype)
    gene1_ct <- ct_gene1_exp[2, ] / (ct_gene1_exp[1, ] + ct_gene1_exp[2, ]) >= 0.1
    gene2_ct <- ct_gene2_exp[2, ] / (ct_gene2_exp[1, ] + ct_gene2_exp[2, ]) >= 0.1
    goi_ct <- colnames(ct_gene1_exp)[gene1_ct * gene2_ct == 1]

    ct_aoi <- subset(aoi_seurat, celltype == aoi_celltype)
        # fraction
    frac_1 <- fraction_param_table(ct_aoi,
            threshold = 0.5,
            param = "gene1_pos",
            goi = gene2)
    head(frac_1)
    # print(head(frac_table))
    ps1 <- corr_fraction_plot(frac_1, goi = gene2, param = "param")
    
    p1 <- VlnPlot(ct_aoi,
                features = gene2,
                group.by = "groups",
                split.by = "gene1_pos",
                # filp = TRUE,
                add.noise = TRUE,
                pt.size = 0.05,
                split.plot = TRUE,
                same.y.lims = TRUE,
                ) +
                labs(
                    # title = paste0(gene2, " in ", gene1, " highest expressed cell type: ", gene1_most_ct)
                    title = paste0("Cell type: ", aoi_celltype)
                    ) +
                theme(legend.text = element_text(colour = 'black',
                        size = 16, angle = 0),
                        axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_text(colour = 'black',
                        size = 20, angle = 0))
    
    # ct_aoi <- subset(aoi_seurat, celltype == gene2_most_ct)
    p2 <- VlnPlot(ct_aoi,
                features = gene1,
                group.by = "groups",
                split.by = "gene2_pos",
                # filp = TRUE,
                add.noise = TRUE,
                pt.size = 0.05,
                split.plot = TRUE,
                same.y.lims = TRUE,
                ) +
                labs(
                    # title = paste0(gene1, " in ", gene2, " highest expressed cell type: ", gene2_most_ct)
                    title = paste0("Cell type: ", aoi_celltype)
                    ) +
                    theme(legend.text = element_text(colour = 'black',
                        size = 16, angle = 0),
                    axis.text.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.title = element_blank(),
                    axis.text.x = element_text(colour = 'black',
                        size = 20, angle = 0))
    p_title <- paste0(output_dir, groupname, "_CorrViolinplot.png")
    ggsave(filename = p_title,
        plot = p1 + p2, width = 12, height = 10, units = "in")

    # Corr_fraction
    frac_2 <- fraction_param_table(ct_aoi,
        threshold = 0.5,
        param = "gene2_pos",
        goi = gene1)
    # print(head(frac_table))
    ps2 <- corr_fraction_plot(frac_2, goi = gene1, param = "param")
    p_title <- paste0(output_dir, groupname,
                    "_CorrFraction.png")
    ggsave(filename = p_title,
            plot = wrap_plots(c(ps1, ps2)), width = 12, height = 10, units = "in")

    cat("FeatureScatter here\n")
    goi_dat <- FetchData(aoi_seurat, vars = c(gene1, gene2))
    cell_idx <- rownames(goi_dat[goi_dat[, 1] > 1 & goi_dat[, 2] > 1, ])
    group_list <- unique(aoi_seurat$groups)
    cat("Groups of current aoi_seurat are: ", group_list, "\n")
    ps <- list()
    for (j in 1:2) {
        param_idx <- names(aoi_seurat$groups[aoi_seurat$groups == group_list[j]])
        ps[[j]] <- FeatureScatter(aoi_seurat,
                    feature1 = gene1,
                    feature2 = gene2,
                    group.by = "groups",
                    cols = colour[j],
                    cells = intersect(cell_idx, param_idx)) +
                theme(axis.title = element_text(size = 12),
                axis.text = element_text(size = 12),
                # legend.key.size = unit(4, "cm"),
                legend.text = element_text(size = 16),
                plot.title = element_text(size = 20, hjust = 0.5))
    }
    p_title <- paste0(output_dir, groupname,
                        "_FeatureScatterPlot.png")
    ggsave(filename = p_title,
            plot = wrap_plots(ps), width = 12, height = 10, units = "in")
}

suppressWarnings(Correlated_GEC(aoi_seurat, gene1,
                gene2, aoi_celltype = aoi_celltype,
                output_dir = output_dir))