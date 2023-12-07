library(monocle)
library(ggplot2)
library(ggsci)
library(patchwork)
library(tidyverse)
output_dir <- "Monocle/"
monocle_param <- c("groups", "DISEASE", "TISSUE", "aging", "gender",
        "gse_alias", "macro_celltype", "Pseudotime", "State")

colour <- c("#DC143C", "#0000FF", "#20B2AA",
        "#FFA500", "#9370DB", "#98FB98",
        "#F08080", "#1E90FF", "#7CFC00",
        "#FFFF00", "#808000", "#FF00FF",
        "#FA8072", "#7B68EE", "#9400D3")
# suppressWarnings(dir.create(output_dir, recursive = TRUE))

deg_filter <- function(deg_table) {
        ribo_mt <- "^MT-|^RPL|^RPS|^ATP|^RP.*\\-|^A.*\\.[0-9]|^RNU|^RN7S|^RNVU|^CT.*\\-|^LINC|ENSG|^C.*orf|^LOC|^RF00|^C.*ORF|^XXbac|^Metazoa|[a-z]|^Y-RNA|^CH17-|.*\\..*|^MIR[0-9]"
        deg_table <- deg_table[grep(ribo_mt, deg_table$gene, invert = TRUE), ]
        deg_table <- deg_table[deg_table$p_val_adj < 0.05 &
                        deg_table$avg_log2FC > 0, ]
        deg_table <- deg_table[grep("\\.|LINC",
                        deg_table$gene, invert = TRUE), ]
        return(deg_table)
}

order_monocle_by_markers <- function(aoi_seurat = NULL,
                            cds_ct_aoi,
                            celltype_options,
                            output_dir) {
        library(dplyr)
        cat("Order Monocle with GroupDEG.\n")
        groupwiseDEG_path <- paste0(gsub("Monocle/", "Step1/", output_dir), celltype_options, "_groupwise_DEG_table.csv")
        if (file.exists(groupwiseDEG_path)) {
                deg_table <- read.csv(groupwiseDEG_path)
                deg_table <- subset(deg_table, celltype %in% celltype_options)
        } else {
                deg_table <- Groupwise_DEG(aoi_seurat = aoi_seurat,
                        output_dir = gsub("Monocle/", "Step1/", output_dir),
                        celltype_options = celltype_options,
                        log2fc = 0.5, pval = 0.05,
                        top_num = 500)
        }
        deg_table <- deg_filter(deg_table)
        ct_markers <- deg_table %>%
                group_by(cluster) %>% slice_min(p_val_adj, n = 500)
        # print(ct_markers$gene)
        cds_ct_aoi <- setOrderingFilter(cds_ct_aoi, ct_markers$gene)
        print(cds_ct_aoi[[1]])
        if (nrow(ct_markers) < 50 || length(cds_ct_aoi[[1]]) < 1500) {
            cat("Not enough DEG for auto param selection \n")
            cds_ct_aoi <- reduceDimension(cds_ct_aoi,
                                norm_method = "none",
                                reduction_method = "DDRTree",
                                max_components = 3,
                                scaling = TRUE,
                                verbose = TRUE,
                                auto_param_selection = FALSE,
                                pseudo_exp = 0)
        } else {
            cds_ct_aoi <- reduceDimension(cds_ct_aoi,
                                norm_method = "none",
                                reduction_method = "DDRTree",
                                max_components = 3,
                                scaling = TRUE,
                                verbose = TRUE,
                                # auto_param_selection = FALSE,
                                pseudo_exp = 0)
        }
        
        #psudotime value
        cds_ct_aoi <- orderCells(cds_ct_aoi)
        return(cds_ct_aoi)
}


monocle_preparation <- function(aoi_seurat, celltype_options, output_dir = "Monocle/") {
        cat("Simple Monocle Preparation without DEG specified \n")
        ct_aoi_counts <- aoi_seurat[["RNA"]]@counts
        pheno_data <- new("AnnotatedDataFrame", data = aoi_seurat@meta.data)
        ##feature data
        genes <- data.frame(gene_short_name = rownames(ct_aoi_counts))
        rownames(genes) <- rownames(ct_aoi_counts)
        genes <- new("AnnotatedDataFrame", data = genes)

        cds_ct_aoi <- newCellDataSet(ct_aoi_counts,
                        phenoData = pheno_data, featureData = genes,
                        expressionFamily = uninormal())
        cds_ct_aoi <- detectGenes(cds_ct_aoi, min_expr = 1)
        cds_ct_aoi <- cds_ct_aoi[fData(cds_ct_aoi)$num_cells_expressed > 50, ]
        cds_ct_aoi <- estimateSizeFactors(cds_ct_aoi)
        cds_ct_aoi <- order_monocle_by_markers(aoi_seurat = aoi_seurat,
                        cds_ct_aoi = cds_ct_aoi,
                        celltype_options = celltype_options,
                        output_dir = output_dir)
        save(cds_ct_aoi, file = paste0(output_dir, celltype_options, "_cds.RData"))
        return(cds_ct_aoi)
}



monocle_BPCell_preparation <- function(aoi_counts, aoi_meta) {
        # phenodata
        pheno_data <- new("AnnotatedDataFrame", data = aoi_meta)
        ##feature data
        genes <- data.frame(gene_short_name = rownames(aoi_counts))
        rownames(genes) <- rownames(aoi_counts)
        genes <- new("AnnotatedDataFrame", data = genes)

        cds_ct_aoi <- newCellDataSet(aoi_counts,
                        phenoData = pheno_data, featureData = genes,
                        expressionFamily = uninormal())
        cds_ct_aoi <- detectGenes(cds_ct_aoi, min_expr = 1)
        cds_ct_aoi <- cds_ct_aoi[fData(cds_ct_aoi)$num_cells_expressed > 50, ]
        cds_ct_aoi <- estimateSizeFactors(cds_ct_aoi)
        return(cds_ct_aoi)
}



cell_trajectory_paramPlot <- function(cds_ct_aoi, param = "groups",
                        output_dir = output_dir) {
        p1 <- plot_cell_trajectory(cds_ct_aoi,
                color_by = param,
                cell_name_size = 0.5,
                legend.key.size = unit(2.5, "cm"),
                legend.text = element_text(size = 3),
                cell_size = 0.1) +
                facet_wrap(~groups, nrow = 1)
        p2 <- plot_cell_trajectory(cds_ct_aoi,
                        color_by = "Pseudotime",
                        cell_size = 0.1,
                        cell_link_size = 0.25,
                        cell_name_size = 0.5)
        p3 <- plot_cell_trajectory(cds_ct_aoi,
                color_by = "State",
                cell_size = 0.1,
                cell_link_size = 0.25,
                cell_name_size = 0.5) +
                scale_color_manual(values = colour)
        p_title <- paste0(output_dir, param, "_trajectory.png")
        p <- p1 + p2 + p3 + plot_layout(ncol = length(p1$layers) + 1)
        ggsave(filename = p_title,
                plot = p1 + p2 + p3, width = 14, height = 8, units = "in")
}

branching_aoi <- function(cds_ct_aoi, branch_options = 1,
                        param = "groups",
                    topnum = 100,
                    output_dir = "Monocle/") {
    cat("Calculate BEAM on Current User's Selected Branch Point\n")
    branch_beam <- BEAM(cds = cds_ct_aoi,
                progenitor_method = "sequential_split",
                branch_point = branch_options, cores = 8)

    sig_branch_beam <- subset(branch_beam, qval < 1e-4)

    write.csv(sig_branch_beam, paste0(output_dir,
            "beam_branch_point_branchrow", branch_options,
            ".csv"))
#     branch_beam <- deg_filter(branch_beam)
    print(head(branch_beam))
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
            "beam_branch_point", branch_options, "heatmap.png")
    ggsave(p_title, p$ph_res,
            width = 16, height = 12, units = "in")

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
        plot = p, width = 16, height = 12, units = "in")
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
#      p$ph_res <- p$ph_res + theme(axis.text = element_text(size = 30, angle = 45),
#                         legend.key.size = unit(1.5, "cm"),
#                         legend.text = element_text(size = 10),
#                         plot.title = element_text(size = 20, hjust = 0.5),
#                         panel.grid.major = element_line(colour = "grey"),
#                         panel.grid.minor = element_line(colour = "grey"),
#                         panel.background = element_rect(fill = "white"),
#                         plot.subtitle = element_text(size = 14))
    p_title <- paste0(output_dir, 
            "beam_branch_point", branch_options, "_heatmap.png")
    ggsave(p_title, p$ph_res,
            width = 8, height = 6, units = "in")

    cat("Plot Branching Gene Heatmap \n")
    # brow_path <- paste0(output_dir,
    #         beam_header,
    #          branch_options, "_branchrow.csv")
    # if (!file.exists(brow_path)) {
    #         branched_row <- p$annotation_row
    #         branched_row <- data.frame(cluster = branched_row$Cluster,
    #                                 gene = row.names(branched_row),
    #                                 stringsAsFactors = FALSE)
    #         write.csv(branched_row,
    #                 paste0(output_dir,
    #                 beam_header,
    #                 branch_options, "_branchrow.csv"))
    # } else {
    #     branched_row <- read.csv(brow_path)
    # }
    # # branching cluster of current selected branch point
    # pdf(paste0(output_dir,
    #     "beam_branch_genes_branchpoint",
    #     branch_options, "_heatmap.pdf"))
    # branched_gene <- branched_row$gene
    # p <- plot_pseudotime_heatmap(cds_ct_aoi[branched_gene, ],
    #                     num_clusters = 4,
    #                     show_rownames = TRUE)
    # p
    # dev.off()
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
