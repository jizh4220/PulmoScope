featureDim_plot <- function(aoi_seurat, goi,
                            reduction = "umap", param = "disease") {
  #Featureplot by param
  # Idents(aoi_seurat) <- aoi_seurat[[param]]
  p1 <- FeaturePlot(aoi_seurat,
            features = goi,
            split.by = param,
            reduction = reduction, pt.size = 0.2, repel = TRUE,
            label = TRUE, label.size = 4)
  p2 <- DimPlot(object = aoi_seurat,
              group.by = "celltype",
              reduction = reduction, pt.size = 0.2, repel = TRUE,
              label = TRUE, label.size = 4)
  p <- p1 + p2 + plot_layout(ncol = length(unique(aoi_seurat[[param]])[, 1]) + 1)
  return(p)
}

ViolinComparePlot <- function(aoi_seurat, goi,
                                 reduction = "umap", param = "disease") {
  p <- VlnPlot(aoi_seurat,
                features = goi,
                group.by = "celltype",
                split.by = param,
                # filp = TRUE,
                add.noise = TRUE,
                pt.size = 0.05,
                split.plot = TRUE,
                same.y.lims = TRUE,
                ) + theme(axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.title = element_blank(),
                          axis.text = element_text(size = 20),
                          axis.text.x = element_text(vjust = 0.9, angle = 90))
  return(p)
}



fraction_param_table <- function(seurat,
                                threshold = 0.5, param = param,
                                goi = goi) {
        # preprocess meta data names
        seurat$param <- seurat[[param]][, 1]
        meta <- seurat@meta.data
        meta_param <- c("disease", "tissue", "aging", "gender",
                        "experiment", "celltype", "groups", "accession", "param")
        frac_meta <- as.data.frame(meta[, intersect(meta_param, colnames(meta))])
        frac_meta <- distinct(frac_meta, accession, param, .keep_all = TRUE)
        seurat@meta.data <- meta

        # Every celltype cell counts in every accession
        df <- FetchData(seurat,
                vars = goi)
        df$accession <- seurat$accession
        df$param <- seurat$param
        acc_table <- as.data.frame(table(df$accession, df$param))
        frac_meta$idx <- paste0(frac_meta$param, frac_meta$accession)
        acc_table <- subset(acc_table, Freq > 0)
        acc_table$idx <- paste0(acc_table$Var2, acc_table$Var1)
        acc_table <- acc_table[match(frac_meta$idx, acc_table$idx), ]
        frac_meta$cell_counts <- acc_table$Freq
        cluster_frac <- list()
        # Every gene expressed counts
        for (i in 1:length(goi)) {
                #current gene cluster-wise expression
                tmp <- df[, i]
                exp_cells <- df[tmp >= threshold, ]
                exp_table <- as.data.frame(table(exp_cells$accession, exp_cells$param))
                exp_table <- subset(exp_table, Freq > 0)
                exp_table$idx <- paste0(exp_table$Var2, exp_table$Var1)
                goi_frac <- frac_meta
                goi_frac$gene <- goi[i]
                goi_frac$expressed_counts <- 0
                goi_frac[match(exp_table$idx, goi_frac$idx), ]$expressed_counts <- exp_table$Freq
                goi_frac$expressed_ratio <- goi_frac$expressed_counts / goi_frac$cell_counts
                cluster_frac[[i]] <- goi_frac
        }

        frac_table <- do.call(rbind, cluster_frac)
        frac_table$idx <- NULL
        return(frac_table)
}

fraction_plot <- function(frac_table, goi = goi,
                          threshold = 0.5,
                          param = "disease",
                          floor = 50
                          ) {
    frac_table$param <- frac_table[, grep(param, colnames(frac_table))]
    sample_names <- names(table(frac_table$param))
    # frac_table$ct_disease <- paste0(frac_table$celltype,
    #               "+", frac_table$disease)
    ps <- list()
    for (i in 1:length(goi)) {
        cur_frac <- frac_table[frac_table$gene == goi[i], ]
        all_sample_table <- c(length(unique(cur_frac$accession[cur_frac$param == sample_names[1]])),
                          length(unique(cur_frac$accession[cur_frac$param == sample_names[2]])))
        cur_frac <- cur_frac[cur_frac$cell_counts >= floor &
                                cur_frac$expressed_ratio > 0, ]
        if (nrow(cur_frac) == 0) {
          cat("Current gene:",
            goi[i],
            "is not found in Fraction Table \n")
          next
        }
        cur_frac$param_ct <- paste0(cur_frac$celltype, "_", cur_frac$param)
        sample_table <- c(length(unique(cur_frac$accession[cur_frac$param == sample_names[1]])),
                          length(unique(cur_frac$accession[cur_frac$param == sample_names[2]])))
        p <- ggplot(cur_frac,
                    aes(x = param_ct, y = expressed_ratio,
                    fill = param)) +
                    geom_boxplot() +
                    # geom_violin(trim = FALSE) +
                    geom_point(
                          position = position_jitter(
                          seed = 1, width = 0.2, height = 0),
                          alpha = 2,
                          size = 1.5) +
                    # geom_boxplot(width = 0.05, fill = "white") +
                    # facet_wrap(~gene) +
                    labs(x = "Genes of Interest",
                        y = "Fraction of Cells Expressing Target Genes",
                        size = 100) +
                theme(axis.text = element_text(size = 12, angle = 45),
                      axis.text.x = element_text(vjust = 0.5, angle = 45),
                      legend.key.size = unit(1.5, "cm"),
                      legend.text = element_text(size = 10),
                      plot.title = element_text(size = 20, hjust = 0.5),
                      panel.grid.major = element_line(colour = "grey"),
                      panel.grid.minor = element_line(colour = "grey"),
                      panel.background = element_rect(fill = "white"),
                      plot.subtitle = element_text(size = 14))
                # theme(panel.background = element_rect(fill = "white"))
        p <- p + labs(title = paste0("Fraction plot of GOI: ",
                                    goi[i]),
                      subtitle = paste0(
                                "Total ",
                                sample_names[1],
                                " Sample Number: ",
                                all_sample_table[1],
                                "\nTotal ",
                                sample_names[1],
                                " Sample Number Pass Cutoff: ",
                                sample_table[1],
                                "\nTotal ",
                                sample_names[2],
                                " Sample Number: ",
                                all_sample_table[2],
                                "\nTotal ",
                                sample_names[2],
                                " Sample Number Pass Cutoff: ",
                                sample_table[2],
                                # length(unique(frac_table$cell_counts)),
                                "\nTotal Cell Counts: ",
                                sum(unique(frac_table$cell_counts)),
                                "\nPass Cutoff Cell Counts: ",
                                round(sum(cur_frac$cell_counts *
                                  cur_frac$expressed_ratio), digits = 0)
                                ,
                                "\nExpression Cutoff: ", threshold,
                                "\nMinimum Cell Number: ", floor
                                ),
                     # tag = paste0("Normal Sample Numbe: ")
                    )
        ps[[i]] <- p
    }
    return(ps)
    cat("Successfully Plot Fraction Plot\n")
}
