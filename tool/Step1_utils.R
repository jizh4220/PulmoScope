
top_markers_Dotplot_heatmap <- function(output_dir = output_dir,
            aoi_seurat, top_num = 5, pval = 0.05) {
    # Marker Genes Representation via Violinplot
    Idents(aoi_seurat) <- aoi_seurat$macro_celltype
    aoi_seurat <- subset(aoi_seurat, downsample = 2000)
    if (!file.exists(
        paste0(output_dir,
        "TopMarkers_table.csv"))) {
        deg_table <- FindAllMarkers(aoi_seurat,
          logfc.threshold = 0, test.use = "wilcox", min.pct = 0.1,
          min.diff.pct = 0,
          only.pos = TRUE,
          return.thresh = 0.05,
          assay = "RNA")
        # if (deg_table)
        deg_table <- subset(deg_table,
            p_val_adj < 0.05)
        write.csv(deg_table,
            paste0(output_dir,
            "TopMarkers_table.csv"))
    } else {
        cat("Found TopMarkers Table, Plot Heatmap \n")
        deg_table <- read.csv(paste0(output_dir,
             "TopMarkers_table.csv"))
    }
    deg_table <- deg_filter(deg_table)
    deg_table <- subset(deg_table,
            p_val_adj <= pval)
    top_markers <- deg_table %>% group_by(cluster) %>%
        slice_max(avg_log2FC, n = top_num)

    p <- DotPlot(aoi_seurat,
                    group.by = "macro_celltype",
                    features = unique(top_markers$gene)) +
                    coord_flip() +
                    theme(axis.title = element_text(size = 15),
                    axis.text = element_text(size = 18, angle = 45),
                    legend.key.size = unit(4, "cm"),
                    legend.text = element_text(size = 26),
                    plot.title = element_text(size = 20, hjust = 0.5))
    p_title <- paste0(output_dir, "TopMarkers_Dotplot.png")
    ggsave(filename = p_title,
        plot = p, width = 25, height = 30, units = "in")
    #Marker Genes Representation via Heatmap
    
    p <- DoHeatmap(aoi_seurat,
                features = top_markers$gene,
                group.by = "macro_celltype",
                size = 8,
                assay = "RNA") +
                scale_fill_gradientn(colors = c("blue", "white", "red")) +
                theme(text = element_text(size = 25))
    p_title <- paste0(output_dir, "TopMarkers_heatmap.png")
    ggsave(filename = p_title,
        plot = p, width = 25, height = 30, units = "in")
    return(deg_table)
}

acc_aoi_seurat <- function(all_folders, sample_info, threshold = 10000) {
    aoi_acc <- unique(sample_info[, "accession"])
    aoi_group <- sample_info[, "groups"]
    ps <- list()
    threshold <- threshold / nrow(sample_info)
    for (i in 1:length(aoi_acc)) {
        ps[[i]] <- locate_seurat_per_acc(aoi_acc[i],
                    seurat_folders = all_folders,
                    group_info = aoi_group[i],
                    threshold)
    }
    ps <- ps[!sapply(ps, is.null)]
    if (is.null(ps)) {
        cat("Cannot Retrieve Accessions For Now \n")
        return(NULL)
    }
    cat("Successfully Retrieve All Accession Seurat \n")
    seurat <- merge(ps[[1]], ps[-1])
    seurat <- harmony_merge(seurat)
    seurat <- seurat %>%
            FindNeighbors(reduction = "harmony", dims = 1:20) %>%
            RunUMAP(reduction = "harmony", dims = 1:20)
    # clean meta data
    meta <- seurat@meta.data
    meta <- meta_alignment(meta)
    seurat@meta.data <- meta
    suppressWarnings(dir.create("data", recursive = TRUE))
    readr::write_rds(seurat, "data/aoi_seurat.rds")
    aoi_table <- ct_disease_acc_table(seurat, 0)
    aoi_table <- as.data.frame(aoi_table)
    aoi_table$celltype <- rownames(aoi_table)
    # jsonData <- toJSON(aoi_table)
    write.csv(aoi_table,
            "Step1/celltype_disease_acctable.csv")
    write.table(rownames(seurat),
            "data/gene_names.txt",
            row.names = FALSE, col.names = FALSE)
    cat("Successfully generate AOI seurat \n")
    cat("############Generating Basic Plots############\n")
    return(seurat)
}



BPCell_Integrate <- function(all_folders, sample_info, threshold = 10000) {
    aoi_acc <- unique(sample_info[, "accession"])
    aoi_group <- sample_info[, "groups"]
    ps <- list()
    threshold <- threshold / nrow(sample_info)
    for (i in 1:length(aoi_acc)) {
        ps[[i]] <- locate_seurat_per_acc(aoi_acc[i],
                    seurat_folders = all_folders,
                    group_info = aoi_group[i],
                    threshold)
    }
    ps <- ps[!sapply(ps, is.null)]
    if (is.null(ps)) {
        cat("Cannot Retrieve Accessions For Now \n")
        return(NULL)
    }
    cat("Successfully Retrieve All Accession Seurat \n")
    seurat <- merge(ps[[1]], ps[-1])
    bp_counts <- bp_seurat[["RNA"]]$counts
    bp_counts <- as(as.matrix(bp_counts), "dgCMatrix")
    write_matrix_dir(bp_counts, bp_counts_path)
    seurat <- harmony_merge(seurat)
    seurat <- seurat %>%
            FindNeighbors(reduction = "harmony", dims = 1:20) %>%
            RunUMAP(reduction = "harmony", dims = 1:20)
    # clean meta data
    meta <- seurat@meta.data
    meta <- meta_alignment(meta)
    seurat@meta.data <- meta
    suppressWarnings(dir.create("data", recursive = TRUE))
    readr::write_rds(seurat, "data/aoi_seurat.rds")
    aoi_table <- ct_disease_acc_table(seurat, 0)
    aoi_table <- as.data.frame(aoi_table)
    aoi_table$celltype <- rownames(aoi_table)
    # jsonData <- toJSON(aoi_table)
    write.csv(aoi_table,
            "Step1/celltype_disease_acctable.csv")
    write.table(rownames(seurat),
            "data/gene_names.txt",
            row.names = FALSE, col.names = FALSE)
    cat("Successfully generate AOI seurat \n")
    cat("############Generating Basic Plots############\n")
    return(seurat)
}

Step1_Run <- function(output_dir = "Step1/", aoi_seurat) {
    suppressWarnings(dir.create(output_dir, recursive = TRUE))
    # disease_table <- table(aoi_seurat$DISEASE, aoi_seurat$orig.ident)
    # disease_table <- as.data.frame(disease_table)
    # colnames(disease_table) <- c("group", "accession", "counts")
    # write.csv(disease_table,
    #         "Step1/Disease_acctable.csv", row.names = TRUE)
    # tissue_table <- table(aoi_seurat$TISSUE, aoi_seurat$orig.ident)
    # tissue_table <- as.data.frame(tissue_table)
    # colnames(tissue_table) <- c("group", "accession", "counts")
    # write.csv(tissue_table,
    #         "Step1/Tissue_acctable.csv", row.names = TRUE)
    # gender_table <- table(aoi_seurat$gender, aoi_seurat$orig.ident)
    # gender_table <- as.data.frame(gender_table)
    # colnames(gender_table) <- c("group", "accession", "counts")
    # write.csv(gender_table,
    #         "Step1/Gender_acctable.csv", row.names = TRUE)
    # aging_table <- table(aoi_seurat$age, aoi_seurat$orig.ident)
    # aging_table <- as.data.frame(aging_table)
    # colnames(aging_table) <- c("group", "accession", "counts")
    # write.csv(aging_table,
    #         "Step1/Age_acctable.csv", row.names = TRUE)
    # accession, disease, tissue, gender, aging, celltype, experiment, group
    meta <- aoi_seurat@meta.data
    meta <- meta_alignment(meta)
    if ("orig.ident" %in% colnames(meta)) {
        # colnames(meta)[colnames(meta) == "orig.ident"] <- "accession"
        # colnames(meta)[colnames(meta) == "gse_alias"] <- "experiment"
        # colnames(meta)[colnames(meta) == "DISEASE"] <- "disease"
        # colnames(meta)[colnames(meta) == "TISSUE"] <- "tissue"
        # colnames(meta)[colnames(meta) == "celltype"] <- "detail_celltype"
        # colnames(meta)[colnames(meta) == "macro_celltype"] <- "celltype"
        aoi_seurat@meta.data <- clean_metanoise(meta)
    }
    # meta_param <- c("disease", "tissue", "aging", "gender",
    #             "experiment", "celltype", "groups", "accession")
    # meta_param <- meta_param[meta_param %in% colnames(meta)]
    meta_param <- colnames(aoi_seurat@meta.data)
    # meta dimplot
    # if_meta_dim <- metadata_dimplot(output_dir = output_dir,
    #     aoi_seurat = aoi_seurat, dimplot_param = meta_param)
    cat("Plotting Celltype Compo and Meta Plot", "\n")
    Compo_Meta_plot(output_dir = output_dir, aoi_seurat = aoi_seurat,
            param = meta_param)
    
    # compo_plot(aoi_seurat = aoi_seurat,
    #                 output_dir = output_dir)
    # if (!file.exists(
    #         paste0(output_dir,
    #         "Celltype_marker_table.csv"))) {
    #     cat("Calculating Top Marker and Plotting Dotplot and Heatmap ", "\n")
    #     deg_table <- top_markers_Dotplot_heatmap(output_dir = output_dir,
    #                     aoi_seurat = aoi_seurat,
    #                     top_num = 10,
    #                     deg_table = FALSE)
    # }
    cat("Finish Calculating Top Marker and Plotting Dotplot and Heatmap \n")
}

# Step1_Module <- function(aoi_seurat, output_dir) {
#     plot_aoi_seurat(output_dir = "Step1/", aoi_seurat)
# }