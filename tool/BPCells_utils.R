library(BPCells)
library(Seurat)
library(SeuratObject)
library(Matrix)
library(dplyr)
library(harmony)
library(data.table)
options(Seurat.object.assay.version = "v5")

BPCells_Routine <- function(seurat) {
        seurat$percent.mt <- PercentageFeatureSet(
                    object = seurat, pattern = "^MT-")
        seurat$log10GenesPerUMI <- log10(
        seurat$nFeature_RNA) / log10(seurat$nCount_RNA)

        seurat <- subset(seurat, percent.mt <= 10)

        seurat <- NormalizeData(seurat,
                                normalization.method = "LogNormalize",
                                scale.factor = 1e4)
        seurat <- FindVariableFeatures(seurat)
        seurat <- ScaleData(seurat)
        return(seurat)
}

#' BPCells_Preparation
#'
#' @counts c(IterableMatrix, dgCMatrix)
#'
#' @meta metadata
#' @storeDir dir to store all necessary files
#'

BPCells_Preparation <- function(counts, meta, storeDir) {
        cat("Generate New Seurat Object\n")
        # append calculation into a seurat object
        seurat <- CreateSeuratObject(counts = counts, meta.data = meta)
        seurat <- BPCells_Routine(seurat)
        # PCA
        pca_path <- file.path(storeDir, "PCA_Embeddings.csv")
        if (!file.exists(pca_path)) {
                cat("Calculate PCA\n")
                # mat_norm <- counts_mat[variable_genes, ]
                # mat_norm <- mat_norm %>% write_matrix_dir(tempfile("mat"))
                # gene_means <- stats$row_stats["mean", variable_genes]
                # gene_vars <- stats$row_stats["variance", variable_genes]
                # mat_norm <- (mat_norm - gene_means) / gene_vars
                # svd <- irlba::irlba(mat_norm, nv = 50)
                # pca <- multiply_cols(svd$v, svd$d)
                # rownames(pca) <- rownames(meta)
                # colnames(pca) <- paste0("PC_", seq(1, ncol(pca)))
                seurat <- RunPCA(seurat, features = VariableFeatures(seurat))
                pca <- Embeddings(seurat, "pca")
                fwrite(pca, pca_path)
        } else {
                pca <- as.matrix(fread(pca_path))
                rownames(pca) <- rownames(meta)
                seurat[["pca"]] <- CreateDimReducObject(embeddings = pca,
                                global = TRUE, assay = "RNA")
        }
        # Harmony Integration
        harmony_path <- file.path(storeDir, "Harmony_Embeddings.csv")
        if (!file.exists(harmony_path)) {
                cat("Performing Harmony Integration via BPCells PCA\n")
                harmony_embeddings <- HarmonyMatrix(pca, meta, "orig.ident", do_pca = FALSE)
                rownames(harmony_embeddings) <- rownames(pca)
                colnames(harmony_embeddings) <- paste0("harmony_",
                        seq(1, ncol(harmony_embeddings)))
                # Could Directly Append Harmony Embeddings into whatever template you want, seurat/scanpy
                fwrite(harmony_embeddings, harmony_path)
        } else {
                harmony_embeddings <- as.matrix(fread(harmony_path))
                rownames(harmony_embeddings) <- rownames(pca)
        }
        seurat[["harmony"]] <- CreateDimReducObject(embeddings = harmony_embeddings, global = T, assay = "RNA")
        umap_path <- file.path(storeDir, "Harmony_UMAP.csv")
        if (!file.exists(umap_path)) {
                # UMAP on harmony corrected embeddings
                cat("UMAP via Harmony Embeddings\n")
                seurat <- seurat %>%
                        FindNeighbors(reduction = "harmony", dims = 1:20) %>%
                        RunUMAP(reduction = "harmony", dims = 1:20)
                umap <- Embeddings(seurat, "umap")
                # rownames(umap) <- rownames(pca)
                # colnames(umap) <- paste0("UMAP_", seq(1, ncol(umap)))
                fwrite(umap, umap_path)
        } else {
                umap <- as.matrix(fread(umap_path))
                rownames(umap) <- rownames(pca)
                seurat[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", global = T, assay = "RNA")
        }
        return(seurat)
}

#' Seurat2BPCell
#' Specifically used to update Seurat object to Seurat V5 object
#' @bp_path input to read in Seurat object, should be .rds type
#'

Seurat2BPCell <- function(bp_path, bp_seurat = NULL) {
        library(data.table)
        bp_dir_path <- fs::path_dir(bp_path)
        if (is.null(bp_seurat)) {
                bp_seurat <- readRDS(bp_path)
        }
        if (class(bp_seurat[["RNA"]]) == "Assay5") {
                cat("Already with BPCell Counts\n")
                        bp_counts <- bp_seurat[["RNA"]]@layers$counts
        } else {
                cat("Prepare BPCell Counts\n")
                bp_counts_path <- gsub(".rds", "BPCell_counts", bp_path)
                if (!file.exists(bp_counts_path)) {
                        bp_counts <- bp_seurat[["RNA"]]$counts
                        bp_counts <- as(as.matrix(bp_counts), "dgCMatrix")
                        write_matrix_dir(bp_counts, bp_counts_path)
                }
                bp_counts <- open_matrix_dir(bp_counts_path)
        }
        bp_meta <- bp_seurat@meta.data
        fwrite(bp_meta, gsub(".rds", "BPCell_metadata.csv", bp_path), row.names = TRUE)
        cat("Integration and UMAP time cost:", system.time({bp_seurat <- BPCells_Preparation(counts = bp_counts, meta = bp_meta, storeDir = bp_dir_path) }))
        bp_seurat_path <- gsub(".rds", "BPCell_seurat.rds", bp_path)
        saveRDS(
                object = bp_seurat,
                file = bp_seurat_path
        )
        cat("Successfully Run Seurat2BPCell\n")
        return(bp_seurat)
}

#' Path2BPCell
#' An integrated function to convert .csv/.h5ad/.rds objects to BPCell object
#' and finally to Seurat V5 object or Scanpy object.
#' 
#' @bp_path path to load necessary files
#'
#'

Path2BPCell <- function(bp_path) {
        library(data.table)
        count_types <- tools::file_ext(bp_path)
        bp_dir_path <- fs::path_dir(bp_path)
        if (count_types == "h5ad") {
          cat("Reading an h5ad count\n")
          bp_counts <- open_matrix_anndata_hdf5(bp_path)
          bp_meta <- as.data.frame(fread(gsub(".h5ad", "_meta.csv", bp_path)))
          rownames(bp_meta) <- bp_meta[, 1]
          bp_meta <- bp_meta[, -1]
          cat( system.time({bp_seurat <- BPCells_Preparation(counts = bp_counts,
                        meta = bp_meta, storeDir = bp_dir_path) }), "\n")
          bp_seurat_path <- gsub("scanpy.h5ad", "BPCell_seurat.rds", bp_path)
          saveRDS(
                  object = bp_seurat,
                  file = bp_seurat_path
          )
        } else if (count_types == "rds") {
          try(Seurat2BPCell(bp_path))
          return(NULL)
        } else if (count_types == "csv") {
          # TODO:Finish Counts Scenario
          # bp_counts <- fread(bp_counts)
          # bp_counts <- as(as.matrix(bp_counts), "dgCMatrix")
        } else {
          stop("Unknown file type, cannot read counts, please double check your input file\n")
        }
}

#' Counts_Data_Extract
#' Specifically used to extract Seurat V3/V4/V5 objects to proper format
#' @bp_path input to read in Seurat object, should be .rds type
#'

Counts_Data_Extract <- function(aoi_seurat, module, storeDir) {
        # align meta data first
        

        # CellChat Requires Normalized Data
        if (module != "CellChat") {
                if (class(aoi_seurat[["RNA"]]) == "Assay5") {
                        cat("Converting Seurat V5 to h5ad preserving all information \n")
                        counts <- aoi_seurat[["RNA"]]@layers$counts     
                } else {
                        cat("Converting Seurat V3 to h5ad preserving all information \n")
                        counts <- aoi_seurat[["RNA"]]@counts
                }
        } else {
                if (class(aoi_seurat[["RNA"]]) == "Assay5") {
                        cat("Converting Seurat V5 to h5ad preserving all information \n")
                        counts <- aoi_seurat[["RNA"]]@layers$counts     
                } else {
                        cat("Converting Seurat V3 to h5ad preserving all information \n")
                        counts <- aoi_seurat[["RNA"]]@counts
                }
        }
        # output
        
}


# https://bnprks.github.io/BPCells/reference/index.html
bpcells_DEG <- function(tp_data, cur_meta) {
    deg_table <- marker_features(tp_data,
            groups = cur_meta$groups,
            method = "wilcoxon")
    clean_table <- deg_table %>%
                        # select(foreground_mean, background_mean) %>%
                        filter(p_val_raw <= 0.05) %>%
                        mutate(gene = feature,
                                cluster = foreground,
                                avg_log2FC = (foreground_mean - background_mean)/log(2),
                                # pct.1 = foreground_mean,
                                # pct.2 = apply(data.group2, 1, function(x) sum(x > 0) / length(group2)),
                                p_val = p_val_raw,
                                p_val_adj = p.adjust(p_val_raw, method = "BH", n = length(p_val_raw))) %>%
                        # we don't need duplicates
                        filter(avg_log2FC >= 0) %>%
                        select(gene, p_val_adj, avg_log2FC, cluster, p_val)
    return(clean_table)
}


GroupBPDEG_Calculation <- function(aoi_seurat = NULL,
                    aoi_counts = NULL, top_num = 20,
                    log2fc = 0.5, pval = 0.05,
                    celltype_options = "all", output_dir) {
    cur_meta <- aoi_seurat@meta.data    
    if (is.null(aoi_counts)) {
        stop("Double check your cur_bpcounts!\n")
        return(NULL)
    }
    counts <- aoi_counts
    cat("Normalize current BPCounts\n")
    # count matrix normalization
    counts_mat <- multiply_cols(counts, 1 / colSums(counts))
    counts_mat <- log1p(counts_mat * 10000) # Log normalization
    all_flags <- FALSE
    if (celltype_options == "all") {
        celltype_options <- unique(aoi_seurat$celltype)
        all_flags <- TRUE
    }

    ps <- list()
    for (i in 1:length(celltype_options)) {
        cat("Calculate BPCell DEG for celltype:", celltype_options[i], "\n")
        sub_meta <- subset(cur_meta, celltype == celltype_options[i])
        sub_tp <- counts_mat[, sub_meta$idx]
        clean_table <- bpcells_DEG(sub_tp, sub_meta)
        clean_table$celltype <- celltype_options[i]
        idx1 <- which(sub_meta$groups == unique(sub_meta$groups)[1])
        idx2 <- which(sub_meta$groups != unique(sub_meta$groups)[1])
        sub_tp <- sub_tp[unique(clean_table$gene), ]
        group1 <- sub_tp[, idx1]
        group1 <- drop0(as(group1, "dgCMatrix"))
        group2 <- sub_tp[, idx2]
        group2 <- drop0(as(group2, "dgCMatrix"))
        pct.1 <- Matrix::rowSums(sign(group1)) / ncol(group1)
        pct.1 <- pct.1[match(clean_table$gene, names(pct.1))]
        pct.2 <- Matrix::rowSums(sign(group2)) / ncol(group2)
        pct.2 <- pct.2[match(clean_table$gene, names(pct.2))]
        # to fit seurat deg pct naming
        clean_table$pct.1 <- pct.2
        clean_table$pct.2 <- pct.1
        ps[[i]] <- clean_table
        write.csv(clean_table, paste0(output_dir, celltype_options[i], "_groupwise_DEG_table.csv"))
        cat("Plot Heatmap with BPCells DEG\n")
        clean_table <- deg_filter(clean_table)
        if (nrow(clean_table) == 0) {
                png(paste0(output_dir, "sig_", celltype_options[i], "_DEG_heatmap.png"))
                cat("Warning: No DEG identified \n")
                dev.off()
                return(NULL)
        }
        clean_table <- subset(clean_table, pct.1 >= 0.05 & pct.2 >= 0.05 & avg_log2FC >= log2fc & p_val_adj <= pval)
        top_markers <- clean_table %>% group_by(cluster) %>%
                slice_min(p_val_adj, n = 2*top_num)
    }
    if (all_flags == TRUE) {
        deg_table <- do.call(rbind, ps)
        deg_table <- as.data.frame(deg_table)

        write.csv(deg_table, paste0(output_dir, "sig_groupwise_DEG_table.csv"))
        return(deg_table)
    }

    return(clean_table)
}


Groupwise_BPCells_DEG <- function(aoi_seurat = NULL,
                    aoi_counts = NULL, top_num = 20,
                    log2fc = 0.5, pval = 0.05,
                    celltype_options = "all", output_dir) {
    cur_meta <- aoi_seurat@meta.data    
    if (is.null(aoi_counts)) {
        stop("Double check your cur_bpcounts!\n")
        return(NULL)
    }
    counts <- aoi_counts
    cat("Normalize current BPCounts\n")
    # count matrix normalization
    counts_mat <- multiply_cols(counts, 1 / colSums(counts))
    counts_mat <- log1p(counts_mat * 10000) # Log normalization
    all_flags <- FALSE
    if (celltype_options == "all") {
        celltype_options <- unique(aoi_seurat$celltype)
        all_flags <- TRUE
    }

    ps <- list()
    for (i in 1:length(celltype_options)) {
        cat("Calculate BPCell DEG for celltype:", celltype_options[i], "\n")
        sub_meta <- subset(cur_meta, celltype == celltype_options[i])
        sub_tp <- counts_mat[, sub_meta$idx]
        clean_table <- bpcells_DEG(sub_tp, sub_meta)
        clean_table$celltype <- celltype_options[i]
        idx1 <- which(sub_meta$groups == unique(sub_meta$groups)[1])
        idx2 <- which(sub_meta$groups != unique(sub_meta$groups)[1])
        sub_tp <- sub_tp[unique(clean_table$gene), ]
        group1 <- sub_tp[, idx1]
        group1 <- drop0(as(group1, "dgCMatrix"))
        group2 <- sub_tp[, idx2]
        group2 <- drop0(as(group2, "dgCMatrix"))
        pct.1 <- Matrix::rowSums(sign(group1)) / ncol(group1)
        pct.1 <- pct.1[match(clean_table$gene, names(pct.1))]
        pct.2 <- Matrix::rowSums(sign(group2)) / ncol(group2)
        pct.2 <- pct.2[match(clean_table$gene, names(pct.2))]
        # to fit seurat deg pct naming
        clean_table$pct.1 <- pct.2
        clean_table$pct.2 <- pct.1
        ps[[i]] <- clean_table
        write.csv(clean_table, paste0(output_dir, celltype_options[i], "_groupwise_DEG_table.csv"))
        cat("Plot Heatmap with BPCells DEG\n")
        clean_table <- deg_filter(clean_table)
        if (nrow(clean_table) == 0) {
                png(paste0(output_dir, "sig_", celltype_options[i], "_DEG_heatmap.png"))
                cat("Warning: No DEG identified \n")
                dev.off()
                return(NULL)
        }
        clean_table <- subset(clean_table, pct.1 >= 0.05 & pct.2 >= 0.05 & avg_log2FC >= log2fc & p_val_adj <= pval)
        top_markers <- clean_table %>% group_by(cluster) %>%
                slice_min(p_val_adj, n = 2*top_num)
        ct_aoi <- subset(aoi_seurat, celltype == celltype_options[i])
        ct_aoi <- ScaleData(object = ct_aoi, features = rownames(ct_aoi))
        # ct_aoi <- ScaleData(ct_aoi)
        p <- DoHeatmap(ct_aoi,
                    features = unique(top_markers$gene),
                    group.by = "groups",
                    size = 8,
                    disp.max = 2,
                    # slot = "data",
                    assay = "RNA") +
                    scale_fill_gradientn(colors = c("blue", "white", "red")) +
                theme(text = element_text(size = 25))
        p_title <- paste0(output_dir, "sig_",
                celltype_options[i], "_DEG_heatmap.png")
        ggsave(filename = p_title,
            plot = p, width = 16, height = 12, units = "in")
    }

    if (all_flags == TRUE) {
        deg_table <- do.call(rbind, ps)
        deg_table <- as.data.frame(deg_table)

        write.csv(deg_table, paste0(output_dir, "sig_groupwise_DEG_table.csv"))
        return(deg_table)
    }

    return(clean_table)
}

BPCells_Topmarker <- function(aoi_seurat = NULL,
                    aoi_counts = NULL, top_num = 5,
                    log2fc = log2fc, pval = pval,
                    output_dir) {
    cur_meta <- aoi_seurat@meta.data
    if (is.null(aoi_counts)) {
        stop("Double check your cur_bpcounts!\n")
        return(NULL)
    }
    counts <- aoi_counts
    cat("Normalize current BPCounts\n")
    # count matrix normalization
    counts_mat <- multiply_cols(counts, 1 / colSums(counts))
    counts_mat <- log1p(counts_mat * 10000) # Log normalization
    cur_meta$groups <- cur_meta$celltype
    clean_table <- bpcells_DEG(counts_mat, cur_meta)
    idx1 <- which(cur_meta$groups == unique(cur_meta$groups)[1])
    idx2 <- which(cur_meta$groups != unique(cur_meta$groups)[1])
    counts_mat <- counts_mat[unique(clean_table$gene), ]
    group1 <- counts_mat[, idx1]
    group1 <- drop0(as(group1, "dgCMatrix"))
    group2 <- counts_mat[, idx2]
    group2 <- drop0(as(group2, "dgCMatrix"))
    pct.1 <- Matrix::rowSums(sign(group1)) / ncol(group1)
    pct.1 <- pct.1[match(clean_table$gene, names(pct.1))]
    pct.2 <- Matrix::rowSums(sign(group2)) / ncol(group2)
    pct.2 <- pct.2[match(clean_table$gene, names(pct.2))]
    # to fit seurat deg pct naming
    clean_table$pct.1 <- pct.2
    clean_table$pct.2 <- pct.1
    write.csv(clean_table, paste0(output_dir,
            "TopMarkers_table.csv"))
    cat("Plot Heatmap with BPCells TopMarkers\n")
    clean_table <- deg_filter(clean_table)
    if (nrow(clean_table) == 0) {
            pdf(paste0(output_dir, "TopMarkers_heatmap.png"))
            cat("Warning: No DEG identified \n")
            dev.off()
            return(NULL)
    }
    clean_table <- subset(clean_table, pct.1 >= 0.05 & pct.2 >= 0.05 & avg_log2FC >= log2fc & p_val_adj <= pval)
    top_markers <- clean_table %>% group_by(cluster) %>%
            slice_max(avg_log2FC, n = top_num)
    if (ncol(aoi_seurat) > 5000) {
      cat("Downsample current seurat object\n")
      aoi_seurat <- subset(aoi_seurat, downsample = 1500)
    }
    p <- DoHeatmap(aoi_seurat,
                features = top_markers$gene,
                group.by = "celltype",
                size = 8,
                assay = "RNA") +
                scale_fill_gradientn(colors = c("blue", "white", "red")) +
                theme(text = element_text(size = 16),
                       axis.text.y = element_text(size = 20)) # Increase this value to increase y-axis text size)
    p_title <- paste0(output_dir, "TopMarkers_heatmap.png")
    ggsave(filename = p_title,
        plot = p, width = 16, height = 16, units = "in")

    p <- DotPlot(aoi_seurat,
                    group.by = "celltype",
                    features = unique(top_markers$gene)) +
                    coord_flip() +
                    theme(axis.title = element_text(size = 12),
                    axis.text = element_text(size = 14),
                    axis.text.x = element_text(size = 24, vjust = 0.5, angle = 45),
                    axis.text.y = element_text(size = 12, vjust = 0.5, angle = 0),
                    legend.key.size = unit(2, "cm"),
                    legend.text = element_text(size = 24),
                    plot.title = element_text(size = 20, hjust = 0.5))
    p_title <- paste0(output_dir, "TopMarkers_Dotplot.png")
    ggsave(filename = p_title,
        plot = p, width = 20, height = 20, units = "in")
}