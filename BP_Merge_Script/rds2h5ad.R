#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(reticulate)
library(Seurat)
library(BPCells)
library(qs)
rds2h5ad <- function(rds_path, aoi_seurat = NULL, obsm = FALSE, if_norm = FALSE, overwrite = FALSE) {
    if (!exists("sc")) {
        sc <- import("scanpy")
    }
    fname <- gsub(".rds|.Rdata|.qs", ".h5ad", rds_path, ignore.case = TRUE)
     if (!grepl("\\.h5ad$", fname)) {
        stop("Error: The file name must end with .h5ad")
    }
    if (file.exists(fname) && overwrite == FALSE) {
        cat("Current RDS object has already been converted to Scanpy object\n")
        return(NULL)
    }
    if (is.null(aoi_seurat)) {
        if (grepl("\\.qs$", fname)) {
            aoi_seurat <- qread(rds_path)
        } else {
            aoi_seurat <- readRDS(rds_path)
        }
    }
    cat("Current RDS object:", rds_path, "\n")
    if (class(aoi_seurat[["RNA"]]) == "Assay5") {
        cat("Converting Seurat V5 to h5ad preserving all information \n")
        if (if_norm == TRUE) {
            counts <- aoi_seurat[["RNA"]]@layers$scale.data
        } else {
            counts <- aoi_seurat[["RNA"]]@layers$counts
        }
        row.names <- counts@features
        counts <- as(counts, "dgCMatrix")
    } else {
        cat("Converting Seurat V3 to h5ad preserving all information \n")
        if (if_norm == TRUE) {
            idx <- rownames(aoi_seurat) %in% VariableFeatures(aoi_seurat)
            counts <- aoi_seurat[["RNA"]]@scale.data[idx, ]
        } else {
            counts <- aoi_seurat[["RNA"]]@counts
        }
        row.names <- rownames(counts)
    }
    adata_seurat <- sc$AnnData(
                    X   = t(as.sparse(counts)),
                    obs = aoi_seurat[[]],
                    var = data.frame(row.names = row.names))
    if (obsm == TRUE) {
        names(aoi_seurat@reductions) <- tolower(names(aoi_seurat@reductions))
        redc_names <- names(aoi_seurat@reductions)
        if ("pca" %in% redc_names) {
            adata_seurat$obsm$update(X_pca = Embeddings(aoi_seurat, "pca"))
        }
        if ("umap" %in% redc_names) {
            adata_seurat$obsm$update(X_umap = Embeddings(aoi_seurat, "umap"))
        }
        if ("tsne" %in% redc_names) {
            adata_seurat$obsm$update(X_tsne = Embeddings(aoi_seurat, "tsne"))
        }
        if ("cca" %in% redc_names) {
            adata_seurat$obsm$update(X_cca = Embeddings(aoi_seurat, "cca"))
        }
        if ("cca_aligned" %in% redc_names) {
            adata_seurat$obsm$update(X_cca_aligned = Embeddings(aoi_seurat, "cca.aligned"))
        }
        if ("harmony" %in% redc_names) {
            adata_seurat$obsm$update(X_harmony = Embeddings(aoi_seurat, "harmony"))
        }
         if ("umap-mvg2000" %in% redc_names) {
            adata_seurat$obsm$update(X_umap_mvg2000 = Embeddings(aoi_seurat, "umap-mvg2000"))
        }
    }
    adata_seurat$write_h5ad(fname, compression = "gzip", compression_opts = 5)
    cat("Successfully convert Seurat to Scanpy\n")
    return(NULL)
}

rds_list <- args[1]
rds_list <- data.table::fread(rds_list, header = FALSE)[, 1]
rds_list <- grep("BPCell|.h5ad", rds_list, invert = TRUE, value = TRUE)
if (length(rds_list) == 0) {
    stop("Current pattern cannot be retrieved from stored rds data\n")
}

lapply(rds_list, rds2h5ad, adata_flag = FALSE)
