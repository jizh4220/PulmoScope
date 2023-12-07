#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
library(BPCells)
library(Seurat)
library(Matrix)
library(dplyr)
library(harmony)

options(Seurat.object.assay.version = "v5")

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
        seurat <- NormalizeData(seurat)
        seurat <- FindVariableFeatures(seurat)
        cat("Normalization and Variable Genes\n")
        counts_mat <- multiply_cols(counts, 1 / colSums(counts))
        counts_mat <- log1p(counts_mat * 10000) # Log normalization
        stats <- matrix_stats(counts_mat, row_stats = "variance")
        variable_genes <- order(stats$row_stats["variance", ],
                                decreasing = TRUE) %>%
                                head(2000) %>%
                                sort()

        # PCA
        pca_path <- file.path(storeDir, "PCA_Embeddings.csv")
        if (!file.exists(pca_path)) {
                cat("Calculate PCA\n")
                mat_norm <- counts_mat[variable_genes, ]
                mat_norm <- mat_norm %>% write_matrix_dir(tempfile("mat"))
                gene_means <- stats$row_stats["mean", variable_genes]
                gene_vars <- stats$row_stats["variance", variable_genes]
                mat_norm <- (mat_norm - gene_means) / gene_vars
                svd <- irlba::irlba(mat_norm, nv = 50)
                pca <- multiply_cols(svd$v, svd$d)
                rownames(pca) <- rownames(meta)
                colnames(pca) <- paste0("PC_", seq(1, ncol(pca)))
                fwrite(pca, pca_path)
        } else {
                pca <- as.matrix(fread(pca_path))
                rownames(pca) <- rownames(meta)
        }
        seurat[["pca"]] <- CreateDimReducObject(embeddings = pca,
            global = TRUE, assay = "RNA")
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
                umap <- uwot::umap(harmony_embeddings)
                rownames(umap) <- rownames(pca)
                colnames(umap) <- paste0("UMAP_", seq(1, ncol(umap)))
                fwrite(umap, umap_path)
        } else {
                umap <- as.matrix(fread(umap_path))
                rownames(umap) <- rownames(pca)
        }
        seurat[["umap"]] <- CreateDimReducObject(embeddings = umap, key = "UMAP_", global = T, assay = "RNA")
        return(seurat)
}

#' Seurat2BPCell
#' Specifically used to update Seurat object to Seurat V5 object
#' @bp_path input to read in Seurat object, should be .rds type
#'

Seurat2BPCell <- function(bp_path) {
        library(data.table)
        bp_dir_path <- fs::path_dir(bp_path)
        bp_seurat <- readRDS(bp_path)
        if (class(bp_seurat[["RNA"]]) == "Assay5") {
                cat("Already with BPCell Counts\n")
                        bp_counts <- bp_seurat[["RNA"]]@layers$counts
        } else {
                cat("Prepare BPCell Counts\n")
                bp_counts_path <- gsub(".rds", "BPCell_counts", bp_path)
                if (!file.exists(bp_counts_path)) {
                        bp_counts <- bp_seurat[["RNA"]]$counts
                        bp_counts <- as(as.matrix(bp_counts), "dgCMatrix")
                        bp_counts <- write_matrix_dir(bp_counts, bp_counts_path)
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

path <- args[1]
Path2BPCell(path)
