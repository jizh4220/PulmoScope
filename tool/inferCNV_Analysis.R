#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 4) {
  stop("Four arguments must be supplied \n", call. = FALSE)
}

print("=====================================")
print("=====        inferCNV Module        =====")
print("=====================================")


run_infercnv <- function(aoi_seurat = aoi_seurat,
                        max_cells = 100,
                        ref_group = "Normal",
                        gene_order = gene_order,
                        output_dir = "inferCNV/",
                        rawcounts = NULL) {
        library(infercnv)
        anno <- as.data.frame(aoi_seurat$groups)
        write.table(anno, paste0("data/aoi_cellAnnotations.txt"),
                        sep = " ", quote = FALSE,
                        col.names = FALSE, row.names = TRUE)
        # if (is.null(rawcounts)) {
        #         rawcounts <- fread_counts("data/aoi_rawcounts.txt")
        #         rawgenes <- read.table("data/rawgenes.txt")
        #         rownames(rawcounts) <- rawgenes[, 1]
        # }
        rawcounts <- aoi_seurat[["RNA"]]@counts
        
        infercnv_obj <- CreateInfercnvObject(
                        raw_counts_matrix = as.matrix(rawcounts),
                        max_cells_per_group = as.numeric(max_cells),
                        annotations_file = "data/aoi_cellAnnotations.txt",
                        delim = " ",
                        gene_order_file = gene_order,
                        ref_group_names = c(ref_group)
                        )
        infercnv_obj <- infercnv::run(
                        infercnv_obj,
                        cutoff = 0.1, #1 for smart-seq, 0.1 for 10x-genomics,
                        out_dir = output_dir,
                        cluster_by_groups = TRUE,
                        analysis_mode = "subclusters",
                        tumor_subcluster_partition_method = "random_trees",
                        denoise = TRUE,
                        HMM = TRUE,
                        no_plot = FALSE,
                        num_threads = 20,
                        plot_steps = TRUE
                        )
}

cluster_cnv_trees <- function(output_dir = "inferCNV/", cnv_group_path) {
        input_dir <- paste0(output_dir, "Inputs/")
        suppressWarnings(dir.create(input_dir, recursive = TRUE))
        file.copy(paste0(output_dir, cnv_group_path), input_dir)
        library(reticulate)
        library(stringr)
        use_condaenv("~/miniconda3/envs/infercnv/bin/python3")
        source_python("~/lungDB_backend/tool/uphyloplot2.py")
        cluster_cell_path <- list.files("CNV_files/", full.names = TRUE)
        #make sure there is a rand_trees cell grouings
        cluster_cell_path <- grep("rand_trees",
                                cluster_cell_path, value = TRUE)
        if (length(cluster_cell_path) < 0) {
                cat("Cannot find cluster tree in current infercnv objects! \n")
        } else {
                cluster_cell_groupings <- read.csv(cluster_cell_path)
        }
        
        inputs_cellgrouings <- read.csv(paste0("Inputs/",
                                        cnv_group_path),
                                        sep = "\t")
        new_groupnames <- gsub("[A-Za-z]|_", "",
                        inputs_cellgrouings[, 1])
        inputs_cellgrouings$groups <- gsub(paste0(unique(aoi_seurat$groups)[1], "."), "", inputs_cellgrouings$cell_group_name)
        inputs_cellgrouings$groups <- gsub(paste0(unique(aoi_seurat$groups)[2], "."), "", inputs_cellgrouings$groups)
        inputs_cellgrouings$cluster_tree <- NA
        inputs_cellgrouings$cluster_val <- NA
        for (i in 1:nrow(cluster_cell_groupings)) {
                idx <- inputs_cellgrouings$groups ==
                        rownames(cluster_cell_groupings)[i]
                if (all(idx == FALSE)) {
                        next
                } else{
                        inputs_cellgrouings$cluster_val[idx] <-
                                cluster_cell_groupings[i, 1]
                        inputs_cellgrouings$cluster_tree[idx] <-
                                cluster_cell_groupings[i, 2]
                }
        }
        file.rename(paste0("CNV_files/", cnv_group_path, ".csv"),
                paste0("CNV_files/cnvtree_celltype_groups_table.csv"))
        write.csv(inputs_cellgrouings,
                paste0(
                "CNV_files/groupings_with_cluster_tree.csv"))
        return(inputs_cellgrouings)
}

project_cnvtree <- function(aoi_seurat, inputs_cellgrouings) {
        # print(head(aoi_seurat))
        aoi_seurat$cluster_tree <- "unknown"
        ds_seurat <- aoi_seurat
        # aoi_seurat$cellid <- colnames(aoi_seurat)
        # ds_seurat <- subset(aoi_seurat, idx %in% inputs_cellgrouings$cell)
        print(inputs_cellgrouings)
        ds_groupings <- inputs_cellgrouings[match(
                        ds_seurat$idx, inputs_cellgrouings$cell), ]
        #print(head(ds_groupings))
        ds_seurat$cluster_tree <- ds_groupings$cluster_tree
        print(table(ds_seurat$cluster_tree))
        ct_clustree <- table(ds_seurat$cluster_tree, ds_seurat$celltype)
        ds_clustree <- table(ds_seurat$cluster_tree, ds_seurat$groups)

        write.csv(cbind(ct_clustree, ds_clustree),
                paste0("CNV_files/cnvtree_strength_table.csv"), row.names = TRUE)
        p1 <- DimPlot(ds_seurat,
                        group.by = "groups",
                        reduction = "umap", pt.size = 0.2, repel = TRUE,
                        label = TRUE, label.size = 6)
        p2 <- DimPlot(ds_seurat,
                        group.by = "cluster_tree",
                        reduction = "umap", pt.size = 0.2, repel = TRUE,
                        label = TRUE, label.size = 6)
        p_title <- paste0("CNV_files/", "cnvtree_dimplot.png")
        ggplot2::ggsave(filename = p_title,
                plot = p1 + p2, width = 16, height = 12, units = "in")
        qsave(ds_seurat, paste0("CNV_files/", "cnvtree_seurat.qs"))
}

# if (length(system.file(package = "infercnv")) == 0) {
#         if (!require("BiocManager", quietly = TRUE))
#                 install.packages("BiocManager")
#         BiocManager::install("infercnv")
# }

library(infercnv)
library(qs)
# library(Seurat)
user_id <- args[4]
celltype_options <- args[1]
# param <- "group"
ref_group <- args[2]
max_cells <- as.numeric(args[3])
ds_flags <- TRUE
gene_order <- "../../hg38_gene_pos.txt"
print(getwd())
setwd(file.path(user_id))
u_dir <- getwd()
library(Seurat)

system.time({
        aoi_seurat <- qread("data/aoi_seurat.qs")
        aoi_seurat <- subset(aoi_seurat, celltype %in% celltype_options)
        cat("By default we downsample max cells per group to 100, time costs more as max cells increase \n")

        Idents(aoi_seurat) <- aoi_seurat$groups
        aoi_seurat <- subset(aoi_seurat, downsample = max_cells)
        aoi_seurat$idx <- colnames(aoi_seurat)

        output_dir <- "inferCNV/"
        suppressWarnings(dir.create(output_dir, recursive = TRUE))


        cnv_group_path <- "17_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.cell_groupings"

        run_infercnv(aoi_seurat = aoi_seurat,
                        max_cells = max_cells,
                        ref_group = ref_group,
                        gene_order = gene_order,
                        output_dir = output_dir)
        cat("Project CNV tree\n")
        inputs_cellgrouings <- cluster_cnv_trees(output_dir = output_dir, cnv_group_path = cnv_group_path)
        # print(head(inputs_cellgrouings))
        cat("Project CNV tree prediction on UMAP \n")
        project_cnvtree(aoi_seurat = aoi_seurat, inputs_cellgrouings = inputs_cellgrouings)
        library(reticulate)
        setwd(u_dir)
        tmp <- read.csv("inferCNV/CNV_files/cnvtree_celltype_groups_table.csv", header = FALSE)
        colnames(tmp) <- c("tree_groups","value","tree_ID")
        write.csv(tmp, "inferCNV/CNV_files/cnvtree_celltype_groups_table.csv")
        source_python("../../csv2json.py")
        path <- 'inferCNV/CNV_files/cnvtree_celltype_groups_table.csv'
        bulk_csv2json(path)
})



