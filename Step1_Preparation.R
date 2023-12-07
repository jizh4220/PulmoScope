#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 2) {
  stop("Two arguments must be supplied \n", call. = FALSE)
}

cat("############ Current Module is Integration Explorer ############\n")
library(data.table)

`%!in%` <- Negate(`%in%`)

library(future)
plan()
plan("multicore", workers = 8)
options(future.globals.maxSize = 51539607552)

aoi_seurat_preparation <- function(all_folders, sample_info, threshold = 10000) {
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
    suppressWarnings(dir.create("data", recursive = TRUE))
    seurat <- merge(ps[[1]], ps[-1])
    cat("Seurat Integration time cost: \n",
        system.time({seurat <- harmony_merge(seurat)
                    seurat <- seurat %>%
                            FindNeighbors(reduction = "harmony", dims = 1:20) %>%
                            RunUMAP(reduction = "harmony", dims = 1:20)
                    # clean meta data
                    meta <- seurat@meta.data
                    meta <- meta_alignment(meta)
                    meta_path <- "data/meta_Info.csv"
                    fwrite(meta, meta_path)
                    seurat@meta.data <- meta
        }), "\n")
    
    saveRDS(seurat, "data/aoi_seurat.rds")
    write.table(rownames(seurat),
            "data/gene_names.txt",
            row.names = FALSE, col.names = FALSE)
    write.table(colnames(seurat),
            "data/cell_idx.txt",
            row.names = FALSE, col.names = FALSE)
    cat("Successfully generate AOI seurat \n")
    cat("############Generating Basic Plots############\n")
    return(seurat)
}


save_embeddings <- function(seurat, dir_path) {
    library(data.table)
    reduct <- names(seurat@reductions)
    for (i in 1:length(reduct)) {
        embed_path <- file.path(dir_path, paste0(reduct[i],
                    "_Embeddings.csv"))
        embed <- Embeddings(seurat, reduct[i])
        fwrite(embed, embed_path)
    }
    cat("All embeddings have been successfully prepared")
}

save_seurat_data <- function(seurat, output_dir) {
    library(data.table)
    meta <- seurat@meta.data
    fwrite(meta, file.path(output_dir, "metadata.csv"), row.names = TRUE)
    gene_name <- rownames(seurat)
    writeLines(gene_name, file.path(output_dir, "gene_names.txt"))
    cell_idx <- colnames(seurat)
    writeLines(cell_idx, file.path(output_dir, "cell_index.txt"))

    ct_table <- table(seurat$celltype)
    sig_ct <- names(ct_table[ct_table >= 100])
    writeLines(sig_ct, file.path(output_dir, "celltype_options.txt"))
}

prepare_BPCounts <- function(cur_MG_path, located_meta) {
    if (length(cur_MG_path) == 1) {
        counts <- open_matrix_dir(cur_MG_path)
        counts <- counts[, intersect(located_meta$idx, colnames(counts))]
        return(as(counts, "dgCMatrix"))
    }
    rn <- list()
    ps <- list()
    # TODO: Finish different MacroGroup situations
    for (i in 1:length(cur_MG_path)) {
        # now we will merge different counts in real-time
        counts <- open_matrix_dir(cur_MG_path[i])
        rn[[i]] <- rownames(counts)
        ps[[i]] <- counts[, intersect(located_meta$idx, colnames(counts))]
    }
    cat("Successfully fetch all counts with cell index.\n")
    shared_genes <- Reduce(intersect, rn)
    for (i in 1:length(cur_MG_path)) {
        # now we will merge different counts in real-time
        counts <- ps[[i]]
        ps[[i]] <- as(counts[shared_genes, ], "dgCMatrix")
    }
    counts <- do.call(cbind, ps)
    # counts <- BPCounts_routine(counts)
    return(counts)
}


locate_sample_info_global <- function(sample_info, global_meta, threshold = threshold) {
    if (length(grep("[A-Z]", rownames(global_meta))) == 0) {
        stop("Please double check your row name of global_meta. \n")
    }
    aoi_acc <- unique(sample_info[, "accession"])
    aoi_group <- sample_info[, "groups"]
    acc_meta <- global_meta[global_meta$accession %in% aoi_acc, ]
    for (i in 1:length(aoi_acc)) {
        acc_meta$groups[acc_meta$accession == aoi_acc[i]] <- aoi_group[i]
    }
    cat("Exclude OOD celltype\n")
    acc_meta <- acc_meta[!duplicated(acc_meta$idx), ]
    df <- table(acc_meta$groups, acc_meta$celltype)
    ct_fix_df <- df[, !apply(df, 2, function(x) any(x < 100))]
    # Find columns that contain values smaller than the threshold
    acc_meta <- subset(acc_meta, celltype %in% colnames(ct_fix_df))
    # i <- 1
    cat("Downsampling based on fixed cell type distribution\n")
    # Calculate proportions
    proportions <- prop.table(ct_fix_df)
    # Apply square root function to proportions
    smoothed_proportions <- sqrt(proportions)
    # Normalize smoothed proportions so they sum up to 1
    smoothed_proportions <- smoothed_proportions / sum(smoothed_proportions)
    # Multiply normalized smoothed proportions by desired total sum
    sampling_size_df <- round(smoothed_proportions * threshold)
    # Downsample the acc_meta data frame for each disease type
    sampled_acc_meta <- acc_meta %>%
            group_by(groups, celltype) %>%
            group_modify(~ slice_sample(.x, n = round(sampling_size_df[unique(.y$groups), unique(.y$celltype)]), replace = FALSE))
    # print(head(sampled_acc_meta))
    # locate specific cell idx
    located_meta <- as.data.frame(sampled_acc_meta)
    df <- table(located_meta$groups, located_meta$celltype)
    ct_propo <- df[1, ] / df[2, ]
    for (i in 1:ncol(df)) {
        if (ct_propo[i] > 10 || ct_propo[i] < 0.1) {
            cat("Reset outlier proportion to maximum 5 folds\n")
            df[which.max(df[, i]), i] <- min(df[, i]) * 5
        }
    }
    located_meta <- located_meta %>%
            group_by(groups, celltype) %>%
            group_modify(~ slice_sample(.x, n = round(df[unique(.y$groups), unique(.y$celltype)]), replace = FALSE))
    located_meta <- as.data.frame(located_meta)
    print(table(located_meta$groups, located_meta$celltype))
    return(located_meta)
}



locate_BPCells_seurat <- function(MG_folders, sample_info, global_meta, threshold = 20000) {
    library(stringr)
    library(BPCells)
    library(dplyr)
    dir.create("data/")
    located_meta <- locate_sample_info_global(sample_info = sample_info,
            global_meta = global_meta, threshold = threshold)
    table(located_meta$groups, located_meta$celltype)
    fwrite(located_meta, "data/metadata.csv")
    # go back to Seurat objects
    mg <- unique(located_meta$MacroGroup)
    # print(mg)
    mg_idx <- str_extract(MG_folders, "MacroGroup[0-9]+")
    cur_MG_path <- MG_folders[mg_idx %in% mg]
    cat("Fetch from MG:", mg, "\n")
    counts <- prepare_BPCounts(cur_MG_path, located_meta = located_meta)
    located_meta <- located_meta[match(colnames(counts), located_meta$idx), ]
    located_meta <- located_meta[!duplicated(located_meta$idx), ]
    rownames(located_meta) <- located_meta$idx
    located_meta <- located_meta %>%
                replace(. == "" | is.na(.), "unknown") 
    counts <- counts[, rownames(located_meta)]

    seurat <- CreateSeuratObject(counts, meta.data = located_meta, min.cells = 10,
                            min.features = 200, names.delim = "_")
    # print(head(seurat))
    seurat$orig.ident <- seurat$accession
    # seurat$accession <- NULL
    seurat <- harmony_merge(seurat)
    save_seurat_data(seurat = seurat, output_dir = "data/")
    qsave(seurat, file.path("data/aoi_seurat.qs"))
    return(seurat)
}

# Celltype Composition
Compo_Meta_plot <- function(
              aoi_seurat,
              param = "groups",
              output_dir = output_dir
              ) {
    aoi_meta <- aoi_seurat@meta.data
    group_table <- table(aoi_seurat$groups)
    group_name <- names(group_table)
    for (i in 1:length(param)) {
      cur_param <- param[i]
      if (all(is.na(aoi_seurat[[cur_param]][1]))) {
        next
      }
      # print(table(aoi_meta[, cur_param], aoi_meta[, "accession"]))
      # Acctable
      if (cur_param != "accession") {
        acc_param <- table(aoi_meta[, cur_param], aoi_meta[, "accession"])
        acc_param <- as.data.frame(acc_param)
        colnames(acc_param) <- c(cur_param, "accession", "counts")
        write.csv(acc_param, paste0(output_dir, cur_param, 
                "_acctable.csv"), row.names = TRUE)
      }
      cat("Plotting Dimplot and CompoPlot by Param:", cur_param, "\n")
      
      if (cur_param == "accession") {
        # Dimplot
        p <- DimPlot(aoi_seurat,
              group.by = cur_param,
              reduction = "umap",
              pt.size = 0.2, repel = FALSE,
              label = TRUE, label.size = 4) +
              # theme(legend.key.size = unit(0, "cm"),
              #       legend.text = element_text(size = 0))
              NoLegend()
        p_title <- paste0(output_dir, "meta_", cur_param, "_DimPlot.png")
        ggsave(filename = p_title,
            plot = p, width = 12, height = 9, units = "in")

        param_comp <- table(aoi_meta[, cur_param], aoi_meta[, "celltype"])
        param_comp <- as.data.frame(param_comp)
        colnames(param_comp) <- c("param", "celltype", "counts")

        # Barplot
        p <- ggplot(param_comp,
                  aes(x = celltype,
                      y = counts,
                      fill = param)) +
          theme_classic(base_size = 15) +
          geom_col(position = "fill", width = 0.5) +
          labs(x = "Param",
              y = "Param Composition",
              title = "Param Composition Per Celltype",
              subtitle = paste0("Cell Number of ", group_name[1],
                      " is: ", group_table[1],
                      "\n",
                      "Cell Number of ", group_name[2],
                      " is: ", group_table[2])) +
          theme(axis.title = element_text(size = 15),
                axis.text = element_text(size = 15, angle = 45),
                legend.position = 'none',
                # legend.key.size = unit(1.2, "cm"),
                # legend.text = element_text(size = 10),
                plot.title = element_text(size = 20, hjust = 0.5))
        p_title <- paste0(output_dir, cur_param, 
            "_Compoplot.png")
        ggsave(filename = p_title,
              plot = p, width = 12, height = 8, units = "in")
        write.csv(param_comp, paste0(output_dir, cur_param, 
            "_CompoTable.csv"))
        next
      } else {
        #Dimplot
        p <- DimPlot(aoi_seurat,
                group.by = cur_param,
                reduction = "umap", pt.size = 0.2, repel = TRUE,
                label = TRUE, label.size = 8)
        p_title <- paste0(output_dir, "meta_", cur_param, "_DimPlot.png")
        ggsave(filename = p_title,
            plot = p, width = 12, height = 8, units = "in")
      }
      if (cur_param != "groups") {
        #Compo Plot
        param_comp <- table(aoi_meta[, cur_param], aoi_meta[, "groups"])
        param_comp <- as.data.frame(param_comp)
        colnames(param_comp) <- c("param", "groups", "counts")
        p <- ggplot(param_comp,
                  aes(x = groups,
                      y = counts,
                      fill = param)) +
          theme_classic(base_size = 15) +
          geom_col(position = "fill", width = 0.5) +
          labs(x = "Param",
              y = "Group Composition",
              title = "Param Composition Per Group",
              subtitle = paste0("Cell Number of ", group_name[1],
                      " is: ", group_table[1],
                      "\n",
                      "Cell Number of ", group_name[2],
                      " is: ", group_table[2])) +
              theme(axis.title = element_text(size = 15),
                    axis.text = element_text(size = 20),
                    axis.text.x = element_text(vjust = 0.5, angle = 45),
                    legend.key.size = unit(1.2, "cm"),
                    legend.text = element_text(size = 10),
                    plot.title = element_text(size = 20, hjust = 0.5))
        p_title <- paste0(output_dir, cur_param, 
            "_Compoplot.png")
        ggsave(filename = p_title,
              plot = p, width = 12, height = 8, units = "in")
        write.csv(param_comp, paste0(output_dir, cur_param, 
            "_CompoTable.csv"))
      } else {
        # Groups Compo Plot
        param_comp <- table(aoi_meta[, cur_param], aoi_meta[, "celltype"])
        param_comp <- as.data.frame(param_comp)
        colnames(param_comp) <- c("param", "celltype", "counts")
        p <- ggplot(param_comp,
                  aes(x = celltype,
                      y = counts,
                      fill = param)) +
          theme_classic(base_size = 15) +
          geom_col(position = "fill", width = 0.5) +
          labs(x = "Param",
              y = "Group Composition",
              title = "Group Composition Per Celltype",
              subtitle = paste0("Cell Number of ", group_name[1],
                      " is: ", group_table[1],
                      "\n",
                      "Cell Number of ", group_name[2],
                      " is: ", group_table[2])) +
          theme(axis.title = element_text(size = 15),
                    axis.text = element_text(size = 20),
                    axis.text.x = element_text(vjust = 0.5, angle = 45),
                    legend.key.size = unit(1.2, "cm"),
                    legend.text = element_text(size = 10),
                    plot.title = element_text(size = 20, hjust = 0.5))
        p_title <- paste0(output_dir, cur_param, 
            "_Compoplot.png")
        ggsave(filename = p_title,
              plot = p, width = 12, height = 8, units = "in")

        param_comp <- table(aoi_meta[, cur_param], aoi_meta[, "celltype"])
        param_comp <- t(param_comp)
        write.csv(param_comp, paste0(output_dir, cur_param, 
            "_CompoTable.csv"))
        compo_table <- read.csv(paste0(output_dir, cur_param, 
                    "_CompoTable.csv"))
        colnames(compo_table)[1] <- "celltype"
        write.csv(compo_table, paste0(output_dir, cur_param, 
                    "_CompoTable.csv"))
      }
    }
    cat("Finish Compo Meta Table and Plot per Accession, Experiment, and Group \n")
    return(TRUE)
}

Routine_Plot <- function(output_dir = "Step1/", aoi_seurat) {
    meta_param <- c("disease", "tissue", "aging", "gender",
                "experiment", "celltype", "groups", "accession")
    cat("Plotting Celltype Compo and Meta Plot", "\n")
    Compo_Meta_plot(output_dir = output_dir, aoi_seurat = aoi_seurat,
            param = meta_param)
    
    cat("Finish Calculating Top Marker and Plotting Dotplot and Heatmap \n")
}

fread_counts <- function(path, header = FALSE) {
    library(data.table)
    if (header == TRUE) {
        counts <- as.data.frame(fread(path, header = TRUE))
    } else {
        counts <- as.data.frame(fread(path))
    }
    return(counts)
}


if (!exists("aoi_seurat")) {
    print(getwd())
    library(data.table)

    MG_folders <- list.files("/mnt/root/lungDB_backend/global/MGDB",
                pattern = "val",
                recursive = TRUE, full.names = TRUE)
    MG_folders <- fs::path_dir(MG_folders)
    global_meta <- fread_counts("/mnt/root/lungDB_backend/global/global_meta_fourth_clean.csv")
    rownames(global_meta) <- global_meta$idx
    user_id <- args[2]
    setwd(file.path(user_id))
    sample_info <- suppressWarnings(read.table(args[1],
            header = TRUE))
    output_dir <- "Step1/"
    
    suppressWarnings(dir.create(output_dir, recursive = TRUE))

    # Regular R Integration
    cat("Overall Seurat Integration Time Cost is:", system.time({
        if (!file.exists("data/aoi_seurat.qs")) {
            cat("Locate Samples of Interest from MacroGroup\n")
            system.time({aoi_seurat <- locate_BPCells_seurat(MG_folders = MG_folders,
                        global_meta = global_meta,
                        sample_info = sample_info, threshold = 20000)})
        } else{
            aoi_seurat <- qread("data/aoi_seurat.qs")
        }
        cat("Succesfully generate aoi seurat.\n")
        Routine_Plot(output_dir = "Step1/", aoi_seurat)
        print(getwd())
        library(reticulate)
        source_python("../../csv2json.py")
        path <- 'Step1/*.csv'
        bulk_csv2json(path)
        cat("Successfully finish Step1 preparation \n")
    }), "\n")
    
}