#!/usr/bin/env Rscript

# Rscript fraction_plot.R seurat_path output_dir gene_list
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 3) {
  stop("Three arguments must be supplied \n", call. = FALSE)
}
seurat_path <- args[1]
output_dir <- args[2]
goi <- read.csv(args[3])
aoi_seurat <- readRDS(seurat_path)

refine_metanames <- function(meta) {
        colnames(meta)[colnames(meta) == "orig.ident"] <- "accession"
        colnames(meta)[colnames(meta) == "gse_alias"] <- "experiment"
        colnames(meta)[colnames(meta) == "DISEASE"] <- "disease"
        colnames(meta)[colnames(meta) == "TISSUE"] <- "tissue"
        colnames(meta)[colnames(meta) == "celltype"] <- "detail_celltype"
        colnames(meta)[colnames(meta) == "macro_celltype"] <- "celltype"
        return(meta)
}

# Fraction_plot necessary
fraction_param_table <- function(seurat,
                                threshold = 0.5, param = param,
                                goi = goi) {
        # preprocess meta data names
        seurat$param <- seurat[[param]][, 1]
        meta <- seurat@meta.data
        if ("orig.ident" %in% colnames(meta)) {
                meta <- refine_metanames(meta)
        }

        meta_param <- c("disease", "tissue", "aging", "gender",
                        "experiment", "celltype", "groups", "accession", "param")
        frac_meta <- as.data.frame(meta[, meta_param])
        frac_meta <- distinct(frac_meta, accession, celltype, .keep_all= TRUE)
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
                          floor = 50
                          ) {
    frac_table <- frac_table[frac_table$cell_counts >= floor &
        frac_table$expressed_ratio > 0, ]
    # frac_table$ct_disease <- paste0(frac_table$celltype,
    #               "+", frac_table$disease)
    ps <- list()
    for (i in 1:length(goi)) {
        cur_frac <- frac_table[frac_table$gene == goi[i], ]
        if (nrow(cur_frac) == 0) {
          cat("Current gene:",
            goi[i],
            "is not found in Fraction Table \n")
          next
        }
        cur_frac$disease_ct <- paste0(cur_frac$celltype, "_", cur_frac$disease)
        sample_table <- table(cur_frac$disease)
        p <- ggplot(cur_frac,
                    aes(x = disease_ct, y = expressed_ratio,
                    fill = disease)) +
                    geom_boxplot() +
                    # geom_violin(trim = FALSE) +
                    geom_point(position = position_jitter(
                          seed = 1, width = 0.2)) +
                    # geom_boxplot(width = 0.05, fill = "white") +
                    # facet_wrap(~gene) +
                    labs(x = "Genes of Interest",
                        y = "Fraction of Cells Expressing Target Genes",
                        size = 100) +
                theme(axis.text = element_text(size = 12, angle = 45),
                        legend.key.size = unit(1.5, "cm"),
                        legend.text = element_text(size = 10),
                        plot.title = element_text(size = 20, hjust = 0.5),
                        plot.subtitle = element_text(size = 14))
        p <- p + labs(title = paste0("Fraction plot of GOI: ",
                                    goi[i]),
                      subtitle = paste0(
                                  "Normal Sample Number: ",
                                  sample_table["Normal"],
                                  "\nDisease Sample Number: ",
                                  sample_table[names(sample_table) != "Normal"],
                                  # length(unique(frac_table$cell_counts)),
                                  "\nTotal Cell Counts: ",
                                  sum(unique(frac_table$cell_counts)),
                                  "\nPass Cutoff Cell Counts: ",
                                  round(sum(cur_frac[cur_frac$gene
                                    == goi[i], ]$cell_counts *
                                    cur_frac[cur_frac$gene
                                    == goi[i], ]$expressed_ratio), digits = 0)
                                  ,
                                  "\nCutoff: ", threshold,
                                  "\nCell Number Floor: ", floor
                                  ),
                     # tag = paste0("Normal Sample Numbe: ")
                    )
        ps[[i]] <- p
    }
    return(ps)
    cat("Successfully Plot Fraction Plot\n")
}

frac_table_path <- paste0(output_dir, paste(unlist(goi), collapse = "_"),
            "_fractionTable.csv")
if (file.exists(frac_table_path)) {
    cat("Previous Fraction Table Calculation has already been stored \n")
    frac_table <- read.csv(frac_table_path)
} else {
    frac_table <- fraction_param_table(seurat = aoi_seurat,
                                threshold = 0.5,
                                param = "celltype", #will not change
                                goi = goi)
    write.csv(frac_table, frac_table_path)
}
ps <- fraction_plot(frac_table = frac_table, goi = goi)
for (i in 1:length(ps)) {
    p_title <- paste0(output_dir,
            goi[i], "_fractionPlot.pdf")
    ggsave(filename = p_title,
            plot = ps[[i]], width = 10, height = 12, units = "in")
}