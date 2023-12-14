# source ~/conda/bin/activate BPCells
# Rscript global_Feature.R goi.txt user_id
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 2) {
  stop("Two arguments must be supplied \n", call. = FALSE)
}
# img_id <- "LUAD;whole_lung;GSM5702473"

user_id <- args[2]
library(ggplot2)
library(ggdendro)
library(cowplot)
library(ggtree)
library(aplot)
library(tidyverse)
library(patchwork)
library(Seurat)
library(BPCells)
`%!in%` <- Negate(`%in%`)

print("=====================================")
print("=====         Global Feature Explorer         =====")
print("=====================================")

setwd(file.path(user_id))
goi <- suppressWarnings(as.character(read.table(args[1], header = FALSE, sep = ";")))

output_dir <- "global_Feature/"
suppressWarnings(dir.create(output_dir, recursive = TRUE))

# load necessary libraries

# goi <- c("FOXP3", "NAT10", "SAMD11", "NOC2L")

# should be one single giotto object
global_Feature <- function(seurat, goi, output_dir = "global_Feature/") {
    data.features <- FetchData(object = seurat, vars = goi)
    meta <- seurat@meta.data
    data.features$id <- paste(meta$groups, meta$celltype, sep = ";")
    
    for (i in 1:length(goi)) {
        result <- lapply(X = unique(x = data.features$id), FUN = function(y) {
                data.use <- data.features[data.features$id == y, i,
                                        drop = FALSE]
                avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
                        return(mean(x = expm1(x = x)))
                })
                pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                        threshold = 0)
                result <- cbind(avg.exp = log2(avg.exp), pct.exp = pct.exp)
                rownames(result) <- y

                return(result)
        })
        result <- do.call(rbind, result)
        gene_cluster <- as.data.frame(result)
        gene_cluster$groups <- gsub(";.*", "", rownames(gene_cluster))
        gene_cluster$celltype <- gsub(".*;", "", rownames(gene_cluster))
        gene_cluster$avg.exp <- as.numeric(gene_cluster$avg.exp)
        gene_cluster$pct.exp <- as.numeric(gene_cluster$pct.exp)
        # ps[[i]] <- gene_cluster
        dotplot <- gene_cluster %>% 
                mutate(`% Expressing` = pct.exp * 100) %>% 
                filter(avg.exp != 0, `% Expressing` > 0) %>% 
                ggplot(aes(x = groups, y = celltype, color = avg.exp, size = `% Expressing`)) + 
                geom_point() + 
                cowplot::theme_cowplot() + 
                # theme(axis.line  = element_blank()) +
                theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5),
                        panel.grid.major = element_line(colour = "grey"),
                        panel.grid.minor = element_line(colour = "grey"),
                        panel.background = element_rect(fill = "white"),) +
                xlab('') +
                ylab('') +
                # theme(axis.ticks = element_blank()) +
                scale_color_gradientn(colours = viridis::viridis(20), limits = c(0, 5), oob = scales::squish, name = 'log2 (avgExp)')
                # scale_y_discrete(position = "right")

        mat <- gene_cluster %>% 
                select(-pct.exp) %>%  # drop unused columns to faciliate widening
                pivot_wider(names_from = groups, values_from = avg.exp) %>% 
                mutate_all(~ifelse(is.infinite(.), 500, .)) %>%  # replace Inf with largest finite number
                data.frame()
        # mat[is.infinite(mat)] <- .Machine$double.xmax
        mat_dist <- dist((mat[, -1]  %>% as.matrix() %>% t()))
        # mat_dist[is.infinite(mat_dist)] <- ""
        groups_clust <- hclust(mat_dist)
        ddgram_col <- as.dendrogram(groups_clust)
        ggtree_plot_col <- ggtree(ddgram_col) +
                        # coord_flip() +
                        layout_dendrogram() +
                        # coord_flip(theta='y') +
                        # theme_tree2() +
                        theme(legend.position='none') +
                        ylim2(dotplot) + 
                        scale_x_reverse()
                        # scale_y_continuous()
        # revts(p)
        # p <- ggtree_plot_col %>% insert_bottom(dotplot)
        p <-    ggtree_plot_col +
                # coord_flip() +
                dotplot +
                plot_layout(ncol = 1, heights = c(6, 15, 0, 0, 0)) +
                plot_annotation(title = paste0("Global Feature Dot Plot: ",
                                        goi[i]),
                                theme = theme(plot.title = element_text(size = 20, hjust = 0.5)))
        # p <- plot_grid(dotplot, nrow = 1, rel_widths = c(0.5,2), align = 'h')
        p_title <- paste0(output_dir, goi[i], "_global_feature.png")
        ggsave(filename = p_title,
                plot = p, width = 12, height = 8, units = "in")
        # break
    }
    cat("All genes of interest have been plotted\n")
}

aoi_seurat <- qread("/mnt/root/lungDB_backend/global/global_smaller_sampler_harmony_seurat.qs")
global_Feature(aoi_seurat, goi, output_dir = "global_Feature/")
