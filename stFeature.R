# source ~/conda/bin/activate giotto
# Rscript stFeature.R img_id goi.txt user_id
args <- commandArgs(trailingOnly = TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) != 3) {
  stop("Three arguments must be supplied \n", call. = FALSE)
}
# img_id <- "LUAD;whole_lung;GSM5702473"
img_id <- args[1]
user_id <- args[3]
output_dir <- "stFeature/"
suppressWarnings(dir.create(output_dir, recursive = TRUE))

library(ggplot2)
library(ggpubr)
library(patchwork)
`%!in%` <- Negate(`%in%`)

print("=====================================")
print("=====         Spatial Gene Explorer         =====")
print("=====================================")

my_theme <- theme(axis.title = element_text(size = 10),
                axis.text = element_text(size = 16, angle = 45),
                legend.key.size = unit(1, "cm"),
                legend.text = element_text(size = 8),
                plot.subtitle = element_text(size = 14),
                plot.title = element_text(size = 18, hjust = 0.5))


# goi <- c("FOXP3", "NAT10", "SAMD11", "NOC2L")

# should be one single giotto object
st_Features <- function(giotto_samples, img_id, giotto_meta, goi, output_dir = "stFeature/") {
    # locate the specific metadata information
    cur_meta <- giotto_meta[giotto_meta$img_name %in% img_id, ]
    cur_disease <- unique(cur_meta$disease)
    cur_tissue <- unique(cur_meta$tissue)
    cat( paste0(cur_tissue, " in condition: ", cur_disease, "\n"))
    for (j in 1:length(goi)) {
        # replace 'your_gene' with the name of your gene
        if (goi[j] %!in% giotto_samples@feat_ID[[1]]) {
            # Create a blank plot
            png(filename = paste0(output_dir, goi[j], "_stFeatureplot.png"))
            plot(1, 1, type="n", ann=FALSE, axes=FALSE)
            # Add your text
            text(x = 1, y = 1, paste0( "Current gene:", goi[j], "cannot be found in this giotto object\n"))
            # Close the PNG device
            dev.off()
            next
        }
        p <- Giotto::spatFeatPlot2D(giotto_samples,
                    # spat_unit = "1",
                    expression_values = 'normalized',
                    cell_color_gradient = c("blue","white","red"),
                    gradient_limits = c(0, 2.5),
                    gradient_midpoint = 0.5,
                    feats = goi[j],
                    cow_rel_h = 0.5,
                    cow_rel_w = 2,
                    # show_image = TRUE,
                    show_plot = FALSE,
                    save_plot = FALSE,
                    return_plot = TRUE,
                    cow_n_col = 1,
                    point_size = 1
                    ) +
                    plot_annotation(title = paste0(cur_tissue,
                                    " in condition: ",
                                    cur_disease),
                        subtitle = unique(cur_meta$orig.ident)
                    ) +
                    my_theme
        ggsave(filename = paste0(output_dir, goi[j], "_stFeatureplot.png"),
                plot = p, width = 5, height = 4)
    }
    cat("All genes of interest have been plotted\n")
}


all_obj <- list.files("/mnt/root/lungDB_backend/global/spatial/data", pattern = ".qs", full.names = TRUE)
img_id <- gsub(".*\\;", "", img_id)
giotto_samples <- qs::qread(grep(img_id, all_obj, value = TRUE))
giotto_meta <- read.csv("/mnt/root/lungDB_backend/global/spatial/data/global_27imgs_seurat_0711_meta.csv")

setwd(file.path(user_id))
goi <- suppressWarnings(as.character(read.table(args[2], header = FALSE, sep = ";")))
st_Features(giotto_samples, img_id, giotto_meta, goi, output_dir = "stFeature/")
