mostexp_DEG_table <- function(deg_path,
                output_dir = output_dir,
                ct = "Macrophages", aoi_groups = "IPF") {
    all_f <- list.files(deg_path,
        pattern = "_DEG.csv", full.names = TRUE)
    f_ct <- gsub(".*sig_|_DEG.*|_mostexp.*", "", all_f)
    deg_f <- all_f[f_ct %in% ct]
    deg_f <- read.csv(paste0(deg_path, "sig_groupwise_DEG_table.csv"))
    deg_table <- subset(deg_f, celltype %in% ct)
    go_enrichment_goi(deg_table = deg_table, ct = ct,
                output_dir = output_dir, aoi_groups = "IPF")
    return(deg_table)
}

enrichment_go_KEGG <- function(enrich, ct = ct,
                output_dir = output_dir, up_down = "up") {
        groupname <- paste0(ct, "_", up_down)
        cat("Running", up_down, "GO enrichment \n")
        # GO Enrichment
        enrich_go <- enrichGO(gene = enrich$ENTREZID,
                        OrgDb = "org.Hs.eg.db",
                        ont = "ALL",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 0.05, # p value threshold
                        qvalueCutoff = 1,
                        readable = TRUE)
        # print(head(enrich_go))
        if (nrow(enrich_go) == 0) {
                cat("Cannot find significant GO enrichment")
        } else {
                write.csv(as.data.frame(enrich_go),
                paste0(output_dir,
                groupname,
                "DEG_GO.csv"))
        }
        p <- barplot(enrich_go, showCategory = 40,
                main = paste0("EnrichmentGO_",
                        up_down, "DEG", ct))
        p_title <- paste0(output_dir,
                groupname,
                "DEG_GO.pdf")
        ggsave(filename = p_title,
                plot = p, width = 6, height = 14, units = "in")
        # KEGG
        cat("Running", up_down, "KEGG enrichment \n")
        R.utils::setOption("clusterProfiler.download.method", "auto")
        # KEGG
        enrich_kegg <- enrichKEGG(gene = enrich$ENTREZID,
                        organism = "hsa",
                        pvalueCutoff = 0.05)
        if (nrow(as.data.frame(enrich_kegg)) < 10) {
                enrich_kegg <- enrichKEGG(gene = enrich$ENTREZID,
                        organism = "hsa",
                        pvalueCutoff = 0.05)
        }
        # print(head(enrich_kegg))
        write.csv(as.data.frame(enrich_kegg),
                paste0(output_dir, groupname,
                "DEG_KEGG.csv"))
        p <- barplot(enrich_kegg, showCategory = 40,
                title = paste0("EnrichmentKEGG_",
                        up_down, "DEG_", ct))
        p_title <- paste0(output_dir,
                groupname, "DEG_KEGG.pdf")
        ggsave(filename = p_title,
                plot = p, width = 6, height = 14, units = "in")
}

go_enrichment_goi <- function(deg_table, ct,
        # kegg_go = "GO",
        output_dir = "Enrichment/", aoi_groups = "IPF") {
        library(clusterProfiler)
        library(dplyr)
        library(ggplot2)
        groupname <- "disease"
        ribo_mt <- "^MT-|^RPL|^RPS|^ATP"
        deg_table <- deg_table[grep(ribo_mt, deg_table$gene, invert = TRUE), ]
        if (ct != "all") {
                deg_table <- subset(deg_table, celltype %in% ct)
        }
        cat("Calculate UP DEG enrichment of Cluster:",
                ct, "Group:", aoi_groups, "\n")
        if (length(grep(aoi_groups, deg_table$cluster)) == 0) {
            pdf(paste0(output_dir, groupname, "_warning.pdf"))
            cat("Groups of Interest not found in current DEG table \n")
            dev.off()
            stop("Groups of Interest not found in current DEG table \n")
        }
        up_deg <- deg_table[deg_table$avg_log2FC > 0
                & deg_table$cluster == aoi_groups, ]
        up_goi <- unique(up_deg$gene)
        up_enrich <- bitr(up_goi, fromType = "SYMBOL",
                        toType = c("ENSEMBL", "ENTREZID"),
                        OrgDb = "org.Hs.eg.db")
        enrichment_go_KEGG(enrich = up_enrich,
                ct = ct,
                output_dir = output_dir, up_down = "up")

        # down_deg <- deg_table[deg_table$avg_log2FC < 0 &
        #         deg_table$cluster == aoi_groups, ]
        # down_goi <- unique(down_deg$gene)
        # down_enrich <- bitr(down_goi, fromType = "SYMBOL",
        #                 toType = c("ENSEMBL", "ENTREZID"),
        #                 OrgDb = "org.Hs.eg.db")
        # enrichment_go_KEGG(enrich = down_enrich,
        #         ct = ct,
        #         output_dir = output_dir, up_down = "down")
}

GSEA_allplots <- function(deg_table, celltype_options, output_dir = "GSEA/", kegg_go) {
        library(msigdbr)
        library(DOSE)
        library(clusterProfiler)
        library(enrichplot)
        library(ggplot2)
        library(forcats)
        library(dplyr)
        library(stringr)
        # library(org.Hs.eg.db)
        deg_table <- deg_filter(deg_table)
        deg_table <- deg_table %>% group_by(cluster) %>%
                slice_min(p_val_adj, n = 1000)
        deg_table[deg_table$cluster == unique(deg_table$cluster)[2], ]$avg_log2FC <- -deg_table[deg_table$cluster == unique(deg_table$cluster)[2], ]$avg_log2FC 
        ct_deg <- subset(deg_table, celltype == celltype_options)
        gene_list <- ct_deg$avg_log2FC
        # name the vector
        names(gene_list) <- ct_deg$gene
        # omit any NA values
        gene_list <- na.omit(gene_list)
        # sort the list in decreasing order (required for clusterProfiler)
        gene_list <- sort(gene_list, decreasing = TRUE)
        # print(gene_list)
        # print(gene_list)
        if (kegg_go == "KEGG") {
                # KEGG
                genesets <- msigdbr(species = "Homo sapiens", category = "C2")
                # genesetsuse <- subset(genesets, gs_subcat %in%  "CP:KEGG", select = c("gs_name", "gene_symbol"))
                genesetsuse <- subset(genesets, gs_subcat != "HPO", select = c("gs_name", "gene_symbol"))
                fun_key <- "enrichKEGG"
        } else {
                # GO
                c5_genesets = msigdbr(species = "Homo sapiens", category = "C5")
                # genesets <- rbind(c5_genesets, h_genesets, c2_genesets)
                genesetsuse <- subset(c5_genesets, gs_subcat != "HPO", select = c("gs_name", "gene_symbol"))
                fun_key <- "enrichGO"
        }

        groupname <- paste0(celltype_options, "_", kegg_go)

        res <- GSEA(gene_list, TERM2GENE = genesetsuse, minGSSize = 5,
                        pvalueCutoff = 1, nPermSimple = 5000, verbose = F)
        if (nrow(res) == 0) {
                cat("No significant GSEA terms\n")
                return(NULL)
        }
 
        # dotplot
        p <- enrichplot::dotplot(res, showCategory = 10) + facet_grid(~.sign) +
                                theme(axis.text = element_text(size = 14, angle = 0),
                        #legend.key.size = unit(1.5, "cm"),
                        #legend.text = element_text(size = 5),
                        plot.title = element_text(size = 20, hjust = 0.5),
                        panel.grid.major = element_line(colour = "grey"),
                        panel.grid.minor = element_line(colour = "grey"),
                        panel.background = element_rect(fill = "white"),
                        plot.subtitle = element_text(size = 14))
        # p_title <- paste0(output_dir,
        #                 groupname,
        #                 "Enrich_Dotplot.pdf")
        # ggsave(filename = p_title,
        #         plot = p, width = 8, height = 8, units = "in")
        p_title <- paste0(output_dir,
                        groupname,
                        "_Enrich_Dotplot.png")
        ggsave(filename = p_title,
                plot = p, width = 12, height = 10, units = "in")


        # barplot
        df_res <- data.frame(res)
        df_res$NES_sign <- ifelse(df_res$NES > 0, "Positive", "Negative")
        write.csv(df_res, file = paste0(output_dir,
                        groupname,
                        "GSEA.csv"), row.names = TRUE)
        df_res_filtered <- df_res %>%
                group_by(NES_sign) %>%
                top_n(5, abs(NES))
        library(stringr)

        p <- ggplot(df_res_filtered,
                aes(NES, fct_reorder(Description, NES),
                fill=-log10(p.adjust))) + 
                geom_bar(stat='identity') + 
                scale_fill_continuous(low='red', high='blue', guide=guide_colorbar(reverse=TRUE)) + 
                theme_minimal() + ylab(NULL) +
                scale_y_discrete(labels=function(x) str_wrap(str_replace_all(x,"_"," "), width=30)) +
                theme(
                        axis.text = element_text(size = 14, angle = 0),
                        plot.title = element_text(size = 20, hjust = 0.5),
                        panel.grid.major = element_line(colour = "grey"),
                        panel.grid.minor = element_line(colour = "grey"),
                        panel.background = element_rect(fill = "white"),
                        plot.subtitle = element_text(size = 14)
                )

        p_title <- paste0(output_dir,
                        groupname,
                        "_Enrich_Barplot.png")
        ggsave(filename = p_title,
                plot = p, width = 12, height = 10, units = "in")
        # head(res)
        heat_res <- GSEA(gene_list[1:50], TERM2GENE = genesetsuse, minGSSize = 5,
                        pvalueCutoff = 1, nPermSimple = 10000, verbose = F)
        #heatplot
        p <- heatplot(heat_res, showCategory = 10, foldChange = gene_list[1:50]) +
                theme(axis.text = element_text(size = 14, angle = 0),
                        legend.key.size = unit(1.5, "cm"),
                        legend.text = element_text(size = 5),
                        plot.title = element_text(size = 20, hjust = 0.5),
                        panel.grid.major = element_line(colour = "grey"),
                        panel.grid.minor = element_line(colour = "grey"),
                        panel.background = element_rect(fill = "white"),
                        plot.subtitle = element_text(size = 14))
        p_title <- paste0(output_dir,
                        groupname,
                        "_Enrich_heatplot.png")
        ggsave(filename = p_title,
                plot = p, width = 12, height = 10, units = "in")
        # Emapplot
        pw_res <- pairwise_termsim(res)
        p <- enrichplot::emapplot(pw_res, showCategory = 30,
                        cex.params = list(category_node = 1))
        p_title <- paste0(output_dir,
                        groupname,
                        "_Enrich_Emaplot.png")
        ggsave(filename = p_title,
                plot = p, width = 12, height = 10, units = "in")

        return(res)
}
