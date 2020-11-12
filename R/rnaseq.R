
# ----------
# PART 1 - Plotting Functions
# ----------


#' PCA plot from vsd object
#'
#' @param vsd A vsd object
#' @param var A string - the name of the variable matching metadata columns
#' @param pal A string - palette name of \code{RColorBrewer}
#'
#' @return A ggplot2 plot
#' @importFrom ggplot2 ggplot geom_point aes scale_color_brewer labs theme_bw theme
#' @export

plot_pca_vsd <- function(vsd, var, pal) {

    pca <- DESeq2::plotPCA(vsd, intgroup = var, returnData = TRUE)

    ggplot(pca) +
        geom_point(aes(x = .data$PC1,
                       y = .data$PC2,
                       color = .data$group), size = 2) +
        scale_color_brewer(palette = pal) +
        labs(color = var) +
        theme_bw() +
        theme(aspect.ratio = 1)
}

# ----------

#' MA plot from DESeq2Results object
#'
#' @param res A DESeq2Results object
#' @param p_co A double - the cutoff of adjusted p-value
#' @param lfc_co A double - the cutoff of log2 fold change
#' @param lfc_plot_lim A double - the y-limit of log2 fold change plot. Default value is \code{5}.
#'
#' @return A ggplot2 plot
#' @importFrom dplyr arrange
#' @importFrom ggplot2 scale_color_manual scale_shape_manual element_text geom_hline
#' @importFrom grid unit
#' @export

plot_deseq_ma <- function(res, p_co, lfc_co, lfc_plot_lim = 5) {

    res <- res %>%
        res_to_tibble(p_co, lfc_co)

    res %>%
        mutate(shape = ifelse(.data$log2FoldChange > lfc_plot_lim | .data$log2FoldChange < -lfc_plot_lim, TRUE, FALSE),
               log2FoldChange = replace(.data$log2FoldChange, .data$log2FoldChange > lfc_plot_lim, lfc_plot_lim),
               log2FoldChange = replace(.data$log2FoldChange, .data$log2FoldChange < -lfc_plot_lim, -lfc_plot_lim)) %>%
        arrange(factor(.data$significant, levels = c("Not Sig","Down","Up"))) %>%
        ggplot() +
        geom_point(aes(x = log10(.data$baseMean),
                       y = .data$log2FoldChange,
                       color = .data$significant,
                       shape = .data$shape),
                   size = 2) +
        scale_color_manual(values = xmc_constants()$tricolor) +
        scale_shape_manual(values = c(16, 17)) +
        theme_bw() +
        labs(y = expression(Log[2]~Fold~Change), x = expression(Log[10]~Mean~Normalized~Count)) +
        theme(legend.position = "none",
              plot.margin = unit(rep(1,4), "cm"),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14)) +
        geom_hline(yintercept = 0, color = "#984ea3", size = 1.5, alpha = 0.5)

}

# ----------

#' Volcano plot from DESeq2Results object
#'
#' @param res A DESeq2Results object
#' @param p_co A double - the cutoff of adjusted p-value
#' @param lfc_co A double - the cutoff of log2 fold change
#' @param p_plot_lim A double - the y-limit of -log10 adjusted p-value. Default value is \code{5}.
#' @param lfc_plot_lim A double - the x-limit of log2 fold change plot. Default value is \code{5}.
#'
#' @return A ggplot2 plot
#' @export

plot_deseq_volcano <- function(res, p_co, lfc_co,
                          p_plot_lim = 5,
                          lfc_plot_lim = 5) {

    res <- res %>%
        res_to_tibble(p_co, lfc_co)

    res %>%
        mutate(shape1 = ifelse(.data$log2FoldChange > lfc_plot_lim | .data$log2FoldChange < -lfc_plot_lim, TRUE, FALSE),
               log2FoldChange = replace(.data$log2FoldChange, .data$log2FoldChange > lfc_plot_lim, lfc_plot_lim),
               log2FoldChange = replace(.data$log2FoldChange, .data$log2FoldChange < -lfc_plot_lim, -lfc_plot_lim)) %>%
        mutate(shape2 = ifelse(-log10(.data$padj) > p_plot_lim, TRUE, FALSE),
               padj = ifelse(-log10(.data$padj) > p_plot_lim, 10^(-p_plot_lim), .data$padj),
               shape = .data$shape1 | .data$shape2) %>%
        arrange(factor(.data$significant, levels = c("Not Sig","Down","Up"))) %>%
        ggplot() +
        geom_point(aes(x = .data$log2FoldChange,
                       y = -log10(.data$padj),
                       color = .data$significant,
                       shape = factor(.data$shape)),
                   size = 2) +
        scale_color_manual(values = xmc_constants()$tricolor) +
        scale_shape_manual(values = c(16, 17), guide = FALSE) +
        theme_bw() +
        labs(y = expression(-Log[10]~Adjusted~p-value), x = expression(Log[2]~Fold~Change)) +
        theme(aspect.ratio = 1,
              legend.position = "none",
              plot.margin = unit(rep(1,4), "cm"),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))

}

# ----------

#' Barplot from GSEA results
#'
#' @param gsea A tibble of GSEA results
#' @param pattern A string - the pattern to remove in the plot. Default value is \code{"HALLMARK_"}.
#'
#' @return A ggplot2 plot
#' @importFrom stringr str_remove
#' @importFrom stats reorder
#' @importFrom ggplot2 geom_bar scale_fill_gradient2 coord_flip
#' @export

plot_deseq_gsea <- function(gsea, pattern = "HALLMARK_") {

    color_set <- xmc_constants()$tricolor

    gsea %>%
        gsea_rm_pattern() %>%
        mutate(color = -log10(.data$padj) * ifelse(.data$padj <= 1, 1, 0) * ifelse(.data$NES > 0, 1, -1)) %>%
        ggplot() +
        geom_bar(aes(x = reorder(.data$pathway, .data$NES),
                     y = .data$NES, fill = .data$color), stat = "identity") +
        scale_fill_gradient2(low = color_set[1],
                             mid = color_set[2],
                             high = color_set[3],
                             midpoint = 0) +
        coord_flip() +
        labs(x = "Pathway",
             y = "Normalized Enrichment Score (NES)",
             fill = expression(-Log[10]~Adjusted~P-value)) +
        theme_bw() +
        theme(legend.position = "bottom",
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 12))
}

# ----------

#' Dotplot from a list of GSEA results
#'
#' @param gsea_list A list of tibble of GSEA results
#' @param p_co A double - the cutoff of adjusted p-value
#' @param pattern A string - the pattern to remove in the plot. Default value is \code{"HALLMARK_"}.
#'
#' @return A ggplot2 plot
#' @importFrom ggplot2 scale_color_gradient2
#' @export

plot_deseq_gsea_list <- function(gsea_list, p_co, pattern = "HALLMARK_") {

    gsea_list_new <- list()

    for (i in 1:length(gsea_list)) {
        gsea_list_new[[i]] <- gsea_list[[i]] %>%
            gsea_rm_pattern(pattern = pattern) %>%
            add_column(Comparison = as.character(names(gsea_list)[i]))
    }

    gsea_df <- do.call("rbind", gsea_list_new)

    color_set <- xmc_constants()$tricolor

    gsea_df %>%
        mutate(color = -log10(.data$padj) * ifelse(.data$padj <= p_co, 1, 0) * ifelse(.data$NES > 0, 1, -1)) %>%
        filter(.data$color != 0) %>%
        ggplot() +
        geom_point(aes(x = .data$pathway,
                       y = factor(.data$Comparison),
                       size = abs(.data$NES),
                       color = .data$color)) +
        scale_color_gradient2(low = color_set[1],
                              mid = color_set[2],
                              high = color_set[3],
                              midpoint = 0) +
        labs(color = bquote("-log"[10] ~ p-value %*% direction),
             size = "abs(Normalized enrichment score)",
             x = "Pathway",
             y = "Comparison") +
        coord_flip() +
        theme_bw() +
        theme(legend.position = "bottom",
              axis.text = element_text(size = 10),
              axis.title = element_text(size = 12),
              axis.text.x = element_text(angle = 45, vjust = 0.5))
}


# ----------
# PART 2 - Helper Functions
# ----------

#' Generate vsd object from counts and metadata
#'
#' @param counts A matrix - the RNA-seq count matrix
#' @param metadata A metadata - the metadata matrix
#'
#' @return A vsd object
#' @export

cts_to_vsd <- function(counts, metadata) {

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                          colData = metadata,
                                          design= ~ 1)
    dds <- DESeq2::DESeq(dds)
    vsd <- DESeq2::vst(dds, blind = FALSE)

    return(vsd)
}

# ----------

#' Transform DESeq2Results object to tibble
#'
#' @param res A DESeq2Results object
#' @param p_co A double - the cutoff of adjusted p-value
#' @param lfc_co A double - the cutoff of log2 fold change
#'
#' @return A tibble
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr filter mutate
#' @export

res_to_tibble <- function(res, p_co, lfc_co) {

    res <- res %>%
        as.data.frame() %>%
        rownames_to_column(var = "symbol") %>%
        as_tibble() %>%
        filter(!is.na(.data$padj)) %>%
        mutate(significant = ifelse(.data$padj <= p_co & .data$log2FoldChange >= lfc_co,
                                    "Up",
                                    ifelse(.data$padj <= p_co & .data$log2FoldChange <= -lfc_co,
                                           "Down",
                                           "Not Sig")))
    return(res)
}

# ----------

#' GSEA from DESeq2Results object
#'
#' @param res A DESeq2Results object
#' @param pathways A list - the list of pathway genes
#'
#' @return A tibble
#' @export

res_to_gsea <- function(res, pathways) {

    stat <- res$stat
    names(stat) <- rownames(res)

    gsea <- fgsea::fgsea(pathways, stat, eps = 0) %>% as_tibble()

    return(gsea)
}

# ----------

#' Remove repeated string pattern from the pathway names of GSEA results
#'
#' @param gsea A tibble of GSEA results
#' @param pattern A string - the pattern to remove in the plot. Default value is \code{"HALLMARK_"}
#'
#' @return A tibble of GSEA results
#' @importFrom stringr str_remove
#' @export

gsea_rm_pattern <- function(gsea, pattern = "HALLMARK_") {

    gsea_new <- gsea %>%
        mutate(pathway = str_remove(string = .data$pathway, pattern = pattern))

    return(gsea_new)
}
