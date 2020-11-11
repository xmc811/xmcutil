
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
        scale_color_manual(values = c("#1f78b4","#d9d9d9", "#e31a1c")) +
        scale_shape_manual(values = c(16, 17)) +
        theme_bw() +
        labs(y = expression(Log[2]~Fold~Change), x = expression(Log[10]~Mean~Normalized~Count)) +
        theme(legend.position = "none",
              plot.margin = unit(rep(1,4), "cm"),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14)) +
        geom_hline(yintercept = 0, color = "#984ea3", size = 1.5, alpha = 0.5)

}


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
        scale_color_manual(values = c("#1f78b4","#d9d9d9", "#e31a1c")) +
        scale_shape_manual(values = c(16, 17), guide = FALSE) +
        theme_bw() +
        labs(y = expression(-Log[10]~Adjusted~p-value), x = expression(Log[2]~Fold~Change)) +
        theme(aspect.ratio = 1,
              legend.position = "none",
              plot.margin = unit(rep(1,4), "cm"),
              axis.text = element_text(size = 12),
              axis.title = element_text(size = 14))

}

