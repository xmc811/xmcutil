
# ----------
# PART 1 - Plotting Functions
# ----------

#' PCA plot from vsd object
#'
#' @param vsd A vsd object
#' @param var A string - the name of the variable matching metadata columns
#' @param pal A string - palette name of \code{RColorBrewer}. Default value is \code{NULL}, taking the palettes set up in \code{xmc_constants()}.
#' @param dir An integer - \code{1} or \code{-1}, to adjust the direction of colors. Default value is \code{1}.
#'
#' @importFrom ggplot2 ggplot geom_point aes scale_color_brewer scale_color_distiller labs theme_bw theme
#'
#' @return A ggplot2 plot
#'
#' @export

plot_pca_vsd <- function(vsd, var, pal = NULL, dir = 1) {

    pca <- DESeq2::plotPCA(vsd, intgroup = var, returnData = TRUE)

    if (is.null(pal)) {
        pal_set <- xmc_constants()$palette
        if (is.numeric(pca[[var]])) {
            pal <- pal_set[2]
        } else {
            pal <- pal_set[1]
        }
    }

    if (is.numeric(pca[[var]])) {
        ggplot(pca) +
            geom_point(aes(x = .data$PC1,
                           y = .data$PC2,
                           color = .data$group), size = 2) +
            scale_color_distiller(palette = pal, direction = dir) +
            labs(color = var) +
            theme_bw() +
            theme(aspect.ratio = 1)
    } else {
        ggplot(pca) +
            geom_point(aes(x = .data$PC1,
                           y = .data$PC2,
                           color = .data$group), size = 2) +
            scale_color_brewer(palette = pal) +
            labs(color = var) +
            theme_bw() +
            theme(aspect.ratio = 1)
    }
}

# ----------

#' Heatmap from vsd object
#'
#' @param vsd A vsd object
#' @param var A string - the name of the variable matching metadata columns
#' @param pal A string - palette name of \code{RColorBrewer}. Default value is \code{NULL},taking the palettes set up in \code{xmc_constants()}.
#' @param dir An integer - \code{0} or \code{1}, to adjust the direction of colors. Default value is \code{1}.
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom circlize colorRamp2
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom stats dist quantile
#' @importFrom grid gpar
#'
#' @return A heatmap
#'
#' @export

plot_heatmap_vsd <- function(vsd, var, pal = NULL, dir = 1) {

    if (is.null(pal)) pal <- xmc_constants()$palette[2]

    sampleDistMatrix <- as.matrix(dist(t(SummarizedExperiment::assay(vsd))))

    rownames(sampleDistMatrix) <- vsd[[var]]
    colnames(sampleDistMatrix) <- vsd[[var]]

    sampleDistMatrix <- log2(sampleDistMatrix + 1)

    colors <- get_all_colors(pal)

    if (dir) colors <- rev(colors)

    col_fun <- colorRamp2(seq(from = quantile(sampleDistMatrix,
                                              1/nrow(sampleDistMatrix)),
                              to = max(sampleDistMatrix),
                              length.out = length(colors)), colors)

    Heatmap(sampleDistMatrix,
            col = col_fun,
            rect_gp = gpar(col = "white", lwd = 2),
            heatmap_width = unit(1, "npc"),
            heatmap_height = unit(1, "npc"),
            row_names_gp = gpar(fontsize = 12),
            column_names_gp = gpar(fontsize = 12),
            heatmap_legend_param = list(title = "Distance",
                                        border = "black"))
}

# ----------

#' MA plot from DESeq2Results object
#'
#' @param res A DESeq2Results object
#' @param p_co A double - the cutoff of adjusted p-value. Default values is \code{0.05}.
#' @param lfc_co A double - the cutoff of log2 fold change. Default values is \code{2}.
#' @param lfc_plot_lim A double - the y-limit of log2 fold change plot. Default value is \code{5}.
#'
#' @importFrom dplyr arrange
#' @importFrom ggplot2 scale_color_manual scale_shape_manual element_text geom_hline
#' @importFrom grid unit
#'
#' @return A ggplot2 plot
#'
#' @export

plot_deseq_ma <- function(res, p_co = 0.05, lfc_co = 2, lfc_plot_lim = 5) {

    res <- res %>%
        res_to_tibble(p_co, lfc_co)

    xmc_const <- xmc_constants()

    res %>%
        res_add_shape(lfc_plot_lim) %>%
        arrange(factor(.data$significant, levels = xmc_const$deg_levels)) %>%
        ggplot() +
        geom_point(aes(x = log10(.data$baseMean),
                       y = .data$log2FoldChange,
                       color = .data$significant,
                       shape = .data$shape1),
                   size = 2) +
        scale_color_manual(values = xmc_const$tricolor) +
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
#' @param p_co A double - the cutoff of adjusted p-value. Default values is \code{0.05}.
#' @param lfc_co A double - the cutoff of log2 fold change. Default values is \code{2}.
#' @param p_plot_lim A double - the y-limit of -log10 adjusted p-value. Default value is \code{5}.
#' @param lfc_plot_lim A double - the x-limit of log2 fold change plot. Default value is \code{5}.
#'
#' @return A ggplot2 plot
#'
#' @export

plot_deseq_volcano <- function(res,
                               p_co = 0.05,
                               lfc_co = 2,
                               p_plot_lim = 5,
                               lfc_plot_lim = 5) {

    res <- res %>%
        res_to_tibble(p_co, lfc_co)

    xmc_const <- xmc_constants()

    res %>%
        res_add_shape(lfc_plot_lim) %>%
        mutate(shape2 = ifelse(-log10(.data$padj) > p_plot_lim, TRUE, FALSE),
               padj = ifelse(-log10(.data$padj) > p_plot_lim, 10^(-p_plot_lim), .data$padj),
               shape = .data$shape1 | .data$shape2) %>%
        arrange(factor(.data$significant, levels = xmc_const$deg_levels)) %>%
        ggplot() +
        geom_point(aes(x = .data$log2FoldChange,
                       y = -log10(.data$padj),
                       color = .data$significant,
                       shape = factor(.data$shape)),
                   size = 2) +
        scale_color_manual(values = xmc_const$tricolor) +
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

#' Sample-gene matrix heatmap from DESeqDataSet object
#'
#' @param dds A DESeqDataSet object
#' @param genes A string vector - list of genes
#' @param pal A string - palette name of \code{RColorBrewer}
#' @param dir An integer - \code{1} or \code{-1}, to adjust the direction of colors.
#'
#' @return A heatmap
#'
#' @export

plot_sample_gene_mtx <- function(dds, genes, pal, dir) {

    mtx <- get_mtx_dds(dds, genes) %>% mtx_rescale()

    colors <- get_all_colors(pal)

    if (dir) colors <- rev(colors)

    col_fun <- colorRamp2(seq(from = -1,
                              to = 1,
                              length.out = length(colors)), colors)

    Heatmap(mtx,
            col = col_fun,
            rect_gp = gpar(col = "white", lwd = 2))
}

# ----------

#' Boxplot from DESeqDataSet object
#'
#' @param dds A DESeqDataSet object
#' @param genes A string vector - list of genes
#' @param var A string - the name of the variable matching metadata columns
#' @param pal A string - palette name of \code{RColorBrewer}
#'
#' @importFrom ggplot2 facet_wrap geom_boxplot scale_fill_brewer
#' @importFrom rlang sym
#'
#' @return A heatmap
#'
#' @export

plot_gene_boxplot <- function(dds, genes, var, pal) {

    df <- get_nm_count_dds(dds, genes, var)

    ggplot(df) +
        geom_boxplot(aes(x = !!sym(var),
                         y = log10(.data$count),
                         fill = !!sym(var))) +
        geom_point(aes(x = !!sym(var),
                       y = log10(.data$count))) +
        facet_wrap(~symbol) +
        scale_fill_brewer(palette = pal)
}

# ----------

#' Barplot from GSEA results
#'
#' @param gsea A tibble of GSEA results
#' @param pattern A string - the pattern to remove in the plot. Default value is \code{"HALLMARK_"}.
#'
#' @importFrom stringr str_remove
#' @importFrom stats reorder
#' @importFrom ggplot2 geom_bar scale_fill_gradient2 coord_flip
#'
#' @return A ggplot2 plot
#'
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
#' @param p_co A double - the cutoff of adjusted p-value. Default value is \code{0.05}.
#' @param pattern A string - the pattern to remove in the plot. Default value is \code{"HALLMARK_"}.
#' @param dropNonSig A logical - whether to drop non-significant pathways in the plot. Default value is \code{TRUE}.
#'
#' @importFrom ggplot2 scale_color_gradient2
#'
#' @return A ggplot2 plot
#'
#' @export

plot_deseq_gsea_list <- function(gsea_list, p_co = 0.05, pattern = "HALLMARK_", dropNonSig = TRUE) {

    gsea_list_new <- list()

    for (i in 1:length(gsea_list)) {
        gsea_list_new[[i]] <- gsea_list[[i]] %>%
            gsea_rm_pattern(pattern = pattern) %>%
            add_column(Comparison = as.character(names(gsea_list)[i]))
    }

    gsea_df <- do.call("rbind", gsea_list_new)

    color_set <- xmc_constants()$tricolor

    gsea_df <- gsea_df %>%
        mutate(color = -log10(.data$padj) * ifelse(.data$padj <= p_co, 1, 0) * ifelse(.data$NES > 0, 1, -1))

    if (dropNonSig) gsea_df <- gsea_df %>% filter(.data$color != 0)

    ggplot(gsea_df) +
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
# PART 2 - Pipeline Functions
# ----------

#' Generate vsd object from counts and metadata
#'
#' @param counts A matrix - the RNA-seq count matrix
#' @param metadata A metadata - the metadata matrix
#'
#' @return A vsd object
#'
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

#' GSEA from DESeq2Results object
#'
#' @param res A DESeq2Results object
#' @param pathways A list - the list of pathway genes. Default value is \code{hmks_hs}.
#'
#' @return A tibble
#'
#' @export

res_to_gsea <- function(res, pathways = hmks_hs) {

    stat <- res$stat
    names(stat) <- rownames(res)

    gsea <- fgsea::fgsea(pathways, stat, eps = 0) %>% as_tibble()

    return(gsea)
}


# ----------
# PART 3 - Helper Functions
# ----------

#' Add shape column for log2 fold change limit
#'
#' @param res A tibble
#' @param lfc_plot_lim A double - the x-limit of log2 fold change plot
#'
#' @return A tibble

res_add_shape <- function(res, lfc_plot_lim) {

    res_new <- res %>%
        mutate(shape1 = ifelse(.data$log2FoldChange > lfc_plot_lim | .data$log2FoldChange < -lfc_plot_lim, TRUE, FALSE),
               log2FoldChange = replace(.data$log2FoldChange, .data$log2FoldChange > lfc_plot_lim, lfc_plot_lim),
               log2FoldChange = replace(.data$log2FoldChange, .data$log2FoldChange < -lfc_plot_lim, -lfc_plot_lim))

    return(res_new)
}

# ----------

#' Transform DESeq2Results object to tibble
#'
#' @param res A DESeq2Results object
#' @param p_co A double - the cutoff of adjusted p-value
#' @param lfc_co A double - the cutoff of log2 fold change
#'
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom dplyr filter mutate
#'
#' @return A tibble

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

#' Remove repeated string pattern from the pathway names of GSEA results
#'
#' @param gsea A tibble of GSEA results
#' @param pattern A string - the pattern to remove in the plot. Default value is \code{"HALLMARK_"}
#'
#' @importFrom stringr str_remove
#'
#' @return A tibble of GSEA results

gsea_rm_pattern <- function(gsea, pattern = "HALLMARK_") {

    gsea_new <- gsea %>%
        mutate(pathway = str_remove(string = .data$pathway, pattern = pattern))

    return(gsea_new)
}

# ----------

#' Rescale RNA-seq sample-gene matrix
#'
#' @param mtx A matrix - the sample-gene matrix
#'
#' @return A matrix

mtx_rescale <- function(mtx) {

    mtx2 <- mtx

    for (i in seq_along(1:nrow(mtx))) {
        mtx2[i,] <- (mtx[i,] - min(mtx[i,]))/(max(mtx[i,]) - min(mtx[i,])) * 2 - 1
    }
    return(mtx2)
}

# ----------

#' Get sample-gene matrix from DESeqDataSet object
#'
#' @param dds A DESeqDataSet object
#' @param genes A string vector - list of genes
#' @param raw A logical - whether to get raw counts data. Default value is \code{FALSE}.
#'
#' @return A matrix

get_mtx_dds <- function(dds, genes, raw = F) {

    if (raw) {
        vsd <- DESeq2::counts(dds)
    } else {
        vsd <- DESeq2::vst(dds, blind = FALSE)
        vsd <- as.matrix(vsd@assays@data[[1]])
    }

    mtx <- vsd[genes,,drop = FALSE]

    return(mtx)
}

# ----------

#' Get sample-gene matrix from DESeqDataSet object
#'
#' @param dds A DESeqDataSet object
#' @param genes A string vector - list of genes
#' @param var A string - the name of the variable matching metadata columns
#'
#' @return A tibble

get_nm_count_dds <- function(dds, genes, var) {

    df_list <- list()

    genes <- genes[genes %in% rownames(dds)]

    for (i in seq_along(1:length(genes))) {

        d <- DESeq2::plotCounts(dds,
                                gene = genes[i],
                                intgroup = var,
                                returnData = TRUE)

        d <- d %>%
            rownames_to_column() %>%
            as_tibble() %>%
            mutate(symbol = genes[i])

        df_list[[i]] <- d
    }

    df <- do.call("rbind", df_list)

    return(df)
}


