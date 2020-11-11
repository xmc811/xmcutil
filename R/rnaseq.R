
#' PCA plot from PCA coordinates
#'
#' @param df A dataframe - the PCA coordinates of RNA-seq samples
#' @param var A string - the name of the variable matching metadata columns
#'
#' @return A ggplot2 plot
#' @export

plot_pca_df <- function(df, var) {

    ggplot2::ggplot(df) +
        ggplot2::geom_point(ggplot2::aes(x = .data$PC1,
                                         y = .data$PC2,
                                         color = .data$group), size = 2) +
        ggplot2::labs(color = var) +
        ggplot2::theme_bw() +
        ggplot2::theme(aspect.ratio = 1)
}

#' PCA plot from PCA coordinates
#'
#' @param counts A matrix - the RNA-seq count matrix
#' @param metadata A metadata - the metadata matrix
#' @param vars A string vector - the vector of varibles to plot
#'
#' @return A list ggplot2 plots
#' @export

plot_rna_pca <- function(counts, metadata, vars) {

    dds <- DESeq2::DESeqDataSetFromMatrix(countData = counts,
                                          colData = metadata,
                                          design= ~ 1)
    dds <- DESeq2::DESeq(dds)
    vsd <- DESeq2::vst(dds, blind = FALSE)

    df_list <- purrr::map(.x = vars, .f = DESeq2::plotPCA, object = vsd, returnData = TRUE)

    plot_list <- purrr::map2(.x = df_list, .y = vars, .f = plot_pca_df)
    names(plot_list) <- vars

    return(plot_list)

}

