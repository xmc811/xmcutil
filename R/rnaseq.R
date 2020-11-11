
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

