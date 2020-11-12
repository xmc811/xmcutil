
# General Helper Functions

#' Number of colors in the RColorBrewer palette
#'
#' @param pal A string - the name of \code{RColorBrewer} palette
#'
#' @return An integer
#' @export

num_colors <- function(pal) {

    df <- RColorBrewer::brewer.pal.info

    a <- rownames(df) == pal
    return(df$maxcolors[a])
}
