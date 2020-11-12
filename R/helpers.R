
# General Helper Functions

#' Get all colors in an RColorBrewer palette
#'
#' @param pal A string - the name of \code{RColorBrewer} palette
#'
#' @return An string vector

get_all_colors <- function(pal) {

    colors <- suppressWarnings({RColorBrewer::brewer.pal(12, pal)})
    return(colors)
}
