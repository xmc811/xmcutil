
#' Helper function to generate constants used by other functions
#'
#' @export

xmc_constants <- function() {

    constants <- list()

    constants$tricolor <- c("#1f78b4","#d9d9d9","#e31a1c")
    constants$palette <- c("Set1", "BuPu")
    constants$deg_levels <- c("Down","Not Sig","Up")

    return(constants)
}
