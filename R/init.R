
#' Helper function to generate constants used by other functions
#'
#' @export

xmc_constants <- function() {

    constants <- list()

    constants$tricolor <- c("#1f78b4","#d9d9d9","#e31a1c")
    constants$tripalette <- c("Set1", "BuPu", "Spectral")
    constants$deg_levels <- c("Down","Not Sig","Up")

    return(constants)
}

# Setup default package options

.set_xmcutil_opts <- function(pkgname) {

    opts <- xmc_constants()

    options(xmcutil.tricolor = opts$tricolor)
    options(xmcutil.tripalette = opts$tripalette)
    options(xmcutil.deg_levels = opts$deg_levels)

}


