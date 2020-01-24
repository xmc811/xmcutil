
# Tools for processing CNV data

#' Look up corresponding values from a lookup table
#'
#' @param values A numeric vector - the log2Ratio values for copy number changes.
#' @param cutoff A numeric vector of length 2 - the cutoff for calling amplifications or deletions
#'
#' @return A numeric vector of the same length with \code{values}.
#' @export

test_cnv <- function(values,
                     cutoff = c(-0.3, 0.3)) {

    if (cutoff[2] <= cutoff[1]) {
        stop("Cutoff values not accepted")
    }

    test_1 <- (values > cutoff[2])
    test_2 <- (values >= cutoff[1])

    test <- test_1 + test_2 - 1

    return(test)
}
