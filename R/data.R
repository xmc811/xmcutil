
#' Hg19 gene start-end coordinates for CNV analysis
#'
#' A tibble of hg19 gene coordinates
#'
#' @format A tibble
#'
#' \describe{
#' \item{chrom}{Chromosome number.}
#' \item{start}{The starting coordinate.}
#' \item{end}{The ending coordinate}
#' \item{geneid}{Gene ID}
#' \item{genename}{Gene Symbol}
#' }
#' @examples
#' head(geneInfo)

"geneInfo"


#' Hallmark gene lists of homo sapiens
#'
#' A list of gene lists
#'
#' @format A list
#'
#' @examples
#' hmks_hs[[1]]

"hmks_hs"

