% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq.R
\name{gsea_rm_pattern}
\alias{gsea_rm_pattern}
\title{Remove repeated string pattern from the pathway names of GSEA results}
\usage{
gsea_rm_pattern(gsea, pattern = "HALLMARK_")
}
\arguments{
\item{gsea}{A tibble of GSEA results}

\item{pattern}{A string - the pattern to remove in the plot. Default value is \code{"HALLMARK_"}}
}
\value{
A tibble of GSEA results
}
\description{
Remove repeated string pattern from the pathway names of GSEA results
}
