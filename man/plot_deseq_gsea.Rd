% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq.R
\name{plot_deseq_gsea}
\alias{plot_deseq_gsea}
\title{Barplot from GSEA results}
\usage{
plot_deseq_gsea(gsea, pattern = "HALLMARK_")
}
\arguments{
\item{gsea}{A tibble of GSEA results}

\item{pattern}{A string - the pattern to remove in the plot. Default value is \code{"HALLMARK_"}.}
}
\value{
A ggplot2 plot
}
\description{
Barplot from GSEA results
}
