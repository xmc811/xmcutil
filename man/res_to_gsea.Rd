% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnaseq.R
\name{res_to_gsea}
\alias{res_to_gsea}
\title{GSEA from DESeq2Results object}
\usage{
res_to_gsea(res, pathways = hmks_hs)
}
\arguments{
\item{res}{A DESeq2Results object}

\item{pathways}{A list - the list of pathway genes. Default value is \code{hmks_hs}.}
}
\value{
A tibble
}
\description{
GSEA from DESeq2Results object
}
