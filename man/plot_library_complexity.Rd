% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/counts_statistics.R
\name{plot_library_complexity}
\alias{plot_library_complexity}
\title{Plot the library complexity}
\usage{
plot_library_complexity(dds, show_progress = TRUE)
}
\arguments{
\item{dds}{A DESeqDataSet}

\item{show_progress}{Whether to show a progress bar of the computation.}
}
\value{
A ggplot object of the ggplot2 package that contains the library complexity plot.
}
\description{
Plot per sample the fraction of genes, versus the fraction of total counts.
}
\examples{
library("DESeq2")
set.seed(1)
dds <- makeExampleDESeqDataSet()
plot_library_complexity(dds)

}
