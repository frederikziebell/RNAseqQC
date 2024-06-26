% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/diff_expression.R
\name{plot_ma}
\alias{plot_ma}
\title{MA-plot of a differential testing result}
\usage{
plot_ma(de_res, dds, annotate_top_n = 5, highlight_genes = NULL)
}
\arguments{
\item{de_res}{An object returned by \code{DESeq2::results()} or \code{DESeq2::lfcShrink()}}

\item{dds}{The DESeqDataSet that was used to build the ´de_res´ object. This is needed for gene name annotation.}

\item{annotate_top_n}{Annotate the top n significant genes by fold change (up- and down-regulated)}

\item{highlight_genes}{Vector of gene names or gene IDs to highlight on the plot (overwrites top_n annotation)}
}
\value{
A ggplot object of the ggplot2 package that contains the MA-plot. The plot shows three classes of points:
Light gray points are genes with low counts that are removed from the analysis by independent filtering. Darker gray points
are not significant genes that show a density map to visualize where the majority of non-significant points are located. Finally,
red point show significant genes.
}
\description{
MA-plot of a differential testing result
}
\examples{
library("DESeq2")
set.seed(1)
dds <- makeExampleDESeqDataSet(n=1500, m=6, betaSD=.3, interceptMean=6)
rowData(dds)$gene_name <- rownames(dds)
dds <- DESeq(dds)
de_res <- results(dds)
plot_ma(de_res, dds)

}
