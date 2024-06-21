#' Plot number of detected genes for each sample
#'
#' For specified thresholds, the number of detected genes is shown for each sample.
#' @param dds A DESeqDataSet
#' @param thresholds Vector of thresholds for which the number of genes with
#' counts greater or equal than the thresholds is plotted
#'
#' @return A ggplot object of the ggplot2 package that contains the gene detection plot.
#'
#' @examples
#' library("DESeq2")
#' set.seed(1)
#' dds <- makeExampleDESeqDataSet()
#' plot_gene_detection(dds)
#'
#' @export
plot_gene_detection <- function(dds, thresholds = c(3, 10, 20, 50)) {

  # prevent integer overflow when adding very large counts
  counts <- assay(dds)
  mode(counts) <- "double"

  thresholds %>%
    map_dfr(function(threshold){
      colSums(counts >= threshold) %>%
        enframe(name = "sample_id", value = "number_of_genes") %>%
        mutate(detection_threshold = threshold)
    }) %>%
    mutate(detection_threshold = factor(detection_threshold, levels = sort(thresholds))) %>%
    ggplot(aes(detection_threshold, number_of_genes, group = sample_id, label = sample_id)) +
      geom_path(alpha = 0.4) +
      labs(x = "detection threshold", y = "number of genes") +
      cowplot::theme_cowplot()
}

#' Plot total counts per sample
#'
#' Plot the distribution of the total number of counts per sample as histogram.
#' @param dds A DESeqDataSet
#' @param n_bins Number of histogram bins
#'
#' @return A ggplot object of the ggplot2 package that contains the histogram of total counts per sample.
#'
#' @examples
#' library("DESeq2")
#' set.seed(1)
#' dds <- makeExampleDESeqDataSet(m=30)
#' plot_total_counts(dds)
#'
#' @export
plot_total_counts <- function(dds, n_bins = 50) {

  # prevent integer overflow when adding very large counts
  counts <- assay(dds)
  mode(counts) <- "double"

  data.frame(col_sum = colSums(counts)) %>%
    ggplot(aes(x = col_sum)) +
    geom_histogram(bins = n_bins) +
    labs(x = "total count", y = "number of samples") +
    cowplot::theme_cowplot()
}

# returns for each element in y the position of the
# largest element in x that is less or equal to y
pos_largest_lower_bound <- function(x, y){
  map_dbl(y, function(y){
    res <- which.max(x[x <= y])
    ifelse(length(res)==0,NA,res)
  })
}

#' Plot the library complexity
#'
#' Plot per sample the fraction of genes, versus the fraction of total counts.
#' @param dds A DESeqDataSet
#' @param show_progress Whether to show a progress bar of the computation.
#'
#' @return A ggplot object of the ggplot2 package that contains the library complexity plot.
#'
#' @examples
#' library("DESeq2")
#' set.seed(1)
#' dds <- makeExampleDESeqDataSet()
#' plot_library_complexity(dds)
#'
#' @export
plot_library_complexity <- function(dds, show_progress = TRUE) {

  # prevent integer overflow when adding very large counts
  counts <- assay(dds)
  mode(counts) <- "double"

  df <- map_dfr(1:ncol(dds), function(i) {
    cts <- sort(counts[, i], decreasing = T)
    fraction_of_genes <- seq_along(cts) / length(cts)
    fraction_of_counts <- cumsum(cts) / sum(cts)
    gene_frac_subset_pos <- c(which(fraction_of_genes < 0.01), pos_largest_lower_bound(fraction_of_genes, seq(.01, 1, .01)))

    tibble(
      sample_id = SummarizedExperiment::colnames(dds)[i],
      fraction_of_genes = fraction_of_genes[gene_frac_subset_pos],
      fraction_of_counts = fraction_of_counts[gene_frac_subset_pos]
    )
  }, .progress = show_progress) %>%
    filter(stats::complete.cases(.))


  df %>%
    ggplot(aes(fraction_of_genes, fraction_of_counts, group = sample_id, label = sample_id)) +
    geom_line(alpha = 0.2) +
    theme_minimal() +
    scale_x_sqrt(breaks = c(0.01, 0.05, 0.1, 0.2, 0.5, 1)) +
    scale_y_continuous(breaks = seq(0, 1, 0.2)) +
    labs(x = "fraction of genes", y = "fraction of counts") +
    theme(legend.position = "none") +
    cowplot::theme_cowplot()
}
