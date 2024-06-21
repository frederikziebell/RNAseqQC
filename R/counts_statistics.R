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

  # prevent 'no visible binding for global variable' package warnings
  gene_id <- count <- sample_id <- detection_threshold <- number_of_genes <- fraction_of_genes <- fraction_of_reads <- NULL

  assay(dds) %>%
    as_tibble(rownames = "gene_id") %>%
    pivot_longer(-gene_id, names_to = "sample_id", values_to = "count") %>%
    arrange(-count) %>%
    group_by(sample_id) %>%
    summarize(
      detection_threshold = count,
      number_of_genes = row_number(),
      .groups = "drop"
    ) %>%
    group_by(sample_id, detection_threshold) %>%
    summarize(number_of_genes = max(number_of_genes), .groups = "drop") %>%
    filter(detection_threshold %in% thresholds) %>%
    mutate(detection_threshold = factor(detection_threshold)) %>%
    ggplot(aes(detection_threshold, number_of_genes, group = sample_id, label = sample_id)) +
    geom_path(alpha = .4) +
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

  # prevent 'no visible binding for global variable' package warnings
  x <- NULL

  data.frame(x = colSums(assay(dds))) %>%
    ggplot(aes(x = x)) +
    geom_histogram(bins = n_bins) +
    labs(x = "total count", y = "number of samples") +
    cowplot::theme_cowplot()
}

#' Plot the library complexity
#'
#' Plot per sample the fraction of genes, versus the fraction of total counts.
#' @param dds A DESeqDataSet
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
plot_library_complexity <- function(dds) {

  # prevent 'no visible binding for global variable' package warnings
  gene_id <- count <- sample_id <- fraction_of_genes <- fraction_of_counts <- NULL

  purrr::map_dfr(1:ncol(dds), function(i) {
    cts <- sort(assay(dds)[,i], decreasing=T)
    tibble(
      sample_id = SummarizedExperiment::colnames(dds)[i],
      fraction_of_counts = cumsum(cts)/sum(cts),
      fraction_of_genes = (1:length(cts))/length(cts)
    )
  }) %>%
    ggplot(aes(fraction_of_genes, fraction_of_counts, group = sample_id, label = sample_id)) +
      geom_line(alpha = .2) +
      theme_minimal() +
      scale_x_sqrt(breaks = c(.01, .05, .1, .2, .5, 1)) +
      scale_y_continuous(breaks = seq(0, 1, .2)) +
      labs(x = "fraction of genes", y = "fraction of counts") +
      theme(legend.position = "none") +
      cowplot::theme_cowplot()
}
