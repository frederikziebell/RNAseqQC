#' Make DESeqDataSet from counts matrix and metadata
#'
#' @param counts The genes x samples counts matrix. The row names must be ENSEMBL gene IDs.
#' @param metadata data.frame of sample information. Order of rows corresponds to the order of columns in the counts matrix.
#' @param ah_record ID of AnnotationHub record used to retrieve an `EnsDb` object.
#' @param design The design formula specified in `DESeqDataSet()`
#' To view all valid record IDs, run
#' ```
#' library(AnnotationHub)
#' mcols(AnnotationHub()) %>%
#'  as_tibble(rownames="ah_record") %>%
#'  filter(rdataclass=="EnsDb")
#'  ```
#'
#' @return A `DESeq2::DESeqDataSet` object containing the counts matrix and metadata.
#'
#' @examples
#' \donttest{
#' library("DESeq2")
#' count_mat <- counts(T47D)
#' meta <- data.frame(colData(T47D))
#' dds <- make_dds(counts = count_mat, metadata = meta, ah_record = "AH89426")
#' }
#'
#' @export
make_dds <- function(counts, metadata, ah_record, design = ~1) {

  # prevent 'no visible binding for global variable' package warnings
  seqnames <- gene_id <- gene_name <- tx_id <- gc_content <- gene_biotype <- chromosome <- start <- NULL

  if (ncol(counts) != nrow(metadata)) {
    stop("The number columns in 'counts' does not equal the number of metadata rows.")
  }

  if (is.null(colnames(counts))) {
    stop("Matrix 'counts' must have column names.")
  }

  if (is.null(rownames(counts))) {
    stop("Matrix 'counts' must have row names.")
  } else if (!all(str_detect(rownames(counts), "^ENS"))) {
    stop("Rownames of matrix 'counts' must be ENSEMBL gene IDs.")
  }

  # load gene information, that will be stored in rowData slot
  suppressMessages({
    ah <- AnnotationHub::AnnotationHub()
    edb <- ah[[ah_record]]
  })
  ensembl_info <- ensembldb::genes(edb, filter = AnnotationFilter::GeneIdFilter(rownames(counts))) %>%
    as_tibble() %>%
    mutate(chromosome = as.character(seqnames)) %>%
    select(gene_id, gene_name, gene_biotype, chromosome, chr_start = start)

  gc_info <- ensembldb::transcripts(edb, filter = AnnotationFilter::GeneIdFilter(rownames(counts))) %>%
    as_tibble() %>%
    dplyr::select(tx_id, gene_id, gc_content) %>%
    group_by(gene_id) %>%
    summarize(gc_content = stats::median(gc_content), .groups = "drop")

  ensembl_info <- left_join(ensembl_info, gc_info, by="gene_id")

  not_annotated_genes <- setdiff(rownames(counts), ensembl_info$gene_id)
  if (length(not_annotated_genes)) {
    message(paste0("There are ", length(not_annotated_genes), " gene IDs that could not be annotated."))
  }

  row_data <- data.frame(gene_id = rownames(counts)) %>%
    left_join(ensembl_info, by = "gene_id") %>%
    mutate(gene_name = ifelse(is.na(gene_name), gene_id, gene_name)) %>%
    column_to_rownames("gene_id")

  dds <- suppressMessages(DESeq2::DESeqDataSetFromMatrix(counts, colData=metadata, rowData=row_data, design = design))
  # somehow in the design, there is also an environment stored, which bloats the whole dds object
  environment(DESeq2::design(dds)) <- NULL
  dds
}

#' Filter genes with low counts
#'
#' @param dds A DESeqDataSet
#' @param min_count,min_rep keep genes with at least min_count counts in at least min_rep replicates
#'
#' @return A `DESeq2::DESeqDataSet` object with only those genes that meet the filter criteria.
#'
#' @examples
#' library("DESeq2")
#' dds <- makeExampleDESeqDataSet()
#' filter_genes(dds)
#'
#' @importFrom DESeq2 counts
#' @export
filter_genes <- function(dds, min_count = 5, min_rep = 3) {
  dds[rowSums(counts(dds) >= min_count) >= min_rep, ]
}


#' Create a mean-sd plot
#' Make a scatterplot that shows for each gene its standard deviation versus mean.
#' @param vsd A DESeqTransform object
#'
#' @return A ggplot object of the ggplot2 package that contains the mean-sd plot.
#' @examples
#' \donttest{
#' library("DESeq2")
#' dds <- makeExampleDESeqDataSet(interceptMean=10, n=5000)
#' vsd <- vst(dds)
#' mean_sd_plot(vsd)
#' }
#' @export
mean_sd_plot <- function(vsd) {
  vsn::meanSdPlot(assay(vsd), plot = FALSE)$gg + cowplot::theme_cowplot()
}
