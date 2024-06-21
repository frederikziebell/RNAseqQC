#' for a character vector x,
#' check if all non-NA elements of x
#' can be converted to numeric
#' @param x A character vector
all_numeric <- function(x){

  x_no_na <- x[!is.na(x)]

  if(length(x_no_na) == 0){
    FALSE
  } else {

    suppressWarnings(
      x_numeric <- as.numeric(x_no_na)
    )

    all(!is.na(x_numeric))
  }
}

.compute_pca <- function(obj, n_feats, scale_feats, na_frac, metadata, color_by){

  if (inherits(obj, "SummarizedExperiment")) {
    mat <- tryCatch(as.matrix(assay(obj)), error = function(e) {
      stop("Assay cannot be converted to a matrix.")
    })
    col_data <- as_tibble(colData(obj))
    row_data <- as_tibble(rowData(obj), rownames = "rowname")
  } else if (inherits(obj, "matrix")) {
    mat <- obj
    col_data <- NULL
    row_data <- NULL
  } else {
    stop("Object is not a matrix or SummarizedExperiment")
  }

  mat[is.nan(mat)] <- NA

  # subset to features with maximum specified fraction of NAs
  nrow_mat_original <- nrow(mat)
  mat_NA_filtered <- mat[rowSums(is.na(mat))/ncol(mat) <= na_frac, ]

  if(nrow(mat_NA_filtered) == 0) {
    stop("There are no features passing the 'na_frac' filtering criterion.")
  } else if (nrow(mat_NA_filtered) / nrow_mat_original < 1/100) {
    warning("Less than 1% of the features pass the 'na_frac' filtering criterion.")
  }

  # center rows
  mat_centered <- mat_NA_filtered %>%
    t() %>%
    scale(scale = scale_feats) %>%
    t()
  # mean impute NAs
  mat_centered[is.na(mat_centered)] <- 0
  # subset to top-variable features
  mat_centered <- mat_centered[matrixStats::rowVars(mat_centered, na.rm = T) %>% order(decreasing = T) %>% head(n_feats), ]
  # do PCA
  pca <- prcomp(t(mat_centered))
  pca_data <- as_tibble(pca$x, rownames = "colname")
  var_exp <- 100 * pca$sdev^2 / sum(pca$sdev^2)

  # join annotation
  if (!is.null(col_data)) {
    pca_data <- bind_cols(pca_data, suppressWarnings(select(col_data, -tidyselect::any_of("colname"))))
  }
  if (!is.null(metadata)) {
    pca_data <- bind_cols(pca_data, suppressWarnings(select(metadata, -tidyselect::any_of("colname"))))
  }

  # first annotation, then PCs
  pca_data <- pca_data %>%
    select(colname, !tidyselect::starts_with("PC"), tidyselect::starts_with("PC"))

  loadings_data <- as_tibble(pca$rotation, rownames = "rowname")
  # add row_data columns from SE object
  if (!is.null(row_data)) {
    loadings_data <- left_join(loadings_data, row_data, by = "rowname") %>%
      select(rowname, !tidyselect::starts_with("PC"), tidyselect::starts_with("PC"))
  }

  # if color_by specifies a gene, i.e. a rowname of obj
  # or an element of rowData(obj)$gene_name, extract this
  # column and cbind it to pca_data
  if(!is.null(color_by)){
    color_gene <- NULL
    if(color_by %in% rownames(mat)){
      color_gene <- color_by
    } else if(!is.null(row_data) & "gene_name" %in% colnames(row_data)){
      if(color_by %in% row_data[["gene_name"]]){
        color_gene <- rownames(mat)[color_by == row_data[["gene_name"]]][1]
      }
    }
    if(!is.null(color_gene)){
      df <- data.frame(color_by = mat[color_gene,])
      colnames(df) <- color_by
      pca_data <- bind_cols(pca_data, df)
    }
    if(is.null(color_gene) & !color_by %in% colnames(pca_data)){
      stop("Could not find metadata column or gene specified by 'color_by'.")
    }
  }

  numeric_color_by <- NULL
  if(!is.null(color_by)){
    numeric_color_by <- FALSE
    if(is.numeric(pca_data[[color_by]])){
      numeric_color_by <- TRUE
    } else if(all_numeric(pca_data[[color_by]])){
      numeric_color_by <- TRUE
      pca_data[[color_by]] <- as.numeric(pca_data[[color_by]])
    }
  }

  list(
    "pca_data" = pca_data,
    "var_exp" = var_exp,
    "loadings_data" = loadings_data,
    "numeric_color_by" = numeric_color_by
  )

}

.make_pca_plot <- function(pca_res, PC_x, PC_y, color_by, shape_by, point_alpha, point_rel_size, show_plot, rasterise, ...){

  pca_data <- pca_res$pca_data
  var_exp <- pca_res$var_exp
  numeric_color_by <- pca_res$numeric_color_by

  if(rasterise == TRUE){
    raster_fun <- function(x) {ggrastr::rasterise(x, ...)}
  } else {
    raster_fun <- identity
  }

  plot <- pca_data %>%
    ggplot(aes(
      .data[[paste0("PC", PC_x)]],
      .data[[paste0("PC", PC_y)]],
      label = colname,
      color = if (!is.null(color_by)){.data[[color_by]]} else {NULL},
      shape = if (!is.null(shape_by)){.data[[shape_by]]} else {NULL},
    )) +
    raster_fun(geom_point(size = rel(point_rel_size), alpha = point_alpha)) +
    coord_fixed() +
    labs(
      x = paste0("PC", PC_x, " (", round(var_exp[PC_x],1), "%)"),
      y = paste0("PC", PC_y, " (", round(var_exp[PC_y],1), "%)")
    ) +
    scale_shape_discrete(name = shape_by) +
    cowplot::theme_cowplot()

  # choose colors depending based on whether the variable to color by
  # is numeric, discrete with <= 10 levels or discrete with > 10 levels
  if(!is.null(color_by)){
    if(numeric_color_by){
      plot <- plot + ggplot2::scale_color_viridis_c(name = color_by)
      # at most 10 discrete levels
    } else if(length(unique(pca_data[[color_by]])) <= 10){
      pal_jco <- c(
        "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF",
        "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#4A6990FF"
      )
      plot <- plot + scale_color_manual(name = color_by, values = pal_jco)
    } else {
      plot <- plot + scale_color_viridis_d(name = color_by)
    }
  }

  plot
}

#' Plot results of a principal component analysis
#'
#' @param obj A (features x samples) matrix or SummarizedExperiment object
#' @param PC_x The PC to show on the x-axis.
#' @param PC_y The PC to show on the y-axis.
#' @param n_feats Number of top-variable features to include.
#' @param scale_feats Whether to scale the features.
#' @param na_frac Only consider features with the stated maximum fraction of NAs or NaNs. NA/NaNs will be mean-imputed for PCA.
#' @param metadata A data.frame used for annotating samples. `rownames(metadata)` must match `colnames(obj)`.
#' @param color_by Variable by which to color points. Must be a column in metadata or in `colData(obj)`.
#' Alternatively, it can be the name of a feature (a rowname of obj) or a gene name (an element of rowData(obj)$gene_name).
#' @param shape_by Variable by which to color points. Must be a column in metadata or in `colData(obj)`.
#' @param point_alpha alpha value of `geom_point()`
#' @param point_rel_size relative size of `geom_point()`
#' @param show_plot Whether to show the plot or not
#' @param rasterise Whether to rasterise the point, using ggrastr.
#' @param ... Other parameters passed on to ggrastr::rasterise
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(100*6), ncol=6)
#' data <- t(t(data)+c(-1, -1.1, -1.2, 1, 1.1, 1.2))
#' plot_pca(data)
#'
#' @return The function displays the plot and returns invisible a list of the plot,
#' the data.frame to make the plot, the vector of percentages of variance explained and the loadings matrix.
#' @details If the `metadata` or `colData` of `obj` contain a column `colname`, this colum will be removed in the `$pca_data` slot,
#' because this column contains the colnames of the data matrix. Similarly, for the `$loadings` slot, the column `rowname` is
#' reserved for the rownames of the data matrix.
#' @export
plot_pca <- function(obj,
  PC_x = 1, PC_y = 2,
  n_feats = 500, scale_feats = FALSE, na_frac = .3,
  metadata = NULL, color_by = NULL, shape_by = NULL,
  point_alpha = .7, point_rel_size = 2,
  show_plot = TRUE, rasterise = FALSE, ...
  ) {

  pca_res <- .compute_pca(
    obj = obj, n_feats = n_feats, scale_feats = scale_feats,
    na_frac = na_frac, metadata = metadata, color_by = color_by
  )

  plot <- .make_pca_plot(
    pca_res = pca_res, PC_x = PC_x, PC_y = PC_y, color_by = color_by,
    shape_by = shape_by, point_alpha = point_alpha, point_rel_size = point_rel_size,
    show_plot = show_plot, rasterise = rasterise, ...
  )

  if (show_plot) {
    print(plot)
  }

  invisible(list("plot" = plot, "data" = pca_res$pca_data, "var_exp" = pca_res$var_exp, "loadings" = pca_res$loadings_data))
}

#' Plot loadings of a principal component
#' @param pca_res A result returned from `plot_pca()`
#' @param PC Number of the principal component to plot
#' @param square Whether to plot squared loadings. The squared loading is equal to the fraction of variance explained by the respective feature in the given principal component.
#' @param color_by Variable (column in `pca_res$loadings`) to color points by. Can also be 'pc_sign' to color by the sign of the loading (useful in combination with the square = TRUE parameter).
#' @param annotate_top_n Annotate the top n features with positive or negative loading
#' @param highlight_genes Vector of gene names or gene IDs to highlight on the plot (overwrites top_n annotation)
#' @param show_plot Whether to show the plot
#'
#' @return The function displays the loadings plot and returns invisible a list of the plot,
#' the data.frame of the PCA loadings.
#'
#' @examples
#' set.seed(1)
#' data <- matrix(rnorm(100*6), ncol=6)
#' data <- t(t(data)+c(-1, -1.1, -1.2, 1, 1.1, 1.2))
#' pca_res <- plot_pca(data)
#' plot_loadings(pca_res)
#' @export
plot_loadings <- function(pca_res, PC = 1, square = FALSE, color_by = NULL, annotate_top_n = 0, highlight_genes = NULL, show_plot = TRUE) {

  loadings_df <- pca_res$loadings %>%
    select(!tidyselect::starts_with("PC"), PC = paste0("PC", PC))

  loadings_df$pc_sign <- as.character(sign(loadings_df$PC))
  if(square == TRUE){
    loadings_df$PC = (loadings_df$PC)^2
    loadings_df$PC_rank = rank(-loadings_df$PC)
  } else{
    loadings_df$PC_rank = rank(loadings_df$PC)
  }

  # top n features by loading
  annotate_top_n_features <- loadings_df %>%
    filter(
      (rank(PC) <= annotate_top_n) |
        (rank(-PC) <= annotate_top_n)
    ) %>%
    pull(rowname)

  # if pca_res was generated from a matrix,
  # there is no column gene_name, so we use
  # the rowname as gene_name to be able to plot
  # feature names
  if(! "gene_name" %in% colnames(loadings_df)){
    loadings_df$gene_name <- loadings_df$rowname
  }

  if(!is.null(color_by)){
    numeric_color_by <- FALSE
    if(is.numeric(loadings_df[[color_by]])){
      numeric_color_by <- TRUE
    } else if(all_numeric(loadings_df[[color_by]])){
      numeric_color_by <- TRUE
      loadings_df[[color_by]] <- as.numeric(loadings_df[[color_by]])
    }
  }

  plot <- loadings_df %>%
    ggplot(aes(x = PC_rank, y = PC, label = gene_name, color = if (!is.null(color_by)){.data[[color_by]]} else {NULL})) +
    geom_point() +
    cowplot::theme_cowplot()

  # choose colors depending based on whether the variable to color by
  # is numeric, discrete with <= 10 levels or discrete with > 10 levels
  if(!is.null(color_by)){
    if(numeric_color_by){
      plot <- plot + ggplot2::scale_color_viridis_c(name = color_by)
      # at most 10 discrete levels
    } else if(length(unique(loadings_df[[color_by]])) <= 10){

      pal_jco <- c(
        "#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF", "#7AA6DCFF",
        "#003C67FF", "#8F7700FF", "#3B3B3BFF", "#A73030FF", "#4A6990FF"
      )
      plot <- plot + scale_color_manual(name = color_by, values = pal_jco)

    } else {
      plot <- plot + scale_color_viridis_d(name = color_by)
    }
  }

  if (is.null(highlight_genes)) {
    plot <- plot +
      # highlight top n
      ggrepel::geom_text_repel(data = loadings_df %>% filter(
        rowname %in% annotate_top_n_features
      ))
  } else {
    plot <- plot +
      gghighlight::gghighlight(gene_name %in% highlight_genes | rowname %in% highlight_genes,
        label_key = gene_name, use_group_by = F,
        label_params = list(fill = NA, label.size = NA, fontface = "bold")
      )
  }

  if(square == TRUE){
    plot <- plot +
      ylim(c(0, max(loadings_df$PC))) +
      labs(x = paste0("rank(squared PC", PC, " loading)"), y = paste0("squared PC", PC, " loading"), color = color_by)
  } else {
    plot <- plot +
      ylim(c(-max(abs(loadings_df$PC)), max(abs(loadings_df$PC)))) +
      labs(x = paste0("rank(PC", PC, " loading)"), y = paste0("PC", PC, " loading"), color = color_by)
  }

  if (show_plot) {
    print(plot)
  }

  invisible(list("plot" = plot, "data" = loadings_df))
}

#' Plot matrix of PCA scatter plots
#'
#' @param obj A (features x samples) matrix or SummarizedExperiment object
#' @param n_PCs Number of principal components to plot
#' @param show_var_exp Whether to show a plot of the percentage of variance explained by each PC in the bottom left corner.
#' @param n_feats Number of top-variable features to include.
#' @param scale_feats Whether to scale the features.
#' @param na_frac Only consider features with the stated maximum fraction of NAs or NaNs. NA/NaNs will be mean-imputed for PCA.
#' @param metadata A data.frame used for annotating samples. `rownames(metadata)` must match `colnames(obj)`.
#' @param color_by Variable by which to color points. Must be a column in metadata or in `colData(obj)`.
#' Alternatively, it can be the name of a feature (a rowname of obj) or a gene name (an element of rowData(obj)$gene_name).
#' @param shape_by Variable by which to color points. Must be a column in metadata or in `colData(obj)`.
#' @param point_alpha alpha value of `geom_point()`
#' @param point_rel_size relative size of `geom_point()`
#' @param transpose Wheter to transpose the whole matrix of scatter plots
#' @param rasterise Whether to rasterise the points using ggrastr.
#' @param ... Other parameters passed on to ggrastr::rasterise
#'
#' @examples
#' \donttest{
#' set.seed(1)
#' data <- matrix(rnorm(100*6), ncol=6)
#' data <- t(t(data)+c(-1, -1.1, -1.2, 1, 1.1, 1.2))
#' plot_pca_scatters(data)
#' }
#' @return The function displays the scatter plots of the PCs
#' @export
plot_pca_scatters <- function(obj,
                              n_PCs = min(10,nrow(obj),ncol(obj)), show_var_exp = T, n_feats = 500, scale_feats = FALSE, na_frac = 0.3,
                              metadata = NULL, color_by = NULL, shape_by = NULL, point_alpha = 0.7, point_rel_size = 2, transpose = FALSE, rasterise = FALSE, ...){

  pca_res <- .compute_pca(
    obj = obj, n_feats = n_feats, scale_feats = scale_feats,
    na_frac = na_frac, metadata = metadata, color_by = color_by
  )

  # all pairwise combinations
  plots <- lapply(1:(n_PCs-1), function(PC_x){

    lapply(2:n_PCs, function(PC_y){

      if(PC_x >= PC_y){
        return(cowplot::ggdraw())
      }

      if(transpose){
        tmp <- PC_x
        PC_x <- PC_y
        PC_y <- tmp
      }

      plot <- .make_pca_plot(
        pca_res = pca_res, PC_x = PC_x, PC_y = PC_y,
        color_by = color_by, shape_by = shape_by,
        point_alpha = point_alpha, point_rel_size = point_rel_size,
        show_plot = F, rasterise = rasterise, ...
      )

      suppressMessages({
        plot <- plot +
          labs(x = paste0("PC",PC_x), y = paste0("PC",PC_y)) +
          coord_fixed(
            ratio =
              (max(pca_res$pca_data[[paste0("PC",PC_x)]])-min(pca_res$pca_data[[paste0("PC",PC_x)]])) /
              (max(pca_res$pca_data[[paste0("PC",PC_y)]])-min(pca_res$pca_data[[paste0("PC",PC_y)]]))
          )
      })

      plot

    })
  })

  plots <- purrr::list_flatten(plots)

  if(show_var_exp == TRUE & n_PCs >= 3){
    p_var_exp <- data.frame(var_exp = pca_res$var_exp) %>%
      ggplot(aes(rank(-var_exp), var_exp)) +
      geom_point(size = rel(point_rel_size), alpha = point_alpha) +
      scale_y_sqrt(limits=c(0,max(pca_res$var_exp))) +
      labs(x = "PC", y="% variance explained") +
      cowplot::theme_cowplot()

    plots[[(n_PCs-1)*(n_PCs-2)+1]] <- p_var_exp
  }

  patchwork::wrap_plots(plots, nrow=n_PCs-1, ncol=n_PCs-1, byrow = transpose, guides="collect")
}
