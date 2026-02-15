#' Plot one ECDF per sample
#' @export
plot_sample_ecdfs <- function(se, se_assay = 1, n_steps = 100){
  mat <- assay(se, i = se_assay)
  
  min_val <- min(mat, na.rm = T)
  max_val <- max(mat, na.rm = T)
  
  plot <- ggplot(
    data = data.frame(x = seq(min_val, max_val, length.out = n_steps)),
    mapping = aes(x)
  ) +
    cowplot::theme_cowplot()
  
  walk(colnames(mat), function(col){
    ecdf_fun = ecdf(mat[,col])
    plot <<- plot + geom_function(fun = ecdf_fun, alpha = .2)
  })
  
  plot
}

#' Plot one boxplot per sample
#' @export
plot_sample_boxplots <- function(se, se_assay = 1, ...){
  mat <- assay(se, i = se_assay)
  
  
  plot <- ggplot(mapping = aes(x,y)) +
    cowplot::theme_cowplot()
  
  walk(colnames(mat), function(col){
    plot <<- plot + 
      geom_boxplot(
        data = data.frame(
          x = col,
          y = mat[,col]
        ),
        ...
      )
  })
  
  plot
}