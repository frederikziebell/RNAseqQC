plot_sample_distributions <- function(se, se_assay = 1, n_steps = 100){
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