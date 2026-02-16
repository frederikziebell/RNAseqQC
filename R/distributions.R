#' Plot one ECDF per sample
#' @param se A SummarizedExperiment
#' @param se_assay Name or the number of the assay to plot
#' @param n_steps Number of points at which the ECDF is evaluated and plotted.
#' @param line_alpha Alpha value of the plotted lines.
#' @param group A variable (column name of \code{colData(se)}) by which to group samples.
#' @param range A numeric vector of length two specifying the interval of values to consider for the ECDFs. 
#' Defaults to the maximal range of the values in the assay.
#' @details
#' The function plots for each sample the ECDF of its values. This allows to identify outlier samples.
#' @return A ggplot object of the ggplot2 package.
#' @examples
#' \donttest{
#' library("DESeq2")
#' vsd <- vst(T47D)
#' plot_sample_ecdfs(vsd, group="treatment")
#' }
#' @export
plot_sample_ecdfs <- function(se, se_assay = 1, n_steps = 100, line_alpha = .2, group = NULL, range = NULL){
  
  mat <- assay(se, i = se_assay)
  
  if(!is.null(group)){
    if(length(group)>1){
      stop("Variable group must have length 1.")
    }
    if(!group %in% colnames(colData(se))){
      stop("'",group, "' is not a colum of the colData.")
    }
  }
  
  if(!is.null(range)){
    if(!length(range)==2 & is.numeric(range)){
      stop("Variable range must be a numeric vector of length 2.")
    }
    if(range[1] >= range[2]){
      stop("range[1] must be less than range[2]")
    }
    min_val <- range[1]
    max_val <- range[2]
  } else {
    min_val <- min(mat, na.rm = T)
    max_val <- max(mat, na.rm = T)    
  }
  
  plot <- ggplot(
    data = data.frame(x = seq(min_val, max_val, length.out = n_steps)),
    mapping = aes(x = x)
  ) +
    cowplot::theme_cowplot()
  
  if(!is.null(group)){
    plot <- plot + labs(color = group)
  }
  
  walk(colnames(mat), function(col){
    
    # subset to values in the interval [min_val, max_val]
    values <- mat[,col]
    values <- values[min_val <= values & values <= max_val & !is.na(values)]
    ecdf_fun = ecdf(values)
    
    if(is.null(group)){
      plot <<- plot + geom_function(
        fun = ecdf_fun, 
        alpha = line_alpha
      )  
    } else {
      plot <<- plot + geom_function(
        fun = ecdf_fun, 
        alpha = line_alpha,
        mapping = aes(color = colData(se)[col, group])
      )
    }
    
  })
  
  plot
}

#' Plot one boxplot per sample
#' @param se A SummarizedExperiment
#' @param se_assay Name or the number of the assay to plot
#' @param group Either a variable (column name of \code{colData(se)}) by which to group samples or a
#' single number indicating how many samples to plot per group. Each group will be placed in a separate plot. 
#' This is useful if the SummarizedExperiment contains more than 100 samples.
#' @param ... Other arguments passed on to \code{aes()} inside \code{ggplot()}.
#' @details
#' The function plots for each sample a bloxplot of its values. For more than 100 samples, individual boxplots
#' are hard to see and it's better to use the \code{group} parameter in that case, e.g. \code{group=100}.
#' @return A ggplot object of the ggplot2 package.
#' @examples
#' \donttest{
#' library("DESeq2")
#' vsd <- vst(T47D)
#' plot_sample_boxplots(vsd, group = "treatment")
#' plot_sample_boxplots(vsd, group = "mutation", color = treatment)
#' }
#' @return Either a list of ggplot objects of the ggplot2 package, or, if the list is length 1, its first element.
#' @export
plot_sample_boxplots <- function(se, se_assay = 1, group = NULL, ...){
  
  if(!is.null(group)){
    if(length(group)>1){
      stop("Variable group must have length 1.")
    }
    if(is.numeric(group)){
      # list of column indices to appear in the same plot
      idx <- seq_len(ncol(se))
      idx_list <- split(idx, (idx-1) %/% group)
    }
    if(is.character(group)){
      if(!group %in% colnames(colData(se))){
        stop("'",group, "' is not a colum of the colData.")
      }
      lvl <- unique(se[[group]])
      idx_list <- map(lvl, function(l){
        which(se[[group]]==l)
      }) %>% 
        `names<-`(lvl)  
    }
  } else {
    idx_list <- list(seq_len(ncol(se)))
  }
  
  plots <- imap(idx_list, function(idx, name){
    
    plot <- ggplot(mapping = aes(x = x, y = y, ...)) +
      cowplot::theme_cowplot() +
      labs(x = NULL) +
      theme(
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
      )
    
    walk(idx, function(col){
      col_name <- colnames(se)[col]
      
      plot <<- plot + 
        geom_boxplot(
          data = data.frame(
            x = col_name,
            y = assay(se, i = se_assay)[,col]
          ) %>% 
            left_join(
              as_tibble(colData(se), rownames = "col_name"),
              by = c("x" = "col_name")
            )
        )    
    })
    
    plot
  })

  if(length(plots)==1){
    plots[[1]]
  } else {
    plots
  }
}