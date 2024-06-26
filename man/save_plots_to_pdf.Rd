% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/other.R
\name{save_plots_to_pdf}
\alias{save_plots_to_pdf}
\title{Save list of plots to PDF}
\usage{
save_plots_to_pdf(
  plots,
  file = "plots.pdf",
  ncol,
  nrow,
  subfig_width = subfig_height * 16/9,
  subfig_height = 2.5,
  legend_position = "original"
)
}
\arguments{
\item{plots}{List of plots that is passed to the \code{plotlist} argument of \code{cowplot::plot_grid}}

\item{file}{file where the plots are saved}

\item{ncol}{number of columns per page for the grid of plots}

\item{nrow}{number of rows per page for the grid of plots}

\item{subfig_width}{width of a plot of the grid in inches}

\item{subfig_height}{height of a plot of the grid in inches}

\item{legend_position}{either 'original' if the original legend of each sub-plot is shown, 'none', if no legend should be shown in any of the sub-plots, 'bottom',
if no legend should be shown in the sub plots and one shared legend at the bottom or 'right', which is same as 'bottom', but shown on the right}
}
\value{
The function returns nothing but is called for it's side effect, which is to save a pdf of plots to the filesystem.
}
\description{
This function takes a list of plots as input and makes a pdf with \code{ncol} x \code{nrow} plots per page.
}
\examples{
\donttest{
library("ggplot2")
manuf <- unique(mpg$manufacturer)
plots <- lapply(manuf, function(x){
  df <- mpg[mpg$manufacturer==x,]
  ggplot(df, aes(cty, hwy)) +
    geom_point() +
    labs(title=x)
})
save_plots_to_pdf(plots, ncol=3, nrow=2)
}

}
