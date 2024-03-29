% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/genomePlot.R
\name{genomePlot}
\alias{genomePlot}
\title{A Function for Plotting Mean Copy Number Values and Differences Across the Genome}
\usage{
genomePlot(
  inputList,
  lwdVec = rep(1, 3),
  ltyVec = c(1:3),
  lineColorVec = c("red", "blue", "black"),
  ylimLow = -1,
  ylimHigh = 1,
  chrLabel = TRUE,
  xaxisLabel = "Chromosome",
  yaxisLabel = NULL,
  mainLabel = NULL,
  rectColors = c("light gray", "gray"),
  axisCex = 1,
  labelCex = 1,
  xaxisLine = 2.5,
  yaxisLine = 2.5,
  mainLine = 0,
  marginVec = c(4, 4, 3, 3),
  legendText = NULL,
  highThreshold = NULL,
  lowThreshold = NULL,
  showLegend = FALSE,
  legendXQuantile = 0.55,
  legendYCoord = 1
)
}
\arguments{
\item{inputList}{A list produced by dataPrep.}

\item{lwdVec}{A vector of line widths.  Default = rep(1, 3).  See \code{\link{par}}.}

\item{ltyVec}{A vector of line types.  Default = c(1:3).  See \code{\link{par}}.}

\item{lineColorVec}{A vector of line colors.  Default = c("red", "blue", "black").  See \code{\link{par}}.}

\item{ylimLow}{The lower limit of the y-values in the plot.  Default = -1.  See \code{\link{plot}}.}

\item{ylimHigh}{The upper limit of the y-values in the plot.  Default = 1.  See \code{\link{plot}}.}

\item{chrLabel}{Binary value determining whether or not chromosomes are labeled.  Default = TRUE.}

\item{xaxisLabel}{Label for the x-axis in the plot.  Default = "Chromosome".  See \code{\link{plot}}.}

\item{yaxisLabel}{Label for the y-axis in the plot.  Default = NULL.  See \code{\link{plot}}.}

\item{mainLabel}{Main label in the plot.  Default = NULL.  See \code{\link{plot}}.}

\item{rectColors}{Background colors for different chromosomes.  Default = c("light gray", "gray").}

\item{axisCex}{Point size for the scale on the axis.  Default = 1.  See \code{\link{par}}.}

\item{labelCex}{Point size for the axis label.  Default = 1.  See \code{\link{par}}.}

\item{xaxisLine}{Numerical value used to specify the location (line) of the x-axis label.  Default = 2.5.  See \code{\link{mtext}}.}

\item{yaxisLine}{Numerical value used to specify the location (line) of the y-axis label.  Default = 2.5.  See \code{\link{mtext}}.}

\item{mainLine}{Numerical value used to specify the location (line) of the main.label.  Default = 0.  See \code{\link{mtext}}.}

\item{marginVec}{Numerical vector specifying margin sizes.  Default = c(4, 4, 3, 3).  See \code{\link{par}}.}

\item{legendText}{Character vector used to legend.  Only shown if showLegend = TRUE.  Default = NULL.  See \code{\link{legend}}.}

\item{highThreshold}{Numerical value representing the position of the upper horizontal line, e.g., a threshold for
assessing statistical significance.  Default = NULL.}

\item{lowThreshold}{Numerical value representing the position of the lower horizontal line, e.g., a threshold for
assessing statistical significance.  Default = NULL.}

\item{showLegend}{Binary value determining whether or not the legend is shown.  Default = FALSE.  See \code{\link{legend}}.}

\item{legendXQuantile}{Quantile to specify the "x" location of the legend.  Only relevant if showLegend = TRUE  Default = 0.55.  See \code{\link{legend}}.}

\item{legendYCoord}{Numerical value to specify the "y"location of the legend.  Only relevant if showLegend = TRUE.  Default = 1.  See \code{\link{legend}}.}
}
\value{
Creates a genomewide plot of mean copy number values and differences.
}
\description{
This function plots mean copy number values from one or two cohorts at a common set of markers across the genome.
}
\details{
This function is used to visualize copy number values and copy number alterations
across the genome.  If Y = NULL in the input list, then the plot shows a single line
corresponding to the mean DNA copy number values based on the entries in X.
If both X and Y are specified, the plot shows three lines corresponding to the mean
DNA copy number values in X, the mean DNA copy number values in Y, and the difference
of the mean DNA copy number values.
}
\examples{
genomeChrPlot(inputList = pD, ylimLow = -1.4, ylimHigh = 1.4)

}
