% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PlotClusters.R
\name{PlotClusters}
\alias{PlotClusters}
\title{Plot clustering result from the BayesNonparCluster estimation.}
\usage{
PlotClusters(
  df = NULL,
  celltype = NULL,
  Assign_obj = NULL,
  mode = "segment",
  consensus = F,
  maxv = 2,
  fontsize = 10,
  fontsize_row = 10,
  fontsize_col = 10,
  allele = F,
  show_annot = T
)
}
\arguments{
\item{df}{A matrix/data.frame with each row being a cell and each column being a segment (from the WES/WGS). The values are the estimated fold change for each cell in each reion (from EstRegionCov).}

\item{celltype}{A matrix with 2+ columns: COL1- cell barcodes; COL2+- cell types and other information for plotting.}

\item{mode}{Character for the mode of heatmap to plot. The value should be one of \link{"segment","genome"}."segment": plot each segment with equal length; "genome": genome view for each segment.}

\item{maxv}{Numeric. Set the ceiling number for plotting.}

\item{Zest}{Estimated cell cluster identity from the BayesNonparCluster estimation.The order should be the same as that of the input matrix.}
}
\value{
A heatmap showing the result from the Bayesian non-parametric clustering. Each row is a cell and each column is a region. The values are the coverage change.
}
\description{
Plot clustering result from the BayesNonparCluster estimation.
}
