% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MCMCtrim.R
\name{MCMCtrim}
\alias{MCMCtrim}
\title{Trimming and thinning the MCMC interations.}
\usage{
MCMCtrim(clusteringObj = NULL, burnin = NULL, thinning = 1, allele = F)
}
\arguments{
\item{clusteringObj}{A list object. Output of the "BayesNonparCluster" function.}

\item{burnin}{Integer. Number of iterations to throw away for the MCMC process.}

\item{thinning}{Integer. Subsampling every integer number of iterations.}

\item{allele}{Logical(TRUE/FALSE). Whether the data is allele-specific or not.}
}
\value{
A list with trimmed "results" and "priors"/
"results" contains Zall: estimated subclone identity of each cell (column) for each iteration (row); Uall: estimated mean values for each subclone and each iteration;
sigma_all: estimated standard deviation of each region (column) for each iteration (row); Likelihood: total data likelihood for each iteration.
"priors" stores the prior values used in the process.
}
\description{
Trimming and thinning the MCMC interations.
}
