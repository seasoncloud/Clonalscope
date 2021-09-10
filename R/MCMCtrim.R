#' Trimming and thinning the MCMC interations.
#'
#' @param clusteringObj A list object. Output of the "BayesNonparCluster" function.
#' @param burnin Integer. Number of iterations to throw away for the MCMC process.
#' @param thinning Integer. Subsampling every integer number of iterations.
#'
#' @return A list with trimmed "results" and "priors"/
#' "results" contains Zall: estimated subclone identity of each cell (column) for each iteration (row); Uall: estimated mean values for each subclone and each iteration;
#' sigma_all: estimated standard deviation of each region (column) for each iteration (row); Likelihood: total data likelihood for each iteration.
#' "priors" stores the prior values used in the process.
#'
#' @export
MCMCtrim=function(clusteringObj=NULL, burnin=NULL, thinning=1){

  if(is.null(burnin)){
    burnin=min(100, nrow(clusteringObj$results$Zall)/2)
  }

  if(burnin>=nrow(clusteringObj$results$Zall)){
    stop("burnin values is too large! Please specify a number < number of iterations.")
  }

  # remove burn-in
  results=clusteringObj$results
  results$Zall=results$Zall[-c(1:burnin),]
  results$Uall=results$Uall[-c(1:burnin)]
  results$sigma_all=results$sigma_all[-c(1:burnin)]
  results$Likelihood=results$Likelihood[-c(1:burnin)]

   # remove thinning
  results$Zall=results$Zall[((0:(ceiling(nrow(results$Zall)/thinning)-1))*thinning+1),]
  results$Uall=results$Uall[(0:(ceiling(nrow(results$Zall)/thinning)-1))*thinning+1]
  results$sigma_all=results$sigma_all[(0:(ceiling(nrow(results$Zall)/thinning)-1))*thinning+1]
  results$Likelihood=results$Likelihood[((0:(ceiling(nrow(results$Zall)/thinning)-1))*thinning+1)]

  clusteringObj$results=results
  clusteringObj$trim=list(burnin=burnin, thinning=thinning)

  message("Object was succefully trimmed and thinned!")

  return(clusteringObj)
}
