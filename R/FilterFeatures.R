#' Filter out cycle genes and HLA genes (Gao et al., 2021)
#'
#' @param mtx A list object. Output of the "BayesNonparCluster" function.
#' @param barcodes Integer. Number of iterations to throw away for the MCMC process.
#' @param features Integer. Subsampling every integer number of iterations.
#' @param rcycle Lpgical. Remove cycle genes or not.
#' @param rHLA Logical. Remove HLA genes or not.
#' @param rMT Logical. Remove MT genes.
#'
#' @return An object with the filtered input files.
#'
#' @export
FilterFeatures=function(mtx=NULL, barcodes=NULL, features=NULL, rcycle=TRUE, rHLA=TRUE,rMT=TRUE, cyclegenes=NULL){

  if(is.null(mtx) | is.null(barcodes) | is.null(features)){
    stop("Please provide the three inputs: mtx, barcodes, and features!")
  }

## filter out genes
  if(rcycle==TRUE){
    ind_rm1=which(features[,2] %in% cyclegenes[,1])
  }else{
    ind_rm1=c()
  }

  if(rHLA==TRUE){
    ind_rm2=which(grepl("^HLA-", features[,2]))
  }else{
    ind_rm2=c()
  }

  if(rMT==TRUE){
    ind_rm3=which(grepl("^MT-", features[,2]))
  }else{
    ind_rm3=c()
  }

  ind_rm=sort(unique(c(ind_rm1, ind_rm2, ind_rm3)))
  mtx=mtx[-ind_rm,]
  features=features[-ind_rm,]


  Obj_filtered=list('mtx'=mtx, 'barcodes'=barcodes, 'features'=features, 'rcycle'=rcycle, 'rHLA'=rHLA, 'rMT'=rMT)
  message("Input files were succefully filtered!")

  return(Obj_filtered)
}
