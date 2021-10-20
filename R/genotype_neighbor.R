#' Genotype each cell based on coverage and major haplotype proportion (theta_hat) in a region.
#'
#' @param X: A ncell by 2 dataframe. Column 1: normalized coverage (rho_hat); Column 2: theta_hat
#' @noRd
#' @keywords internal
#' @export
genotype_neighbor=function(X, cluster_name=F, maxcp=3){

  mu=NULL
  for(ii in 1:maxcp){
    mu_tmp=c(rep(0.5*ii,ii+1),c(0,( (1:(ii))/ii)))
    mu=rbind(mu, matrix(mu_tmp, byrow=F, ncol=2))
  }

  center_dict=1:nrow(mu)
  names(center_dict)=paste0("rho",mu[,1],"_theta", round(mu[,2],2))

  d_matrix=matrix(nrow=dim(X)[1], ncol=dim(mu)[1])
  for(ii in 1:dim(mu)[1]){
    dd=apply(X,1,function(ss) sqrt(sum((ss-mu[ii,])^2)))
    d_matrix[,ii]=dd}

  cluster=apply(d_matrix,1, function(ss) which(ss==min(ss))[1])
  sel_cluster=as.numeric(names(table(cluster))[which(table(cluster)>0)])

  ## update
  mu=mu[sel_cluster,, drop=F]
  d_matrix=matrix(nrow=dim(X)[1], ncol=dim(mu)[1])
  for(ii in 1:dim(mu)[1]){
    dd=apply(X,1,function(ss) sqrt(sum((ss-mu[ii,])^2)))
    d_matrix[,ii]=dd}
  cluster=apply(d_matrix,1, function(ss) which(ss==min(ss))[1])
  names(cluster)=paste0("rho",mu[cluster,1],"_theta", round(mu[cluster,2],2))
  if(cluster_name==FALSE){
  cluster_ind=center_dict[names(cluster)]
  }else{
  cluster_ind=names(cluster)
  }
  return(cluster_ind)
}
