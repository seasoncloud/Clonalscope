#' Cluster cells based on a nested Bayesian non-parametric model (on the estimated coverage change for each region) using Gibb's sampling.
#'
#' @param Xir A matrix/data.frame with each row being a cell and each column being a segment (from teh WES/WGS). The values are the estimated fold change for each cell in each reion.
#' @param cna_states_WGS A numeric vector with the same length as the number of rows of Xir. Each value indicates the estimated states on WGS/WES data (1 means no coverage change).
#' @param alpha The concentration parameter used in the nested Chinese restaurant process controlling the probability that a new CNV state is sampled other than the existing ones for each region.
#' @param beta The concentration parameter used in the nested Chinese restaurant process controlling the probability that a new subclone is sampled for each cell.
#' @param niter Integer. Number of iteration for the Gibb's sampling.
#' @param sigmas0 A numeric vector with the same length as the number of columns of Xir. The prior standard deviation for the normal distributoion.
#'
#' @return A list with "results" and "priors"/
#' "results" contains Zall: estimated subclone identity of each cell (column) for each iteration (row); Uall: estimated mean values for each subclone and each iteration;
#' sigma_all: estimated standard deviation of each region (column) for each iteration (row); Likelihood: total data likelihood for each iteration.
#' "priors" stores the prior values used in the process.
#'
#' @export
BayesNonparCluster=function(Xir=NULL,cna_states_WGS=NULL,alpha=0.1, beta=0.1, niter=100, sigmas0=NULL){

  #cna_states_WGS=U[2,]
  # set values
  N=nrow(Xir)
  R=ncol(Xir)
  Xir=as.matrix(Xir)

  # priors
  U0=matrix(c(rep(1,R), cna_states_WGS), byrow=T, ncol=R)
  P0=apply(Xir,1, function(x) sum(dnorm(x, U0[1,],0.3, log = T)))
  P1=apply(Xir,1, function(x) sum(dnorm(x, U0[2,],0.3, log = T)))
  Z0=rep(1, N)
  Z0[which(P0<P1)]=2

  if(length(sigmas0) !=R){
  sigmas0=rep(0.5,R)
  }
  colnames(U0)=names(sigmas0)=paste0("R", 1:R)
  rownames(U0)=paste0("Cluster", 1:nrow(U0))


  # initialization
  LL=numeric(niter)
  Uall=vector("list",niter)
  sigma_all=vector("list", niter)
  Zall=matrix(nrow=niter, ncol=N)
  Ut=U0
  Zt=Z0
  sigmas=sigmas0
  kt=nrow(Ut)


  npoints_percluster=as.numeric(table(Zt))
  #names(npoints_percluster)=paste0("c",1:2)

  #set.seed(100)
  # Gibbs sampling
  for(tt in 1:100){
    for(ii in 1:N){
      mu_new=numeric(R)
      Zi=Zt[ii]
      again=TRUE
      while(again==T){
        for(rr in 1:R){
          tmpr=runif(1, min=0.3, max=3)
          mu_rk=c(Ut[,rr],tmpr)
          indr=sample.int(kt+1,1, replace = F, prob=c(npoints_percluster, alpha))
          mu_r_new=mu_rk[indr]
          mu_new[rr]=mu_r_new
        }
        if(!any(apply(Ut,1, function(x) identical(x,mu_new)))){again=FALSE}
      }
      npoints_percluster[Zt[ii]]=npoints_percluster[Zt[ii]]-1
      #p_k= apply(Ut,1, function(x) prod(dnorm(Xir[ii,],x , log = F)))*npoints_percluster
      p_k= sapply(1:nrow(Ut), function(x) prod(dnorm(Xir[ii,],mean=Ut[x,], sd=sigmas , log = F)))*npoints_percluster
      p_new=prod(dnorm(Xir[ii,], mean=mu_new,sd=sigmas, log=F))*beta
      Zt[ii]=sample.int(kt+1, 1, replace = F, prob=c(p_k, p_new))
      #Zt[ii]=sample.int(kt, 1, replace = F, prob=c(p_k))
      if(Zt[ii]>kt){
        cname= paste0('c',kt+1)
        npoints_percluster=c(npoints_percluster, 1)
        Ut=rbind(Ut, mu_new)
      }else{
        npoints_percluster[Zt[ii]]=npoints_percluster[Zt[ii]]+1
        }


      kt=length(npoints_percluster)
    #print(ii)
      }
    Ut=matrix(nrow=kt, ncol=R)
    for(kk in 1:kt){
      meanX=colMeans(Xir[which(Zt==kk),, drop=F])
      meanX[is.na(meanX)]=0
      if(!all(meanX==0)){
      muX=rnorm(n = R,meanX, sigmas/sqrt(sum(Zt==kk)))
      muX[which(muX<0)]=0
      }else{
        muX=rep(0, R)
      }
      Ut[kk,]=muX
    }
    sigmas=sqrt(1/(N)*colSums((Xir-Ut[Zt,])^2))
    colnames(Ut)= names(sigmas)=paste0("R", 1:R)
    rownames(Ut)=paste0('Cluster',1:nrow(Ut))
    Uall[[tt]]=Ut
    Zall[tt,]=Zt
    sigma_all[[tt]]=sigmas


    LL[tt]=sum(sapply(1:N, function(x) sum(dnorm(Xir[x,], Ut[Zt, , drop=F][x,],sd=sigmas, log = T))))
    print(paste0("Iteration:", tt))
    print(paste0("nclusters=",nrow(Ut)))
    print(table(Zt))
  }

  rownames(Zall)=paste0("Iter",1:niter)
  colnames(Zall)=paste0("Cell", 1:N)

  results=list(Zall=Zall, Uall=Uall, sigma_all=sigma_all, Likelihood=LL)
  priors=list(Z0=Z0, U0=U0, sigmas0=sigmas0, alpha=alpha, beta=beta)
  return(list(results=results, priors=priors))

}
