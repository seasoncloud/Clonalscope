#' Cluster cells based on a nested Bayesian non-parametric model (on the estimated coverage change for each region) using Gibb's sampling.
#'
#' @param Xir A matrix/data.frame with each row being a cell and each column being a segment (from the WES/WGS). The values are the estimated fold change for each cell in each reion (from EstRegionCov).
#' @param cna_states_WGS A numeric vector with the same length as the number of rows of Xir. Each value indicates the estimated states on WGS/WES data (1 means no coverage change).
#' @param alpha The concentration parameter used in the nested Chinese restaurant process controlling the probability that a new CNV state is sampled other than the existing ones for each region.
#' @param beta The concentration parameter used in the nested Chinese restaurant process controlling the probability that a new subclone is sampled for each cell.
#' @param niter Integer. Number of iteration for the Gibb's sampling.
#' @param sigmas0 A numeric vector with the same length as the number of columns of Xir (or length 1). The prior standard deviation for the normal distributoion.
#' @param U0 A matrix/data.frame with each row being a cluster and each column being a region (ncol(U0) should be the same as ncol(Xir)). The values represent prior mean of relatvie coverage for each cluster.
#' @param Z0 A numeric vector indicating the cluster identity for each cell. The order show be the same as that of Xir's row. The values should represent the order of the cluster the same as that of U0's row.
#' @param seed Integer. Set the seed for the stochastic process.
#'
#' @return A list with "results" and "priors"/
#' "results" contains Zall: estimated subclone identity of each cell (column) for each iteration (row); Uall: estimated mean values for each subclone and each iteration;
#' sigma_all: estimated standard deviation of each region (column) for each iteration (row); Likelihood: total data likelihood for each iteration.
#' "priors" stores the prior values used in the process.
#'
#' @import amap
#' @export
BayesNonparCluster=function(Xir=NULL,cna_states_WGS=NULL,alpha=1, beta=1, niter=200, sigmas0=NULL, U0=NULL, Z0=NULL, seed=200){
  #library(amap)
  #cna_states_WGS=U[2,]
  # set values
  N=nrow(Xir)
  R=ncol(Xir)
  Xir=as.matrix(Xir)
  cna_states_WGS=as.numeric(cna_states_WGS)

  # priors
  if(is.null(U0)){
    wU0=FALSE
    if(!is.null(cna_states_WGS)){
  U0=matrix(c(rep(1,R), cna_states_WGS), byrow=T, ncol=R)
    }else{
      U0=matrix(c(rep(1,R), byrow=T, ncol=R))
      if(nrow(U0)!=1){
      cna_states_WGS=U0[which(!is.na(result$Uest[,1]))[2],]
      }else{
        cna_states_WGS=U0[1,]
      }
    }
  }else if(ncol(U0)==R){
    wU0=TRUE
    U0=U0
    U0=U0[1:max(which(!is.na(U0[,1]))),, drop=F]
    if(nrow(U0)!=1){
      cna_states_WGS=U0[which(!is.na(result$Uest[,1]))[2],]
    }else{
      cna_states_WGS=U0[1,]
    }
  }else{
    stop("Please specify a U0 matrix with ncol the same as that of Xir!")
  }


  if(!is.null(sigmas0)){
    if(length(sigmas0)==1){
      sigmas0=rep(sigmas0, R)
    }else if(length(sigmas0)==R){
      sigmas0=sigmas0
    }else{
      stop("Please specify sigmas0 with length R or length 1.")
    }
  }else{
    sigmas0=rep(0.25,R)
  }

  #weights=abs(as.numeric(apply(U0==1, 2, all))-1)+1
  weights=rep(1, ncol(U0))
  # weights_mtx=matrix(rep(weights, nrow(Xir)), byrow=T, ncol=ncol(Xir))
  # Xir=Xir*weights_mtx
  #
  # weights_mtx2=matrix(rep(weights, nrow(U0)), byrow=T, ncol=ncol(U0))
  # U0=U0*weights_mtx2


  # chr=sapply(strsplit(colnames(df2),"-"),'[',1)
  # start=as.numeric(sapply(strsplit(colnames(df2),"-"),'[',2))
  # end=as.numeric(sapply(strsplit(colnames(df2),"-"),'[',3))
  # seglen=end-start
  # P0=apply(Xir,1, function(x) sum(dnorm(x, U0[1,],0.3, log = T)))
  # P1=apply(Xir,1, function(x) sum(dnorm(x, U0[2,],0.3, log = T)))
  # Z0=rep(1, N)
  # Z0[which(P0<P1)]=2

  if(is.null(Z0)){
  PP=matrix(0, ncol=nrow(U0), nrow=nrow(Xir))
  for(ii in 1:nrow(U0)){
    pp=apply(Xir,1, function(x) sum(dnorm(x, U0[ii,],sigmas0, log = T)))
    PP[,ii]=pp
  }
  Z0=apply(PP, 1, function(x) which.max(x))
  }else if(length(Z0)!=N){
    stop("Please specify a vector with length the same as the number of cells!")
  }else{
    Z0[is.na(Z0)]=0
  }

  #U0=U0[1:max(Z0),, drop=F]


  # extreme condition
  # if(length(table(Z0))==1){
  #   dd=Dist(Xir, method='correlation')
  #   hc=hclust(dd)
  #   ct=cutree(hc, k = 2)
  #   Z0=ct
  #   U0=matrix(c(colMeans(Xir[which(ct==1),]), colMeans(Xir[which(ct==2),])), byrow=T, ncol=R)
  #   #dd=as.matrix(dd)
  # }


  colnames(U0)=names(sigmas0)=paste0("R", 1:R)
  #rownames(U0)=paste0("Cluster", 1:nrow(U0))

  # initialization
  L0=sum(sapply(1:N, function(x) sum(dnorm(Xir[x,], U0[Z0, , drop=F][x,],sd=sigmas0, log = T))))
  LL=numeric(niter)
  Uall=vector("list",niter)
  sigma_all=vector("list", niter)
  Zall=matrix(nrow=niter, ncol=N)
  Ut=U0
  Zt=Z0
  sigmas=sigmas0
  kt=nrow(Ut)


  #npoints_percluster=as.numeric(table(Zt))
  npoints_percluster=table(Zt)[as.character(1:max(Zt))]
  npoints_percluster[is.na(npoints_percluster)]=0

  #names(npoints_percluster)=paste0("c",1:2)


set.seed(seed)
  #set.seed(100)
  # Gibbs sampling
  for(tt in 1:niter){
    for(ii in 1:N){
      mu_new=numeric(R)
      Zi=Zt[ii]

      again=TRUE
      while(again==T){
        for(rr in 1:R){
          tmpr=runif(1, min=0.3, max=2.5)
          #tmpr=sample(c(0.5,1,1.5,2, 2.5),1)
          mu_rk=c(Ut[,rr],tmpr)
          indr=sample.int(kt+1,1, replace = F, prob=c(npoints_percluster, alpha))
          mu_r_new=mu_rk[indr]
          mu_new[rr]=mu_r_new
        }
        if(!any(apply(Ut,1, function(x) identical(x,mu_new)))){again=FALSE}
      }
      #mu_new=Xir[sample(1:nrow(Xir),1),, drop=F]

      # again=TRUE
      # while(again==T){
      #   indr=sample(1:R, 1, replace = F)
      #   indc=sample(which(!is.na(Ut[,1])), 1, replace = F)
      #   mu_new=Ut[indc,, drop=T]
      #   mu_new[indr]=runif(1, min=0.3, max=2.5)
      #   if(!any(apply(Ut,1, function(x) identical(x,mu_new)))){again=FALSE}
      # }


      npoints_percluster[Zt[ii]]=npoints_percluster[Zt[ii]]-1
      #p_k= apply(Ut,1, function(x) prod(dnorm(Xir[ii,],x , log = F)))*npoints_percluster
      #p_k= sapply(1:nrow(Ut), function(x) prod(dnorm(Xir[ii,],mean=Ut[x,], sd=sigmas , log = F)))*npoints_percluster
      #p_new=prod(dnorm(Xir[ii,], mean=mu_new,sd=sigmas, log=F))*beta

      p_k0= sapply((1:nrow(Ut))[which(npoints_percluster!=0)], function(x) sum(weights*dnorm(Xir[ii,],mean=Ut[x,], sd=sigmas , log = T)))+
        #sapply(1:nrow(Ut), function(x) sum(dnorm(Xir_allele[ii,],mean=mu[Ut[x,],2], sd=sigmas_allele , log = T)))+
        log(npoints_percluster[which(npoints_percluster!=0)])
      p_new0=sum(weights*dnorm(Xir[ii,], mean=mu_new,sd=sigmas, log=T))+
        #sum(dnorm(Xir_allele[ii,], mean=mu[mu_new,2],sd=sigmas_allele, log=T))+
        log(beta)

      p_k1=sapply(p_k0, function(x) exp(x-max(c(p_k0,p_new0))))
      p_k=rep(0, nrow(Ut));
      p_k[which(npoints_percluster!=0)]=p_k1

      #p_k=sapply(p_k0, function(x) exp(x-max(c(p_k0,p_new0))))
      p_new=exp(p_new0-max(c(p_k0,p_new0)))

      Zt[ii]=sample.int(kt+1, 1, replace = F, prob=c(p_k, p_new))
      #Zt[ii]=sample.int(kt, 1, replace = F, prob=c(p_k))
      if(Zt[ii]>kt){
        #cname= paste0('c',kt+1)
        npoints_percluster=c(npoints_percluster, 1)
        Ut=rbind(Ut, mu_new)
      }else{
        npoints_percluster[Zt[ii]]=npoints_percluster[Zt[ii]]+1
        }


      kt=length(npoints_percluster)
    #print(ii)
      }
    Ut=matrix(nrow=kt, ncol=R)
    for(kk in (1:kt)[which(npoints_percluster!=0)]){
      meanX=colMeans(Xir[which(Zt==kk),, drop=F],  na.rm = T)
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
    #sigmas=rep(0.3, length(sigmas))
    #sigmas=pmax(sigmas, 0.2)
    colnames(Ut)= names(sigmas)=paste0("R", 1:R)
    rownames(Ut)=paste0('Cluster',1:nrow(Ut))
    Uall[[tt]]=Ut
    Zall[tt,]=Zt
    sigma_all[[tt]]=sigmas


    LL[tt]=sum(sapply(1:N, function(x) sum(weights*dnorm(Xir[x,], Ut[Zt, , drop=F][x,],sd=sigmas, log = T))))
    print(paste0("Iteration:", tt))
    print(paste0("nclusters=",nrow(Ut)))
    print(table(Zt))
    #print(Ut[which(!is.na(Ut[,1])),c(7,8,19), drop=F])

  }

  ## update the results for the last iteration for stability
  Ut=matrix(nrow=kt, ncol=R)
  for(kk in (1:kt)[which(npoints_percluster!=0)]){
    meanX=colMeans(Xir[which(Zt==kk),, drop=F],  na.rm = T)
    meanX[is.na(meanX)]=0
    Ut[kk,]=meanX
  }

  sigmas=sqrt(1/(N)*colSums((Xir-Ut[Zt,])^2))
  colnames(Ut)= names(sigmas)=paste0("R", 1:R)
  rownames(Ut)=paste0('Cluster',1:nrow(Ut))
  Uall[[tt]]=Ut
  Zall[tt,]=Zt
  sigma_all[[tt]]=sigmas


  LL[tt]=sum(sapply(1:N, function(x) sum(weights*dnorm(Xir[x,], Ut[Zt, , drop=F][x,],sd=sigmas, log = T))))



  rownames(Zall)=paste0("Iter",1:niter)
  colnames(Zall)=rownames(Xir)

  results=list(Zall=Zall, Uall=Uall, sigma_all=sigma_all, Likelihood=LL)
  priors=list(Z0=Z0, U0=U0, sigmas0=sigmas0, alpha=alpha, beta=beta, cna_states_WGS=cna_states_WGS, L0=L0, wU0=wU0, seed=seed)
  return(list(results=results, priors=priors, data=Xir))

}
