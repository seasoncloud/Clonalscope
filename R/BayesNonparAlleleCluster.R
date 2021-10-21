#' Cluster cells based on a nested Bayesian non-parametric model (on the estimated coverage change for each region) using Gibb's sampling.
#'
#' @param Xir A matrix/data.frame with each row being a cell and each column being a segment (from the WES/WGS). The values are the estimated fold change (cov) and allelic ratio (allele) for each cell in each region.
#' @param cna_states_WGS A numeric vector with the same length as the number of rows of Xir. Each value indicates the estimated allele specific copy number states from WGS/WES data (4 means no coverage or allelic change).
#' @param alpha The concentration parameter used in the nested Chinese restaurant process controlling the probability that a new CNV state is sampled other than the existing ones for each region.
#' @param beta The concentration parameter used in the nested Chinese restaurant process controlling the probability that a new subclone is sampled for each cell.
#' @param niter Integer. Number of iteration for the Gibb's sampling.
#' @param sigmas0 A numeric vector with the same length as the number of columns of Xir (or length 1). The prior standard deviation for the normal distribution for both coverage change and allelic ratio.
#' @param U0 A matrix/data.frame with each row being a cluster and each column being a region (ncol(U0) should be the same as ncol(Xir)). The values represent prior mean of relatvie coverage for each cluster.
#' @param Z0 A numeric vector indicating the cluster identity for each cell. The order show be the same as that of Xir's row. The values should represent the order of the cluster the same as that of U0's row.
#' @param maxcp Integer. Setting the maximum number of copies for the analysis.
#' @param seed Integer. Set the seed for the stochastic process.
#'
#' @return A list with "results" and "priors"/
#' "results" contains Zall: estimated subclone identity of each cell (column) for each iteration (row); Uall: estimated mean values for each subclone and each iteration;
#' sigma_all: estimated standard deviation of each region (column) for each iteration (row); Likelihood: total data likelihood for each iteration.
#' "priors" stores the prior values used in the process.
#'
#' @import amap
#' @export
BayesNonparAlleleCluster=function(Xir_cov=NULL,Xir_allele=NULL, cna_states_WGS=NULL,alpha=1, beta=1, niter=200, sigmas0_cov=NULL, sigmas0_allele=NULL, U0=NULL, Z0=NULL, maxcp=6,seed=200){
  library(amap)
  # cna_states_WGS=U[2,]
  # set values
  #mu=matrix(c(c(0.5,0), c(0.5,1),c(1,0),c(1,0.5),c(1,1),c(1.5,0),c(1.5,0.33),c(1.5, 0.66),c(1.5, 1),c(2, 0), c(2, 0.25), c(2,0.5),c(2,0.75),c(2,1),c(2.5, 0),c(2.5,0.2),c(2.5,0.4), c(2.5, 0.6), c(2.5, 0.8), c(2.5, 1), c(3,0),c(3,1/6),c(3,2/6),c(3,3/6),c(3,4/6),c(3,5/6),c(3,1))
  #          , byrow=T, ncol=2)
  mu=NULL
  for(ii in 1:maxcp){
    mu_tmp=c(rep(0.5*ii,ii+1),c(0,( (1:(ii))/ii)))
    mu=rbind(mu, matrix(mu_tmp, byrow=F, ncol=2))
  }

  N=nrow(Xir_cov)
  R=ncol(Xir_cov)
  Xir_cov=as.matrix(Xir_cov)
  Xir_allele=as.matrix(Xir_allele)
  cna_states_WGS=as.numeric(cna_states_WGS)


  # priors
  if(is.null(U0)){
    if(!is.null(cna_states_WGS)){
      U0=matrix(c(rep(4,R), cna_states_WGS), byrow=T, ncol=R)
    }else{
      U0=matrix(c(rep(4,R), byrow=T, ncol=R))
      cna_states_WGS=U0[min(2, nrow(U0)),]
    }
  }else if(ncol(U0)==R){
    U0=U0
    U0=U0[1:max(which(!is.na(U0[,1]))),, drop=F]
    cna_states_WGS=U0[min(2, nrow(U0)),]
  }else{
    stop("Please specify a U0 matrix with ncol the same as that of Xir!")
  }

  U0c=matrix(mu[U0,1],byrow=F, ncol=ncol(U0))
  U0a=matrix(mu[U0,2],byrow=F, ncol=ncol(U0))
  U0a[U0a==0]=0.05; U0a[U0a==1]=0.95

  if(!is.null(sigmas0_cov)){
    if(length(sigmas0_cov)==1){
      sigmas0_cov=rep(sigmas0_cov, R)
    }else if(length(sigmas0_cov)==R){
      sigmas0_cov=sigmas0_cov
    }else{
      stop("Please specify sigmas0_cov with length R or length 1.")
    }
  }else{
    sigmas0_cov=rep(0.3,R)
  }

  if(!is.null(sigmas0_allele)){
    if(length(sigmas0_allele)==1){
      sigmas0_allele=rep(sigmas0_allele, R)
    }else if(length(sigmas0_allele)==R){
      sigmas0_allele=sigmas0_allele
    }else{
      stop("Please specify sigmas0_allele with length R or length 1.")
    }
  }else{
    sigmas0_allele=rep(0.15,R)
    # sdtheta=sqrt(theta.kjs[k,notNA]*(1-theta.kjs[k,notNA])/snpCoverages[i,notNA])
    # cloneLogLikelihoods[i,k]=sum(-log(sqrt(2*pi)*sigma.rho)-(rhohats[i,notNA]-rho.kjs[k,notNA])^2/(2*sigma.rho^2))+
    #   sum(-log(sqrt(2*pi)*sdtheta)-(thetahats[i,notNA]-theta.kjs[k,notNA])^2/(2*sdtheta^2))
  }

  if(is.null(Z0)){
    PP=matrix(0, ncol=nrow(U0), nrow=nrow(Xir_cov))
    for(ii in 1:nrow(U0)){
      pp=apply(Xir_cov,1, function(x) sum(dnorm(x[!is.na(x)], mu[U0[ii,],1][!is.na(x)],sigmas0_cov, log = T)))+apply(Xir_allele,1, function(x) sum(dnorm(x[!is.na(x)], mu[U0[ii,],2][!is.na(x)],sigmas0_allele, log = T)))
      PP[,ii]=pp
    }
    Z0=apply(PP, 1, function(x) which.max(x))
  }else if(length(Z0)!=N){
    stop("Please specify a vector with length the same as the number of cells!")
  }else{
    Z0[is.na(Z0)]=0
  }


  # # priors
  # U0=matrix(c(rep(4,R), cna_states_WGS), byrow=T, ncol=R)
  # P0=apply(Xir_cov,1, function(x) sum(dnorm(x[!is.na(x)], mu[U0[1,],1][!is.na(x)],0.3, log = T)))+apply(Xir_allele,1, function(x) sum(dnorm(x[!is.na(x)], mu[U0[1,],2][!is.na(x)],0.3, log = T)))
  # P1=apply(Xir_cov,1, function(x) sum(dnorm(x[!is.na(x)], mu[U0[2,],1][!is.na(x)],0.3, log = T)))+apply(Xir_allele,1, function(x) sum(dnorm(x[!is.na(x)], mu[U0[2,],2][!is.na(x)],0.3, log = T)))
  # Z0=rep(1, N)
  # Z0[which(P0<P1)]=2





  colnames(U0)=names(sigmas0_cov)=names(sigmas0_allele)=paste0("R", 1:R)
  rownames(U0)=paste0("Cluster", 1:nrow(U0))


  # initialization
  L0=sum(sapply(1:N, function(x) sum(dnorm(Xir_cov[x,][!is.na(Xir_cov[x,])], U0c[Z0, , drop=F][x,which(!is.na(Xir_cov[x,]))],sd=sigmas0_cov[!is.na(Xir_cov[x,])], log = T))))+
    sum(sapply(1:N, function(x) sum(dnorm(Xir_allele[x,][!is.na(Xir_allele[x,])], U0a[Z0, , drop=F][x,which(!is.na(Xir_allele[x,]))],sd=sigmas0_allele[!is.na(Xir_allele[x,])], log = T))))

  LL=numeric(niter)
  Uall=vector("list",niter)
  Ucov=vector("list",niter)
  Uallele=vector("list", niter)
  sigma_cov_all=vector("list", niter)
  sigma_allele_all=vector("list", niter)
  Zall=matrix(nrow=niter, ncol=N)
  Ut=U0
  Zt=Z0
  sigmas_cov=sigmas0_cov
  sigmas_allele=sigmas0_allele
  kt=nrow(Ut)


  npoints_percluster=as.numeric(table(Zt))
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
          #tmpr=runif(1, min=0.3, max=3)
          tmpr=sample((1:nrow(mu))[which(! (1:nrow(mu)) %in% Ut[,rr])],1)
          mu_rk=c(Ut[,rr],tmpr)
          indr=sample.int(kt+1,1, replace = F, prob=c(npoints_percluster, alpha))
          mu_r_new=mu_rk[indr]
          mu_new[rr]=mu_r_new
        }
        if(!any(apply(Ut,1, function(x) identical(x,mu_new)))){again=FALSE}
      }
      npoints_percluster[Zt[ii]]=npoints_percluster[Zt[ii]]-1
      #p_k= apply(Ut,1, function(x) prod(dnorm(Xir[ii,],x , log = F)))*npoints_percluster

      # p_k= sapply(1:nrow(Ut), function(x) prod(dnorm(Xir_cov[ii,],mean=mu[Ut[x,],1], sd=sigmas_cov , log = F)))*
      #    sapply(1:nrow(Ut), function(x) prod(dnorm(Xir_allele[ii,],mean=mu[Ut[x,],2], sd=sigmas_allele , log = F)))*
      #    npoints_percluster
      # p_new=prod(dnorm(Xir_cov[ii,], mean=mu[mu_new,1],sd=sigmas_cov, log=F))*
      #  prod(dnorm(Xir_allele[ii,], mean=mu[mu_new,2],sd=sigmas_allele, log=F))*beta


      p_k0= sapply((1:nrow(Ut))[which(npoints_percluster!=0)], function(x) sum(dnorm(Xir_cov[ii,][!is.na(Xir_cov[ii,])],mean=mu[Ut[x,],1][!is.na(Xir_cov[ii,])], sd=sigmas_cov[!is.na(Xir_cov[ii,])] , log = T)))+
        sapply((1:nrow(Ut))[which(npoints_percluster!=0)], function(x) sum(dnorm(Xir_allele[ii,][!is.na(Xir_allele[ii,])],mean=mu[Ut[x,],2][!is.na(Xir_allele[ii,])], sd=sigmas_allele[!is.na(Xir_allele[ii,])] , log = T)))+
        log(npoints_percluster[which(npoints_percluster!=0)])
      p_new0=sum(dnorm(Xir_cov[ii,][!is.na(Xir_cov[ii,])], mean=mu[mu_new,1][!is.na(Xir_cov[ii,])],sd=sigmas_cov[!is.na(Xir_cov[ii,])], log=T))+
        sum(dnorm(Xir_allele[ii,][!is.na(Xir_allele[ii,])], mean=mu[mu_new,2][!is.na(Xir_allele[ii,])],sd=sigmas_allele[!is.na(Xir_allele[ii,])], log=T))+log(beta)


      p_k1=sapply(p_k0, function(x) exp(x-max(c(p_k0,p_new0))))
      p_k=rep(0, nrow(Ut));
      p_k[which(npoints_percluster!=0)]=p_k1

      p_new=exp(p_new0-max(c(p_k0,p_new0)))


      Zt[ii]=sample.int(kt+1, 1, replace = F, prob=c(p_k, p_new))
      #Zt[ii]=sample.int(kt, 1, replace = F, prob=c(p_k))
      if(Zt[ii]>kt){
        npoints_percluster=c(npoints_percluster, 1)
        Ut=rbind(Ut, mu_new)
      }else{
        npoints_percluster[Zt[ii]]=npoints_percluster[Zt[ii]]+1
      }


      kt=length(npoints_percluster)
      #print(ii)
    }
    Ut=matrix(0,nrow=kt, ncol=R)
    Utc=matrix(0,nrow=kt, ncol=R)
    Uta=matrix(0,nrow=kt, ncol=R)
    for(kk in (1:kt)[which(npoints_percluster!=0)]){
      meanXcov=colMeans(Xir_cov[which(Zt==kk),, drop=F], na.rm = T)
      meanXcov[is.na(meanXcov)]=0
      meanXallele=colMeans(Xir_allele[which(Zt==kk),, drop=F], na.rm = T)
      meanXallele[is.na(meanXallele)]=0
      if(!all(meanXcov==0)){
        muXcov=rnorm(n = R,meanXcov, sigmas_cov/sqrt(sum(Zt==kk)))
        muXcov[which(muXcov<0)]=0
        muXallele=rnorm(n = R,meanXallele, sigmas_allele/sqrt(sum(Zt==kk)))
        muXallele[which(muXallele<0)]=0
      }else{
        muXcov=rep(0, R)
        muXallele=rep(0, R)
      }
      gt=sapply(1:ncol(Xir_cov), function(x) genotype_neighbor(cbind(muXcov[x], muXallele[x]),maxcp = 6 ))
      Ut[kk,]=gt
      Utc[kk,]=muXcov
      Uta[kk,]=muXallele

    }
    #sigmas_cov=sqrt(1/(N)*colSums(((Xir_cov-Utc[Zt,])^2)))
    #sigmas_allele=sqrt(1/(N)*colSums((Xir_allele-Uta[Zt,])^2))
    ncov=colSums(!is.na(Xir_cov))
    nallele=colSums(!is.na(Xir_allele))
    sigmas_cov=sqrt(1/(ncov)*colSums(((Xir_cov-Utc[Zt,])^2), na.rm = T))
    sigmas_allele=sqrt(1/(nallele)*colSums((Xir_allele-Uta[Zt,])^2, na.rm = T))

    sigmas_cov[which(sigmas_cov==0)]=0.01
    sigmas_allele[which(sigmas_allele==0)]=0.01
    colnames(Ut)=colnames(Utc)=colnames(Uta)=names(sigmas_cov)=names(sigmas_allele)=paste0("R", 1:R)
    rownames(Ut)=rownames(Utc)=rownames(Uta)=paste0('Cluster',1:nrow(Ut))
    Uall[[tt]]=Ut
    Uallele[[tt]]=Uta
    Ucov[[tt]]=Utc
    Zall[tt,]=Zt
    sigma_allele_all[[tt]]=sigmas_allele
    sigma_cov_all[[tt]]=sigmas_cov



    LL[tt]=sum(sapply(1:N, function(x) sum(dnorm(Xir_cov[x,][!is.na(Xir_cov[x,])], Utc[Zt, , drop=F][x,which(!is.na(Xir_cov[x,]))],sd=sigmas_cov[!is.na(Xir_cov[x,])], log = T))))+
      sum(sapply(1:N, function(x) sum(dnorm(Xir_allele[x,][!is.na(Xir_allele[x,])], Uta[Zt, , drop=F][x,which(!is.na(Xir_allele[x,]))],sd=sigmas_allele[!is.na(Xir_allele[x,])], log = T))))
    print(paste0("Iteration:", tt))
    print(paste0("nclusters=",nrow(Ut)))
    print(table(Zt))
  }

  ## update the results for the last iteration for stability
  Ut=matrix(0,nrow=kt, ncol=R)
  Utc=matrix(0,nrow=kt, ncol=R)
  Uta=matrix(0,nrow=kt, ncol=R)
  for(kk in (1:kt)[which(npoints_percluster!=0)]){
    meanXcov=colMeans(Xir_cov[which(Zt==kk),, drop=F], na.rm = T)
    meanXcov[is.na(meanXcov)]=0
    meanXallele=colMeans(Xir_allele[which(Zt==kk),, drop=F], na.rm = T)
    meanXallele[is.na(meanXallele)]=0

    gt=sapply(1:ncol(Xir_cov), function(x) genotype_neighbor(cbind(meanXcov[x], meanXallele[x]),maxcp = 6 ))
    Ut[kk,]=gt
    Utc[kk,]=meanXcov
    Uta[kk,]=meanXallele

  }

  ncov=colSums(!is.na(Xir_cov))
  nallele=colSums(!is.na(Xir_allele))
  sigmas_cov=sqrt(1/(ncov)*colSums(((Xir_cov-Utc[Zt,])^2), na.rm = T))
  sigmas_allele=sqrt(1/(nallele)*colSums((Xir_allele-Uta[Zt,])^2, na.rm = T))

  sigmas_cov[which(sigmas_cov==0)]=0.01
  sigmas_allele[which(sigmas_allele==0)]=0.01
  colnames(Ut)=colnames(Utc)=colnames(Uta)=names(sigmas_cov)=names(sigmas_allele)=paste0("R", 1:R)
  rownames(Ut)=rownames(Utc)=rownames(Uta)=paste0('Cluster',1:nrow(Ut))
  Uall[[tt]]=Ut
  Uallele[[tt]]=Uta
  Ucov[[tt]]=Utc
  Zall[tt,]=Zt
  sigma_allele_all[[tt]]=sigmas_allele
  sigma_cov_all[[tt]]=sigmas_cov



  LL[tt]=sum(sapply(1:N, function(x) sum(dnorm(Xir_cov[x,][!is.na(Xir_cov[x,])], Utc[Zt, , drop=F][x,which(!is.na(Xir_cov[x,]))],sd=sigmas_cov[!is.na(Xir_cov[x,])], log = T))))+
    sum(sapply(1:N, function(x) sum(dnorm(Xir_allele[x,][!is.na(Xir_allele[x,])], Uta[Zt, , drop=F][x,which(!is.na(Xir_allele[x,]))],sd=sigmas_allele[!is.na(Xir_allele[x,])], log = T))))



  rownames(Zall)=paste0("Iter",1:niter)
  colnames(Zall)=paste0("Cell", 1:N)


  results=list(Zall=Zall, Uall=Uall,Ucov=Ucov,Uallele=Uallele,
               sigma_cov_all=sigma_cov_all,sigma_allele_all=sigma_allele_all,
               Likelihood=LL)
  priors=list(Z0=Z0, U0=U0, sigmas0_cov=sigmas0_cov, sigmas0_allele=sigmas0_allele, alpha=alpha, beta=beta, cna_states_WGS=cna_states_WGS, L0=L0)
  data=list(Xir_cov=Xir_cov, Xir_allele=Xir_allele)
  return(list(results=results, priors=priors, data=data))

}

