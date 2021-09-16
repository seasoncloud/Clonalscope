AssignCluster=function(cluster_obj=NULL, mincell=20, cuttoff=0.4){

  ## set values
  R=ncol(cluster_obj$priors$U0)
  Zest=cluster_obj$results$Zall
  Uall=cluster_obj$results$Uall
  Xir=cluster_obj$data
  N=nrow(Xir)

  ## take majority vote
  set.seed(100)
  maj_vote=apply(Zest, 2, function(x) which.max(table(x)[as.character(1:max(as.numeric(cluster_obj$results$Zall)))]))
  Zest=maj_vote

  ## filter out clusters < mincells
  ind_clusters=as.numeric(names(table(Zest))[which(table(Zest)>mincell)])
  ind_cell=which(!Zest %in% ind_clusters)
  Zest[ind_cell]=NA

  ## compute new means
  Usub=matrix(nrow=nrow(Uall[[length(Uall)]]), ncol=R)
  for(kk in ind_clusters){
    meanX=colMeans(Xir[which(Zest==kk),, drop=F],  na.rm = T)
    meanX[is.na(meanX)]=0
    Usub[kk,]=meanX
  }

  sigmas=sqrt(1/(N-length(ind_cell))*colSums((Xir[-ind_cell,]-Usub[Zest[-ind_cell],])^2))

  ## Reassign cells not in the current clusters to one of the current cluster
  for(ii in ind_cell){
    P=rep(NA, nrow(Usub))
    P[ind_clusters]=apply(Usub[ind_clusters,],1, function(x) sum(dnorm(Xir[ii,], x, sigmas, log=T)))
    P[ind_clusters]=P[ind_clusters]-max(P[ind_clusters])
    #Zest[ii]=which.max(P)
    P[is.na(P)]=min(P, na.rm = T)-100
    P=exp(P)
    Zest[ii]=sample.int(nrow(Usub), 1, replace = F, prob=P)
  }

  ## Recompute new means and variance
  Usub=matrix(nrow=nrow(Uall[[length(Uall)]]), ncol=R)
  for(kk in ind_clusters){
    meanX=colMeans(Xir[which(Zest==kk),, drop=F],  na.rm = T)
    meanX[is.na(meanX)]=0
    Usub[kk,]=meanX
  }
  colnames(Usub)=paste0("R",1:ncol(Usub))
  rownames(Usub)=paste0("Cluster", 1:nrow(Usub))

  sigmas=sqrt(1/(N)*colSums((Xir-Usub[Zest,])^2))

  ## order clusters and classify normal/tumor cells
  corrs=apply(Usub, 1, function(x) cor(as.numeric(cluster_obj$priors$U0[2,]),x))
  corrs[is.na(corrs)]=0
  annot=as.numeric(corrs)
  annot[annot<=cutoff]="N"
  annot[annot!="N"]=paste0("T")

  message("Succeed!")

  return(list(Zest=Zest, corrs=corrs, annot=annot,Uest=Usub, sigmas_est=sigmas, Zmaj= maj_vote))


}
