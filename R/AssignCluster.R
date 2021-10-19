AssignCluster=function(cluster_obj=NULL, mincell=20, cutoff=0.5){

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
    P[ind_clusters]=apply(Usub[ind_clusters,, drop=F],1, function(x) sum(dnorm(Xir[ii,], x, sigmas, log=T)))
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
  priors=cluster_obj$priors$cna_states_WGS#cluster_obj$priors$U0[2,]

  if(all(priors==1)){
    priors=rep(1.5, ncol(Usub))
    #priors=Usub[as.numeric(names(table(Zest)))[which.max(table(Zest))],]
  }



  # correlation
  # while(sd(priors)==0){
  #   priors=sample(c(0.001, -0.001), length(priors), replace = T)+priors
  # }

  #ind_cna=which(priors != 1)
  corrs=apply(Usub[,,drop=F], 1, function(x)
  sum((as.numeric(x)-1)*(as.numeric(priors[])-1)))#/(sqrt(sum((as.numeric(x)-1)^2))*sqrt(sum((as.numeric(priors[])-1)^2))))
  #cor(as.numeric(priors),x)) ###c
  corrs[is.na(corrs)]=0
  corrs=corrs/max(corrs)

  #if(dist=='correlation'){
  # if(sd(priors)==0){
  #   Manhattan <- function(a, b) (sum(abs(a - b)))
  #   corrs=apply(Usub, 1, function(x) Manhattan(as.numeric(priors),x))/length(priors)
  # }else{

  ## cross-product
  # library(lsa)
  # crossproduct=function(x){
  #   if(all(is.na(x))==F){
  #     cos=cosine(priors, x)
  #     sin=sqrt(1-cos^2)
  #     #abscp=sqrt(sum(x^2))*sqrt(sum(priors2^2))*sin/sqrt(sum(priors2^2))
  #     x2=abs(x-priors)+1
  #     abscp=sqrt(sum(x2^2))/sqrt(sum(1*R))*sin
  #     # x2=x#/(sqrt(sum(x^2)))
  #     #   return(sqrt(sum((crossn(as.numeric(priors2),x2)))^2))
  #   }else{
  #     return(100)
  #   }}
  #
  # corrs=apply(Usub, 1, crossproduct)
  # corrs[is.na(corrs)]=100

  ## eucledian
 # corrs=apply(Usub, 1, function(x) sqrt(sum((priors-x)^2))/length(priors))
  #corrs[is.na(corrs)]=0
  #}
  # }else{
  #   Manhattan <- function(a, b) (sum(abs(a - b)))
  #   corrs=apply(Usub, 1, function(x) Manhattan(as.numeric(priors),x))/length(priors)
  # }
  annot=as.numeric(corrs)

  if(cutoff=='km'){
  km=kmeans(annot[which(annot!=100)],2)
  cutoff=mean(km$centers)
  }

  annot[annot>cutoff]="T"
  annot[annot!="T"]=paste0("N")

  ## label clusters with priors
  # Zest0=Zest
  # if(!is.null(rownames(cluster_obj$priors$U0))){
  # tmp=1:max(Zest)
  # tmp[1:nrow(cluster_obj$priors$U0)]=rownames(cluster_obj$priors$U0)
  # Zest=tmp[Zest]
  # }

  message("Succeed!")

  return(list(Zest=Zest, corrs=corrs[Zest], annot=annot[Zest],Uest=Usub, sigmas_est=sigmas, Zmaj= maj_vote, cutoff=cutoff, U0=cluster_obj$priors$U0))


}
