AssignCluster=function(cluster_obj=NULL, mincell=20, cutoff=0.4, allele=FALSE, cna_states_WGS=NULL, rm_extreme=F, cap_corr=1.5){

  if(allele==FALSE){
    ## set values
    R=ncol(cluster_obj$priors$U0)
    Zest=cluster_obj$results$Zall
    Uall=cluster_obj$results$Uall
    Xir=cluster_obj$data
    N=nrow(Xir)

    #if(is.null(corr_region_ind)){
    corr_region_ind=1:R
    #}
    #print(corr_region_ind)

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

    sigmas=sqrt(1/(N-length(ind_cell))*colSums((Xir[-ind_cell,]-Usub[Zest[-ind_cell],])^2, na.rm = T))

    ## Reassign cells not in the current clusters to one of the current cluster
    for(ii in ind_cell){
      P=rep(NA, nrow(Usub))
      P[ind_clusters]=apply(Usub[ind_clusters,, drop=F],1, function(x) sum(dnorm(as.numeric(Xir[ii,]), x, sigmas, log=T), na.rm = T))
      P[ind_clusters]=P[ind_clusters]-max(P[ind_clusters], na.rm = T)
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
    #colnames(Usub)=paste0("R",1:ncol(Usub))
    colnames(Usub)=colnames(Xir)
    rownames(Usub)=paste0("Cluster", 1:nrow(Usub))

    sigmas=sqrt(1/(N)*colSums((Xir-Usub[Zest,])^2, na.rm = T))

    ## order clusters and classify normal/tumor cells
    if(is.null(cna_states_WGS)){
      priors=cluster_obj$priors$cna_states_WGS#cluster_obj$priors$U0[2,]
    }else{
      priors=cna_states_WGS
    }
    if(all(priors==1)){
      priors=rep(1.5, ncol(Usub))
      #priors=Usub[as.numeric(names(table(Zest)))[which.max(table(Zest))],]
    }

    # mu=NULL
    # for(ii in 1:6){
    #   mu_tmp=c(rep(0.5*ii,ii+1),c(0,( (1:(ii))/ii)))
    #   mu=rbind(mu, matrix(mu_tmp, byrow=F, ncol=2))
    # }





    # correlation
    # while(sd(priors)==0){
    #   priors=sample(c(0.001, -0.001), length(priors), replace = T)+priors
    # }

    #ind_cna=which(priors != 1)
    Usub2=apply(Usub,2, function(x) pmin(x, cap_corr))
    #Usub2=apply(Usub,2, function(x) pmax(x,0.8))
    # Usub2=apply(Usub, c(1,2), function(x){
    #   if(is.na(x)){
    #     return(NA)
    #   }else if(x>1){return(1.5)
    #   }else if(x<1){
    #       return(0.5)
    #   }else if(x==1){
    #       return(1)
    #   }
    # })
    #Usub2=apply(Usub,2, function(x) pmax(x,0.5))
    corrs=apply(Usub2[,corr_region_ind,drop=F], 1, function(x){
      tmp=sort((as.numeric(x)-1)*(as.numeric(priors[corr_region_ind])-1))
      if(length(tmp)==0){
        return(NA)
      }else{
        if(rm_extreme==T){
          # if(sum(tmp[2:(length(tmp)-1)])/(length(which(as.numeric(priors[corr_region_ind])!=1))-2)<0.1){
          #   warning("WGS and scRNA do not match well. annot might not be accurate!")
          # }
          # return(sum(tmp[2:(length(tmp)-1)]))
        }else{
          # if(sum(tmp)/(length(which(as.numeric(priors[corr_region_ind])!=1)))<0.1){
          #   warning("WGS and scRNA do not match well. annot might not be accurate!")
          # }
          return(sum(tmp))}}
    }
    )#/(sqrt(sum((as.numeric(x)-1)^2))*sqrt(sum((as.numeric(priors[])-1)^2))))
    #cor(as.numeric(priors),x)) ###c

    corrs[is.na(corrs)]=0
    if(max(corrs, na.rm = T)/(length(which(as.numeric(priors[corr_region_ind])!=1)))<0.1){
      warning("WGS and scRNA do not match well. annot might not be accurate!")
    }
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

    return(list(Zest=Zest, corrs=corrs[Zest], annot=annot[Zest],Uest=Usub, sigmas_est=sigmas, Zmaj= maj_vote, cutoff=cutoff, U0=cluster_obj$priors$U0,  wU0=cluster_obj$priors$wU0))


  }else{  ### allele

    R=ncol(cluster_obj$priors$U0)
    Zest=cluster_obj$results$Zall
    Uall=cluster_obj$results$Uall
    Ucov=cluster_obj$results$Ucov
    Uallele=cluster_obj$results$Uallele
    Xir_cov=cluster_obj$data$Xir_cov
    Xir_allele=cluster_obj$data$Xir_allele
    N=nrow(Xir_cov)

    #if(is.null(corr_region_ind)){
    corr_region_ind=1:R
    #}
    ## take majority vote
    set.seed(100)
    maj_vote=apply(Zest, 2, function(x) which.max(table(x)[as.character(1:max(as.numeric(cluster_obj$results$Zall)))]))
    Zest=maj_vote

    ## filter out clusters < mincells
    ind_clusters=as.numeric(names(table(Zest))[which(table(Zest)>mincell)])
    ind_cell=which(!Zest %in% ind_clusters)
    Zest[ind_cell]=NA

    ## compute new means
    Usubc=matrix(nrow=nrow(Uall[[length(Uall)]]), ncol=R)
    Usuba=matrix(nrow=nrow(Uall[[length(Uall)]]), ncol=R)
    for(kk in ind_clusters){
      meanXcov=colMeans(Xir_cov[which(Zest==kk),, drop=F],  na.rm = T)
      meanXcov[is.na(meanXcov)]=0
      Usubc[kk,]=meanXcov

      meanXallele=colMeans(Xir_allele[which(Zest==kk),, drop=F],  na.rm = T)
      meanXallele[is.na(meanXallele)]=0
      Usuba[kk,]=meanXallele
    }

    Xir_cov_tmp=t(sapply(1:nrow(Xir_cov), function(x) if(x %in% ind_cell){return(rep(FALSE, ncol(Xir_cov)))}else{return(rep(TRUE, ncol(Xir_cov)))}))
    Xir_allele_tmp=t(sapply(1:nrow(Xir_allele), function(x) if(x %in% ind_cell){return(rep(FALSE, ncol(Xir_allele)))}else{return(rep(TRUE, ncol(Xir_allele)))}))

    ncov=colSums(!is.na(Xir_cov) & Xir_cov_tmp)
    nallele=colSums(!is.na(Xir_allele) & Xir_allele_tmp)

    sigmas_cov=sqrt(1/(ncov)*colSums(((Xir_cov[-ind_cell,, drop=F]-Usubc[Zest[-ind_cell],])^2), na.rm = T))
    sigmas_allele=sqrt(1/(nallele)*colSums((Xir_allele[-ind_cell,, drop=F]-Usuba[Zest[-ind_cell],])^2, na.rm = T))




    ## Reassign cells not in the current clusters to one of the current cluster
    for(ii in ind_cell){
      P=rep(NA, nrow(Usubc))
      P[ind_clusters]=apply(Usubc[ind_clusters,, drop=F],1, function(x) sum(dnorm(Xir_cov[ii,], x, sigmas_cov, log=T), na.rm = T))+
        apply(Usuba[ind_clusters,, drop=F],1, function(x) sum(dnorm(Xir_allele[ii,], x, sigmas_allele, log=T), na.rm = T))
      P[ind_clusters]=P[ind_clusters]-max(P[ind_clusters])
      #Zest[ii]=which.max(P)
      P[is.na(P)]=min(P, na.rm = T)-100
      P=exp(P)
      Zest[ii]=sample.int(nrow(Usuba), 1, replace = F, prob=P)
    }

    ## Recompute new means and variance
    Usubc=matrix(nrow=nrow(Uall[[length(Uall)]]), ncol=R)
    Usuba=matrix(nrow=nrow(Uall[[length(Uall)]]), ncol=R)
    Usub=matrix(nrow=nrow(Uall[[length(Uall)]]), ncol=R)
    for(kk in ind_clusters){
      meanXcov=colMeans(Xir_cov[which(Zest==kk),, drop=F],  na.rm = T)
      meanXcov[is.na(meanXcov)]=0
      Usubc[kk,]=meanXcov

      meanXallele=colMeans(Xir_allele[which(Zest==kk),, drop=F],  na.rm = T)
      meanXallele[is.na(meanXallele)]=0
      Usuba[kk,]=meanXallele


      gt=sapply(1:ncol(Xir_cov), function(x) genotype_neighbor(cbind(meanXcov[x], meanXallele[x]),maxcp = 6 ))
      Usub[kk,]=gt
    }
    #colnames(Usubc)=colnames(Usuba)=colnames(Usub)=paste0("R",1:ncol(Usubc))
    colnames(Usubc)=colnames(Usuba)=colnames(Usub)=colnames(Xir)
    rownames(Usubc)=rownames(Usuba)=rownames(Usub)=paste0("Cluster", 1:nrow(Usubc))




    ncov=colSums(!is.na(Xir_cov))
    nallele=colSums(!is.na(Xir_allele))

    sigmas_cov=sqrt(1/(ncov)*colSums(((Xir_cov[-ind_cell,, drop=F]-Usubc[Zest[-ind_cell],])^2), na.rm = T))
    sigmas_allele=sqrt(1/(nallele)*colSums((Xir_allele[-ind_cell,, drop=F]-Usuba[Zest[-ind_cell],])^2, na.rm = T))


    ## order clusters and classify normal/tumor cells
    #priors=cluster_obj$priors$cna_states_WGS#cluster_obj$priors$U0[2,]
    if(is.null(cna_states_WGS)){
      priors=cluster_obj$priors$cna_states_WGS#cluster_obj$priors$U0[2,]
    }else{
      priors=cna_states_WGS
    }

    if(all(priors==4)){
      priors=rep(8, ncol(Usubc))
      #priors=Usub[as.numeric(names(table(Zest)))[which.max(table(Zest))],]
    }

    mu=NULL
    for(ii in 1:6){
      mu_tmp=c(rep(0.5*ii,ii+1),c(0,( (1:(ii))/ii)))
      mu=rbind(mu, matrix(mu_tmp, byrow=F, ncol=2))
    }

    # correlation
    # while(sd(priors)==0){
    #   priors=sample(c(0.001, -0.001), length(priors), replace = T)+priors
    # }
    priorsc=mu[priors,1]

    #ind_cna=which(priors != 1)
    #corrs=apply(Usubc[,corr_region_ind,drop=F], 1, function(x)
    #  sum((as.numeric(x)-1)*(as.numeric(priorsc[corr_region_ind])-1)))#/(sqrt(sum((as.numeric(x)-1)^2))*sqrt(sum((as.numeric(priors[])-1)^2))))

    corrs=apply(Usub2[,corr_region_ind,drop=F], 1, function(x){
      tmp=sort((as.numeric(x)-1)*(as.numeric(priorsc[corr_region_ind])-1))
      if(length(tmp)==0){
        return(NA)
      }else{
        if(rm_extreme==T){
          return(sum(tmp[2:(length(tmp)-1)]))
        }else{
          return(sum(tmp))}}
    }
    )

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

    return(list(Zest=Zest, corrs=corrs[Zest], annot=annot[Zest],Uest=Usub,Ucest=Usubc,Uaest=Usuba, sigmas_cov_est=sigmas_cov,sigmas_allele_est=sigmas_allele,  Zmaj= maj_vote, cutoff=cutoff, U0=cluster_obj$priors$U0, wU0=cluster_obj$priors$wU0))


  }


}
