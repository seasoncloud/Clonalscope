RunCovCluster=function(mtx=NULL, barcodes=NULL, features=NULL, bed=NULL, celltype0=NULL, var_pt=0.99, var_pt_ctrl=0.99, include='tumor',
                       alpha_source='all', ctrl_region=NULL, seg_table_filtered=NULL,size=NULL, dir_path=NULL, breaks=50,
                       ngene_filter=150, est_cap=3, alpha=2, beta=2, niter=200,
                       sigmas0=NULL, U0=NULL, Z0=NULL, seed=200, clustering0=NULL,
                       threshold_2nd=-0.4, burnin=NULL, thinning=1 ,mincell = 20, cutoff = 0.3,
                       Est_read1=FALSE,Est_read2=FALSE, Clust_read1=FALSE, Clust_read2=FALSE){

  plot_path=paste0(dir_path,"/cov_hist.pdf")
  plot_path2=paste0(dir_path,"/cov_hist_updated.pdf")
  result_all=list()

  if(is.null(celltype0)){
    celltype0=cbind(barcodes[,1], rep("normal", nrow(barcodes)))
  }

  if( Est_read1==TRUE){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
    deltas_all=readRDS(paste0(dir_path, "/deltas_allseg.rds"))
  }else{
    deltas_all=EstRegionCov(mtx=mtx, barcodes=barcodes, features=features, bed=bed, celltype0=celltype0, var_pt=var_pt, var_pt_ctrl=var_pt_ctrl, include=include,
                            alpha_source=alpha_source, ctrl_region=ctrl_region, seg_table_filtered=Obj_filtered$seg_table_filtered, size=size,
                            plot_path=plot_path, breaks=breaks)

    saveRDS(deltas_all, paste0(dir_path, "/deltas_allseg.rds"))
  }


  df_obj=PrepCovMatrix(deltas_all=deltas_all, ngene_filter=ngene_filter)
  df=df_obj$df
  print(paste0("Dimension after filtering:",dim(df)))

  seg_table_filtered=df_obj$seg_table_filtered

  cna_states_WGS=seg_table_filtered$states
  cna_states_WGS[which(cna_states_WGS>1)]= 1.5
  cna_states_WGS[which(cna_states_WGS<1)]= 0.5
  cna_states_WGS=as.numeric(cna_states_WGS)

  df2=apply(df, c(2), function(x) pmin(x, est_cap))

  if(!is.null(clustering0)){
    clustering02=MCMCtrim(clustering0)
    result0=AssignCluster(clustering02)

    ## subset the regions
    df2=df2[,which(colnames(df2) %in% colnames(clustering0$data)), drop=F]
    U0=result0$Uest
    U0=U0[,which(colnames(clustering0$data) %in% colnames(df2)), drop=F]
    cna_states_WGS=clustering0$priors$cna_states_WGS[which(colnames(clustering0$data) %in% colnames(df2))]
  }

  if( Clust_read1==TRUE){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
    clustering=readRDS(paste0(dir_path, "/nonpara_clustering.rds"))

  }else{

    clustering=BayesNonparCluster(Xir=df2, cna_states_WGS =cna_states_WGS , alpha=alpha, beta=beta, niter = niter , sigmas0 =  sigmas0, U0 = U0 , Z0 = Z0 , seed = seed)
    saveRDS(clustering,paste0(dir_path, "/nonpara_clustering.rds"))
  }
  clustering2=MCMCtrim(clustering)

  result=AssignCluster(clustering2, mincell = mincell, cutoff = cutoff)
  Zest=result$Zest
  print(paste0("Zest: ",table(Zest)))


  result1=list( deltas_all= deltas_all, cna_states_WGS=cna_states_WGS, df_obj=df_obj, clustering=clustering, clustering2=clustering2, result=result)
  result_all[['result1']]=result1

  # update coverage estimation if needed
  result2=NULL
  tmp=1
  while(any(result$corrs< threshold_2nd)){
    tmp=tmp+1
    message("Estimated clusters with negative correlation to the WGS data.")
    message(paste0("Start estimation ",tmp))

    new_normal=rownames(clustering2$data)[which(result$corrs %in% unique(result$corrs)[which(unique(result$corrs)< threshold_2nd)])]
    celltype0[which(celltype0[,2]=='normal'),2]='normal_pre'
    celltype0[which(celltype0[,1] %in% new_normal),2]='normal'

    #   print(paste0('tmp:', tmp))
    # print(paste0('Est_read2:', Est_read2))
    # print(paste0('Clust_read2:', Clust_read2))

    if( Est_read2==TRUE & paste0("deltas_allseg_updated",tmp,".rds") %in% list.files(dir_path)){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
      deltas_all=readRDS(paste0(dir_path, "/deltas_allseg_updated",tmp,".rds"))
    }else{
      deltas_all=EstRegionCov(mtx=mtx, barcodes=barcodes, features=features, bed=bed, celltype0=celltype0, var_pt=var_pt, var_pt_ctrl=var_pt_ctrl, include=include,
                              alpha_source=alpha_source, ctrl_region=ctrl_region, seg_table_filtered=Obj_filtered$seg_table_filtered, size=size,
                              plot_path=plot_path2, breaks=breaks)

      saveRDS(deltas_all, paste0(dir_path, "/deltas_allseg_updated",tmp,".rds"))
    }
    df_obj=PrepCovMatrix(deltas_all=deltas_all, ngene_filter=ngene_filter)
    df=df_obj$df
    print(paste0("Dimension after filtering:",dim(df)))
    seg_table_filtered=df_obj$seg_table_filtered

    cna_states_WGS=seg_table_filtered$states
    cna_states_WGS[which(cna_states_WGS>1)]= 1.5
    cna_states_WGS[which(cna_states_WGS<1)]= 0.5
    cna_states_WGS=as.numeric(cna_states_WGS)

    df2=apply(df, c(2), function(x) pmin(x, est_cap))

    if(!is.null(clustering0)){
      clustering02=MCMCtrim(clustering0)
      result0=AssignCluster(clustering02)

      ## subset the regions
      df2=df2[,which(colnames(df2) %in% colnames(clustering0$data)), drop=F]
      U0=result0$Uest
      U0=U0[,which(colnames(clustering0$data) %in% colnames(df2)), drop=F]
      cna_states_WGS=clustering0$priors$cna_states_WGS[which(colnames(clustering0$data) %in% colnames(df2))]
    }

    if( Clust_read2==TRUE & paste0("nonpara_clustering_updated",tmp,".rds") %in% list.files(dir_path)){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
      clustering=readRDS(paste0(dir_path, "/nonpara_clustering_updated",tmp,".rds"))

    }else{


      clustering=BayesNonparCluster(Xir=df2, cna_states_WGS =cna_states_WGS , alpha=alpha, beta=beta, niter = niter , sigmas0 =  sigmas0, U0 = U0 , Z0 = Z0 , seed = seed)
      saveRDS(clustering,paste0(dir_path, "/nonpara_clustering_updated",tmp,".rds"))
    }
    clustering2=MCMCtrim(clustering, burnin = burnin, thinning = thinning)

    result=AssignCluster(clustering2, mincell = mincell, cutoff = cutoff)
    table(result$Zest)
    Zest=result$Zest

    print(paste0("Zest: ",table(Zest)))

    result2=list( deltas_all= deltas_all, cna_states_WGS=cna_states_WGS, df_obj=df_obj, clustering=clustering, clustering2=clustering2, result=result, new_normal=new_normal)
    result_all[[paste0('result',tmp)]]=result2

    #Est_read2=FALSE
    #Clust_read2=FALSE
  }

  result_all[['celltype0']]=celltype0
  return(result_all)
}
