#' A wrapper function to run "EstRegionCov", "PrepCovMatrix", "BayesNonparCluster", "MCMCtrim", and "AssignCluster".
#'
#' Please see the documentation of each function for more details.
#'
#' @param mtx A (sparse) matrix with each row being a feature and each column being a cell barcode. (Output from cellranger)
#' @param barcodes A matrix/data.frame with barcodes for each cell in the first column.
#' @param features A matrix/data.frame with the columns- 1st: gendID; 2nd: gene name. Each row is a feature whose order should correspond to the rows of "mtx".
#' @param bed A matrix with the bed file format. COL1- chromosome; COL2- start coordinate; COL3- end coordinate; COL4- the gene ID matches the one in "features" (COL1).
#' @param celltype0 A matrix with two columns: COL1- cell barcodes; COL2- cell types (Tumor cells should be labeled as "tumor" and the set of cells used as control should be labeled as "normal").
#' @param var_pt A decimal between [0,1]. Exclude the extremely variable genes with variance > the var_pt percentile for the tested cells.
#' @param var_pt_ctrl A decimal between [0,1]. Exclude the extremely variable genes with variance > the var_pt percentile for the control cells.
#' @param include The cells included for testing. The value should be one of ["tumor","all","control"].
#' @param alpha_source The regions used to compute normalizing factors. The value should be one of ["all", "control"]. If "all", total read counts are used; if "control", a subset of chromosomes are used (ctrl_region needs to be specified).
#' @param ctrl_region A vector indicating the chromosomes used to compute normalizing factors. The format should be c("chr1", "chr3", ...)
#' @param seg_table_filtered A segmentation table from the "Segmentation_bulk" function.
#' @param size A matrix with two columns: col1: different chromosome; col2: for the size (bp) of different chromosomes.
#' @param dir_path A path for saving the objects and plots.
#' @param breaks A parameter used in the "hist" function.
#' @param ngene_filter Integer. Include only the segments/regions covering at lest ngene_filter genes.
#' @param prep_mode Character. Please specify either "intersect" or "union". "Intersect" selects cells with non-zero read counts for all the selected regions. "Union" includes all cells in the datasets.
#' @param est_cap Integer. Set the maximum values for the coverage change across all cells and regions.
#' @param alpha The concentration parameter used in the nested Chinese restaurant process controlling the probability that a new CNV state is sampled other than the existing ones for each region.
#' @param beta The concentration parameter used in the nested Chinese restaurant process controlling the probability that a new subclone is sampled for each cell.
#' @param niter Integer. Number of iteration for the Gibb's sampling.
#' @param sigmas0 A numeric vector with the same length as the number of columns of Xir (or length 1). The prior standard deviation for the normal distributoion.
#' @param U0 A matrix/data.frame with each row being a cluster and each column being a region (ncol(U0) should be the same as ncol(Xir)). The values represent prior mean of relatvie coverage for each cluster.
#' @param Z0 A numeric vector indicating the cluster identity for each cell. The order show be the same as that of Xir's row. The values should represent the order of the cluster the same as that of U0's row.
#' @param clust_mode Character. Cluster the cells based on either "all" regions or "cna_only" (regions where CNAs are observed in the matched WGS/WES data).
#' @param clustering_barcodes A vector specifying the cells (barcodes) that are used in the clustering.
#' @param seed Integer. Set the seed for the stochastic process.
#' @param clustering0 A list object which is the output of BayesNonparCluster from another associated dataset (for subclone tracing).
#' @param result0 A list object which is the output of AssignCluster from another associated dataset (for subclone tracing).
#' @param verbose Logical. Whether or not to print out the results for each iteration for the BayesNonparCluster function.
#' @param burnin Integer. Number of iterations to throw away for the MCMC process.
#' @param thinning Integer. Subsampling every integer number of iterations.
#' @param mincell Integer. A cluster will be considered if it contains at least mincell cells. Otherwise, cells in the invalid clusters will be re-assigned to the existing clusters.
#' @param cutoff Integer. Set the threshold for classify clusters as malignant or not by comparing to the WGS copy number profiles. This threshould may depend on the datasets, but the relative values of the malignancy index are still informative.
#' @param re_est. Logical. Whether or not to initiate the next round estimation.
#' @param threshold_2nd. Numeric usually negative values. If any of the clusters with malignancy index < threshold_2nd, initiate the next round estimation. Very negative values comparing to the WGS copy number profiles suggest usage of the control/reference cells in coverage change estimation is unreliable. The parameter is effective only when "re_est" is 'TRUE'.
#' @param Est_read Logical. Whether or not to read the output objects from the EstRegionCov function from the specified path. (So no need to rerun the function)
#' @param Clust_read Logical. Whether or not to read the output objects from the BayesNonparCluster function from the specified path. (So no need to rerun the function)
#'
#' @return A list with the output results from "EstRegionCov" (deltas_all), "PrepCovMatrix" (df_obj), "BayesNonparCluster" (clustering), "MCMCtrim" (clustering2), and "AssignCluster" (result).
#' "result"$: A list with the output files from the $ iteration.
#' "result_final": A list with the output files from the last iteration.
#' "celltype0": cell type used in the analysis.
#'
#' @export
RunCovCluster=function(mtx=NULL, barcodes=NULL, features=NULL, bed=NULL, celltype0=NULL, var_pt=0.99, var_pt_ctrl=0.99, include='tumor',
                       alpha_source='all', ctrl_region=NULL, seg_table_filtered=NULL,size=NULL, dir_path=NULL, breaks=50,
                       ngene_filter=150, prep_mode='intersect',
                       est_cap=3, alpha=2, beta=2, niter=200,
                       sigmas0=NULL, U0=NULL, Z0=NULL,clust_mode='all',clustering_barcodes=NULL, seed=200, clustering0=NULL,result0=NULL,verbose=FALSE,
                       burnin=NULL, thinning=1 ,mincell = NULL, cutoff = 0.2,
                       threshold_2nd=-0.2, re_est=NULL,Est_read1=FALSE,Est_read2=FALSE, Clust_read1=FALSE, Clust_read2=FALSE,dense_mat=T){

  plot_path=paste0(dir_path,"/cov_hist.pdf")
  plot_path2=paste0(dir_path,"/cov_hist_updated.pdf")
  result_all=list()
  
  # check for matrix sparsity, convert sparse matrix to dense matrix
  if(dense_mat==T){
    mtx = as.matrix(mtx)
  }
  
  if(is.null(re_est)){
    if(is.null(celltype0)){
      re_est=TRUE
    }else if(!is.null(celltype0)){
      re_est=FALSE
    }
  }

  if(is.null(clustering_barcodes)){
    clustering_barcodes=barcodes[,1]
  }

  if(is.null(celltype0)){
    celltype0=cbind(barcodes[,1], rep("normal", nrow(barcodes)))
  }

  if( Est_read1==TRUE){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
    deltas_all=readRDS(paste0(dir_path, "/deltas_allseg.rds"))
  }else{
    deltas_all=EstRegionCov(mtx=mtx, barcodes=barcodes, features=features, bed=bed, celltype0=celltype0, var_pt=var_pt, var_pt_ctrl=var_pt_ctrl, include=include,
                            alpha_source=alpha_source, ctrl_region=ctrl_region, seg_table_filtered=seg_table_filtered, size=size,
                            plot_path=plot_path, breaks=breaks)

    saveRDS(deltas_all, paste0(dir_path, "/deltas_allseg.rds"))
  }


  df_obj=PrepCovMatrix(deltas_all=deltas_all, ngene_filter=ngene_filter, prep_mode=prep_mode)
  df=df_obj$df
  #rownames(df)=gsub("[.]","-", rownames(df))
  print(paste0("Dimension after filtering:",dim(df)))


  # if(is.null(mincell)){
  #   if(nrow(df)<500){
  #     mincell=5
  #   }else if(nrow(df)<1500){
  #     mincell=10
  #   }else{
  #     mincell=20
  #   }
  # }

  seg_table_filtered=df_obj$seg_table_filtered

  cna_states_WGS=seg_table_filtered$states
  cna_states_WGS[which(cna_states_WGS>1)]= 1.5
  cna_states_WGS[which(cna_states_WGS<1)]= 0.5
  cna_states_WGS=as.numeric(cna_states_WGS)

  df2=apply(df, c(2), function(x) pmin(x, est_cap))

  if(!is.null(clustering0)){
    if(is.null(result0)){
      clustering02=MCMCtrim(clustering0)
      result0=AssignCluster(clustering02)
    }else{
      if(ncol(clustering0$data)==ncol(result0$Uest)){
        result0=result0
      }else{
        stop("clustering0 and result0 do not match!")
      }
    }
    ## subset the regions
    df2=df2[,which(colnames(df2) %in% colnames(clustering0$data)), drop=F]
    U0=result0$Uest
    U0=U0[,which(colnames(clustering0$data) %in% colnames(df2)), drop=F]
    cna_states_WGS=clustering0$priors$cna_states_WGS[which(colnames(clustering0$data) %in% colnames(df2))]
  }

  if( Clust_read1==TRUE){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
    clustering=readRDS(paste0(dir_path, "/nonpara_clustering.rds"))

  }else{
    if(clust_mode=='all'){
      selr=1:length(cna_states_WGS)
     # clustering=BayesNonparCluster(Xir=df2[which(rownames(df2) %in% clustering_barcodes), , drop=F], cna_states_WGS =cna_states_WGS , alpha=alpha, beta=beta, niter = niter , sigmas0 =  sigmas0, U0 = U0 , Z0 = Z0[which(rownames(df2) %in% clustering_barcodes)] , seed = seed)
    }else if(clust_mode=='cna_only'){
      selr=which(as.character(cna_states_WGS)!='1')
    }else{
      stop("Please specify a valid clust_mode.")
    }
    clustering=BayesNonparCluster(Xir=df2[which(rownames(df2) %in% clustering_barcodes),selr], cna_states_WGS =cna_states_WGS[selr] , alpha=alpha, beta=beta, niter = niter , sigmas0 =  sigmas0, U0 = U0 , Z0 = Z0[which(rownames(df2) %in% clustering_barcodes)] , seed = seed, verbose=verbose)

    saveRDS(clustering,paste0(dir_path, "/nonpara_clustering.rds"))
  }
  clustering2=MCMCtrim(clustering, burnin = burnin, thinning = thinning)

  ## sel mincelll
  if(is.null(mincell)){
   mincell=round(length(clustering_barcodes)*0.01)
  }

  result=AssignCluster(clustering2, mincell = mincell, cutoff = cutoff)
  Zest=result$Zest
  print(paste0("Zest: ",table(Zest)))

  df_obj$df=df_obj$df[which(rownames(df_obj$df) %in% clustering_barcodes),, drop=F]

  result1=list( deltas_all= deltas_all, cna_states_WGS=cna_states_WGS, df_obj=df_obj, clustering=clustering, clustering2=clustering2, result=result)
  result_all[['result1']]=result1

  # update coverage estimation if needed
  result2=NULL
  tmp=1
  while((any(result$corrs< threshold_2nd) & tmp<10) & re_est==TRUE){
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
                              alpha_source=alpha_source, ctrl_region=ctrl_region, seg_table_filtered=seg_table_filtered, size=size,
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
      # clustering02=MCMCtrim(clustering0)
      # result0=AssignCluster(clustering02)

      ## subset the regions
      df2=df2[,which(colnames(df2) %in% colnames(clustering0$data)), drop=F]
      U0=result0$Uest
      U0=U0[,which(colnames(clustering0$data) %in% colnames(df2)), drop=F]
      cna_states_WGS=clustering0$priors$cna_states_WGS[which(colnames(clustering0$data) %in% colnames(df2))]
    }

    if( Clust_read2==TRUE & paste0("nonpara_clustering_updated",tmp,".rds") %in% list.files(dir_path)){# | !("deltas_allseg.rds" %in% list.files(paste0(dir_path,"rds")))){
      clustering=readRDS(paste0(dir_path, "/nonpara_clustering_updated",tmp,".rds"))

    }else{
      if(clust_mode=='all'){
        selr=1:length(cna_states_WGS)
        #clustering=BayesNonparCluster(Xir=df2, cna_states_WGS =cna_states_WGS , alpha=alpha, beta=beta, niter = niter , sigmas0 =  sigmas0, U0 = U0 , Z0 = Z0 , seed = seed)
      }else if(clust_mode=='cna_only'){
        selr=which(as.character(cna_states_WGS)!='1')
      }else{
        stop("Please specify a valid clust_mode.")
      }
      clustering=BayesNonparCluster(Xir=df2[which(rownames(df2) %in% clustering_barcodes),selr], cna_states_WGS =cna_states_WGS[selr] , alpha=alpha, beta=beta, niter = niter , sigmas0 =  sigmas0, U0 = U0 , Z0 = Z0[which(rownames(df2) %in% clustering_barcodes)] , seed = seed, verbose=verbose)

      saveRDS(clustering,paste0(dir_path, "/nonpara_clustering_updated",tmp,".rds"))
    }
    clustering2=MCMCtrim(clustering, burnin = burnin, thinning = thinning)

    result=AssignCluster(clustering2, mincell = mincell, cutoff = cutoff)
    table(result$Zest)
    Zest=result$Zest

    print(paste0("Zest: ",table(Zest)))

    df_obj$df=df_obj$df[which(rownames(df_obj$df) %in% clustering_barcodes),, drop=F]

    result2=list( deltas_all= deltas_all, cna_states_WGS=cna_states_WGS, df_obj=df_obj, clustering=clustering, clustering2=clustering2, result=result, new_normal=new_normal)
    result_all[[paste0('result',tmp)]]=result2

    if(tmp==10){
      warning("After 10 rounds of iterations, there is still a cluster with corrs<-0.3.")
    }
    #Est_read2=FALSE
    #Clust_read2=FALSE
  }

  result_all[['celltype0']]=celltype0
  ll=sort(as.numeric(gsub("result","",names(result_all)[grepl("result", names(result_all))])))[length(names(result_all)[grepl("result", names(result_all))])]
  result_all[['result_final']]=result_all[[paste0('result',ll )]]
  return(result_all)
}
