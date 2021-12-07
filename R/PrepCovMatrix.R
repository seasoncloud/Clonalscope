PrepCovMatrix=function(deltas_all=NULL, ngene_filter=200, ncell_filter=NULL, prep_mode='union'){
  ## library required
  library(Matrix)
  df=deltas_all$deltas_all
  ngenes=deltas_all$ngenes
  seg_table_filtered=deltas_all$seg_table_filtered

  ## parameter setting
  #ngene_filter
  if(is.null(ncell_filter)){
    ncells=as.numeric(lapply(df, length))
    ncell_filter=max(ncells)-(max(ncells)*0.1)
  }
  message(paste0("Selecting regions with > ", ngene_filter, " genes \n"))
  message(paste0("Selecting regions with > ", ncell_filter, " cells \n"))

  ## generate matrix
  cell_list<-lapply(df, function(x) {names(x)})

  df=df[unlist(lapply(cell_list, function(x) (length(x)>ncell_filter)))]
  ngenes=ngenes[which(unlist(lapply(cell_list, function(x) (length(x)>ncell_filter))))]
  seg_table_filtered=seg_table_filtered[which(unlist(lapply(cell_list, function(x) (length(x)>ncell_filter)))),, drop=F]

  cell_list<-lapply(df, function(x) {names(x)})
  if(prep_mode=='union'){
  cell_intersect <- Reduce(union, cell_list)
  }else if(prep_mode=='intersect'){
    cell_intersect <- Reduce(intersect, cell_list)
  }else{
    stop("Please specify the prep_mode as either union or intersect!")
  }
  df <- sapply(df,function(x){
    x[match(cell_intersect, names(x))]
  })

  df=data.frame(df, stringsAsFactors = F)
  rownames(df)=cell_intersect
  colnames(df)=sapply(strsplit(colnames(df),"alphaall"),'[',2)

  df=df[,which(ngenes>ngene_filter)]
  seg_table_filtered=seg_table_filtered[which(ngenes>ngene_filter),,drop=F]
  colnames(df)=paste0("chr",seg_table_filtered$chr,"-",seg_table_filtered$start, "-", seg_table_filtered$end )

  message("Matrix generated!")

  return(list(df=df, ngene_filter=ngene_filter, ncell_filter=ncell_filter, seg_table_filtered=seg_table_filtered))

}
