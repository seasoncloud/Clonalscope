#' Filter the output from EstRegionCov and prepare the input matrix for BayesNonparCluster.
#'
#' @param deltas_all A list object. Output of the "EstRegionCov" function storeing values of estimated coverage change for each cell and each region
#' @param ngene_filter Integer. Include only the segments/regions covering at lest ngene_filter genes.
#' @param ncell_filter Integer. Include only the segments/regions with at least ncell_filter cells carrying non-zero read counts.
#' @param prep_mode Character. Please specify either "intersect" or "union". "Intersect" selects cells with non-zero read counts for all the selected regions. "Union" includes all cells in the datasets.
#'
#' @return A list with the filtered data frame (df), the segments used (seg_table_filtered), and the arguments used.
#'
#' @export
PrepCovMatrix=function(deltas_all=NULL, ngene_filter=200, ncell_filter=NULL, prep_mode='intersect'){
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
