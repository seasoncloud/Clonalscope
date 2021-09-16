#' Plot clustering result from the BayesNonparCluster estimation.
#'
#' @param df A matrix/data.frame with each row being a cell and each column being a segment (from the WES/WGS). The values are the estimated fold change for each cell in each reion (from EstRegionCov).
#' @param celltype A matrix with 2+ columns: COL1- cell barcodes; COL2+- cell types and other information for plotting.
#' @param Zest Estimated cell cluster identity from the BayesNonparCluster estimation.The order should be the same as that of the input matrix.
#' @param mode Character indicating cell annotation in the heatmap. The value should be one of ["cluster","celltype"]."cluster": cell annotation with only estimated cluster identity; "celltype": cell annotation with both clsuter identity and cell types.
#' @param maxv Numeric. Set the ceiling number for plotting.
#'
#' @return A heatmap showing the result from the Bayesian non-parametric clustering. Each row is a cell and each column is a region. The values are the coverage change.
#'
#' @import pheatmap
#' @export
PlotClusters=function(df=NULL, celltype=NULL, Assign_obj=NULL, mode="cluster", maxv=3){

  Zest=Assign_obj$Zest
  corrs=Assign_obj$corrs[Zest]
  annot=Assign_obj$annot[Zest]

  #celltype=cbind(celltype[,1:2], as.character(Zest)[match(celltype[,1],rownames(df))])
  celltype=cbind(celltype[,],as.character(annot)[match(celltype[,1],rownames(df))], as.character(Zest)[match(celltype[,1],rownames(df))], as.numeric(corrs)[match(celltype[,1],rownames(df))])
  colnames(celltype)=c(colnames(celltype)[1:(ncol(celltype)-3)],"annot", "Zest","corr")
  df2=apply(df, c(2), function(x) pmin(x, maxv))
  celltype_cluster=celltype[,(ncol(celltype)), drop=F]
  rownames(celltype_cluster)=celltype[,1]
  od=order(as.numeric(celltype_cluster[,1]))
  df2=df2[od,, drop=F]
  celltype_cluster=celltype_cluster[od,, drop=F]
  celltype_cluster=data.frame(celltype_cluster, stringsAsFactors = F)

  if(mode=='celltype'){
    celltype0=celltype[,-c(1,ncol(celltype)), drop=F]
    rownames(celltype0)=celltype[,1]
    celltype0=celltype0[od,, drop=F]
    celltype0=data.frame(celltype0, stringsAsFactors = F)
    pheatmap(df2,cluster_cols = F, cluster_rows = F,clustering_distance_rows = "correlation",clustering_method = "ward.D2", show_rownames = F, show_colnames = T, annotation_row = celltype0)
  }else{
    celltype0=celltype[,(ncol(celltype)-1), drop=F]
    rownames(celltype0)=celltype[,1]
    celltype0=data.frame(celltype0, stringsAsFactors = F)
    pheatmap(df2,cluster_cols = F, cluster_rows = F,clustering_distance_rows = "correlation",clustering_method = "ward.D2", show_rownames = F, show_colnames = T, annotation_row = celltype0)
  }
}
