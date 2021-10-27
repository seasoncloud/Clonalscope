#' Plot clustering result from the BayesNonparCluster estimation.
#'
#' @param df A matrix/data.frame with each row being a cell and each column being a segment (from the WES/WGS). The values are the estimated fold change for each cell in each reion (from EstRegionCov).
#' @param celltype A matrix with 2+ columns: COL1- cell barcodes; COL2+- cell types and other information for plotting.
#' @param Zest Estimated cell cluster identity from the BayesNonparCluster estimation.The order should be the same as that of the input matrix.
#' @param mode Character. Mode used to color the cells in the UMAP. The value should be one of ["celltype","Zest", "annot"].
#' @param maxv Numeric. Set the ceiling number for plotting.
#'
#' @return A UMAP showing the result from the Bayesian non-parametric clustering.
#'
#' @import pheatmap
#' @import uwot
#' @export
PlotUMAP=function(df=NULL, celltype=NULL, Assign_obj=NULL, mode="Zest",  maxv=3, size=1){
  library(pheatmap)
  Zest=Assign_obj$Zest
  #corrs=Assign_obj$corrs
  annot=Assign_obj$annot
  U0=Assign_obj$U0

  #celltype=cbind(celltype[,1:2], as.character(Zest)[match(celltype[,1],rownames(df))])
  celltype=celltype[match(rownames(df), celltype[,1]), , drop=F]

  celltype=cbind(celltype[,],as.character(annot)[match(celltype[,1],rownames(df))], as.character(Zest)[match(celltype[,1],rownames(df))])
  if(is.null(colnames(celltype))){
    colnames(celltype)=paste0('COL', 1:ncol(celltype))
  }

  colnames(celltype)=c(colnames(celltype)[1:(ncol(celltype)-2)],"annot", "Zest")
  df2=apply(df, c(2), function(x) pmin(x, maxv))

  celltype0=celltype[,-c(1), drop=F]
  rownames(celltype0)=celltype[,1]
  #celltype0=celltype0[od,, drop=F]
  celltype0=data.frame(celltype0, stringsAsFactors = F)
  # set colors

  ann_colors=list()
  for(ii in 1:ncol(celltype0)){
    if(colnames(celltype0)[ii]=='Zest'){
      if(is.numeric(Zest)){
        cols0=as.numeric(celltype0[,ii])
      }else{
        cols0=(celltype0[,ii])}
      new_indc=which(! names(table(cols0)) %in% as.character(1:max(which(!is.na(U0[,1])))))
      if((length(new_indc)+length(1:max(which(!is.na(U0[,1])))))>20){
        col_use=c(colors()[c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)],
                  colors()[sample(1:length(colors())[-c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)],(length(new_indc)+length(1:max(which(!is.na(U0[,1])))))-20)])
      }else{
        col_use=c(colors()[c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:(length(new_indc)+length(1:max(which(!is.na(U0[,1])))))]
      }
      names(col_use)=as.character(c(which(!is.na(U0[,1])), names(table(cols0))[new_indc]))
      col_use=col_use[which(names(col_use) %in% names(table(cols0)))]

    }else{
      cols0=(celltype0[,ii])

      if(length(table(cols0))>20){
        #if((length(new_indc)+length(1:max(which(!is.na(U0[,1])))))>20){
        col_use=c(colors()[c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)],
                  colors()[sample(1:length(colors())[-c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)],length(table(cols0))-20)])
      }else{
        col_use=c(colors()[c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:length(table(cols0))]
      }
      names(col_use)=names(table(cols0))
    }
    #table(cols0)

    #tmp=1:length(table(cols0))
    #names(tmp)=paste0('c',names(table(cols0)))
    #celltype_df$celltype_df=as.character(tmp[paste0('c', cols0)])

    ## names(col_use)=names(table(cols0))
    ann_colors[[colnames(celltype0)[ii]]]=col_use
  }

  umap_emb=uwot::umap(df2)#,

  emb=data.frame(#TSNE1=tsne_emb[,1], TSNE2=tsne_emb[,2],
    UMAP1=umap_emb[,1], UMAP2=umap_emb[,2],
    celltype=celltype0[,1], Zest=as.character(celltype0[,3]),annot=celltype0[,2])#,
  #condition=celltype00[match(rownames(df2[,]), celltype00[,1]),3])#,



  library(ggplot2)

  if(mode=='Zest'){
  pp=ggplot(emb, aes(x=UMAP1, y=UMAP2))+geom_point(aes(color=Zest), size=size)+
    scale_color_manual(breaks = names(ann_colors[[colnames(celltype0)[3]]]),
                       values=ann_colors[[colnames(celltype0)[3]]])
  }else if(mode=='celltype'){
    pp=ggplot(emb, aes(x=UMAP1, y=UMAP2))+geom_point(aes(color=celltype), size=size)+
      scale_color_manual(breaks = names(ann_colors[[colnames(celltype0)[1]]]),
                         values=ann_colors[[colnames(celltype0)[1]]])
  }else if(mode=='annot'){
    pp=ggplot(emb, aes(x=UMAP1, y=UMAP2))+geom_point(aes(color=annot), size=size)+
      scale_color_manual(breaks = names(ann_colors[[colnames(celltype0)[2]]]),
                         values=ann_colors[[colnames(celltype0)[2]]])
  }else{
    stop("Please specify a valid mode")
  }

  plot(pp) ###


}
