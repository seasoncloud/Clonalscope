#' Plot clustering result from the BayesNonparCluster estimation.
#'
#' @param df A matrix/data.frame with each row being a cell and each column being a segment (from the WES/WGS). The values are the estimated fold change for each cell in each reion (from EstRegionCov).
#' @param celltype A matrix with 2+ columns: COL1- cell barcodes; COL2 cell types and COL2+ other information for plotting.
#' @param Zest Estimated cell cluster identity from the BayesNonparCluster estimation.The order should be the same as that of the input matrix.
#' @param mode Character for the mode of heatmap to plot. The value should be one of ["segment","genome"]."segment": plot each segment with equal length; "genome": genome view for each segment.
#' @param maxv Numeric. Set the ceiling number for plotting.
#' @param lab_mode Numeric for the annotation mode. 1: T/N; 2: prob; FALSE: ignore.
#'
#' @return A heatmap showing the result from the Bayesian non-parametric clustering. Each row is a cell and each column is a region. The values are the coverage change.
#'
#' @import pheatmap
#' @import RColorBrewer
#' @export
PlotClusters=function(df=NULL, celltype=NULL, Assign_obj=NULL, mode="segment", consensus=F, maxv=2, fontsize=10,  fontsize_row = 10 , fontsize_col = 10, allele=F, lab_mode="annot", od_mode=2 ){
  maxv=pmax(maxv,2)
  if(allele==F){
    library(pheatmap)
    Zest=Assign_obj$Zest
    corrs=Assign_obj$corrs
    annot=Assign_obj$annot
    U0=Assign_obj$U0
    wU0=Assign_obj$wU0
    
    if(wU0==FALSE){
      U0=U0[1,, drop=F]
    }
    
    corrs=round(corrs,6)
    #corrs=(corrs-min(corrs))/(max(corrs)-min(corrs))
    
    #celltype=cbind(celltype[,1:2], as.character(Zest)[match(celltype[,1],rownames(df))])
    celltype=celltype[match(rownames(df), celltype[,1]), , drop=F]
    
    celltype=cbind(celltype[,],as.character(annot)[match(celltype[,1],rownames(df))], as.character(Zest)[match(celltype[,1],rownames(df))], as.numeric(corrs)[match(celltype[,1],rownames(df))])
    if(is.null(colnames(celltype))){
      colnames(celltype)=paste0('COL', 1:ncol(celltype))
    }
    
    colnames(celltype)=c(colnames(celltype)[1:(ncol(celltype)-3)],"annot", "Zest","corr")
    df2=apply(df, c(2), function(x) pmin(x, maxv))
    celltype_cluster=celltype[,c(2,ncol(celltype)), drop=F]
    rownames(celltype_cluster)=celltype[,1]
    if(od_mode==1){
      od=order(as.numeric(celltype_cluster[,2]), decreasing = F)
    }else if(od_mode==2){
      od=order(as.numeric(celltype_cluster[,2]),celltype_cluster[,1], decreasing = F)
    }else{
      stop("od_mode should be 1 or 2.")
    }
    df2=df2[od,, drop=F]
    celltype_cluster=celltype_cluster[od,, drop=F]
    celltype_cluster=data.frame(celltype_cluster, stringsAsFactors = F)
    
    # ## order with hierarchical clustering
    # corr_name=celltype_cluster[which(!duplicated(celltype_cluster[,1])),1]
    # df3=NULL
    # for(ii in 1:length(corr_name)){
    #   df_tmp=df2[which(celltype_cluster[,1] == corr_name[ii]),,drop=F]
    #   #DD=dist(df_tmp)
    #   DD=Dist(df_tmp, method='correlation')
    #   hh=hclust(DD, method=clust_method)
    #   df_tmp=df_tmp[hh$order,,drop=F]
    #   df3=rbind(df3, df_tmp)
    # }
    # colnames(df3)=colnames(df2)
    # df2=df3
    
    if(lab_mode!=FALSE){
      if(lab_mode=='annot'){
        celltype0=celltype[,-c(1,ncol(celltype)), drop=F]
      }else if(lab_mode=='corr'){
        celltype0=celltype[,-c(1,(ncol(celltype)-2)), drop=F]
        celltype0=celltype0[,c(1:(ncol(celltype0)-2), ncol(celltype0), (ncol(celltype0)-1)) , drop=F]
      }else{
        stop("Please specify a valid lab_mode mode!")
      }
    }else{
      celltype0=celltype[,-c(1,(ncol(celltype)-2),ncol(celltype)), drop=F]
    }
    rownames(celltype0)=celltype[,1]
    celltype0=celltype0[od,, drop=F]
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
        if((length(new_indc)+length(which(!is.na(U0[,1]))))>20){  ##
          col_use=c(colors()[c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)],
                    colors()[sample(1:length(colors())[-c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)],(length(new_indc)+length(which(!is.na(U0[,1]))))-20)])
        }else{
          col_use=c(colors()[c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:(length(new_indc)+length(which(!is.na(U0[,1]))))]
        }
        names(col_use)=as.character(c(which(!is.na(U0[,1])), names(table(cols0))[new_indc]))
        col_use=col_use[which(names(col_use) %in% names(table(cols0)))]
      }else{
        #}else if(colnames(celltype0)[ii]!='corr'){
        cols0=(celltype0[,ii])
        
        if(length(table(cols0))>20){
          #if((length(new_indc)+length(1:max(which(!is.na(U0[,1])))))>20){
          col_use=c(colors()[c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)],
                    colors()[sample(1:length(colors())[-c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)],length(table(cols0))-20)])
        }else{
          col_use=c(colors()[c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:length(table(cols0))]
        }
        names(col_use)=names(table(cols0))
        # }else{
        #   col_use=c("black", "firebrick")
      }
      #table(cols0)
      
      #tmp=1:length(table(cols0))
      #names(tmp)=paste0('c',names(table(cols0)))
      #celltype_df$celltype_df=as.character(tmp[paste0('c', cols0)])
      
      ## names(col_use)=names(table(cols0))
      ann_colors[[colnames(celltype0)[ii]]]=col_use
    }
    
    if(consensus==T){
      df2_colnamae=colnames(df2)
      df2_rownames=rownames(df2)
      df2=Assign_obj$Uest[as.numeric(celltype[od,(ncol(celltype)-1), drop=T]),,drop=F]
      df2=apply(df2, c(2), function(x) pmin(x, maxv))
      colnames(df2)=df2_colnamae
      rownames(df2)=df2_rownames
      
    }
    
    breaksList = seq(0-(maxv-2), maxv, by = (maxv+(maxv-2))/100)
    
    if(mode=='segment'){
      pheatmap(df2,cluster_cols = F, cluster_rows = F,clustering_distance_rows = "correlation",clustering_method = "ward.D2",
               show_rownames = F, show_colnames = T, annotation_row = celltype0,annotation_colors = ann_colors,
               fontsize = fontsize,  fontsize_row =fontsize_row , fontsize_col = fontsize_col,
               color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100), breaks = breaksList)
      
    }else{
      
      chr=sapply(strsplit(colnames(df2),"-"),'[',1)
      start=as.numeric(sapply(strsplit(colnames(df2),"-"),'[',2))
      end=as.numeric(sapply(strsplit(colnames(df2),"-"),'[',3))
      seglen=end-start
      repn=round(seglen/min(seglen))
      chr_rep=rep(chr, repn)
      chrn=table(chr_rep)[paste0('chr',1:22)]
      chrn=chrn[!is.na(chrn)]
      df2_plot=NULL
      for(ii in 1:length(repn)){
        mm=matrix(rep(df2[,ii], repn[ii]), ncol=repn[ii], byrow=F)
        df2_plot=cbind(df2_plot, mm)
      }
      rownames(df2_plot)=rownames(df2)
      coln=rep("", sum(repn))
      coln[ceiling(c(0,cumsum(chrn)[1:(length(chrn)-1)])+chrn/2)]=names(chrn)
      colnames(df2_plot)=coln
      
      pheatmap(df2_plot,cluster_cols = F, cluster_rows = F,clustering_distance_rows = "correlation",clustering_method = "ward.D2",
               show_rownames = F, show_colnames = T, annotation_row = celltype0, gaps_col = cumsum(chrn),
               annotation_colors = ann_colors, fontsize = fontsize, fontsize_row =fontsize_row , fontsize_col = fontsize_col,
               color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100),
               breaks = breaksList)
    }
    #   }else{
    #   celltype0=celltype[,(ncol(celltype)-1), drop=F]
    #   rownames(celltype0)=celltype[,1]
    #   celltype0=data.frame(celltype0, stringsAsFactors = F)
    #   pheatmap(df2,cluster_cols = F, cluster_rows = F,clustering_distance_rows = "correlation",clustering_method = "ward.D2", show_rownames = F, show_colnames = T, annotation_row = celltype0)
    # }
  }else{
    Zest=Assign_obj$Zest
    corrs=Assign_obj$corrs
    annot=Assign_obj$annot
    U0=Assign_obj$U0
    df_cov=df$Xir_cov
    df_allele=df$Xir_allele
    wU0=Assign_obj$wU0
    
    if(wU0==FALSE){
      U0=U0[1,, drop=F]
    }
    corrs=round(corrs,3)
    #corrs=(corrs-min(corrs))/(max(corrs)-min(corrs))
    
    
    #celltype=cbind(celltype[,1:2], as.character(Zest)[match(celltype[,1],rownames(df))])
    celltype=celltype[match(rownames(df_cov), celltype[,1]), , drop=F]
    
    celltype=cbind(celltype[,],as.character(annot)[match(celltype[,1],rownames(df_cov))], as.character(Zest)[match(celltype[,1],rownames(df_cov))], as.numeric(corrs)[match(celltype[,1],rownames(df_cov))])
    if(is.null(colnames(celltype))){
      colnames(celltype)=paste0('COL', 1:ncol(celltype))
    }
    
    colnames(celltype)=c(colnames(celltype)[1:(ncol(celltype)-3)],"annot", "Zest","corr")
    df2_cov=apply(df_cov, c(2), function(x) pmin(x, maxv))
    celltype_cluster=celltype[,c(2,ncol(celltype)), drop=F]
    rownames(celltype_cluster)=celltype[,1]
    
    if(od_mode==1){
      od=order(as.numeric(celltype_cluster[,2]), decreasing = F)
    }else if(od_mode==2){
      od=order(as.numeric(celltype_cluster[,2]),celltype_cluster[,1], decreasing = F)
    }else{
      stop("od_mode should be 1 or 2.")
    }
    
    
    df2_cov=df2_cov[od,, drop=F]
    celltype_cluster=celltype_cluster[od,, drop=F]
    celltype_cluster=data.frame(celltype_cluster, stringsAsFactors = F)
    
    if(lab_mode!=FALSE){
      if(lab_mode=='annot'){
        celltype0=celltype[,-c(1,ncol(celltype)), drop=F]
      }else if(lab_mode=='corr'){
        celltype0=celltype[,-c(1,(ncol(celltype)-2)), drop=F]
        celltype0=celltype0[,c(1:(ncol(celltype0)-2), ncol(celltype0), (ncol(celltype0)-1)) , drop=F]
      }else{
        stop("Please specify a valid lab_mode mode!")
      }
    }else{
      celltype0=celltype[,-c(1,(ncol(celltype)-2),ncol(celltype)), drop=F]
    }
    
    rownames(celltype0)=celltype[,1]
    celltype0=celltype0[od,, drop=F]
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
        if((length(new_indc)+(length(which(!is.na(U0[,1])))))>20){
          col_use=c(colors()[c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)],
                    colors()[sample(1:length(colors())[-c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)],(length(new_indc)+length((which(!is.na(U0[,1])))))-20)])
        }else{
          col_use=c(colors()[c(609, 536, 62, 652, 611, 463, 498, 71, 258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])[1:(length(new_indc)+length(which(!is.na(U0[,1]))))]
        }
        names(col_use)=as.character(c(which(!is.na(U0[,1])), names(table(cols0))[new_indc]))
        col_use=col_use[which(names(col_use) %in% names(table(cols0)))]
        
        #}else if(colnames(celltype0)[ii]!='corr'){
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
        # }else{
        #   col_use=c("black", "firebrick")
      }
      #table(cols0)
      
      #tmp=1:length(table(cols0))
      #names(tmp)=paste0('c',names(table(cols0)))
      #celltype_df$celltype_df=as.character(tmp[paste0('c', cols0)])
      
      ## names(col_use)=names(table(cols0))
      ann_colors[[colnames(celltype0)[ii]]]=col_use
    }
    
    if(consensus==T){
      df2_colnames=colnames(df2_cov)
      df2_rownames=rownames(df2_cov)
      df2=Assign_obj$Uest[as.numeric(celltype[od,(ncol(celltype)-1), drop=T]),,drop=F]
      df2=apply(df2, c(2), function(x) pmin(x, maxv))
      colnames(df2)=df2_colnames
      rownames(df2)=df2_rownames
      
    }else{
      df2_colnames=colnames(df2_cov)
      df2_rownames=rownames(df2_cov)
      df2=sapply(1:ncol(Xir_cov), function(x) genotype_neighbor(cbind(Xir_cov[,x], Xir_allele[,x]), F, maxcp = 6))
      df2=df2[od,, drop=F]
      colnames(df2)=df2_colnames
      rownames(df2)=df2_rownames
      
    }
    
    col=c('#4d9efa','#0323a1','#9ecae1','#b0b0b0','#00d9ff',
          "#fff1ba","#ffb521","#DC7633","#BA4A00",
          "#fde0dd","#fcc5c0","#f768a1","#ae017e","#49006a",
          "#c7e9b4","#7fcdbb","#41b6c4","#41ab5d","#006d2c", "#000000",
          "#7B241C", "#7B241C", "#7B241C","#7B241C","#7B241C","#7B241C","#7B241C")
    
    breaksList = seq(0-(maxv-2), maxv, by = (maxv+(maxv-2))/100)
    
    if(mode=='segment'){
      
      pheatmap(df2,cluster_cols = F, cluster_rows = F,clustering_distance_rows = "correlation",clustering_method = "ward.D2",
               show_rownames = F, show_colnames = T, annotation_row = celltype0,annotation_colors = ann_colors, color = col,breaks = 0:26,
               fontsize = fontsize,  fontsize_row =fontsize_row , fontsize_col = fontsize_col,
               color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100),
               breaks=breaksList)
      
    }else{
      
      chr=sapply(strsplit(colnames(df2),"-"),'[',1)
      start=as.numeric(sapply(strsplit(colnames(df2),"-"),'[',2))
      end=as.numeric(sapply(strsplit(colnames(df2),"-"),'[',3))
      seglen=end-start
      repn=round(seglen/min(seglen))
      chr_rep=rep(chr, repn)
      chrn=table(chr_rep)[paste0('chr',1:22)]
      chrn=chrn[!is.na(chrn)]
      df2_plot=NULL
      for(ii in 1:length(repn)){
        mm=matrix(rep(df2[,ii], repn[ii]), ncol=repn[ii], byrow=F)
        df2_plot=cbind(df2_plot, mm)
      }
      rownames(df2_plot)=rownames(df2)
      coln=rep("", sum(repn))
      coln[ceiling(c(0,cumsum(chrn)[1:(length(chrn)-1)])+chrn/2)]=names(chrn)
      colnames(df2_plot)=coln
      
      pheatmap(df2_plot,cluster_cols = F, cluster_rows = F,clustering_distance_rows = "correlation",clustering_method = "ward.D2",
               show_rownames = F, show_colnames = T, annotation_row = celltype0, gaps_col = cumsum(chrn), annotation_colors = ann_colors, color = col,breaks = 0:26,
               fontsize = fontsize, fontsize_row =fontsize_row , fontsize_col = fontsize_col,
               color = colorRampPalette(rev(RColorBrewer::brewer.pal(n = 7, name ="RdBu")))(100),
               breaks=breaksList)
    }
    #   }else{
    #   celltype0=celltype[,(ncol(celltype)-1), drop=F]
    #   rownames(celltype0)=celltype[,1]
    #   celltype0=data.frame(celltype0, stringsAsFactors = F)
    #   pheatmap(df2,cluster_cols = F, cluster_rows = F,clustering_distance_rows = "correlation",clustering_method = "ward.D2", show_rownames = F, show_colnames = T, annotation_row = celltype0)
    # }
  }
  
  
  
}
