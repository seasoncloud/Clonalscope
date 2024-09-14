#library(Seurat)
#library(BiocFileCache)
#library(GSEABase)
#library(AUCell)
#library(pheatmap)
#library(stringr)

#' Automatically obtains an initial estimation of normal cells/spots, when paired DNA is not available.
#'
#' @param counts A count matrix of transcirptomics data (gene x cell/spots)
#' @param gene_symbols Gene symbols, same dimension and order as the row of count matrix.
#' @param method Options from c("marker","pca")"marker" - marker gene based. "pca" - PCA based (similar to STARCH).
#' @param marker_genes Default NULL. A vector of manually input marker genes. If NULL, single cell signatures from MSigDB will be used.
#' @param marker_source Options from c("manual","MSigDB"). 
#' @param min_cells Number of minimum cells requried to filter genes (only used in PCA method).
#' @param min.cells.seurat Minimum cells required for Seurat filtering. Default value is 3.
#' @param min.features.seurat Minimum features required for Seurat preprocessing. Default value is 200.
#' @param dims.seurat Number of PCs for Seurat PCA and UMAP. Default 15.
#' @param resolution.seurat Clustering resolution for Seurat. Default 0.5.
#' @param nfeatures.seurat Number of variable genes for Seurat. Default 2000.
#' @param scale.facto.seurat Number of scaling factor for Seurat. Default 10000.
#' @param tissue_type Tissue type of the input sample. Available options from SCSig are 
#' "Cord_Blood","Esophagus","Stomach","Small_Intestine","Large_Intestine","PFC","Embryonic_CTX","Midbrain|Neuro",
#' "Bone_Marrow","Liver","Fetal_Kidney","Adult_Kidney","Fetal_Retina","Pancreas".
#' @param normal_celltype Keywords of possible normal celltypes, whose single cell marker genes will be 
#' searched in SCSig (for the same tissue type). Keywords shall be seperated by "|".
#' Example: normal_celltype=c("immune|endothelial|stromal|fibroblast"). Case insensitive. 
#' @param plot=F Plot out AUCell assignment score of each celltype.
#' @param save=F Logical Value. Whether to save Seurat object.
#' @param save_path Path to save Seurat object.
#' @return initial_normal_spots, a vector which contains identified barcodes/names of the normal cells/spots.
#'
#' @import Seurat,AUCell,BiocFileCache,GSEABase,pheatmap,stringr
#' @export
FindNormalReference <- function(counts, gene_symbols=NULL, method=c("marker","pca"),
                                marker_genes=NULL,marker_source=c("manual","MSigDB"),min_cells=NULL,
                                min.cells.seurat = 3, min.features.seurat = 200, dims.seurat=15,resolution.seurat=0.5,
                                nfeatures.seurat=2000, scale.facto.seurat=10000,
                                tissue_type=c("Cord_Blood","Esophagus","Stomach","Small_Intestine","Large_Intestine",
                                              "PFC","Embryonic_CTX","Midbrain|Neuro","Bone_Marrow","Liver","Fetal_Kidney",
                                              "Adult_Kidney","Fetal_Retina","Pancreas"),
                                normal_celltype=c("immune|endothelial|stromal|fibroblast"),
                                plot=F,save=F,save_path=NULL){
  message("Finding Normal Celltypes as Reference ...
  For customized input of normal reference cells, view the tutorial here: 
  https://github.com/seasoncloud/Clonalscope/tree/identify_normal_cells_noWGS/samples/V11Y04-378-A1")
  rownames(counts)=as.character(gene_symbols)
  
  if(method=="marker"){
    if(is.null(marker_genes) & (marker_source =="manual")){
      message("No marker genes provided, using default markers for immune,stromal and tumor cells.")
      # marker genes for normal cells
      marker_genes=list(immune=c("CD45","CD34","CD14","CD4","CD8",
                                 "CD19","CD68","SPP1","CD44","TIMP1",
                                 "TIMP3","CALD1"), 
                        stromal=c("FAP","PDGFRA","VIM","COL1A1","COL1A2",
                                  "COL5A1","LUM","TPM2","PDGFA","VWF"),
                        tumor = c("TP53"))
    }
    if(marker_source =="MSigDB"){
      # deciding initial normal clusters based on marker genesets from MSigDB 
      bfc <- BiocFileCache(ask=FALSE)
      scsig.path <- bfcrpath(bfc, file.path("http://software.broadinstitute.org",
                                            "gsea/msigdb/supplemental/scsig.all.v1.0.symbols.gmt"))
      scsigs <- getGmt(scsig.path)
    }
    
    message("Clustering by gene expression with Seurat...")
    # create seurat object
    seurat_obj <- CreateSeuratObject(counts = counts, project = "tutorial_normal",
                                     min.cells = min.cells.seurat, min.features = min.features.seurat)
    # pre-processing
    seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = scale.facto.seurat)
    seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = nfeatures.seurat)
    all.genes <- rownames(seurat_obj)
    seurat_obj <- ScaleData(seurat_obj, features = all.genes)
    # linear reduction
    seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
    seurat_obj <- FindNeighbors(seurat_obj, dims = 1:dims.seurat)
    seurat_obj <- FindClusters(seurat_obj, resolution = resolution.seurat)
    # UMAP
    seurat_obj <- RunUMAP(seurat_obj, dims = 1:dims.seurat)
    if(save){
      saveRDS(seurat_obj,paste0(save_path,"/seurat_obj.rds"))
    }
    
    # AUCell for celltype annotation
    counts.rankings <- AUCell_buildRankings(as.matrix(counts),
                                            plotStats=FALSE, verbose=FALSE)
    markers.aucs <- AUCell_calcAUC(scsigs, counts.rankings)
    # # Applying MsigDB genesets to data:
    # scsig.aucs <- AUCell_calcAUC(scsigs, counts.rankings)
    # scsig.results <- t(scsig.aucs@assays@data$AUC)
    # full.labels <- colnames(scsig.results)[max.col(scsig.results)]
    # names(full.labels) <- rownames(scsig.results)
    # tab <- table(full.labels[names(Idents(seurat_obj))], Idents(seurat_obj))
    # fullheat <- pheatmap(log10(tab+10), color=viridis::viridis(100), silent=TRUE)
    # 
    # # Restricting to the subset of tissue type:
    # scsigs.sub <- scsigs[grep("Immune|Endothelial|Stromal", names(scsigs))]
    scsigs.sub <- scsigs[grep(tissue_type, names(scsigs))]
    sub.aucs <- AUCell_calcAUC(scsigs.sub, counts.rankings)
    sub.results <- t(sub.aucs@assays@data$AUC)
    sub.labels <- colnames(sub.results)[max.col(sub.results)]
    names(sub.labels) <- rownames(sub.results)
    tab <- table(sub.labels[names(Idents(seurat_obj))], Idents(seurat_obj))
    subheat <- pheatmap(log10(tab+10), color=viridis::viridis(100), silent=TRUE)
    if(plot){
      print(subheat)
    }
    
    normal_info = tab[grep(normal_celltype,rownames(tab),ignore.case = T),]
    if(is.null(dim(normal_info))){ # only one celltype found
      normal_cluster = colnames(tab)[which.max(normal_info)]
    }else{
      normal_cluster = unique(colnames(tab)[sapply(1:dim(normal_info)[1],
                                                   function(i){order(normal_info[i,],decreasing = T)[1:2]})])
    }
    
    initial_normal_spots = names(Idents(seurat_obj))[Idents(seurat_obj) %in% normal_cluster]
    return(initial_normal_spots)
    
  }else if(method=="pca"){
    message("Identifying PCs..")
    initial_normal_spots=starch_normal(mtx=counts,min_cells=NULL,min_umi_perspot=10,max_value=3)
  }
  
  return(initial_normal_spots)
}

#' This function clusters cells/spots based on the first principal component of gene expression.
#' It is a R equivalent of the method that STARCH identifies normal spots.
#'
#' @param mtx A count matrix of transcirptomics data (gene x cell/spots)
#' @param min_cells Number of minimum cells required to filter genes.
#' @param min_umi_perspot Number of minimum UMI required per spot/cell. Default 10.
#' @param max_value Capping value for log normalzied gene expression.
#' @return normal_spots, a vector which contains identified barcodes/names of the normal cells/spots.
#'
#' @import
#' 
starch_normal <- function(mtx=NULL,min_cells=NULL,min_umi_perspot=10,max_value=3){
  # minimum cells requried for each gene
  if(is.null(min_cells)){
    min_cells = dim(mtx)[2]/20
  }
  # filter genes
  gene_idx = rowSums(mtx) > min_cells
  # filter spots
  spots_idx = colSums(mtx) > min_umi_perspot
  data= mtx[gene_idx,spots_idx]
  # normalize by library size
  m = median(colSums(data))
  data = t(t(data)/colSums(data)) *m
  data = log(data+1)
  # threshold data
  data= pmin(data,max_value)
  data=pmax(data,-max_value)
  # PCA (SVD) on the filtered data
  pca.res <- svd(t(data))
  cell_pc_1 = pca.res$u[,1]
  # Kmeans clustering with K=2, on data with 1st principal component
  kmeans.res <- kmeans(cell_pc_1,centers=2)
  cluster1_mean = mean(data[,colnames(data)[kmeans.res$cluster == 1]],na.rm=T)
  cluster2_mean = mean(data[,colnames(data)[kmeans.res$cluster == 2]],na.rm=T)
  if(cluster1_mean > cluster2_mean){
    normal_spots = colnames(data)[kmeans.res$cluster == 2] 
  }else{
    normal_spots = colnames(data)[kmeans.res$cluster == 1] 
  }
  
  return(normal_spots)
}

#' This function plots celltype/clustering results for spatial transcriptomics datasets.
#' @param spot_data Standard 10X VISIUM spatial data frame, columns being barcode, include, row, col, image_row, image_col. 
#' @param celltype A dataframe containing two columns - barcode and clusters/celltype
#' @param save Logical. Whether to save the plot.
#' @param output_path Plot saving path.
#' @param title Plot Title.
#' @param cluster Cluster mode. If TRUE, clonalscope cluster will be plotted. If FALSE, plotting celltypes.
#' @param plot_colors A vector of plotting colors. Default NULL, default colors will be used.
#' @param pt_size UMAP point size. Default size 3. 
#' 
#' @import scales, ggplot2
#' @export
SpatialPlot <- function(spot_data,celltype,save=F,output_path=NULL,title="",cluster=F,plot_colors=NULL,pt_size=3){
  celltype_temp=celltype;rownames(celltype_temp) = celltype_temp[,1] # barcodes as rownames
  plot_df= spot_data[spot_data$include == 1, c(3,4)]
  if(cluster){
    plot_df$celltype = "Non-Tumor"
    if(is.null(plot_colors)){
      plot_colors = c(hue_pal()(length(levels(as.factor(celltype[,2])))),"grey")
      names(plot_colors) <- c(levels(as.factor(as.numeric(celltype[,2]))),"Non-Tumor")
    }
  }else{
    plot_df$celltype = "Unknown"
    if(is.null(plot_colors)){
      plot_colors = c(hue_pal()(length(levels(as.factor(celltype[,2])))),"grey")
      names(plot_colors) <- c(levels(as.factor(celltype[,2])),"Unknown")
    }
  }
  plot_df[celltype_temp[,1], "celltype"] = celltype_temp[,2]
  colnames(plot_df) = c("x","y","celltype")
  plot_df$celltype = factor(plot_df$celltype,levels=names(plot_colors))
  
  
  g<- ggplot(plot_df,aes(x=y,y=-x,color=celltype)) + 
    geom_point(size=pt_size) +
    ggtitle(title) +
    #xlim(-max(spot_data[,3]),min(spot_data[,3]))+ylim(-max(spot_data[,4]),-min(spot_data[,4]))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"),
                       plot.title = element_text(hjust = 0.5,size=20),
                       legend.text = element_text(size=15),
                       legend.title = element_text(size=15)) +
    scale_colour_manual(values=plot_colors)
  if(save){
    png(paste0(output_path,"/",title,".png"),height=800,width=800)
    print(g)
    dev.off()
  }
  return(g)
}

#' This function calculates cosine similarity of a CNV profile 1 with CNV profile 2.
#'
cosine_similarity <- function(cnv_1,cnv_2){
  cos_sim = sum((cnv_1 - 1)*(cnv_2-1))/(sqrt(sum((cnv_1-1)^2))*sqrt(sum((cnv_2-1)^2)))
  return(cos_sim)
}

#' This function assigns the final tumor/normal identity.
#' The clusters with smallest CNV loads are chosen as normal reference.
#' All the rest clusters are compared with the 
#' @param Cov_obj A Clonalscope Object that contains clustering results.
#' @param cutoff Numeric value from [-1,1]. Cosine Similarity threshold, if lower, assigned as normal cells.
#' @return An object with CNV loads, cosine similarity and final assignment of tumor/normal status.
#'
#' @import 
#' @export
MalignantAssignment<- function(Cov_obj,cutoff=0.5){
  celltype= Cov_obj$celltype0 # input celltype
  result=Cov_obj$result_final$result
  cnv_df= Cov_obj$result_final$df_obj$df
  # calculate average CNV signals 
  cluster_cnv_load = sapply(sort(unique(Cov_obj$result_final$result$Zest)),function(c){
    cell_bcs = names(Cov_obj$result_final$result$Zest)[Cov_obj$result_final$result$Zest == c]
    median_cnv =  apply(cnv_df[intersect(rownames(cnv_df),cell_bcs),],2,median)
    median_cnv 
  })
  cluster_cnv_load  = t(cluster_cnv_load)
  rownames(cluster_cnv_load) = sort(unique(Cov_obj$result_final$result$Zest))
  
  cnv_load = sapply(1:dim(cluster_cnv_load)[1],function(i){sum((cluster_cnv_load[i,] - 1)^2)/length(cluster_cnv_load[i,])})
  names(cnv_load) = rownames(cluster_cnv_load)
  
  normal_cluster =rownames(cluster_cnv_load)[order(cnv_load)[1:2]]#)[which.min(cnv_load)] 
  
  cluster_cos_sim = sapply(1:dim(cluster_cnv_load)[1],function(i){
    cnv_1=cluster_cnv_load[i,]
    cell_bcs = names(Cov_obj$result_final$result$Zest)[Cov_obj$result_final$result$Zest %in% normal_cluster]
    cnv_2 = apply(cnv_df[intersect(rownames(cnv_df),cell_bcs),],2,median)
    #cnv_2=cluster_cnv_load[normal_cluster,]
    cosine_similarity(cnv_1,cnv_2)
  })
  names(cluster_cos_sim) = rownames(cluster_cnv_load)
  
  final_assignment=Cov_obj$result_final$result$Zest 
  final_assignment[final_assignment %in% names(cluster_cos_sim)[cluster_cos_sim <= cutoff]]= "Tumor" 
  final_assignment[final_assignment %in% names(cluster_cos_sim)[cluster_cos_sim > cutoff]]= "Normal" 
  return(list(cluster_cnv_load=cluster_cnv_load,cnv_load=cnv_load,
              normal_cluster=normal_cluster,
              cluster_cos_sim=cluster_cos_sim,final_assignment=final_assignment))
}
