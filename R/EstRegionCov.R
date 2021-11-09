#' Estimate coverage change for each cell in each segments (the segmentation results from WGS/WES) using a set of normal cells under a Poisson model.
#'
#' @param mtx A (sparse) matrix with each row being a feature and each column being a cell barcode. (Output from cellranger)
#' @param barcodes A matrix/ data.frame with barcodes for each cell in the first column.
#' @param features A matrix/ data.frame with the columns- 1st: gendID; 2nd: gene name. Each row is a feature whose order should correspond to the rows of "mtx".
#' @param bed A matrix with the bed file format. COL1- chromosome; COL2- start coordinate; COL3- end coordinate; COL4- the gene ID matches the one in "features".
#' @param celltype0 A matrix with two columns: COL1- cell barcodes; COL2- cell types (Tumor cells should be labeled as "tumor" and the set of cells used as control should be labeled as "normal").
#' @param var_pt A decimal between [0,1]. Exclude the extremely variable genes with variance > the var_pt percentile for the tested cells.
#' @param var_pt_ctrl A decimal between [0,1]. Exclude the extremely variable genes with variance > the var_pt percentile for the control cells.
#' @param include The cells included for testing. The value should be one of ["tumor","all","control"].
#' @param alpha_source The regions used to compute normalizing factors. The value should be one of ["all", "control"]. If "all", total read counts are used; if "control", a subset of chromosomes are used (ctrl_region needs to be specified).
#' @param ctrl_region A vector indicating the chromosomes used to compute normalizing factors. The format should be c("chr1", "chr3", ...)
#' @param seg_table_filtered A segmentation table from the "Segmentation_bulk" function.
#' @param size A matrix with two columns: col1: different chromosome; col2: for the size (bp) of different chromosomes.
#' @param plot_path A path for plotting histograms.
#' @param breaks A parameter used in the "hist" function.
#'
#' @return A list with "deltas_all" and "ngenes". "deltas all" stores values of estimated coverage change for each cell and each region. "ngenes" shows number of genes in each region.
#'
#' @import rtracklayer
#' @import Matrix
#' @export
EstRegionCov=function(mtx=NULL, barcodes=NULL, features=NULL, bed=NULL, celltype0=NULL, var_pt=0.99, var_pt_ctrl=0.99,ngene_filter=0, include='tumor',
                      alpha_source='all', ctrl_region=NULL, seg_table_filtered=NULL,size=NULL, plot_path=NULL, breaks=30){
  ## test
  ## check bed with "chr"
  if(!grepl("chr",bed[1,1])){
    bed[,1]=paste0('chr', bed[,1])
  }

  if(is.null(celltype0)){
    celltype0=cbind(barcodes[,1], rep("normal", nrow(barcodes)))
  }

  colnames(mtx)=barcodes[,1]
  rownames(mtx)=features[,1]
  celltypes=celltype0

  mtx=mtx[,match(celltypes[,1], colnames(mtx) )]

  rna_var=apply(mtx,1, var)
  ngenes=Matrix::rowSums(mtx)

  #rna=mtx[which(rna_var<quantile(rna_var[which(rna_var!=0)],var_pt) & (rna_var!=0) & ngenes>ngene_filter), ]
  #rna_control=mtx[which(rna_var<quantile(rna_var[which(rna_var!=0)],var_pt_ctrl) & (rna_var!=0) & ngenes>ngene_filter), ]
  rna=mtx[which(rna_var<quantile(rna_var,var_pt) & (rna_var!=0) & ngenes>ngene_filter), ]
  rna_control=mtx[which(rna_var<quantile(rna_var,var_pt_ctrl) & (rna_var!=0) & ngenes>ngene_filter), ]


  bed_sub=bed[match(rownames(rna), bed[,4]),, drop=F]
  bed_sub[is.na(bed_sub)]=0


  #include # 'control', "all" for cells "include"
  #alpha_source='all'#'all' #'control'
  #ctrl_region=paste0('chr', c(1,4,6,9,10,11,12)) # if (alpha_source==control)


  est_name=paste0('pt', var_pt, "_",include,"_alpha",alpha_source)#'deltas_rnavar_alphaall_cov_'

  if(is.null(seg_table_filtered)){
    message("Estimation for each chromosome")
    seg_table_filtered0=data.frame("chr"=gsub("chr","",size[1:22,1]), 'start'=0, 'end'=as.numeric(size[1:22,2]),
                                   'states'=0, 'length'=as.numeric(size[1:22,2]),'mean'=0, 'var'=0, 'Var1'=1:22,'Freq'=50000,'chrr'=size[1:22,1], stringsAsFactors = F)
  }

  deltas_all=list()
  ngenes=c()

  if(!is.null(plot_path)){
    pdf(plot_path, width = 8, height = 6)
    par(mfrow=c(2,3))
  }

  # index for subsetting the seg_table_filtered
  sel_ind=c()

  for(chrr in paste0(seg_table_filtered$chrr)){
    #chrr=c('chr8')
    #gene_ind=which(bed_sub$V1 %in% chrr)
    query=GRanges(seqnames = paste0('chr',seg_table_filtered$chr[which(seg_table_filtered$chrr==chrr)]),
                  ranges = IRanges(as.numeric(seg_table_filtered$start[which(seg_table_filtered$chrr==chrr)]),
                                   as.numeric(seg_table_filtered$end[which(seg_table_filtered$chrr==chrr)])))
    subject=GRanges(seqnames = bed_sub[,1], ranges=IRanges(start=as.numeric(bed_sub[,2]), end=as.numeric(bed_sub[,3])))
    ov=as.matrix(findOverlaps(query=query, subject = subject))
    gene_ind=ov[,2]


    rna_sub=rna[gene_ind,, drop=F]

    sel_cell=which(Matrix::colSums(rna_sub)>0)



    if(length(sel_cell)==0){
      next
    }

    sel_ind=c(sel_ind, chrr)
    rna_sub=rna_sub[,sel_cell, drop=F]


    celltypes=celltype0
    celltypes=celltypes[match(colnames(rna), celltypes[,1]),]
    celltypes=celltypes[sel_cell,, drop=F]

    # subset cells
    # sel_cell=which(colSums(rna_sub)<=quantile(colSums(rna_sub), 0.99) & colSums(rna_sub)>=quantile(colSums(rna_sub), 0.01))
    # rna_sub=rna_sub[,sel_cell]
    # celltypes=celltypes[sel_cell,]

    #
    controlCounts=as.matrix(rna_sub[,which(celltypes[,2]=='normal'), drop=F])
    if(include=='tumor'){
      testCounts=as.matrix(rna_sub[,which(celltypes[,2]=='tumor'), drop=F])
    }else if(include=='control'){
      testCounts=controlCounts
    }else if(include=='all'){
      testCounts=as.matrix(rna_sub)
    }else{
      stop("include variable is not valid!")
    }



    # compute cell size from control region
    if(alpha_source=='all'){
      #alphas=colSums(rna[,sel_cell])
      alphas=Matrix::colSums(rna_control[,sel_cell, drop=F])
    }else if(alpha_source=='control'){
      # from control region
      rna_chr1=rna[which(bed_sub[,1] %in% ctrl_region), , drop=F]
      alphas=Matrix::colSums(rna_chr1[,sel_cell, drop=F])
    }else{
      stop("alpha_source is not valid!")
    }

    alpha_controls=alphas[which(celltypes[,2]=='normal')]/median(alphas)  #mean

    if(include=='tumor'){
      alpha_tests=alphas[which(celltypes[,2]=='tumor')]/median(alphas) #mean
    }else if(include=='control'){
      alpha_tests=alpha_controls
    }else if(include=='all'){
      alpha_tests=alphas/median(alphas)
    }else{
      stop("include variable is not valid!")
    }

    testCounts=testCounts/matrix(rep(alpha_tests, nrow(testCounts)), nrow=nrow(testCounts), byrow = T)
    controlCounts=controlCounts/matrix(rep(alpha_controls, nrow(controlCounts)), nrow=nrow(controlCounts), byrow = T)

    # estimate delta for cell_ir
    # for one region r

    deltas=Matrix::colSums(testCounts)/median(Matrix::colSums(controlCounts))

    names(deltas)=colnames(testCounts)

    if(!is.null(plot_path)){
      celltype_tmp=celltypes[match(names(deltas), celltypes[,1]),2]
      if(length(which(celltype_tmp=='tumor'))>1 & length(which(celltype_tmp=='normal'))>1){
        #c2 <- rgb(173,216,230,max = 255, alpha = 50, names = "lt.blue")
        c2 <- rgb(56,101,195,max = 255, alpha = 50, names = "lt.blue")
        #c1 <- rgb(255,192,203, max = 255, alpha = 50, names = "lt.pink")
        c1 <- rgb(239,90,155, max = 255, alpha = 50, names = "lt.pink")
        ha=hist(deltas[which(celltype_tmp=='tumor')],breaks = breaks, plot=FALSE)
        if(length(unique(celltype_tmp))>1){
          hb=hist(deltas[which(celltype_tmp=='normal')],breaks = breaks, plot = FALSE)
        }
        plot(ha, xlim=c(0,3), main=paste0(chrr, ":", nrow(controlCounts), " genes; ", length(deltas), " cells"), col=c1)
        if(length(unique(celltype_tmp))>1){
          plot(hb, add=TRUE, col=c2)
        }
      }
    }
    ##
    #plot(colSums(rna[,which(celltypes[,2]=='tumor')]), deltas, xlim=c(0,1000))

    #km=kmeans(deltas,2)

    #saveRDS(deltas, "../scRNACNV_analysis/rds/p5931t_chr21_all_deltas.rds")
    #assign(paste0(est_name,chrr),deltas)

    deltas_all[[paste0(est_name,chrr)]]=deltas
    ngene=nrow(controlCounts)
    names(ngene)=paste0(est_name,chrr)
    ngenes=c(ngenes, ngene)
    #assign(paste0("deltas_",chrr),deltas)
    print(chrr)
  }

  if(!is.null(plot_path)){
    dev.off()
  }

  seg_table_filtered=seg_table_filtered[which(seg_table_filtered$chrr %in% sel_ind),, drop=F]

  return(list(deltas_all=deltas_all, ngenes=ngenes, seg_table_filtered=seg_table_filtered, alpha=alphas))
}


