##
# The following code is fully based on and modified from CopyKAT github: https://github.com/navinlabcode/copykat/blob/master/R/copykat.R
# The goal is to identify normal cells from copykat step 1-4.
#
copykat_normal <- function(rawmat=rawdata,id.type="E", cell.line="no",ngene.chr=5, win.size=25, KS.cut=0.1, LOW.DR=0.05, UP.DR=0.1,
                           distance="euclidean",  genome="hg20",norm.cell.names="", n.cores=4, output.seg="FALSE",force_WNS=F){
  set.seed(42)

  print("CopyKAT: step1: read and filter data ...")
  print(paste(nrow(rawmat), " genes, ", ncol(rawmat), " cells in raw data", sep=""))
  
  genes.raw <- apply(rawmat, 2, function(x)(sum(x>0)))
  
  if(sum(genes.raw> 200)==0) stop("none cells have more than 200 genes")
  if(sum(genes.raw<100)>1){
    rawmat <- rawmat[, -which(genes.raw< 200)]
    print(paste("filtered out ", sum(genes.raw<=200), " cells with less than 200 genes; remaining ", ncol(rawmat), " cells", sep=""))
  }
  ##
  der<- apply(rawmat,1,function(x)(sum(x>0)))/ncol(rawmat)
  
  if(sum(der>LOW.DR)>=1){
    rawmat <- rawmat[which(der > LOW.DR), ]; print(paste(nrow(rawmat)," genes past LOW.DR filtering", sep=""))
  }
  
  WNS1 <- "data quality is ok"
  if(nrow(rawmat) < 7000){
    WNS1 <- "low data quality"
    UP.DR<- LOW.DR
    print("WARNING: low data quality; assigned LOW.DR to UP.DR...")
  }
  
  print("step 2: annotations gene coordinates ...")
  if(genome=="hg20"){
    anno.mat <- annotateGenes.hg20(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
  } else if(genome=="mm10"){
    anno.mat <- annotateGenes.mm10(mat = rawmat, ID.type = id.type) #SYMBOL or ENSEMBLE
    dim(rawmat)
  }  
  ### secondary filtering
  ToRemov2 <- NULL
  for(i in 8:ncol(anno.mat)){
    cell <- cbind(anno.mat$chromosome_name, anno.mat[,i])
    cell <- cell[cell[,2]!=0,]
    if(length(as.numeric(cell))< 5){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    } else if(length(rle(cell[,1])$length)<length(unique((anno.mat$chromosome_name)))|min(rle(cell[,1])$length)< ngene.chr){
      rm <- colnames(anno.mat)[i]
      ToRemov2 <- c(ToRemov2, rm)
    }
    i<- i+1
  }
  
  if(length(ToRemov2)==(ncol(anno.mat)-7)) stop("all cells are filtered")
  if(length(ToRemov2)>0){
    anno.mat <-anno.mat[, -which(colnames(anno.mat) %in% ToRemov2)]
  }
  rawmat3 <- data.matrix(anno.mat[, 8:ncol(anno.mat)])
  norm.mat<- log(sqrt(rawmat3)+sqrt(rawmat3+1))
  norm.mat<- apply(norm.mat,2,function(x)(x <- x-mean(x)))
  colnames(norm.mat) <-  colnames(rawmat3)
  
  ##smooth data
  print("step 3: smoothing data with dlm ...")
  dlm.sm <- function(c){
    model <- dlm::dlmModPoly(order=1, dV=0.16, dW=0.001)
    x <- dlm::dlmSmooth(norm.mat[, c], model)$s
    x<- x[2:length(x)]
    x <- x-mean(x)
  }
  test.mc <-parallel::mclapply(1:ncol(norm.mat), dlm.sm, mc.cores = n.cores)
  norm.mat.smooth <- matrix(unlist(test.mc), ncol = ncol(norm.mat), byrow = FALSE)
  # added for debug
  colnames(norm.mat.smooth) = colnames(norm.mat)
  rownames(norm.mat.smooth) = rownames(norm.mat)
  
  if (cell.line=="yes"){
    print("running pure cell line mode")
    relt <- baseline.synthetic(norm.mat=norm.mat.smooth, min.cells=10, n.cores=n.cores)
    norm.mat.relat <- relt$expr.relat
    CL <- relt$cl
    WNS <- "run with cell line mode"
    preN <- NULL
    
  } else if(length(norm.cell.names)>1){
    
    #print(paste(length(norm.cell.names), "normal cells provided", sep=""))
    NNN <- length(colnames(norm.mat.smooth)[which(colnames(norm.mat.smooth) %in% norm.cell.names)])
    print(paste(NNN, " known normal cells found in dataset", sep=""))
    
    if (NNN==0) stop("known normal cells provided; however none existing in testing dataset")
    print("run with known normal...")
    
    basel <- apply(norm.mat.smooth[, which(colnames(norm.mat.smooth) %in% norm.cell.names)],1,median); print("baseline is from known input")
    
    d <- parallelDist::parDist(t(norm.mat.smooth),threads =n.cores, method="euclidean") ##use smooth and segmented data to detect intra-normal cells
    
    km <- 6
    fit <- hclust(d, method="ward.D2")
    CL <- cutree(fit, km)
    
    while(!all(table(CL)>5)){
      km <- km -1
      CL <- cutree(fit, k=km)
      if(km==2){
        break
      }
    } 
    WNS <- "run with known normal"
    preN <- norm.cell.names
    ##relative expression using pred.normal cells
    norm.mat.relat <- norm.mat.smooth-basel
  }else {
    basa <- baseline.norm.cl(norm.mat.smooth=norm.mat.smooth, min.cells=5, n.cores=n.cores)
    basel <- basa$basel
    WNS <- basa$WNS
    preN <- basa$preN
    CL <- basa$cl
    if(force_WNS == T){
      
    }
    if (WNS =="unclassified.prediction"){
      basa <- baseline.GMM(CNA.mat=norm.mat.smooth, max.normal=5, mu.cut=0.05, Nfraq.cut=0.99,RE.before=basa,n.cores=n.cores)
      basel <-basa$basel
      WNS <- basa$WNS
      
      preN <- basa$preN
      
    }
  }
  return(preN)
}

# This function is a R equiavalent of the method that STARCH identifies normal spots 
# (same as filter_spots + get_normal_spots functions in STARCH).
# The filtering parameters are set to be the same as STARCH
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

SpatialPlot <- function(spot_data,celltype,save=F,output_path=NULL,title=""){
  celltype_temp=celltype;rownames(celltype_temp) = celltype_temp[,1] # barcodes as rownames
  plot_df= spot_data[spot_data$include == 1, c(3,4)]
  plot_df$celltype = "Unknown"
  plot_df[celltype_temp[,1], "celltype"] = celltype_temp[,2]
  colnames(plot_df) = c("x","y","celltype")
  g<- ggplot(plot_df,aes(x=y,y=-x,color=celltype)) + 
    geom_point(size=1) +
    ggtitle(title) +
    theme(plot.title = element_text(hjust = 0.5,size=20)) +
    #xlim(-max(spot_data[,3]),min(spot_data[,3]))+ylim(-max(spot_data[,4]),-min(spot_data[,4]))+
    theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
  if(save){
    png(paste0(output_path,"Spatial_Distribution.png"),height=500,width=500)
    print(g)
    dev.off()
  }
  return(g)
}
