#' Generate bin by cell matrix for scATAC-seq data
#'
#' @param bin_bed A matrix of the BED format. The first three columns are 'chr', "start site", "end site". Each row is a bin region.
#' @param barcodes A matrix/ data.frame with barcodes for each cell in the first column.
#' @param path_to_matrix The path to the "matrix.mtx" file with the same format as that from the Cell Ranger software.
#' @param path_to_genes The path to the "genes.tsv" file with the same format as that from the Cell Ranger software.
#' @param path_to_gtf The path to the gtf file with the same format as that from the Cell Ranger software.
#'
#' @import Matrix
#' @import data.table
#' @import rtracklayer
#' @return A Alleloscope object including the necessary information.
#'
#' @export
## generate bed files for the bins
Gen_bin_cell_rna_filtered=function(bin_bed=NULL, barcodes=NULL, gene_matrix=NULL, genes=NULL, gtf=NULL){
  # generate bin by cell matrix from fragment file
  barcodes=barcodes$V1

  cat("Read gene by cell matrix...\n")
  #gene_matrix=readMM(path_to_matrix)
  #genes=read.table(path_to_genes,sep='\t', header=F, stringsAsFactors = F)

  cat(paste0("Total ", length(barcodes)," cells.\n"))
  cat(paste0("Total ", nrow(genes)," genes.\n"))

  ## check bin bed chr column
  if(grepl("chr",bin_bed[1,1])){
    bin_bed[,1]=gsub("chr","",bin_bed[,1])
  }else{
    bin_bed=bin_bed
  }

  chr200k=bin_bed
  chr200k=chr200k[order(as.numeric(chr200k[,1]), as.numeric(chr200k[,2])),]
  bins=paste0('chr',chr200k[,1],':',chr200k[,2],"_", chr200k[,3])
  query=GRanges(paste0('chr',chr200k[,1]), IRanges(chr200k[,2]+1,chr200k[,3]))


  #gtf=read.table(path_to_gtf, sep='\t', stringsAsFactors = F)
  gtf=gtf[which(gtf$V3=='gene'),]
  ensg=sapply(strsplit(gtf$V9,';'),'[',1)
  ensg=sapply(strsplit(ensg,'gene_id '),'[',2)
  chr=paste0(gtf$V1,'-',gtf$V4,'-', gtf$V5)
  names(chr)=ensg
  genes$site=chr[genes[,1]]
  genes[is.na(genes[,4]),4]="0-0-0"

  if(grepl("chr",genes$site[1])){
    subject=GRanges(paste0(sapply(strsplit(genes$site,'-'),'[',1)), IRanges(as.numeric(sapply(strsplit(genes$site,'-'),'[',2)), as.numeric(sapply(strsplit(genes$site,'-'),'[',3))))
  }else{
    subject=GRanges(paste0('chr',sapply(strsplit(genes$site,'-'),'[',1)), IRanges(as.numeric(sapply(strsplit(genes$site,'-'),'[',2)), as.numeric(sapply(strsplit(genes$site,'-'),'[',3))))
  }


  ov=findOverlaps(query, subject )
  ov=as.matrix(ov)
  #tmp=fragments_incell$barcode[ov[,2]]
  #ov=cbind(ov,match(tmp, barcodes))

  cat("Generate bin-by-cell matrix...")
  mtx=as.data.table(as.matrix(gene_matrix))

  chr200_mtx=matrix(nrow=dim(chr200k)[1], ncol=dim(mtx)[2])
  for(ii in 1:dim(chr200k)[1]){
    peak_ind=ov[which(ov[,1]==ii),2]
    gcount=apply(mtx[peak_ind,, drop=FALSE],2, sum)
    chr200_mtx[ii, ]=gcount
    #cat(ii)
  }

  colnames(chr200_mtx)=barcodes
  rownames(chr200_mtx)=paste0("chr",chr200k[,1],"-", chr200k[,2],"-", chr200k[,3])
  chr200_mtx = Matrix(chr200_mtx, sparse = TRUE)

  message("The bin-by-cell matrix has beed successfully generated!")
  return(chr200_mtx)
}
