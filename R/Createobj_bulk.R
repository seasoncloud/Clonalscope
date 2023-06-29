#' Generate Alleloscope object for analysis
#'
#' @param raw_counts A binned coverage matrix (m1 bin by n1 cell) with values being read counts in DNA sequencing data for all chromosomal regions of tumor sample. n1 can be 1 for bulk sample. The rownames should be formatted as: chr1-1-20000.
#' @param ref_counts A binned coverage matrix (m2 bin by n2 cell) with values being read counts in DNA sequencing data for all chromosomal regions of normal sample. n2 can be 1 for bulk sample. The rownames should be formatted as: chr1-1-20000."
#' Numbers of bins (rows) should be the same in the paired chromosomal regions for the paired samples
#' @param samplename Sample name for the data.
#' @param genome_assembly The genome assembly used for sequencing alignment. (ex: "GRCh38" or "GRCh37")
#' @param dir_path Path of the output directory.
#' @param size A matrix with two columns: col1: different chromosome; col2: for the size (bp) of different chromosomes.
#' @param assay A character indicating the type of sequencing data. (ex: "scDNAseq" or "scATACseq")
#'
#' @return A Alleloscope object including the necessary information.
#'
#' @import stringr
#' @export
Createobj_bulk=function(raw_counts=NULL, ref_counts=NULL,samplename='sample',genome_assembly="GRCh38", dir_path='./', size=NULL, assay='WGS'){

  # check parameters
  if(!(nrow(size)>0 & ncol(size)==2)){
    stop("Please provide a matrix/ data.frame with two columns: col1: different chromosome; col2: for the size (bp) of different chromosomes.")
  }else if(length(unlist(strsplit(rownames(raw_counts)[2],'-')))!=3){
    stop("The rownames for the raw_counts matrix should be formatted as: chr1-1-20000.")
  }else if(length(unlist(strsplit(rownames(ref_counts)[2],'-')))!=3){
    stop("The rownames for the ref_counts matrix should be formatted as: chr1-1-20000.")
  }else if(!(nrow(raw_counts)>0 & ncol(raw_counts)>0)){
    stop("raw_counts matrix is not valid.")
  }else if(!(nrow(ref_counts)>0 & ncol(ref_counts)>0)){
    stop("ref_counts matrix is not valid.")
  }
  # check coverage matrix's order, automatically sorting autosomes
  wgst_order = sapply(strsplit(unique(sapply(strsplit(rownames(raw_counts),"-"),"[[",1)),"chr"),"[[",2)
  if(any(wgst_order != wgst_order[order(as.numeric(as.character(wgst_order)))])){
    message("WGS of tumor sample not ordered, reordering numerically...")
    chr_names = sapply(strsplit(rownames(raw_counts),"-"),"[[",1)
    chr_start=sapply(strsplit(rownames(raw_counts),"-"),"[[",2)
    if(grepl('chr', chr_names[1])){
      chr_names=sapply(strsplit(chr_names,'hr'),'[',2)##if chr
    }
    raw_counts=raw_counts[order(as.numeric(chr_names),as.numeric(chr_start)),,drop=F]
  }
  wgsn_order = sapply(strsplit(unique(sapply(strsplit(rownames(ref_counts),"-"),"[[",1)),"chr"),"[[",2)
  if(any(wgsn_order != wgsn_order[order(as.numeric(as.character(wgsn_order)))])){
    message("WGS of normal sample not ordered, reordering numerically...")
    chr_names = sapply(strsplit(rownames(ref_counts),"-"),"[[",1)
    chr_start=sapply(strsplit(rownames(ref_counts),"-"),"[[",2)
    if(grepl('chr', chr_names[1])){
      chr_names=sapply(strsplit(chr_names,'hr'),'[',2)##if chr
    }
    ref_counts=ref_counts[order(as.numeric(chr_names),as.numeric(chr_start)),,drop=F]
  }

  
  
  dir.create(paste0(dir_path,"/plots"))

    chrr=sapply(strsplit(rownames(raw_counts),"-"),'[',1)
    chrr=gsub("chr", "", chrr)

    if(grepl('chr', as.character(size[2,1]))){
      size[,1]=sapply(strsplit(size[,1],'hr'),'[',2)##if chr
    }
    #size=size[which(size[,1] %in% as.character(1:22)),]
    size=size[which(size[,1] %in% chrr),]
    size=size[order(as.numeric(as.character(size[,1]))),]
    size_name=size[,1]
    size=as.numeric(as.character(size[,2]))
    names(size)=size_name

    output=list(raw_counts=raw_counts, ref_counts=ref_counts,"size"=size,"samplename"=samplename,
                "dir_path"=dir_path,"genome_assembly"=genome_assembly,
                #"cell_info" = cell_info,
                "seg_table"=NULL, "seg_table_filtered"=NULL,  "assay"=assay)


  message("Object successfully created!")
  return(output)

}
