#' Generate Alleloscope object for analysis
#'
#' @param samplename Sample name for the data.
#' @param genome_assembly The genome assembly used for sequencing alignment. (ex: "GRCh38" or "GRCh37")
#' @param dir_path Path of the output directory.
#' @param size A matrix with two columns: col1: different chromosome; col2: for the size (bp) of different chromosomes.
#' @param assay A character indicating the type of sequencing data. (ex: "scDNAseq" or "scATACseq")
#'
#' @return A Alleloscope object including the necessary information.
#'
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
