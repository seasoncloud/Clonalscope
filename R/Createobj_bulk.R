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
Createobj_bulk=function(samplename='sample',genome_assembly="GRCh38", dir_path='./', size=NULL, assay='WGS'){

  # check parameters
  if(!(nrow(size)>0 & ncol(size)==2)){
    stop("Please provide a matrix/ data.frame with two columns: col1: different chromosome; col2: for the size (bp) of different chromosomes.")
  }

  dir.create(paste0(dir_path,"/plots"))

    if(grepl('chr', as.character(size[2,1]))){
      size[,1]=sapply(strsplit(size[,1],'hr'),'[',2)##if chr
    }
    size=size[which(size[,1] %in% as.character(1:22)),]
    size=size[order(as.numeric(as.character(size$V1))),]
    size_name=size[,1]
    size=as.numeric(as.character(size[,2]))
    names(size)=size_name

    output=list("size"=size,"samplename"=samplename,
                "dir_path"=dir_path,"genome_assembly"=genome_assembly,
                #"cell_info" = cell_info,
                "seg_table"=NULL, "seg_table_filtered"=NULL,  "assay"=assay)


  message("Object successfully created!")
  return(output)

}
