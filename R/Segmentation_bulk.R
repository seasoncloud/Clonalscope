#' HMM segmentation based on coverage matrix for paired tumor and normal sample.
#'
#' If there is no paired normal, other normal sample with the same genome coordinate also works.
#'
#' @param Obj_filtered An Alleloscope object.
#' @param plot_seg Logical (TRUE/ FALSE). Whether or not to plot the segmentation result.
#' @param hmm_states An ordered vector for the HMM numeric states (deletion, 1-copy gain, 2-copy gains).
#' @param hmm_sd Numeric. Fixed standard deviation for the HMM states.
#' @param hmm_p Numeric. Transition probability for the HMM algorithm.
#' @param nmean Integer. Width of moving window for runmean.
#' @param adj Numeric. Value for tumor coverage adjustment.
#' @param rds_path The path for saving the rds files for the estimated results for each region.
#' @param max_qt Numeric value in [0,1]. Setting the maximum value to the max_qt quantile to avoid extreme values.
#'
#' @return A Alleloscope object with "seg_table" added.
#'
#' @import caTools
#' @import HiddenMarkov
#' @export
Segmentation_bulk=function(Obj_filtered=NULL, hmm_states=c(0.5, 1.5, 1.8), hmm_sd=0.2, hmm_p=0.000001,nmean=100, plot_seg=TRUE,rds_path=NULL, adj=0, max_qt=0.99){
  ## test
  # check parameters
  if(is.null(Obj_filtered)){
    stop("Please provide a valid Alleloscope object for Obj_filtered.")
  }else if(length(hmm_states)!=3){
    stop("hmm_states should be an ordered vector for the three HMM numeric states.")
  }else if(!is.numeric(hmm_sd)){
    stop("hmm_sd should be a numeric value.")
  }else if(!is.numeric(hmm_p)){
    stop("hmn_p should be an ordered vector for the three HMM numeric states.")
  }


  # set values
  assay=Obj_filtered$assay
  dir_path=Obj_filtered$dir_path
  samplename=Obj_filtered$samplename
  genome_assembly=Obj_filtered$genome_assembly
  raw_counts=Obj_filtered$raw_counts
  ref_counts=Obj_filtered$ref_counts

  size_all=Obj_filtered$size
  size_all=size_all[order(as.numeric(as.character(names(size_all))))]


  if(is.null(rds_path)){
    rds_path=paste0(dir_path,"/rds/")
  }
  dir.create(rds_path)

  dir.create(paste0(dir_path,"/plots/"))

  ## check if tumor and normal are from the same genome assembly


  ## selection chr name from the names of "size"
  if(grepl('chr',names(Obj_filtered$size)[1])){
    chr_name=as.character(names(Obj_filtered$size))
  }else{
    chr_name=paste0('chr',names(Obj_filtered$size)) # size does not have 'chr'
  }


  ## raw/ref count matrix info
  raw_chr=sapply(strsplit(rownames(raw_counts),'-'),'[',1)
  raw_start=as.numeric(sapply(strsplit(rownames(raw_counts),'-'),'[',2))
  raw_end=as.numeric(sapply(strsplit(rownames(raw_counts),'-'),'[',3))

  ref_chr=sapply(strsplit(rownames(ref_counts),'-'),'[',1)

  ## check if the tumor and normal are comparable
  if(table(raw_chr)[chr_name[length(chr_name)]] != table(ref_chr)[chr_name[length(chr_name)]]){
    stop("Tumor and control coverage files are with different numbers of bins.")
  }


  ## pool all cells
  chr_sum=c()
  ref_sum=c()
  chromnum=c()
  raw_counts[,1]=as.numeric(raw_counts[,1])
  ref_counts[,1]=as.numeric(ref_counts[,1])

  for(chrr in chr_name){
    chr_mm=raw_counts[which(raw_chr %in% chrr),, drop=F]
    ref_mm=ref_counts[which(ref_chr %in% chrr),, drop=F]
    chr_sub=apply(chr_mm,1, sum)
    ref_sub=apply(ref_mm,1, sum)
    chr_sum=c(chr_sum, chr_sub)
    ref_sum=c(ref_sum, ref_sub)
  }
  chromnum=raw_chr


  cov2=chr_sum/ref_sum
  cov2[is.na(cov2)]=0
  cov2[cov2==Inf]=100

  # if(grepl('chr',chromnum[1])){
  #   chromnum=as.numeric(sapply(strsplit(chromnum,'hr'),'[',2)) ## if chrommnum with chr
  # }else{
  #   chromnum=as.numeric(chromnum)}
  # chrline=cumsum(table(chromnum))
  # maploc=1:length(cov2)


  cov3=cov2
  cov3=pmin(cov3, quantile(cov3, max_qt))
  cov3=cov3
  cov4=cov3/median(cov3)

  # remove low cov bins
  sel_bin=which(cov4>0.4)
  cov4=cov4[sel_bin]
  cov2=cov2[sel_bin]
  cov3=cov3[sel_bin]
  chromnum=chromnum[sel_bin]


  raw_chr=raw_chr[sel_bin]
  raw_start=raw_start[sel_bin]
  raw_end=raw_end[sel_bin]
  ref_chr=ref_chr[sel_bin]


  if(grepl('chr',chromnum[1])){
    chromnum=as.numeric(sapply(strsplit(chromnum,'hr'),'[',2)) ## if chrommnum with chr
  }else{
    chromnum=as.numeric(chromnum)}
  chrline=cumsum(table(chromnum))
  maploc=1:length(cov4)

  ## HMM
  nsnp=c()
  seg_table_all=NULL
  for(ii in sapply(strsplit(chr_name,'hr'),'[',2)){
    cov4=(cov3/median(cov3))[which(chromnum==paste0(ii))]+adj
    cov5=caTools::runmean(cov4, nmean)

    ppa1= hmm_states[3]
    ppa2= hmm_states[2]
    ppn=1
    ppd= hmm_states[1]
    delta <- c(0.1,0.2,0.5,0.2)
    t=hmm_p
    z  <- HiddenMarkov::dthmm(cov5, matrix(c(1-3*t, t, t,t,t, 1-3*t,t,t, t,t,1-3*t,t,t,t,t,1-3*t), byrow=TRUE, nrow=4), delta, "norm", list(mean=c(ppa1, ppa2,ppn, ppd),sd=c(hmm_sd, hmm_sd, hmm_sd, hmm_sd)))


    results <- HiddenMarkov::Viterbi(z)

    #plot(cov2, col=c('red', 'blue')[results], pch=20)
    states=results
    states[results==1]=ppa1
    states[results==2]=ppa2
    states[results==3]=ppn
    states[results==4]=ppd



    ##get seg positions
    var_list_sub=matrix(nrow=length(cov5), ncol=3)
    var_list_sub[,1]=ii

    start=raw_start[which(chromnum==ii)]-1
    end= raw_end[which(chromnum==ii)]
    size=size_all[ii]

    var_list_sub[,2:3]=cbind(start, end)
    colnames(var_list_sub)=c('chr', 'start','end')
    var_list_sub=data.frame(var_list_sub, stringsAsFactors = F)

    changepoints=which(sapply(2:length(states), function(x) (states[x]!=states[x-1])))+1
    changepoints=changepoints[which(changepoints!=2)]

    ff=0
    for(nn in c(changepoints, length(states))){
      nsnp=c(nsnp, nn-ff)
      ff=nn
    }

    seg_table=matrix(0,nrow=(length(changepoints)+1), ncol=7)
    seg_table[,1]=rep(ii, dim(seg_table)[1])
    seg_table[,2]=as.numeric(as.character(c(start[1],(as.numeric(var_list_sub[changepoints,2])))))   ##0 based
    seg_table[,3]=as.numeric(as.character(c((as.numeric(var_list_sub[changepoints,2])), size)))
    seg_table[,4]=c(states[2], states[changepoints])
    seg_table[,5]=as.numeric(seg_table[,3])-as.numeric(seg_table[,2])
    seg_table[,6]=sapply(c(1:(length(changepoints)+1)), function(x) mean(cov4[c(1, changepoints, length(cov4))[x]:(c(1, changepoints, length(cov4))[x+1]-1)]) )
    seg_table[,7]=sapply(c(1:(length(changepoints)+1)), function(x) var(cov4[c(1, changepoints, length(cov4))[x]:(c(1, changepoints, length(cov4))[x+1]-1)]) )
    seg_table_all=rbind(seg_table_all, seg_table)
    cat(paste0('chr',ii," "))

  }

  cat("\n")
  ## combine first points


  colnames(seg_table_all)=c("chr","start", "end", "states","length","mean","var")


  ## visualization
  if(plot_seg==TRUE){
    ylim=c(0,max((cov3/median(cov3))+adj))
    pdf(paste0(dir_path,"/plots/",samplename,"_seg_hmm.pdf"))
    par(mfrow=c(3,1))
    plot(x=maploc, y=(cov3/median(cov3))+adj, ylab="COV",
         main=paste("HMM segmentation across all chroms"),
         xaxt="n", pch=20, cex=0.3, xlab = "chromosome", ylim=ylim)
    points(x=maploc, y=rep(seg_table_all[,4],nsnp), type='l', col='red', lwd=1)
    if(length(unique(chromnum))==1){
      axis(side=1, at=table(chromnum)*0.5, labels = sapply(strsplit(chr_name,'hr'),'[',2), cex.axis=1)
    }else{
      axis(side=1, at=c(0,chrline[1:(length(chrline)-1)])+table(chromnum)*0.5, labels = sapply(strsplit(chr_name,'hr'),'[',2), cex.axis=1)}
    abline(v=chrline, col='blue', lty=1, lwd=1)

    for(ii in sapply(strsplit(chr_name,'hr'),'[',2)){
      plot(x=maploc[which(chromnum==ii)], y=((cov3/median(cov3))+adj)[which(chromnum==ii)], ylab="COV",main=paste("chr",ii), xaxt="n", pch=20, cex=0.3, xlab = "chromosome", ylim=ylim ) # limit to only the "high coverage" SNPs.
      points(x=maploc[which(chromnum==ii)], y=rep(seg_table_all[,4],nsnp)[which(chromnum==ii)], type='l', col='red', lwd=3)
      if(length(unique(chromnum))==1){
        axis(side=1, at=table(chromnum)*0.5, labels = sapply(strsplit(chr_name,'hr'),'[',2), cex.axis=1)
      }else{
        axis(side=1, at=c(0,chrline[1:(length(chrline)-1)])+table(chromnum)*0.5, labels = sapply(strsplit(chr_name,'hr'),'[',2), cex.axis=1)
      }
      abline(v=chrline, col='blue', lty=1, lwd=3)
    }
    dev.off()
  }

  seg_table_all=data.frame(seg_table_all,stringsAsFactors = F)
  seg_table_all$chrr=as.character(paste0(seg_table_all$chr,":", seg_table_all$start))

  Obj_filtered[['seg_table']]=seg_table_all
  saveRDS(seg_table_all, paste0(rds_path, "/seg_table_all_",samplename,".rds"))

  message("Segmentation done!")
  cat("\"seg_table\" was added to the Obj_filtered object.\n")
  cat(paste0("Segmentation plot was saved in the path:", dir_path,"plots/\n"))
  cat(paste0("Segmentation result is stored as seg_table_all.rds in the path:", rds_path,"\n"))


  return(Obj_filtered)
}
