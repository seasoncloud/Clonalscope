Clonalscope (Spatial Transcriptomic subclone detection)
================
Chi-Yun Wu, Zhang Lab, University of Pennsylvania

## Description
Clonalscope enables subclone detection, malignant spot labeling (with matched WGS/WES data), and subclone tracing in spatial transcriptomic data. Here we demonstrated Clonalscope's ability in subclone detection on a spatial transcriptomic dataset.

For more information about the method, please check out the [github](https://github.com/seasoncloud/Clonalscope) and the [manuscript](https://doi.org/10.1101/2022.07.05.498882).
<br/>

## Prepare input files
The following are the input files for different steps.

1. A feature-barcode (sparse) UMI matrix. [EXAMPLE for 1-3](https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.tar.gz) 
* The matrix is the standard output of the Space Ranger software.
* Each element is the number of UMIs associated with a feature (row) and a barcode (column).
<br/>
 
2. A tsv file for all cell barcodes. 
* The "barcodes.tsv" file is the standard output of the Space Ranger software.
* Each row is a barcode indicating cell identity.
<br/>
 
3. A tsv file for feature information. 
* The "features.tsv" file is the standard output of the Space Ranger software.
* Each row is a gene.The first two columns need to be gene ID and gene name.
<br/>

4. A matrix for cell type annotation. 
[EXAMPLE](https://github.com/seasoncloud/Clonalscope/blob/main/data-raw/P5931/scRNA/celltype_all.rds)
* Each row is a cell. The first and second columns are "cell barcodes" and "annotated cell types". 
* This matrix is used to specify the cells to use as the baseline for coverage normalization. Please specify the set of cells to be used as the baseline as "normal". 
<br/>

5. Read count matrices from matched tumor and normal WGS/WES data. (Optional) [EXAMPLE](https://github.com/seasoncloud/Clonalscope/tree/main/data-raw/P5931/WGS) 
* Each row is a bin. The values in the matrices represent total read counts for each bin.
* The matrix has one column showing the summed read counts for each bin.
* The row names should be ordered and in the format of "chr-start-end" such as chr1-1-20000.
* Without paired normal sample, other normal samples aligned to the same reference genome (eg. GRCh38) also work if with matched bins.
* Without matched tumor and normal WGS/WES data, chromosome arms can be used as the segments.
<br/>

## Tutorial for Spatial transcriptomic data
* Here is an example application of subclone detection and lineage tracing to two datasets from two sections of a human breast tumor (https://www.10xgenomics.com/resources/datasets/human-breast-cancer-block-a-section-2-1-standard-1-1-0).
<br/>

#### Step0. Load the input files

* In R, set up the environment and read common files
```
library(Clonalscope) # load the library
setwd("~/Clonalscope/") # set path to the github folder

dir_path <- "./samples/BC_ductal2/ST/output/"; dir.create(dir_path) # set up output directory

# Size of each chromosome (hg19 and GRCh38 are provided.)
size=read.table("data-raw/sizes.cellranger-GRCh38-1.0.0.txt", stringsAsFactors = F)

# List of cyclegenes retrieved from the "CopyKAT"package (https://github.com/navinlabcode/copykat)
cyclegenes=readRDS("data-raw/cyclegenes.rds")

# bed file indicating gene positions (hg19 and GRCh38 are provided.)
bed=read.table("data-raw/hg38.genes.bed", sep='\t', header = T)

```

* Read example files
```
# files for scRNA-seq data
mtx=readMM("data-raw/BC_ductal2/ST/filtered_feature_bc_matrix/matrix.mtx.gz")
barcodes=read.table("data-raw/BC_ductal2/ST/filtered_feature_bc_matrix/barcodes.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
features=read.table("data-raw/BC_ductal2/ST/filtered_feature_bc_matrix/features.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
celltype0=readRDS("data-raw/BC_ductal2/ST/celltype_all.rds")
```
<br/>

#### Step1. Segment the genome using matched WGS/WES data (Optional)

* In this step, prior segments and CNA states are retrieved from matched WGS/WES data for the downstream analysis. 

* Chromosome arms are used as the segments here without mathced WGS/WES data.
```
# Generate segmentation table for each chromosome arm.
chrarm=read.table("data-raw/cytoarm_table_hg38.txt", stringsAsFactors = F, sep='\t', header=T)
chrarm=chrarm[order(as.numeric(chrarm[,1]),as.numeric(chrarm[,3])),]
bin_bed=chrarm[,-2]

seg_table_filtered=data.frame("chr"=bin_bed[,1], 'start'=as.numeric(bin_bed[,2]),
'end'=as.numeric(bin_bed[,3]), 'states'=1, 'length'=as.numeric(bin_bed[,3])-as.numeric(bin_bed[,2]),
'mean'=0, 'var'=0, 'Var1'=1:nrow(bin_bed),'Freq'=50000,
'chrr'=paste0(bin_bed[,1],":", bin_bed[,2]), stringsAsFactors = F)
```

* For segmentation with matched WGS/WES data, please see example codes in  [Step.1](https://github.com/seasoncloud/Clonalscope/tree/main/samples/P5931/scRNA) 

<br/><br/>

#### Step2. Subclone detection in spatial transcriptomic data

* In Step2, copy number states of each region are estimated for each spot in the spatial transcriptomic dataset with the prior information from bulk WGS/WES data. Then subclones based on copy number configurations are detected by a nested Chinese Restaurant Process.

* Filter input files
```
Input_filtered=FilterFeatures(mtx=mtx, barcodes=barcodes, features=features, cyclegenes=cyclegenes)

# Remove raw inputs
rm(mtx); rm(barcodes); rm(features)
```

* Subclone detection in scRNA-seq data
```
set.seed(2022)
Cov_obj=RunCovCluster(mtx=Input_filtered$mtx, barcodes=Input_filtered$barcodes, 
                      features=Input_filtered$features, bed=bed, 
                      celltype0=celltype0, var_pt=0.99, var_pt_ctrl=0.99, include='all',
                      alpha_source='all', ctrl_region=NULL, 
                      seg_table_filtered=seg_table_filtered, size=size,
                      dir_path=dir_path, breaks=50, prep_mode = 'intersect', seed=200) 
# save the object
saveRDS(Cov_obj,paste0(dir_path,"/Cov_obj.rds"))
```
<br/>

#### Step3. Visualization

* Extract clustering results from the object.
```
clustering= Cov_obj$result_final$clustering
clustering2= Cov_obj$result_final$clustering2
result=Cov_obj$result_final$result
Zest=result$Zest
table(result$Zest)
```

* If you want to adjust the resolution (ie. minimum number of cells/spots in each cluster), the following command can be re-run with different "mincell" values. Larger values will result in less clusters.
```
#result=AssignCluster(clustering2, mincell = 100)
#Zest=result$Zest
#table(result$Zest)
```

* Visualize the subclones/clusters in the heatmap.
```
PlotClusters(df = clustering$data, celltype = celltype0, Assign_obj =result, mode = "genome",  fontsize = 7, lab_mode=FALSE)
```

![](../../../inst/plots/BC_ductal2_heatmap.png?raw=true "BC section2 heatmap")

* Visuzlize the subclones/clusters in the UMAP.
```
# Plot the UMAP and save the coordinates
set.seed(2022)
emb=PlotUMAP(df = clustering$data, celltype = celltype0, Assign_obj =result, mode = "Zest")
```

![](../../../inst/plots/BC_ductal2_UMAP.png?raw=true "BC section2 UMAP")
<br/><br/>

#### Step4. Integrate with Spatial information (Optional)

* Load the h5 file (Standard Space Ranger output) as a Seurat object.
* The h5 file for the human breast tumor dataset can be found [here](https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Breast_Cancer_Block_A_Section_2/V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5)
```
library(Seurat)
vs <- Load10X_Spatial("data-raw/BC_ductal2/ST/", 
filename = "V1_Breast_Cancer_Block_A_Section_2_filtered_feature_bc_matrix.h5")
```

* Set the color the same as those in Clonalscope.
```
library(gplots)
cols=col2hex(colors()[c(1,609, 536, 62, 652, 611, 463, 498, 71,
258, 84, 56, 26, 154, 59, 134, 78, 116, 85, 20, 259)])
```

* Visualize subclones in the H&E image.
```
clusters=Zest[match(names(vs$orig.ident), names(Zest))]; clusters[is.na(clusters)]=0
Idents(vs)=factor(clusters, levels=sort(as.numeric(unique(clusters))))
SpatialDimPlot(vs, label = F, label.size = 3, stroke = 0, pt.size.factor = 2, cols=cols)
```

![](../../../inst/plots/BC_ductal2_H&E.png?raw=true "BC section2 H&E")
<br/><br/>


## Citation
Wu, C.-Y. et al. Cancer subclone detection based on DNA copy number in single cell and spatial omic sequencing data. *bioRxiv* (2022): [https://doi.org/10.1101/2022.07.05.498882](https://doi.org/10.1101/2022.07.05.498882)





[Back to the main page](https://github.com/seasoncloud/Clonalscope)
