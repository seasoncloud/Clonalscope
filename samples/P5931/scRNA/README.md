Clonalscope (scRNA-seq)
================
Chi-Yun Wu, Zhang Lab, University of Pennsylvania

## Description
For scRNA-seq data, Clonalscope enables subclone detection and malignant cell labeling.

For more information about the method, please check out the [github](https://github.com/seasoncloud/Clonalscope) and the [manuscript]().
<br/>

## Prepare input files
The following are the input files for different steps.

1. A feature-barcode (sparse) UMI matrix. [EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/alt_all_sub.mtx) 
* The matrix is the standard output of the Cell Ranger software.
* Each element is the number of UMIs associated with a feature (row) and a barcode (column).
<br/>
 
2. A tsv file for all cell barcodes. [EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/barcodes_sub.tsv)
* The "barcodes.tsv" file is the standard output of the Cell Ranger software.
* Each row is a barcode indicating cell identity.
<br/>
 
3. A tsv file for feature information. 
[EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/barcodes_sub.tsv)
* The "features.tsv" file is the standard output of the Cell Ranger software.
* Each row is a gene.The first two columns need to be gene ID and gene name.
<br/>

4. A matrix for cell type annotation. 
[EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/barcodes_sub.tsv)
* Each row is a cell. The first and second columns are "cell barcodes" and "annotated cell types". 
* This matrix is used to specify the cells to use as the baseline for coverage normalization. Please specify the set of cells to be used as the baseline as "normal". 
<br/>

5. Read count matrices from matched tumor and normal WGS/WES data. (Optional) [EXAMPLE](https://github.com/seasoncloud/Alleloscope/blob/main/data-raw/SNU601/scDNA/tumor_sub.txt) 
* Each row is a bin. The values in the matrices represent total read counts for each bin.
* The matrix has four columns: col1: chr1 to 22; col2: start position; col3: end position; col4 summed read counts in the bins
* Without paired normal sample, other normal samples aligned to the same reference genome (eg. GRCh38) also work if with matched bins.
* Without matched tumor and normal WGS/WES data, chromosome arms can be used as the segments.
<br/>

## Tutorial for scRNA-seq data
* Here is an example application to the P5931 scRNA-seq dataset from Sathe et al., 2020.
<br/>

#### Step0. Load the input files

* In R, set up the environment and read common files
```
library(Clonalscope) # load the library
setwd("~/Clonalscope/") # set path to the github folder

dir_path <- "./samples/P5931/scRNA/output/"; dir.create(dir_path) # set up output directory

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
mtx=readMM("data-raw/P5931/scRNA/filtered_feature_bc_matrix/matrix.mtx.gz")
barcodes=read.table("data-raw/P5931/scRNA/filtered_feature_bc_matrix/barcodes.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
features=read.table("data-raw/P5931/scRNA/filtered_feature_bc_matrix/features.tsv.gz", stringsAsFactors = F, sep='\t', header=F)
celltype0=readRDS("data-raw/P5931/scRNA/celltype_all.rds")

# files for matched WGS/WES data
WGSt=readRDS(paste0("data-raw/P5931/WGS/raw_counts.rds"))
WGSn=readRDS(paste0("data-raw/P5931/WGS/ref_counts.rds"))
```
<br/>

#### Step1. Segment the genome using matched WGS/WES data (Optional)

* Chromosome arms can be used as the segments without mathced WGS/WES data. 

* First, create an object for segmentation.
```
Obj_filtered=Createobj_bulk(raw_counts=WGSt, # from matched DNA sequencing (bulk/single)
                            ref_counts=WGSn, # from matched DNA sequencing (bulk/single)
                            samplename= "P5931",
                            genome_assembly="GRCh38", dir_path=dir_path, size=size, assay='WGS')
```

* Segment the genome on matched WGS data.
```
Obj_filtered=Segmentation_bulk(Obj_filtered=Obj_filtered,
                               plot_seg = TRUE, hmm_states = c(0.8, 1.1, 1.2))
```
<br/><br/>

#### Step2. Subclone detection in scRNA-seq data

* Assign cells to be used as baseline for coverage normalization as "normal"
```
celltype=celltype0
celltype[which(!grepl('Epithelial', celltype[,2])),2]='normal'
```

* Select barcodes for clustering only on epithelial cells
```
clustering_barcodes=celltype[which(grepl("Epithelial",celltype[,2])),1]
```

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
                      celltype0=celltype, var_pt=0.99, var_pt_ctrl=0.99, include='all',
                      alpha_source='all', ctrl_region=NULL, 
                      seg_table_filtered=Obj_filtered$seg_table, size=size,
                      dir_path=dir_path, breaks=50, prep_mode = 'intersect', est_cap = 2,
                      clust_mode='cna_only',clustering_barcodes=clustering_barcodes) 
```
<br/>

#### Step3. Visualization

* Extract values from the object.
```
clustering= Cov_obj$result_final$clustering
clustering2= Cov_obj$result_final$clustering2
result=Cov_obj$result_final$result
table(result$Zest)
```

* Visualize the subclones/clusters in the heatmap.
```
PlotClusters(df = Cov_obj$result_final$df_obj$df, celltype = celltype0, Assign_obj =result, mode = "genome",  fontsize = 7, lab_mode='annot')
```

![](../../../inst/plots/P5931_heatmap.png?raw=true "P5931 heatmap")

* Visuzlize the subclones/clusters in the UMAP.
```
* Plot the UMAP and save the coordinates
set.seed(2022)
emb=PlotUMAP(df = Cov_obj$result_final$df_obj$df, celltype = celltype0, Assign_obj =result, mode = "Zest")
```

![](../../../inst/plots/P5931_UMAP.png?raw=true "P5931 UMAP")
<br/><br/>



## Citation
Wu, C.-Y. et al. Cancer subclone detection based on DNA copy number in single cell and spatial omic sequencing data. *bioRxiv* (2022): []()





[Back to the main page](https://github.com/seasoncloud/Clonalscope)
