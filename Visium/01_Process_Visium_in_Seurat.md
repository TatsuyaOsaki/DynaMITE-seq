## February 2025
## Round 1 REVISIONS for SCLC 10X Visium
## Combine Round 1 analysis with Sample A of Round 2


```R
library(Seurat)
library(SeuratExtend)
#library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(limma)
library(openxlsx)

#library(STdeconvolve)

#knitr::opts_chunk$set(fig.width=14)
```

    Loading required package: SeuratObject
    
    Warning message:
    “package ‘SeuratObject’ was built under R version 4.3.3”
    Loading required package: sp
    
    
    Attaching package: ‘SeuratObject’
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, t
    
    
    Loading required package: SeuratExtendData
    
    
    Attaching package: ‘dplyr’
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    



```R
VA2 <- Load10X_Spatial("/immuno/ian/Projects/2024_09_Barbie_Visium_Round2/spaceranger_count/VisA_count/outs/")
VA1 <- Load10X_Spatial("/homes1/iddryg/ian/Projects/2023_02_SCLC_10X_Visium/VA_count/outs/")
VB1 <- Load10X_Spatial("/homes1/iddryg/ian/Projects/2023_02_SCLC_10X_Visium/VB_count/outs/")
VA2$orig.ident <- "VA2"
VA1$orig.ident <- "VA1"
VB1$orig.ident <- "VB1"
# set visium round, will use this for integration / batch correction
VA2$visium_round <- "round2"
VA1$visium_round <- "round1"
VB1$visium_round <- "round1"
```


```R
# Filter out spots that have zero UMIs, otherwise SCTransform won't run
VA2 <- subset(VA2, nCount_Spatial>0)
VA1 <- subset(VA1, nCount_Spatial>0)
VB1 <- subset(VB1, nCount_Spatial>0)
```


```R
VA2
```


    An object of class Seurat 
    17943 features across 4316 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     1 layer present: counts
     1 spatial field of view present: slice1



```R
VA1
```


    An object of class Seurat 
    17943 features across 2806 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     1 layer present: counts
     1 spatial field of view present: slice1



```R
VB1
```


    An object of class Seurat 
    17943 features across 2589 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     1 layer present: counts
     1 spatial field of view present: slice1



```R
Assays(VA1)
```


'Spatial'



```R
Layers(VA1)
```


'counts'



```R
# convert to V5 assays...
VA1[["Spatial"]] <- as(VA1[["Spatial"]], "Assay5")
VB1[["Spatial"]] <- as(VB1[["Spatial"]], "Assay5")
VA2[["Spatial"]] <- as(VA2[["Spatial"]], "Assay5")
```


```R
VA2
```


    An object of class Seurat 
    17943 features across 4316 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     1 layer present: counts
     1 spatial field of view present: slice1



```R
VA1
```


    An object of class Seurat 
    17943 features across 2806 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     1 layer present: counts
     1 spatial field of view present: slice1



```R
VB1
```


    An object of class Seurat 
    17943 features across 2589 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     1 layer present: counts
     1 spatial field of view present: slice1



```R
options(future.globals.maxSize = 3e+09)
```
# perform SCTransform on the slices separately, then merge
VA2 <- SCTransform(VA2, assay = "Spatial", verbose = TRUE, variable.features.n = 10000)
VA2 <- RunPCA(VA2, npcs = 30, verbose = F)# perform SCTransform on the slices separately, then merge
VA1 <- SCTransform(VA1, assay = "Spatial", verbose = TRUE, variable.features.n = 10000)
VA1 <- RunPCA(VA1, npcs = 30, verbose = F)# perform SCTransform on the slices separately, then merge
VB1 <- SCTransform(VB1, assay = "Spatial", verbose = TRUE, variable.features.n = 10000)
VB1 <- RunPCA(VB1, npcs = 30, verbose = F)

```R
options(repr.plot.width=8, repr.plot.height=7)
```


```R
VlnPlot(VA2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend() & geom_hline(yintercept = 300)
```

    Warning message:
    “Default search for "data" layer in "Spatial" assay yielded no results; utilizing "counts" layer instead.”



    
![png](output_18_1.png)
    



```R
VlnPlot(VA1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend() & geom_hline(yintercept = 300)
```

    Warning message:
    “Default search for "data" layer in "Spatial" assay yielded no results; utilizing "counts" layer instead.”



    
![png](output_19_1.png)
    



```R
VlnPlot(VB1, features = "nCount_Spatial", pt.size = 0.1) + NoLegend() & geom_hline(yintercept = 300)
```

    Warning message:
    “Default search for "data" layer in "Spatial" assay yielded no results; utilizing "counts" layer instead.”



    
![png](output_20_1.png)
    



```R
VA2
```


    An object of class Seurat 
    17943 features across 4316 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     1 layer present: counts
     1 spatial field of view present: slice1



```R
VA1
```


    An object of class Seurat 
    17943 features across 2806 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     1 layer present: counts
     1 spatial field of view present: slice1



```R
VB1
```


    An object of class Seurat 
    17943 features across 2589 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     1 layer present: counts
     1 spatial field of view present: slice1


# Integrate datasets, which allows for correction from different batches (round1 and round2)


```R
# merge step here.......
# Combine datasets into a single object, then split into layers
AB1_A2 <- merge(VA2, y = c(VA1, VB1), add.cell.ids = c("VA2", "VA1", "VB1"))
```

    Warning message:
    “Key ‘slice1_’ taken, using ‘slice12_’ instead”
    Warning message:
    “Key ‘slice1_’ taken, using ‘slice13_’ instead”



```R
AB1_A2
Assays(AB1_A2)
Layers(AB1_A2)
```


    An object of class Seurat 
    17943 features across 9711 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     3 layers present: counts.1, counts.2, counts.3
     3 spatial fields of view present: slice1 slice1.2 slice1.3



'Spatial'



<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts.1'</li><li>'counts.2'</li><li>'counts.3'</li></ol>




```R
Layers(AB1_A2[['Spatial']])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts.1'</li><li>'counts.2'</li><li>'counts.3'</li></ol>




```R
# According to this: https://github.com/satijalab/seurat/issues/7316
# we should JoinLayers() and then re-split. 
AB1_A2[['Spatial']]  <- JoinLayers(object = AB1_A2[['Spatial']])
```


```R
AB1_A2[["Spatial"]] <- split(AB1_A2[["Spatial"]], f = AB1_A2$visium_round)
```


```R
AB1_A2
Assays(AB1_A2)
Layers(AB1_A2)
Layers(AB1_A2[['Spatial']])
```


    An object of class Seurat 
    17943 features across 9711 samples within 1 assay 
    Active assay: Spatial (17943 features, 0 variable features)
     2 layers present: counts.round2, counts.round1
     3 spatial fields of view present: slice1 slice1.2 slice1.3



'Spatial'



<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts.round2'</li><li>'counts.round1'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts.round2'</li><li>'counts.round1'</li></ol>




```R
AB1_A2 <- SCTransform(AB1_A2, assay = "Spatial", verbose = TRUE, variable.features.n = 10000)
AB1_A2 <- RunPCA(AB1_A2, npcs = 30, verbose = F)
```

    Running SCTransform on assay: Spatial
    
    Running SCTransform on layer: counts.round2
    
    vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.
    
    Variance stabilizing transformation of count matrix of size 17807 by 4316
    
    Model formula is y ~ log_umi
    
    Get Negative Binomial regression parameters per gene
    
    Using 2000 genes, 4316 cells
    
    Found 136 outliers - those will be ignored in fitting/regularization step
    
    
    Second step: Get residuals using fitted parameters for 17807 genes
    
    Computing corrected count matrix for 17807 genes
    
    Calculating gene attributes
    
    Wall clock passed: Time difference of 34.00261 secs
    
    Determine variable features
    
    Centering data matrix
    
    Running SCTransform on layer: counts.round1
    
    vst.flavor='v2' set. Using model with fixed slope and excluding poisson genes.
    
    Variance stabilizing transformation of count matrix of size 16989 by 5395
    
    Model formula is y ~ log_umi
    
    Get Negative Binomial regression parameters per gene
    
    Using 2000 genes, 5000 cells
    
    Found 246 outliers - those will be ignored in fitting/regularization step
    
    
    Second step: Get residuals using fitted parameters for 16989 genes
    
    Computing corrected count matrix for 16989 genes
    
    Calculating gene attributes
    
    Wall clock passed: Time difference of 32.12096 secs
    
    Determine variable features
    
    Centering data matrix
    
    Getting residuals for block 1(of 2) for round1 dataset
    
    Getting residuals for block 2(of 2) for round1 dataset
    
    Centering data matrix
    
    Finished calculating residuals for round1
    
    Centering data matrix
    
    Centering data matrix
    
    Set default assay to SCT
    



```R
AB1_A2
Assays(AB1_A2)
Layers(AB1_A2)
Layers(AB1_A2[['Spatial']])
```


    An object of class Seurat 
    35796 features across 9711 samples within 2 assays 
    Active assay: SCT (17853 features, 10000 variable features)
     3 layers present: counts, data, scale.data
     1 other assay present: Spatial
     1 dimensional reduction calculated: pca
     3 spatial fields of view present: slice1 slice1.2 slice1.3



<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Spatial'</li><li>'SCT'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts'</li><li>'data'</li><li>'scale.data'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts.round2'</li><li>'counts.round1'</li></ol>




```R
Layers(AB1_A2[['SCT']])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts'</li><li>'data'</li><li>'scale.data'</li></ol>




```R
# Integrate layers
# https://github.com/satijalab/seurat/issues/7672
AB1_A2 <- IntegrateLayers(
  object = AB1_A2,
  method = CCAIntegration,
  orig.reduction = "pca",
  normalization.method = "SCT",  # Now using SCT assay after splitting RNA
  new.reduction = "integrated.cca",
  verbose = TRUE
)
```

    Finding all pairwise anchors
    
    Running CCA
    
    Merging objects
    
    Finding neighborhoods
    
    Finding anchors
    
    	Found 9781 anchors
    
    Merging dataset 1 into 2
    
    Extracting anchors for merged samples
    
    Finding integration vectors
    
    Finding integration vector weights
    
    Integrating data
    



```R
AB1_A2
```


    An object of class Seurat 
    35796 features across 9711 samples within 2 assays 
    Active assay: SCT (17853 features, 10000 variable features)
     3 layers present: counts, data, scale.data
     1 other assay present: Spatial
     2 dimensional reductions calculated: pca, integrated.cca
     3 spatial fields of view present: slice1 slice1.2 slice1.3



```R
# I think we don't need to JoinLayers() if we used SCT normalization. There's no mention of it in the vignette
# AB1_A2 <- JoinLayers(AB1_A2)
```


```R
# Joint analysis workflow
#AB1_A2 <- RunPCA(AB1_A2, npcs = 30, reduction.name = "pca_integrated")
AB1_A2 <- FindNeighbors(AB1_A2, reduction = "integrated.cca", dims = 1:30)
AB1_A2 <- FindClusters(AB1_A2, resolution = 0.8)
AB1_A2 <- RunUMAP(AB1_A2, dims = 1:30, reduction = "integrated.cca")
```

    Computing nearest neighbor graph
    
    Computing SNN
    


    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 9711
    Number of edges: 361141
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.8110
    Number of communities: 12
    Elapsed time: 0 seconds


    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    15:57:16 UMAP embedding parameters a = 0.9922 b = 1.112
    
    15:57:16 Read 9711 rows and found 30 numeric columns
    
    15:57:16 Using Annoy for neighbor search, n_neighbors = 30
    
    15:57:16 Building Annoy index with metric = cosine, n_trees = 50
    
    0%   10   20   30   40   50   60   70   80   90   100%
    
    [----|----|----|----|----|----|----|----|----|----|
    
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    |
    
    15:57:17 Writing NN index file to temp file /tmp/RtmpAHrjea/file1ce2492f9934b9
    
    15:57:17 Searching Annoy index using 1 thread, search_k = 3000
    
    15:57:20 Annoy recall = 100%
    
    15:57:20 Commencing smooth kNN distance calibration using 1 thread
     with target n_neighbors = 30
    
    15:57:21 Initializing from normalized Laplacian + noise (using RSpectra)
    
    15:57:22 Commencing optimization for 500 epochs, with 454462 positive edges
    
    15:57:22 Using rng type: pcg
    
    15:57:33 Optimization finished
    



```R
AB1_A2
```


    An object of class Seurat 
    35796 features across 9711 samples within 2 assays 
    Active assay: SCT (17853 features, 10000 variable features)
     3 layers present: counts, data, scale.data
     1 other assay present: Spatial
     3 dimensional reductions calculated: pca, integrated.cca, umap
     3 spatial fields of view present: slice1 slice1.2 slice1.3



```R
# Visualization
options(repr.plot.width=16, repr.plot.height=8)
DimPlot(AB1_A2, group.by = c("orig.ident", "seurat_clusters"), reduction = "umap")
```


    
![png](output_39_0.png)
    



```R
# Visualization
options(repr.plot.width=16, repr.plot.height=8)
DimPlot(AB1_A2, group.by = c("orig.ident", "seurat_clusters"), reduction = "pca")
```


    
![png](output_40_0.png)
    



```R
# Visualization
options(repr.plot.width=16, repr.plot.height=8)
DimPlot(AB1_A2, group.by = c("orig.ident", "seurat_clusters"), reduction = "integrated.cca")
```


    
![png](output_41_0.png)
    



```R
options(repr.plot.width=28, repr.plot.height=10)
SpatialDimPlot(AB1_A2, crop = TRUE)
```


    
![png](output_42_0.png)
    



```R
AB1_A2@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 8</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA2_AAACAAGTATCTCCCA-1</th><td>VA2</td><td>45867</td><td>8515</td><td>round2</td><td>24897</td><td>7041</td><td>5</td><td>5</td></tr>
	<tr><th scope=row>VA2_AAACACCAATAACTGC-1</th><td>VA2</td><td>11849</td><td>5082</td><td>round2</td><td>22823</td><td>5881</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>VA2_AAACAGAGCGACTCCT-1</th><td>VA2</td><td>43834</td><td>8662</td><td>round2</td><td>24885</td><td>7403</td><td>5</td><td>5</td></tr>
</tbody>
</table>




```R
AB1_A2
```


    An object of class Seurat 
    35796 features across 9711 samples within 2 assays 
    Active assay: SCT (17853 features, 10000 variable features)
     3 layers present: counts, data, scale.data
     1 other assay present: Spatial
     3 dimensional reductions calculated: pca, integrated.cca, umap
     3 spatial fields of view present: slice1 slice1.2 slice1.3



```R
Layers(AB1_A2)
Layers(AB1_A2[["Spatial"]])
Layers(AB1_A2[["SCT"]])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts'</li><li>'data'</li><li>'scale.data'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts.round2'</li><li>'counts.round1'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts'</li><li>'data'</li><li>'scale.data'</li></ol>




```R
# let's add normalized data to the Spatial slot... and run the rest of the standard pipeline separate from SCT
DefaultAssay(AB1_A2) <- "Spatial"
AB1_A2 <- NormalizeData(AB1_A2)
AB1_A2 <- FindVariableFeatures(AB1_A2)
AB1_A2 <- ScaleData(AB1_A2)
AB1_A2 <- RunPCA(AB1_A2)
```

    Normalizing layer: counts.round2
    
    Normalizing layer: counts.round1
    
    Finding variable features for layer counts.round2
    
    Finding variable features for layer counts.round1
    
    Centering and scaling data matrix
    
    PC_ 1 
    Positive:  PGAP3, IFI27, LMO3, ERV3-1, KRT15, SRCIN1, PIK3C2G, RERGL, ZNF117, CP 
    	   PCP4, MED1, MUC4, LENG8, SLITRK6, FTL, PPP1R1B, HIST1H3H, SLFN11, PARP14 
    	   CISD3, FABP5, CKB, IGHA1, LINC00672, TC2N, HIST2H2AB, CACNB1, PVALB, CRABP2 
    Negative:  GRTP1, AKR1C2, SCG3, AZGP1, SCG2, TFDP1, ADPRHL1, GPX2, AKR1C3, GAS6 
    	   CACNA1A, AOC1, AKR1C1, MCF2L, UCHL1, NELL1, MS4A8, CPLX2, SV2B, PLEKHB1 
    	   NKX2-2, SCGN, NEURL1, TNNT1, TMEM255B, IRS2, SGSM1, SUCNR1, PEX5L, RUNDC3A 
    PC_ 2 
    Positive:  DCN, IGKC, GADD45B, SLC34A2, PTGDS, FOSB, COL3A1, COL1A2, CCL18, LYZ 
    	   APOE, CCN1, C3, SPARC, AQP1, IGHG1, THBS1, CAVIN1, SFTPB, CTSD 
    	   IGFBP7, HLA-DMA, C1R, SERPINE1, FBN1, BGN, ACTA2, FN1, IGHJ6, DUSP1 
    Negative:  HIST1H2BH, HIST1H2AG, HIST1H1B, HIST1H4E, MDK, HIST1H3H, EIF4A2, HDAC2, TOP2A, HIST1H2AI 
    	   ALCAM, HIST2H2AB, HIST1H2BG, CKB, HIST1H2BD, NSD2, UBE2C, HIST1H4C, HIST1H3I, SMAD9 
    	   FAM111B, SP3, HIST1H2BL, RBP1, ERV3-1, HIST1H4A, TYMS, ZNF117, PGAP3, PIK3C2G 
    PC_ 3 
    Positive:  NPTX1, HOXB8, PCSK5, ADCY2, HIST1H2BB, TFF3, KLK12, PLLP, TMEM255B, PEG10 
    	   GRP, ERICH5, BDNF, FER1L6, AMPH, RELN, FRMD3, DYNC1I1, GRIN3A, SUCNR1 
    	   MMP9, NOP14, DLGAP1, EYA1, SMC2, CACNA2D1, SMS, PENK, TSPAN8, HHIP 
    Negative:  PDIA2, PHGR1, LIG4, KRT40, F10, PRAME, MLIP, SEZ6L, TAAR1, GRIA2 
    	   GTSF1, MCF2L, MAGEA10, ARX, ZNF208, AL032819.3, LSAMP, ATF5, ZNF676, KIF5C 
    	   MYL9, MAGEA1, POU3F1, GPR158, SERPINI1, AEBP1, FEV, TGFBI, ANKRD30B, MS4A8 
    PC_ 4 
    Positive:  PHGR1, PDIA2, F10, MLIP, PRAME, GRIA2, KRT40, SEZ6L, TAAR1, LAMA1 
    	   SLC7A14, PCSK2, GTSF1, KCNH7, LIG4, HOXD8, HS3ST4, KIF5C, F7, GAS2 
    	   MCF2L, ARX, HGD, CLNK, KCNJ3, AL032819.3, DCX, PTPRN, GPR158, MAGEA10 
    Negative:  S100A6, SCGB1A1, BPIFB1, SCGB3A1, VIM, CAPS, C20orf85, SLPI, MYL9, LBH 
    	   LCN2, SFTPB, DNAAF1, PIGR, TPPP3, HIST1H1B, HLA-DRA, SAA2, CCDC78, TSPAN1 
    	   FTL, HMGB1, CLU, WARS, A2M, FBLN1, SET, IGKC, TMEM190, HIST1H2BD 
    PC_ 5 
    Positive:  C20orf85, SCGB1A1, BPIFB1, SCGB3A1, DNAAF1, TPPP3, TMEM190, CFAP157, CAPS, LCN2 
    	   C11orf88, C9orf24, SAA2, FAM216B, RSPH1, MAPK15, CCDC78, SLPI, CDHR3, TSPAN1 
    	   FAM92B, TUBA4B, CFAP43, CCDC187, ZBBX, DLEC1, ROPN1L, EFCAB1, LRRC46, DNAI1 
    Negative:  COL1A1, SULF1, LUM, TIMP1, COL3A1, POSTN, FN1, COL5A1, BGN, AEBP1 
    	   SFRP4, COL1A2, SPP1, COL6A3, COL11A1, COL10A1, THBS2, IGKC, GREM1, MMP2 
    	   IGHG1, COL5A2, IGKV4-1, VCAN, IGLC1, ISLR, SPARC, INHBA, FTL, CDH11 
    
    Warning message:
    “Number of dimensions changing from 30 to 50”



```R
AB1_A2
```


    An object of class Seurat 
    35796 features across 9711 samples within 2 assays 
    Active assay: Spatial (17943 features, 2000 variable features)
     5 layers present: counts.round2, counts.round1, data.round2, data.round1, scale.data
     1 other assay present: SCT
     3 dimensional reductions calculated: pca, integrated.cca, umap
     3 spatial fields of view present: slice1 slice1.2 slice1.3



```R
AB1_A2@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 8</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA2_AAACAAGTATCTCCCA-1</th><td>VA2</td><td>45867</td><td>8515</td><td>round2</td><td>24897</td><td>7041</td><td>5</td><td>5</td></tr>
	<tr><th scope=row>VA2_AAACACCAATAACTGC-1</th><td>VA2</td><td>11849</td><td>5082</td><td>round2</td><td>22823</td><td>5881</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>VA2_AAACAGAGCGACTCCT-1</th><td>VA2</td><td>43834</td><td>8662</td><td>round2</td><td>24885</td><td>7403</td><td>5</td><td>5</td></tr>
</tbody>
</table>




```R
Layers(AB1_A2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'counts.round2'</li><li>'counts.round1'</li><li>'data.round2'</li><li>'data.round1'</li><li>'scale.data'</li></ol>




```R
# rejoin the layers
AB1_A2 <- JoinLayers(AB1_A2)
Layers(AB1_A2)
AB1_A2
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'data'</li><li>'counts'</li><li>'scale.data'</li></ol>




    An object of class Seurat 
    35796 features across 9711 samples within 2 assays 
    Active assay: Spatial (17943 features, 2000 variable features)
     3 layers present: data, counts, scale.data
     1 other assay present: SCT
     3 dimensional reductions calculated: pca, integrated.cca, umap
     3 spatial fields of view present: slice1 slice1.2 slice1.3



```R
# Save integrated dataset before continuing
saveRDS(AB1_A2, 'Processing/AB1_2_integrated_20250207.rds')
```


```R
options(repr.plot.width=22, repr.plot.height=17)
# set ident to orig.ident
Idents(AB1_A2) <- AB1_A2$orig.ident
plot1 <- VlnPlot(AB1_A2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(AB1_A2, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2, ncol=1)
```


    
![png](output_52_0.png)
    



```R
options(repr.plot.width=22, repr.plot.height=17)
# set ident to orig.ident
Idents(AB1_A2) <- AB1_A2$orig.ident
plot1 <- VlnPlot(AB1_A2, features = "nFeature_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(AB1_A2, features = "nFeature_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2, ncol=1)
```


    
![png](output_53_0.png)
    



```R
png(filename="Output/01_Process_Visium_and_Integrate/Violin_nCount_perSlide.png", width=600, height=400)
print(VlnPlot(AB1_A2, features = "nCount_Spatial", pt.size = 0.1) + NoLegend())
dev.off()
```


<strong>png:</strong> 2



```R
png(filename="Output/01_Process_Visium_and_Integrate/SpatialFeaturePlot_perSlide.png", width=1900, height=700)
print(SpatialFeaturePlot(AB1_A2, features = "nCount_Spatial") + theme(legend.position = "right"))
dev.off()
```


<strong>png:</strong> 2



```R
options(repr.plot.width=19, repr.plot.height=9)
# change ident to seurat cluster
Idents(AB1_A2) <- AB1_A2$seurat_clusters

plot <- DimPlot(AB1_A2, reduction = "umap", group.by = c("seurat_clusters", "orig.ident"),pt.size=1.3)

png(filename="Output/01_Process_Visium_and_Integrate/DimPlot.png", width=1900, height=900)
print(plot)
dev.off()

plot
```


<strong>png:</strong> 2



    
![png](output_56_1.png)
    



```R
options(repr.plot.width=26, repr.plot.height=9)

# reduce background image opacity with image.alpha
plot <- SpatialDimPlot(AB1_A2, image.alpha=0.5)
plot

png(filename="Output/01_Process_Visium_and_Integrate/SpatialDimPlot.png", width=2600, height=900)
print(plot)
dev.off()
```


<strong>png:</strong> 2



    
![png](output_57_1.png)
    



```R
AB1_A2@meta.data %>% head(5)
```


<table class="dataframe">
<caption>A data.frame: 5 × 8</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA2_AAACAAGTATCTCCCA-1</th><td>VA2</td><td>45867</td><td>8515</td><td>round2</td><td>24897</td><td>7041</td><td>5</td><td>5</td></tr>
	<tr><th scope=row>VA2_AAACACCAATAACTGC-1</th><td>VA2</td><td>11849</td><td>5082</td><td>round2</td><td>22823</td><td>5881</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>VA2_AAACAGAGCGACTCCT-1</th><td>VA2</td><td>43834</td><td>8662</td><td>round2</td><td>24885</td><td>7403</td><td>5</td><td>5</td></tr>
	<tr><th scope=row>VA2_AAACAGCTTTCAGAAG-1</th><td>VA2</td><td>46896</td><td>8897</td><td>round2</td><td>24782</td><td>7399</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>VA2_AAACAGGGTCTATATT-1</th><td>VA2</td><td>18610</td><td>6201</td><td>round2</td><td>22740</td><td>6158</td><td>4</td><td>4</td></tr>
</tbody>
</table>




```R
options(repr.plot.width=19, repr.plot.height=9)
plot <- VlnPlot(AB1_A2, features = "nCount_Spatial", group.by=c('seurat_clusters'), pt.size = 0.1) + NoLegend()

png(filename="Output/01_Process_Visium_and_Integrate/Violin_nCount_Spatial.png", width=1900, height=900)
print(plot)
dev.off()

plot
```


<strong>png:</strong> 2



    
![png](output_59_1.png)
    



```R
options(repr.plot.width=30, repr.plot.height=9)
plot <- VlnPlot(AB1_A2, features = "nCount_Spatial", split.by="orig.ident", pt.size = 0.1)

png(filename="Output/01_Process_Visium_and_Integrate/Violin_nCount_Spatial_groupbyslide.png", width=3000, height=900)
print(plot)
dev.off()

plot
```

    The default behaviour of split.by has changed.
    Separate violin plots are now plotted side-by-side.
    To restore the old behaviour of a single split violin,
    set split.plot = TRUE.
          
    This message will be shown once per session.
    



<strong>png:</strong> 2



    
![png](output_60_2.png)
    



```R
options(repr.plot.width=30, repr.plot.height=9)
plot <- VlnPlot(AB1_A2, features = "nCount_SCT", split.by="orig.ident", pt.size = 0.1)

png(filename="Output/01_Process_Visium_and_Integrate/Violin_nCount_SCT_groupbyslide.png", width=3000, height=900)
print(plot)
dev.off()

plot
```


<strong>png:</strong> 2



    
![png](output_61_1.png)
    


# Let's calculate and plot CPM


```R
AB1_A2
```


    An object of class Seurat 
    35796 features across 9711 samples within 2 assays 
    Active assay: Spatial (17943 features, 2000 variable features)
     3 layers present: data, counts, scale.data
     1 other assay present: SCT
     3 dimensional reductions calculated: pca, integrated.cca, umap
     3 spatial fields of view present: slice1 slice1.2 slice1.3



```R
Assays(AB1_A2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Spatial'</li><li>'SCT'</li></ol>




```R
AB1_A2[['Spatial']]
```


    Assay (v5) data with 17943 features for 9711 cells
    Top 10 variable features:
     SCGB1A1, IGLV3-1, IGHG3, IGKV4-1, IGLC1, IGHG4, IGHG2, IGHM, IGHA1,
    IGHG1 
    Layers:
     data, counts, scale.data 



```R
AB1_A2[['SCT']]
```


    SCTAssay data with 17853 features for 9711 cells, and 2 SCTModel(s) 
    Top 10 variable features:
     IGHG1, IGKC, SCGB1A1, IGHG3, IGKV4-1, IGLC1, IGHG2, IGLV3-1, IGHG4,
    IGHA1 



```R
AB1_A2[['SCT']]@counts[1:6,1:6]
AB1_A2[['SCT']]@data[1:6,1:6]
AB1_A2[['SCT']]@scale.data[1:6,1:6]
```


    6 x 6 sparse Matrix of class "dgCMatrix"
            VA2_AAACAAGTATCTCCCA-1 VA2_AAACACCAATAACTGC-1 VA2_AAACAGAGCGACTCCT-1
    SAMD11                       .                      .                      .
    NOC2L                        2                      2                      4
    KLHL17                       .                      .                      .
    PLEKHN1                      .                      .                      .
    PERM1                        .                      .                      .
    HES4                         3                      2                      .
            VA2_AAACAGCTTTCAGAAG-1 VA2_AAACAGGGTCTATATT-1 VA2_AAACAGTGTTCCTGGG-1
    SAMD11                       .                      .                      1
    NOC2L                        2                      4                      1
    KLHL17                       .                      .                      .
    PLEKHN1                      .                      .                      .
    PERM1                        .                      .                      .
    HES4                         1                      1                      1



    6 x 6 sparse Matrix of class "dgCMatrix"
            VA2_AAACAAGTATCTCCCA-1 VA2_AAACACCAATAACTGC-1 VA2_AAACAGAGCGACTCCT-1
    SAMD11                .                      .                      .       
    NOC2L                 1.098612               1.098612               1.609438
    KLHL17                .                      .                      .       
    PLEKHN1               .                      .                      .       
    PERM1                 .                      .                      .       
    HES4                  1.386294               1.098612               .       
            VA2_AAACAGCTTTCAGAAG-1 VA2_AAACAGGGTCTATATT-1 VA2_AAACAGTGTTCCTGGG-1
    SAMD11               .                      .                      0.6931472
    NOC2L                1.0986123              1.6094379              0.6931472
    KLHL17               .                      .                      .        
    PLEKHN1              .                      .                      .        
    PERM1                .                      .                      .        
    HES4                 0.6931472              0.6931472              0.6931472



<table class="dataframe">
<caption>A matrix: 6 × 6 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>VA2_AAACAAGTATCTCCCA-1</th><th scope=col>VA2_AAACACCAATAACTGC-1</th><th scope=col>VA2_AAACAGAGCGACTCCT-1</th><th scope=col>VA2_AAACAGCTTTCAGAAG-1</th><th scope=col>VA2_AAACAGGGTCTATATT-1</th><th scope=col>VA2_AAACAGTGTTCCTGGG-1</th></tr>
</thead>
<tbody>
	<tr><th scope=row>NOC2L</th><td> 0.05645474</td><td> 0.2926448</td><td> 1.7880043</td><td>0.02099814</td><td> 1.6236676</td><td>-0.46664479</td></tr>
	<tr><th scope=row>HES4</th><td> 1.71703982</td><td> 0.8248323</td><td>-0.4751658</td><td>0.19729685</td><td> 0.3509848</td><td>-0.40787556</td></tr>
	<tr><th scope=row>ISG15</th><td>-1.29434794</td><td>-1.2019006</td><td>-1.6410559</td><td>5.25841045</td><td> 0.5305226</td><td>-0.73651935</td></tr>
	<tr><th scope=row>AGRN</th><td> 0.34945732</td><td> 5.7534880</td><td>-0.3886751</td><td>1.86068166</td><td> 0.4418304</td><td> 2.11842016</td></tr>
	<tr><th scope=row>C1orf159</th><td>-0.74435833</td><td>-0.3953798</td><td>-0.7297748</td><td>0.37331865</td><td> 3.2449189</td><td>-0.22182830</td></tr>
	<tr><th scope=row>TNFRSF4</th><td> 1.08404422</td><td>-0.2268084</td><td>-0.5309565</td><td>1.05845497</td><td>-0.3450251</td><td>-0.08540632</td></tr>
</tbody>
</table>




```R
colnames(AB1_A2@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li></ol>




```R
options(repr.plot.width=13, repr.plot.height=11)
Idents(AB1_A2) <- 'orig.ident'
FeatureScatter(object = AB1_A2, feature1 = 'nCount_Spatial', feature2 = 'nFeature_Spatial')
```


    
![png](output_69_0.png)
    



```R
options(repr.plot.width=13, repr.plot.height=11)
Idents(AB1_A2) <- 'seurat_clusters'
FeatureScatter(object = AB1_A2, feature1 = 'nCount_Spatial', feature2 = 'nFeature_Spatial')
```


    
![png](output_70_0.png)
    



```R

```


```R
options(repr.plot.width=28, repr.plot.height=9)
plot <- SpatialFeaturePlot(AB1_A2, features = c("HLA-DRA"))

png(filename="Output/01_Process_Visium_and_Integrate/HLA-DRA.png", width=2600, height=900)
print(plot)
dev.off()

plot
```


<strong>png:</strong> 2



    
![png](output_72_1.png)
    


### Find DE genes


```R
AB1_A2
```


    An object of class Seurat 
    35796 features across 9711 samples within 2 assays 
    Active assay: Spatial (17943 features, 2000 variable features)
     3 layers present: data, counts, scale.data
     1 other assay present: SCT
     3 dimensional reductions calculated: pca, integrated.cca, umap
     3 spatial fields of view present: slice1 slice1.2 slice1.3



```R
AB1_A2@assays$SCT@SCTModel.list
```


    $round2
    An sctransform model.
      Model formula:  y ~ log_umi 
      Parameters stored for 17807 features, 4316 cells.
    
    $round1
    An sctransform model.
      Model formula:  y ~ log_umi 
      Parameters stored for 16989 features, 5395 cells.




```R
# copy the Spatial assay into one named RNA
AB1_A2[["RNA"]] = AB1_A2[["Spatial"]]
AB1_A2@assays
```

    Warning message:
    “Key ‘spatial_’ taken, using ‘rna_’ instead”



    $Spatial
    Assay (v5) data with 17943 features for 9711 cells
    Top 10 variable features:
     SCGB1A1, IGLV3-1, IGHG3, IGKV4-1, IGLC1, IGHG4, IGHG2, IGHM, IGHA1,
    IGHG1 
    Layers:
     data, counts, scale.data 
    
    $SCT
    SCTAssay data with 17853 features for 9711 cells, and 2 SCTModel(s) 
    Top 10 variable features:
     IGHG1, IGKC, SCGB1A1, IGHG3, IGKV4-1, IGLC1, IGHG2, IGLV3-1, IGHG4,
    IGHA1 
    
    $RNA
    Assay (v5) data with 17943 features for 9711 cells
    Top 10 variable features:
     SCGB1A1, IGLV3-1, IGHG3, IGKV4-1, IGLC1, IGHG4, IGHG2, IGHM, IGHA1,
    IGHG1 
    Layers:
     data, counts, scale.data 




```R
# fix to allow for prepsctfindmarkers....
#https://github.com/satijalab/seurat/issues/8235
slot(object = AB1_A2@assays$SCT@SCTModel.list[[2]], name="umi.assay")<-"Spatial"
```


```R
SCTResults(object=AB1_A2, slot="umi.assay")
```


<dl>
	<dt>$round2</dt>
		<dd>'Spatial'</dd>
	<dt>$round1</dt>
		<dd>'Spatial'</dd>
</dl>




```R
#data_AB_prep <- PrepSCTFindMarkers(data_AB, assay = "SCT", verbose = TRUE)
AB1_A2_prep <- PrepSCTFindMarkers(AB1_A2, verbose = TRUE)
```

    Found 2 SCT models. Recorrecting SCT counts using minimum median counts: 3962
    



```R
all_de_markers_AB1_A2 <- FindAllMarkers(AB1_A2_prep)
```

    Calculating cluster 0
    
    Calculating cluster 1
    
    Calculating cluster 2
    
    Calculating cluster 3
    
    Calculating cluster 4
    
    Calculating cluster 5
    
    Calculating cluster 6
    
    Calculating cluster 7
    
    Calculating cluster 8
    
    Calculating cluster 9
    
    Calculating cluster 10
    
    Calculating cluster 11
    



```R
dim(all_de_markers_AB1_A2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>84336</li><li>7</li></ol>




```R
all_de_markers_AB1_A2 %>% head(25)
```


<table class="dataframe">
<caption>A data.frame: 25 × 7</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>cluster</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>IGKV4-1</th><td> 0.000000e+00</td><td>2.7515734</td><td>0.898</td><td>0.654</td><td> 0.000000e+00</td><td>0</td><td>IGKV4-1</td></tr>
	<tr><th scope=row>IGHG1</th><td> 0.000000e+00</td><td>2.6662994</td><td>0.999</td><td>0.928</td><td> 0.000000e+00</td><td>0</td><td>IGHG1  </td></tr>
	<tr><th scope=row>IGKC</th><td> 0.000000e+00</td><td>2.6473265</td><td>1.000</td><td>0.990</td><td> 0.000000e+00</td><td>0</td><td>IGKC   </td></tr>
	<tr><th scope=row>JCHAIN</th><td>2.670260e-292</td><td>1.9611859</td><td>0.840</td><td>0.523</td><td>4.791248e-288</td><td>0</td><td>JCHAIN </td></tr>
	<tr><th scope=row>IGLC1</th><td>2.034847e-266</td><td>2.7933805</td><td>0.932</td><td>0.719</td><td>3.651126e-262</td><td>0</td><td>IGLC1  </td></tr>
	<tr><th scope=row>MZB1</th><td>2.411646e-247</td><td>1.9535183</td><td>0.542</td><td>0.168</td><td>4.327216e-243</td><td>0</td><td>MZB1   </td></tr>
	<tr><th scope=row>IGHG2</th><td>1.387158e-235</td><td>2.5724962</td><td>0.827</td><td>0.537</td><td>2.488978e-231</td><td>0</td><td>IGHG2  </td></tr>
	<tr><th scope=row>COL1A1</th><td>3.138237e-234</td><td>1.1382459</td><td>0.961</td><td>0.833</td><td>5.630938e-230</td><td>0</td><td>COL1A1 </td></tr>
	<tr><th scope=row>IGHG3</th><td>3.219171e-208</td><td>2.8930975</td><td>0.765</td><td>0.543</td><td>5.776159e-204</td><td>0</td><td>IGHG3  </td></tr>
	<tr><th scope=row>COL3A1</th><td>4.336263e-192</td><td>1.0442515</td><td>0.929</td><td>0.770</td><td>7.780558e-188</td><td>0</td><td>COL3A1 </td></tr>
	<tr><th scope=row>COL1A2</th><td>3.056144e-190</td><td>1.0291563</td><td>0.946</td><td>0.805</td><td>5.483639e-186</td><td>0</td><td>COL1A2 </td></tr>
	<tr><th scope=row>SPARC</th><td>4.633183e-187</td><td>0.9466975</td><td>0.921</td><td>0.736</td><td>8.313321e-183</td><td>0</td><td>SPARC  </td></tr>
	<tr><th scope=row>IGLJ1</th><td>5.367360e-187</td><td>2.6147123</td><td>0.283</td><td>0.052</td><td>9.630654e-183</td><td>0</td><td>IGLJ1  </td></tr>
	<tr><th scope=row>PIM2</th><td>9.980481e-178</td><td>1.9814158</td><td>0.494</td><td>0.179</td><td>1.790798e-173</td><td>0</td><td>PIM2   </td></tr>
	<tr><th scope=row>IGHA1</th><td>2.522246e-176</td><td>2.2742247</td><td>0.946</td><td>0.720</td><td>4.525666e-172</td><td>0</td><td>IGHA1  </td></tr>
	<tr><th scope=row>LUM</th><td>1.254825e-172</td><td>1.2202680</td><td>0.843</td><td>0.572</td><td>2.251533e-168</td><td>0</td><td>LUM    </td></tr>
	<tr><th scope=row>IGHM</th><td>4.301719e-172</td><td>1.9596216</td><td>0.727</td><td>0.425</td><td>7.718574e-168</td><td>0</td><td>IGHM   </td></tr>
	<tr><th scope=row>CD79A</th><td>9.457659e-172</td><td>2.1277262</td><td>0.350</td><td>0.089</td><td>1.696988e-167</td><td>0</td><td>CD79A  </td></tr>
	<tr><th scope=row>IGHJ2</th><td>1.295759e-159</td><td>2.8457166</td><td>0.346</td><td>0.095</td><td>2.324981e-155</td><td>0</td><td>IGHJ2  </td></tr>
	<tr><th scope=row>TIMP1</th><td>1.310466e-159</td><td>0.8044986</td><td>0.950</td><td>0.848</td><td>2.351370e-155</td><td>0</td><td>TIMP1  </td></tr>
	<tr><th scope=row>BGN</th><td>5.394207e-156</td><td>0.8125985</td><td>0.929</td><td>0.744</td><td>9.678826e-152</td><td>0</td><td>BGN    </td></tr>
	<tr><th scope=row>FN1</th><td>2.358229e-153</td><td>0.7846811</td><td>0.983</td><td>0.910</td><td>4.231371e-149</td><td>0</td><td>FN1    </td></tr>
	<tr><th scope=row>B2M</th><td>6.994892e-150</td><td>0.4452693</td><td>0.999</td><td>0.987</td><td>1.255093e-145</td><td>0</td><td>B2M    </td></tr>
	<tr><th scope=row>ACTA2</th><td>1.524151e-148</td><td>0.9896355</td><td>0.864</td><td>0.670</td><td>2.734784e-144</td><td>0</td><td>ACTA2  </td></tr>
	<tr><th scope=row>AEBP1</th><td>7.902222e-148</td><td>1.0034582</td><td>0.875</td><td>0.627</td><td>1.417896e-143</td><td>0</td><td>AEBP1  </td></tr>
</tbody>
</table>




```R
require(openxlsx)
AB1_A2_cluster_datasets <- list("cluster0" = subset(all_de_markers_AB1_A2, cluster == 0), 
                         "cluster1" = subset(all_de_markers_AB1_A2, cluster == 1), 
                         "cluster2" = subset(all_de_markers_AB1_A2, cluster == 2), 
                         "cluster3" = subset(all_de_markers_AB1_A2, cluster == 3), 
                         "cluster4" = subset(all_de_markers_AB1_A2, cluster == 4), 
                         "cluster5" = subset(all_de_markers_AB1_A2, cluster == 5), 
                         "cluster6" = subset(all_de_markers_AB1_A2, cluster == 6), 
                         "cluster7" = subset(all_de_markers_AB1_A2, cluster == 7),
                         "cluster8" = subset(all_de_markers_AB1_A2, cluster == 8), 
                           "cluster9" = subset(all_de_markers_AB1_A2, cluster == 9), 
                           "cluster10" = subset(all_de_markers_AB1_A2, cluster == 10), 
                           "cluster11" = subset(all_de_markers_AB1_A2, cluster == 11))
write.xlsx(AB1_A2_cluster_datasets, file = "Output/01_Process_Visium_and_Integrate/AB1_A2_cluster_genes.xlsx")
```


```R
top20 <- all_de_markers_AB1_A2 %>% group_by(cluster) %>% top_n(20, avg_log2FC)
#DoHeatmap(object = obj, genes.use = top20$gene, slim.col.label = TRUE, remove.key = TRUE)
```


```R
options(repr.plot.width=16, repr.plot.height=16)

plot <- DoHeatmap(AB1_A2_prep, features = top20$gene)

png(filename="Output/01_Process_Visium_and_Integrate/Clusters_Heatmap_AB1_A2.png", width=1600, height=1600)
print(plot)
dev.off()

plot
```

    Warning message in DoHeatmap(AB1_A2_prep, features = top20$gene):
    “The following features were omitted as they were not found in the scale.data slot for the Spatial assay: TEX13C, LRIT1, PRLH, NKX1-1, NXPE4, KCNG4, CLC, C16orf92, ANKRD23, CDH20, ZNRF4, OR5M9, PRSS51, SERPINE3, LYZL4, KIAA1549L, SERPINA7, WFDC6, CFAP73, AC113348.1, LRRTM1, BOD1L2, COL9A3, CTNND2, THSD7B, PHOX2A, POU3F3, TERT, CHRNA7, SLURP2, NOS1, EPHA6, CENPVL3, OR1A1, CILP, IGFL3, RSPO3, LMLN2, FCRLA, LDHAL6A, CAPN11, SPATA48, FLRT2, PLPP4, GLI1, C5orf46, ORM1, DUX4, C4orf45, OR5R1, SERPINB10, LILRA4, CYP1A2, TRIM71, SYCE1, AGTR2, MUC21, SLC5A8, ADAMTS9, LRCOL1, PKLR, KLHL33, OTP, AMER2, FAM133A, ANKRD34B, SIX6, CTXND2, HS6ST3, CHRNA4, FRMD1, FBXO47, SPRR2G, FADS6, C3orf22, MS4A3, SNX32, ZNF735, ACAN, APOH, DCAF8L2, AMER3, CNTNAP5, PITX2, KLHL34, GPR143, SLC7A10, MROH2A, VSTM5, RNF175, PTCHD1, DMRT1, AQP6, C1orf94, NGEF, C14orf180, GPR137C, CAPNS2, ZYG11A, SUSD4, AQP5, TRDN, IFNG, FAM71C, ENTHD1, CLEC5A, RETN, ATP6V0D2, PDCD1LG2, TM4SF19, MRO, CCL7, HTRA4, MARCO, AQP9, C1QTNF4, CORO6, RPGRIP1, ZNF781, MCTP1, MN1, EBF1, PDE1A, SKA1, C2orf81, PDE3A, HOXA3, CD69, ZNF396, RASSF7, LRRC37B, ACER2, ABHD17A, LY6G5C, MLNR, CACNG1, BCO1, DOK7, UNC5D, LANCL3, BHLHE23, TIMD4, TNP1, TAC1, AHSP, KIRREL3, AMBN, SBK2, NPY1R, C20orf203, MYO3B, SLC5A1, CYP21A2, PNCK, LOXL4, NME8”



<strong>png:</strong> 2



    
![png](output_85_2.png)
    



```R

```


```R
all_de_markers_AB1_A2 %>%
    group_by(cluster) %>%
    dplyr::filter(avg_log2FC > 1) %>%
    slice_head(n = 20) %>%
    ungroup() -> top20_2

hm2_plot <- DoHeatmap(AB1_A2_prep, features = top20_2$gene) + NoLegend()

png(filename="Output/01_Process_Visium_and_Integrate/Clusters_Heatmap_AB1_A2_2.png", width=1600, height=1600)
print(hm2_plot)
dev.off()

hm2_plot
```

    Warning message in DoHeatmap(AB1_A2_prep, features = top20_2$gene):
    “The following features were omitted as they were not found in the scale.data slot for the Spatial assay: BNIP5, KBTBD7, SERPINA7, ARG2, FGF9, C18orf21, VPS37D, AARSD1, ZMYM3, SLU7, NPAS4, CENPVL3, G0S2, PSAP, UBC, ICAM1, C11orf96, LPCAT1, ZFP36, SFTA3, ARHGAP23, SOCS7, TUBA1B, NKX2-1, PIP4K2B, LASP1, GRB7, STARD3, CDK12, PSMB3, MLLT6, PCGF2, ERBB2, MAPK8IP2, CDK5R2, SCN3A, SMG1, NFIB, CXCL17, TREM1, SPI1, SLC11A1, OGG1, MLLT3, HEATR6, DGLUCY, ZSCAN26, SIN3B, ZNF195, NMD3, FUZ, ZNF212, GATAD2B, RGS3, SMARCA4, MAN2B2, RRP36, INO80E, SLC35F5, UGGT1, INHA, ATP4A, ESRRG, NGF, CIART, GGT1, CYP21A2, FGF11, SLC26A2, PLXDC1, PNCK”



<strong>png:</strong> 2



    
![png](output_87_2.png)
    



```R

```


```R

```


```R

```
