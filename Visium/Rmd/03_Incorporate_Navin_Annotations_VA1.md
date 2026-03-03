## March 2025
## Round 1 REVISIONS for SCLC 10X Visium
## I annotated spots for vasc, Tumor, and MHCI. 
## Add Navin's Annotations for VA1 and see which clusters and groups they represent. 


```R
library(Seurat)
#library(SeuratExtend)
#library(spacexr)
library(Matrix)
library(doParallel)
library(ggplot2)
library(org.Hs.eg.db)
library(tidyverse)
library(dittoSeq)
library(pheatmap)
library(GSVA)
library(ggplot2)
library(scales)
#library('scCustomize')
```

    Loading required package: SeuratObject
    
    Loading required package: sp
    
    
    Attaching package: ‘SeuratObject’
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, t
    
    
    Loading required package: foreach
    
    Loading required package: iterators
    
    Loading required package: parallel
    
    Loading required package: AnnotationDbi
    
    Loading required package: stats4
    
    Loading required package: BiocGenerics
    
    Loading required package: generics
    
    
    Attaching package: ‘generics’
    
    
    The following objects are masked from ‘package:base’:
    
        as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
        setequal, union
    
    
    
    Attaching package: ‘BiocGenerics’
    
    
    The following objects are masked from ‘package:stats’:
    
        IQR, mad, sd, var, xtabs
    
    
    The following objects are masked from ‘package:base’:
    
        anyDuplicated, aperm, append, as.data.frame, basename, cbind,
        colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
        get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
        order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
        rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
        unsplit, which.max, which.min
    
    
    Loading required package: Biobase
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    Loading required package: IRanges
    
    Loading required package: S4Vectors
    
    
    Attaching package: ‘S4Vectors’
    
    
    The following objects are masked from ‘package:Matrix’:
    
        expand, unname
    
    
    The following object is masked from ‘package:utils’:
    
        findMatches
    
    
    The following objects are masked from ‘package:base’:
    
        expand.grid, I, unname
    
    
    
    Attaching package: ‘IRanges’
    
    
    The following object is masked from ‘package:sp’:
    
        %over%
    
    
    
    
    ── [1mAttaching core tidyverse packages[22m ──────────────────────── tidyverse 2.0.0 ──
    [32m✔[39m [34mdplyr    [39m 1.1.4     [32m✔[39m [34mreadr    [39m 2.1.5
    [32m✔[39m [34mforcats  [39m 1.0.0     [32m✔[39m [34mstringr  [39m 1.5.1
    [32m✔[39m [34mlubridate[39m 1.9.4     [32m✔[39m [34mtibble   [39m 3.2.1
    [32m✔[39m [34mpurrr    [39m 1.0.4     [32m✔[39m [34mtidyr    [39m 1.3.1
    ── [1mConflicts[22m ────────────────────────────────────────── tidyverse_conflicts() ──
    [31m✖[39m [34mlubridate[39m::[32m%within%()[39m    masks [34mIRanges[39m::%within%()
    [31m✖[39m [34mpurrr[39m::[32maccumulate()[39m      masks [34mforeach[39m::accumulate()
    [31m✖[39m [34mdplyr[39m::[32mcollapse()[39m        masks [34mIRanges[39m::collapse()
    [31m✖[39m [34mdplyr[39m::[32mcombine()[39m         masks [34mBiobase[39m::combine(), [34mBiocGenerics[39m::combine()
    [31m✖[39m [34mdplyr[39m::[32mdesc()[39m            masks [34mIRanges[39m::desc()
    [31m✖[39m [34mtidyr[39m::[32mexpand()[39m          masks [34mS4Vectors[39m::expand(), [34mMatrix[39m::expand()
    [31m✖[39m [34mdplyr[39m::[32mfilter()[39m          masks [34mstats[39m::filter()
    [31m✖[39m [34mdplyr[39m::[32mfirst()[39m           masks [34mS4Vectors[39m::first()
    [31m✖[39m [34mdplyr[39m::[32mlag()[39m             masks [34mstats[39m::lag()
    [31m✖[39m [34mtidyr[39m::[32mpack()[39m            masks [34mMatrix[39m::pack()
    [31m✖[39m [34mBiocGenerics[39m::[32mPosition()[39m masks [34mggplot2[39m::Position(), [34mbase[39m::Position()
    [31m✖[39m [34mpurrr[39m::[32mreduce()[39m          masks [34mIRanges[39m::reduce()
    [31m✖[39m [34mdplyr[39m::[32mrename()[39m          masks [34mS4Vectors[39m::rename()
    [31m✖[39m [34mlubridate[39m::[32msecond()[39m      masks [34mS4Vectors[39m::second()
    [31m✖[39m [34mlubridate[39m::[32msecond<-()[39m    masks [34mS4Vectors[39m::second<-()
    [31m✖[39m [34mdplyr[39m::[32mselect()[39m          masks [34mAnnotationDbi[39m::select()
    [31m✖[39m [34mdplyr[39m::[32mslice()[39m           masks [34mIRanges[39m::slice()
    [31m✖[39m [34mtidyr[39m::[32munpack()[39m          masks [34mMatrix[39m::unpack()
    [31m✖[39m [34mpurrr[39m::[32mwhen()[39m            masks [34mforeach[39m::when()
    [36mℹ[39m Use the conflicted package ([3m[34m<http://conflicted.r-lib.org/>[39m[23m) to force all conflicts to become errors
    
    Attaching package: ‘scales’
    
    
    The following object is masked from ‘package:purrr’:
    
        discard
    
    
    The following object is masked from ‘package:readr’:
    
        col_factor
    
    


# Load integrated and annotated datasets


```R
# load Seurat Integrated objects
AB1_A2 <- readRDS('Processing/AB1_A2_Annotated_step1_20250326.rds')
```


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
print(AB1_A2@meta.data %>% head(3))
```

                           orig.ident nCount_Spatial nFeature_Spatial visium_round
    VA1_AAACAACGAATAGTTC-1        VA1           3625             2571       round1
    VA1_AAACAAGTATCTCCCA-1        VA1           5804             3329       round1
    VA1_AAACAATCTACTAGCA-1        VA1           4980             3146       round1
                           nCount_SCT nFeature_SCT SCT_snn_res.0.8 seurat_clusters
    VA1_AAACAACGAATAGTTC-1       3649         2571              10              10
    VA1_AAACAAGTATCTCCCA-1       4956         3308               6               6
    VA1_AAACAATCTACTAGCA-1       4711         3146               6               6
                           slice_ident vasc_summed_expr vasc_label Tumor_score
    VA1_AAACAACGAATAGTTC-1         VA1         2.648104  vasc_high          11
    VA1_AAACAAGTATCTCCCA-1         VA1         1.001716   vasc_low          13
    VA1_AAACAATCTACTAGCA-1         VA1         1.101286   vasc_low          19
                           Tumor_label1 Tumor_label MHCI_score MHCI_label1
    VA1_AAACAACGAATAGTTC-1    Tumor_low   Tumor_low          2    MHCI_low
    VA1_AAACAAGTATCTCCCA-1    Tumor_low   Tumor_low          2    MHCI_low
    VA1_AAACAATCTACTAGCA-1   Tumor_high  Tumor_high          3    MHCI_low
                           MHCI_label    Tumor_MHCI_label   Tumor_MHCI_label2
    VA1_AAACAACGAATAGTTC-1   MHCI_low  Tumor_low_MHCI_low           Tumor_low
    VA1_AAACAAGTATCTCCCA-1   MHCI_low  Tumor_low_MHCI_low           Tumor_low
    VA1_AAACAATCTACTAGCA-1   MHCI_low Tumor_high_MHCI_low Tumor_high_MHCI_low
                           vasc_label1
    VA1_AAACAACGAATAGTTC-1   vasc_high
    VA1_AAACAAGTATCTCCCA-1    vasc_low
    VA1_AAACAATCTACTAGCA-1    vasc_low


# Load Navin's manual annotations


```R
navin_annot_VA <- read.csv('/immuno/ian/Projects/2023_02_SCLC_10X_Visium/Annotated_cloupe/Navin_annotations_VA.csv')
```


```R
# Add "VA1_" to Barcode to match Seurat's barcode naming
navin_annot_VA$Barcode <- lapply(navin_annot_VA$Barcode, function(x) paste0('VA1_', trimws(x)))
```


```R
# convert to character
navin_annot_VA$Barcode <- as.character(navin_annot_VA$Barcode)
```


```R
head(navin_annot_VA,3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 2</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Navin_annotations</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>Respiratory epithelium</td></tr>
	<tr><th scope=row>2</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>SCLC                  </td></tr>
	<tr><th scope=row>3</th><td>VA1_AAACAATCTACTAGCA-1</td><td>SCLC                  </td></tr>
</tbody>
</table>



# Try to merge navin's annotations in with the metadata...


```R
dim(navin_annot_VA)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>2806</li><li>2</li></ol>




```R
dim(AB1_A2[[]])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>20</li></ol>




```R
# add barcode as a column of metadata
AB1_A2@meta.data$Barcode <- rownames(AB1_A2@meta.data)
```


```R
'VA_AAACAACGAATAGTTC-1' == 'VA_AAACAACGAATAGTTC-1'
```


TRUE



```R
AB1_A2[[]][0:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 21</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>⋯</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Barcode</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>⋯</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>VA1_AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>⋯</td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>VA1_AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>⋯</td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>VA1_AAACAATCTACTAGCA-1</td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td>4 </td><td>4 </td><td>VA1</td><td>0.907884</td><td>⋯</td><td> 9</td><td>Tumor_low </td><td>Tumor_low </td><td>1</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>VA1_AAACACCAATAACTGC-1</td></tr>
	<tr><th scope=row>VA1_AAACAGCTTTCAGAAG-1</th><td>VA1</td><td>5132</td><td>3219</td><td>round1</td><td>4797</td><td>3218</td><td>4 </td><td>4 </td><td>VA1</td><td>3.751279</td><td>⋯</td><td>14</td><td>Tumor_low </td><td>Tumor_low </td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low          </td><td>vasc_high</td><td>VA1_AAACAGCTTTCAGAAG-1</td></tr>
</tbody>
</table>




```R
# how many barcodes actually overlap?
overlap_list <- Reduce(intersect,list(AB1_A2[[]]$Barcode,navin_annot_VA$Barcode))
print(length(overlap_list))
```

    [1] 2806



```R
AB1_A2@meta.data <- left_join(AB1_A2[[]], navin_annot_VA, by = join_by(Barcode == Barcode))

# add the rownames back in, since they are removed during a merge in R...
rownames(AB1_A2@meta.data) <- AB1_A2@meta.data$Barcode
```


```R
dim(AB1_A2[[]])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>22</li></ol>




```R
AB1_A2[[]][0:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 22</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>⋯</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Barcode</th><th scope=col>Navin_annotations</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>⋯</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>VA1_AAACAACGAATAGTTC-1</td><td>Respiratory epithelium</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>⋯</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>VA1_AAACAAGTATCTCCCA-1</td><td>SCLC                  </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>⋯</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>VA1_AAACAATCTACTAGCA-1</td><td>SCLC                  </td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td>4 </td><td>4 </td><td>VA1</td><td>0.907884</td><td>⋯</td><td>Tumor_low </td><td>Tumor_low </td><td>1</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>VA1_AAACACCAATAACTGC-1</td><td>SCLC+vasculature      </td></tr>
	<tr><th scope=row>VA1_AAACAGCTTTCAGAAG-1</th><td>VA1</td><td>5132</td><td>3219</td><td>round1</td><td>4797</td><td>3218</td><td>4 </td><td>4 </td><td>VA1</td><td>3.751279</td><td>⋯</td><td>Tumor_low </td><td>Tumor_low </td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low          </td><td>vasc_high</td><td>VA1_AAACAGCTTTCAGAAG-1</td><td>SCLC+vasculature      </td></tr>
</tbody>
</table>




```R
table(AB1_A2@meta.data$Navin_annotations, useNA = "always")
```


    
                                                  Alveoli 
                           413                         12 
             Endothelial cells                   Fibrosis 
                            52                        155 
      Immune cells_lymphocytes Immune cells_myeloid cells 
                           192                         29 
      Immune cells+vasculature     Respiratory epithelium 
                            68                         30 
                          SCLC                   SCLC_TSI 
                          1233                        112 
                 SCLC+fibrosis          SCLC+immune cells 
                            76                         69 
              SCLC+vasculature              Smooth muscle 
                           309                         20 
                Tumor necrosis                       <NA> 
                            36                       6905 



```R
colSums(!is.na(AB1_A2@meta.data))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>orig.ident</dt><dd>9711</dd><dt>nCount_Spatial</dt><dd>9711</dd><dt>nFeature_Spatial</dt><dd>9711</dd><dt>visium_round</dt><dd>9711</dd><dt>nCount_SCT</dt><dd>9711</dd><dt>nFeature_SCT</dt><dd>9711</dd><dt>SCT_snn_res.0.8</dt><dd>9711</dd><dt>seurat_clusters</dt><dd>9711</dd><dt>slice_ident</dt><dd>9711</dd><dt>vasc_summed_expr</dt><dd>9711</dd><dt>vasc_label</dt><dd>9711</dd><dt>Tumor_score</dt><dd>9711</dd><dt>Tumor_label1</dt><dd>9711</dd><dt>Tumor_label</dt><dd>9711</dd><dt>MHCI_score</dt><dd>9711</dd><dt>MHCI_label1</dt><dd>9711</dd><dt>MHCI_label</dt><dd>9711</dd><dt>Tumor_MHCI_label</dt><dd>9711</dd><dt>Tumor_MHCI_label2</dt><dd>9711</dd><dt>vasc_label1</dt><dd>9711</dd><dt>Barcode</dt><dd>9711</dd><dt>Navin_annotations</dt><dd>2806</dd></dl>



## Add Tissue_Slice VA and VB into the metadata...


```R
AB1_A2@meta.data$Tissue_Slice <- substr(AB1_A2@meta.data[,'Barcode'], 1, 3)
```


```R
AB1_A2@meta.data[0:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 23</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>⋯</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Barcode</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>⋯</td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>VA1_AAACAACGAATAGTTC-1</td><td>Respiratory epithelium</td><td>VA1</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>⋯</td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>VA1_AAACAAGTATCTCCCA-1</td><td>SCLC                  </td><td>VA1</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>⋯</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>VA1_AAACAATCTACTAGCA-1</td><td>SCLC                  </td><td>VA1</td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td>4 </td><td>4 </td><td>VA1</td><td>0.907884</td><td>⋯</td><td>Tumor_low </td><td>1</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>VA1_AAACACCAATAACTGC-1</td><td>SCLC+vasculature      </td><td>VA1</td></tr>
	<tr><th scope=row>VA1_AAACAGCTTTCAGAAG-1</th><td>VA1</td><td>5132</td><td>3219</td><td>round1</td><td>4797</td><td>3218</td><td>4 </td><td>4 </td><td>VA1</td><td>3.751279</td><td>⋯</td><td>Tumor_low </td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low          </td><td>vasc_high</td><td>VA1_AAACAGCTTTCAGAAG-1</td><td>SCLC+vasculature      </td><td>VA1</td></tr>
</tbody>
</table>




```R
print(unique(AB1_A2@meta.data$Tissue_Slice))
```

    [1] "VA1" "VA2" "VB1"


# Looks like navin's annotations have successfully been integrated with the metadata. Now we can try to plot some things...

# Since we only have navin annotations for VA1, let's subset data to just look at VA for now...


```R
data_VA <- subset(AB1_A2, subset = slice_ident == 'VA1')
print(dim(AB1_A2))
print(dim(data_VA))
```

    Warning message:
    “Not validating Seurat objects”
    Warning message:
    “Not validating Centroids objects”
    Warning message:
    “Not validating Centroids objects”
    Warning message:
    “Not validating FOV objects”
    Warning message:
    “Not validating FOV objects”
    Warning message:
    “Not validating FOV objects”
    Warning message:
    “Not validating FOV objects”
    Warning message:
    “Not validating FOV objects”
    Warning message:
    “Not validating FOV objects”
    Warning message:
    “Not validating Seurat objects”
    Warning message:
    “Not validating Seurat objects”


    [1] 17943  9711
    [1] 17943  2806



```R
print(dim(AB1_A2[[]]))
print(dim(data_VA[[]]))
```

    [1] 9711   23
    [1] 2806   23



```R
options(repr.plot.width=15, repr.plot.height=7)
```


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
data_VA
```


    An object of class Seurat 
    35796 features across 2806 samples within 2 assays 
    Active assay: Spatial (17943 features, 2000 variable features)
     3 layers present: data, counts, scale.data
     1 other assay present: SCT
     3 dimensional reductions calculated: pca, integrated.cca, umap
     1 spatial field of view present: slice1.2



```R
colnames(AB1_A2[[]])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Barcode'</li><li>'Navin_annotations'</li><li>'Tissue_Slice'</li></ol>




```R
colnames(data_VA[[]])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Barcode'</li><li>'Navin_annotations'</li><li>'Tissue_Slice'</li></ol>




```R
Images(AB1_A2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'slice1'</li><li>'slice1.2'</li><li>'slice1.3'</li></ol>




```R
Images(data_VA)
```


'slice1.2'



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
Assays(data_VA)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Spatial'</li><li>'SCT'</li></ol>




```R
dim(data_VA[['Spatial']])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>17943</li><li>2806</li></ol>



# Let's simplify and group some of Navin's annotations together


```R
print(unique(data_VA[[]]$Navin_annotations))
```

     [1] "Respiratory epithelium"     "SCLC"                      
     [3] "SCLC+vasculature"           "Fibrosis"                  
     [5] ""                           "Immune cells+vasculature"  
     [7] "Immune cells_lymphocytes"   "SCLC_TSI"                  
     [9] "Immune cells_myeloid cells" "SCLC+immune cells"         
    [11] "Smooth muscle"              "Endothelial cells"         
    [13] "Tumor necrosis"             "SCLC+fibrosis"             
    [15] "Alveoli"                   



```R
data_VA[[]]$Navin_annotations_simplified <- recode(data_VA[[]]$Navin_annotations, 
                                                   " " = " ",
                                                   `Respiratory epithelium` = "Respiratory epithelium", 
                                                   `SCLC+vasculature` = "SCLC", 
                                                   `SCLC` = "SCLC", 
                                                   `SCLC_TSI` = "SCLC", 
                                                   `SCLC+fibrosis` = "SCLC", 
                                                   `SCLC+immune cells` = "SCLC+Immune", 
                                                   `Immune cells_lymphocytes` = "Immune", 
                                                   `Immune cells+vasculature` = "Immune", 
                                                   `Immune cells_myeloid cells` = "Immune", 
                                                   `Smooth muscle` = "Smooth muscle", 
                                                   `Endothelial cells` = "Endothelial cells", 
                                                   `Tumor necrosis` = "Tumor necrosis", 
                                                   `Alveoli` = "Alveoli", 
                                                   `Fibrosis` = "Fibrosis")
```


```R
print(table(data_VA[[]]$Navin_annotations, useNA = "always"))
print(table(data_VA[[]]$Navin_annotations_simplified, useNA = "always"))
```

    
                                                  Alveoli 
                           413                         12 
             Endothelial cells                   Fibrosis 
                            52                        155 
      Immune cells_lymphocytes Immune cells_myeloid cells 
                           192                         29 
      Immune cells+vasculature     Respiratory epithelium 
                            68                         30 
                          SCLC                   SCLC_TSI 
                          1233                        112 
                 SCLC+fibrosis          SCLC+immune cells 
                            76                         69 
              SCLC+vasculature              Smooth muscle 
                           309                         20 
                Tumor necrosis                       <NA> 
                            36                          0 
    
                                          Alveoli      Endothelial cells 
                       413                     12                     52 
                  Fibrosis                 Immune Respiratory epithelium 
                       155                    289                     30 
                      SCLC            SCLC+Immune          Smooth muscle 
                      1730                     69                     20 
            Tumor necrosis                   <NA> 
                        36                      0 



```R
print(dim(data_VA[[]]))
```

    [1] 2806   24


# do the same for AB1_A2


```R
AB1_A2[[]]$Navin_annotations_simplified <- recode(AB1_A2[[]]$Navin_annotations, 
                                                   " " = " ",
                                                   `Respiratory epithelium` = "Respiratory epithelium", 
                                                   `SCLC+vasculature` = "SCLC", 
                                                   `SCLC` = "SCLC", 
                                                   `SCLC_TSI` = "SCLC", 
                                                   `SCLC+fibrosis` = "SCLC", 
                                                   `SCLC+immune cells` = "SCLC+Immune", 
                                                   `Immune cells_lymphocytes` = "Immune", 
                                                   `Immune cells+vasculature` = "Immune", 
                                                   `Immune cells_myeloid cells` = "Immune", 
                                                   `Smooth muscle` = "Smooth muscle", 
                                                   `Endothelial cells` = "Endothelial cells", 
                                                   `Tumor necrosis` = "Tumor necrosis", 
                                                   `Alveoli` = "Alveoli", 
                                                   `Fibrosis` = "Fibrosis")
```


```R
print(table(AB1_A2[[]]$Navin_annotations, useNA = "always"))
print(table(AB1_A2[[]]$Navin_annotations_simplified, useNA = "always"))
```

    
                                                  Alveoli 
                           413                         12 
             Endothelial cells                   Fibrosis 
                            52                        155 
      Immune cells_lymphocytes Immune cells_myeloid cells 
                           192                         29 
      Immune cells+vasculature     Respiratory epithelium 
                            68                         30 
                          SCLC                   SCLC_TSI 
                          1233                        112 
                 SCLC+fibrosis          SCLC+immune cells 
                            76                         69 
              SCLC+vasculature              Smooth muscle 
                           309                         20 
                Tumor necrosis                       <NA> 
                            36                       6905 
    
                                          Alveoli      Endothelial cells 
                       413                     12                     52 
                  Fibrosis                 Immune Respiratory epithelium 
                       155                    289                     30 
                      SCLC            SCLC+Immune          Smooth muscle 
                      1730                     69                     20 
            Tumor necrosis                   <NA> 
                        36                   6905 



```R
print(dim(AB1_A2[[]]))
```

    [1] 9711   24



```R
colnames(AB1_A2@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Barcode'</li><li>'Navin_annotations'</li><li>'Tissue_Slice'</li><li>'Navin_annotations_simplified'</li></ol>



# Navin also annotated seurat clusters using marker genes and looking at the spots overlapped on the H&E images using the loupe browser. Let's load those annotations in too. 


```R
navin_cluster_annots <- read.csv('Processing/Navin_Annotations_20250331.csv')
```


```R
navin_cluster_annots <- navin_cluster_annots %>% rename('seurat_clusters' = 'Cluster', 'Navin_Cluster_Annotation'='Navin_Annotation')
```


```R
navin_cluster_annots
```


<table class="dataframe">
<caption>A data.frame: 12 × 2</caption>
<thead>
	<tr><th scope=col>seurat_clusters</th><th scope=col>Navin_Cluster_Annotation</th></tr>
	<tr><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td> 0</td><td>Lymphocytes B Cells (and Plasma Cells)</td></tr>
	<tr><td> 1</td><td>SCLC                                  </td></tr>
	<tr><td> 2</td><td>SCLC (TSI)                            </td></tr>
	<tr><td> 3</td><td>SCLC (Stroma)                         </td></tr>
	<tr><td> 4</td><td>SCLC                                  </td></tr>
	<tr><td> 5</td><td>SCLC                                  </td></tr>
	<tr><td> 6</td><td>SCLC + Effector Immune Cells (T+NK)   </td></tr>
	<tr><td> 7</td><td>Normal Lung/Stroma                    </td></tr>
	<tr><td> 8</td><td>Fibroblasts                           </td></tr>
	<tr><td> 9</td><td>SCLC (Apoptotic)                      </td></tr>
	<tr><td>10</td><td>Respiratory Epithelium (Normal Lung)  </td></tr>
	<tr><td>11</td><td>SCLC (Necrotic)                       </td></tr>
</tbody>
</table>


# Example named list
mapping_list <- list("A" = "Group1", "B" = "Group2", "C" = "Group3")

# Example dataframe
df <- data.frame(ID = c("A", "B", "C", "A", "C"), stringsAsFactors = FALSE)

# Use the named list to create a new column
df$new_column <- unlist(mapping_list[df$ID])

# View the updated dataframe
print(df)


```R
AB1_A2@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 24</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>⋯</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Barcode</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>⋯</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>VA1_AAACAACGAATAGTTC-1</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>⋯</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>VA1_AAACAAGTATCTCCCA-1</td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>⋯</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>VA1_AAACAATCTACTAGCA-1</td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td></tr>
</tbody>
</table>




```R
# create a mapping to add these annotations
AB1_A2[[]]$Navin_Clusters1 <- recode(AB1_A2[[]]$seurat_clusters, 
                                                   "0" = "Lymphocytes B Cells (and Plasma Cells)",
                                                   "1" = "SCLC", 
                                                   "2" = "SCLC (TSI)", 
                                                   "3" = "SCLC (Stroma)", 
                                                   "4" = "SCLC", 
                                                   "5" = "SCLC", 
                                                   "6" = "SCLC + Effector Immune Cells (T+NK)", 
                                                   "7" = "Normal Lung/Stroma", 
                                                   "8" = "Fibroblasts", 
                                                   "9" = "SCLC (Apoptotic)", 
                                                   "10" = "Respiratory Epithelium (Normal Lung)", 
                                                   "11" = "SCLC (Necrotic)")
```


```R
# create a mapping to add these annotations
AB1_A2[[]]$Navin_Clusters <- recode(AB1_A2[[]]$seurat_clusters, 
                                                   "0" = "0_Lymphocytes_B_Cells",
                                                   "1" = "1_SCLC", 
                                                   "2" = "2_SCLC_TSI", 
                                                   "3" = "3_SCLC_Stroma", 
                                                   "4" = "4_SCLC", 
                                                   "5" = "5_SCLC", 
                                                   "6" = "6_SCLC_Effector_Immune", 
                                                   "7" = "7_Normal_Lung_Stroma", 
                                                   "8" = "8_Fibroblasts", 
                                                   "9" = "9_SCLC_Apoptotic", 
                                                   "10" = "10_Respiratory_Epithelium_Lung", 
                                                   "11" = "11_SCLC_Necrotic")
```


```R
# create a mapping to add these annotations
AB1_A2[[]]$Navin_Clusters2 <- recode(AB1_A2[[]]$seurat_clusters, 
                                                   "0" = "0_Lymphocytes_B_Cells",
                                                   "1" = "1_4_5_11_SCLC", 
                                                   "2" = "2_SCLC_TSI", 
                                                   "3" = "3_SCLC_Stroma", 
                                                   "4" = "1_4_5_11_SCLC", 
                                                   "5" = "1_4_5_11_SCLC", 
                                                   "6" = "6_SCLC_Effector_Immune", 
                                                   "7" = "7_10_Normal_Lung", 
                                                   "8" = "8_Fibroblasts", 
                                                   "9" = "9_SCLC_Apoptotic", 
                                                   "10" = "7_10_Normal_Lung", 
                                                   "11" = "1_4_5_11_SCLC")
```


```R
AB1_A2@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 27</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>⋯</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Barcode</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>⋯</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>VA1_AAACAACGAATAGTTC-1</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>⋯</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>VA1_AAACAAGTATCTCCCA-1</td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>⋯</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>VA1_AAACAATCTACTAGCA-1</td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td></tr>
</tbody>
</table>




```R
# also add these annotations to data_VA
data_VA[[]]$Navin_Clusters1 <- recode(data_VA[[]]$seurat_clusters, 
                                                   "0" = "Lymphocytes B Cells (and Plasma Cells)",
                                                   "1" = "SCLC", 
                                                   "2" = "SCLC (TSI)", 
                                                   "3" = "SCLC (Stroma)", 
                                                   "4" = "SCLC", 
                                                   "5" = "SCLC", 
                                                   "6" = "SCLC + Effector Immune Cells (T+NK)", 
                                                   "7" = "Normal Lung/Stroma", 
                                                   "8" = "Fibroblasts", 
                                                   "9" = "SCLC (Apoptotic)", 
                                                   "10" = "Respiratory Epithelium (Normal Lung)", 
                                                   "11" = "SCLC (Necrotic)")
data_VA[[]]$Navin_Clusters <- recode(data_VA[[]]$seurat_clusters, 
                                                   "0" = "0_Lymphocytes_B_Cells",
                                                   "1" = "1_SCLC", 
                                                   "2" = "2_SCLC_TSI", 
                                                   "3" = "3_SCLC_Stroma", 
                                                   "4" = "4_SCLC", 
                                                   "5" = "5_SCLC", 
                                                   "6" = "6_SCLC_Effector_Immune", 
                                                   "7" = "7_Normal_Lung_Stroma", 
                                                   "8" = "8_Fibroblasts", 
                                                   "9" = "9_SCLC_Apoptotic", 
                                                   "10" = "10_Respiratory_Epithelium_Lung", 
                                                   "11" = "11_SCLC_Necrotic")
data_VA[[]]$Navin_Clusters2 <- recode(data_VA[[]]$seurat_clusters, 
                                                   "0" = "0_Lymphocytes_B_Cells",
                                                   "1" = "1_4_5_11_SCLC", 
                                                   "2" = "2_SCLC_TSI", 
                                                   "3" = "3_SCLC_Stroma", 
                                                   "4" = "1_4_5_11_SCLC", 
                                                   "5" = "1_4_5_11_SCLC", 
                                                   "6" = "6_SCLC_Effector_Immune", 
                                                   "7" = "7_10_Normal_Lung", 
                                                   "8" = "8_Fibroblasts", 
                                                   "9" = "9_SCLC_Apoptotic", 
                                                   "10" = "7_10_Normal_Lung", 
                                                   "11" = "1_4_5_11_SCLC")
```

# Save data now annotated with Navin's manual annotations


```R
saveRDS(AB1_A2, 'Processing/AB1_A2_Annotated_step2Navin_20250326.rds')
saveRDS(data_VA, 'Processing/A1only_Annotated_step2Navin_20250326.rds')
```

# Dimplots


```R
#options(repr.plot.width=17, repr.plot.height=14)
options(repr.plot.width=20, repr.plot.height=21)

navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Navin_Clusters", "Navin_annotations", 
                                                                "Navin_Clusters2", "Navin_annotations_simplified",
                                                                "Tumor_MHCI_label", "vasc_label"), 
                      ncol = 2,
                      cols = dittoColors()[1:15],
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap

# save
# 1300, 1100 old
png('Output/03_Incorporate_Navin_Annotations_VA1/Navin_Annot_UMAP.png',
    width = 1400,
    height = 1650)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_65_1.png)
    



```R
#options(repr.plot.width=17, repr.plot.height=14)
options(repr.plot.width=20, repr.plot.height=21)

navin_umap2 <- DimPlot(data_VA, reduction = "umap", group.by = c("Navin_Clusters", "Navin_annotations_simplified", 
                                                                 "Navin_Clusters2", 
                                                                 "Tumor_MHCI_label", "vasc_label"), 
                      ncol = 2,
                      cols = dittoColors()[1:12],
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap2

# save
# 1300, 1100 old
png('Output/03_Incorporate_Navin_Annotations_VA1/Navin_Annot_simplified_UMAP.png',
    width = 1400,
    height = 1650)
print(navin_umap2)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_66_1.png)
    


# How much of each seurat cluster is each cell type?


```R
# check ditto color palette...
dittoColors()[1:25]
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'#E69F00'</li><li>'#56B4E9'</li><li>'#009E73'</li><li>'#F0E442'</li><li>'#0072B2'</li><li>'#D55E00'</li><li>'#CC79A7'</li><li>'#666666'</li><li>'#AD7700'</li><li>'#1C91D4'</li><li>'#007756'</li><li>'#D5C711'</li><li>'#005685'</li><li>'#A04700'</li><li>'#B14380'</li><li>'#4D4D4D'</li><li>'#FFBE2D'</li><li>'#80C7EF'</li><li>'#00F6B3'</li><li>'#F4EB71'</li><li>'#06A5FF'</li><li>'#FF8320'</li><li>'#D99BBD'</li><li>'#8C8C8C'</li><li>'#FFCB57'</li></ol>




```R
show_col(dittoColors()[1:25])
```


    
![png](output_69_0.png)
    



```R
options(repr.plot.width=17, repr.plot.height=8)
```


```R
navin_annot_percent_in_seurat_clusters_plot_VA <- dittoBarPlot(
    object = data_VA,
    var = "Navin_annotations",
    group.by = "Navin_Clusters") + 
    theme(text = element_text(size=20))

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/Navin_Annot_Percent_in_Seurat_Clusters_VA1.png",
    width = 1100,
    height = 700)
print(navin_annot_percent_in_seurat_clusters_plot_VA)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/Navin_Annot_Percent_in_Seurat_Clusters_VA1.pdf",
    width = 11,
    height = 7)
print(navin_annot_percent_in_seurat_clusters_plot_VA)
dev.off()

navin_annot_percent_in_seurat_clusters_plot_VA
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_71_2.png)
    



```R
navin_annot_percent_in_seurat_clusters_plot_VA_2 <- dittoBarPlot(
    object = data_VA,
    var = "Navin_annotations_simplified",
    group.by = "Navin_Clusters") + 
    theme(text = element_text(size=20))

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/Navin_Annot_simplified_Percent_in_Seurat_Clusters_VA1.png",
    width = 1100,
    height = 700)
print(navin_annot_percent_in_seurat_clusters_plot_VA_2)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/Navin_Annot_simplified_Percent_in_Seurat_Clusters_VA1.pdf",
    width = 11,
    height = 7)
print(navin_annot_percent_in_seurat_clusters_plot_VA_2)
dev.off()

navin_annot_percent_in_seurat_clusters_plot_VA_2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_72_2.png)
    


# vs. Navin Cluster Annotations


```R
options(repr.plot.width=17, repr.plot.height=10)
```


```R
navin_annot_percent_in_seurat_clusters_plot_VA <- dittoBarPlot(
    object = data_VA,
    var = "Navin_annotations",
    group.by = "Navin_Clusters2") + 
    theme(text = element_text(size=20))

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/Navin_Annot_Percent_in_Navin_Annotated_Clusters_VA1.png",
    width = 1100,
    height = 750)
print(navin_annot_percent_in_seurat_clusters_plot_VA)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/Navin_Annot_Percent_in_Navin_Annotated_Clusters_VA1.pdf",
    width = 11,
    height = 7.5)
print(navin_annot_percent_in_seurat_clusters_plot_VA)
dev.off()

navin_annot_percent_in_seurat_clusters_plot_VA
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_75_2.png)
    



```R
navin_annot_percent_in_seurat_clusters_plot_VA_2 <- dittoBarPlot(
    object = data_VA,
    var = "Navin_annotations_simplified",
    group.by = "Navin_Clusters2") + 
    theme(text = element_text(size=20))

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/Navin_Annot_simplified_Percent_in_Navin_Annotated_Clusters_VA1.png",
    width = 1100,
    height = 750)
print(navin_annot_percent_in_seurat_clusters_plot_VA_2)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/Navin_Annot_simplified_Percent_in_Navin_Annotated_Clusters_VA1.pdf",
    width = 11,
    height = 7.5)
print(navin_annot_percent_in_seurat_clusters_plot_VA_2)
dev.off()

navin_annot_percent_in_seurat_clusters_plot_VA_2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_76_2.png)
    


# how much of each cluster in each tissue slice?

Navin's seurat cluster labels based on top DE genes: 
still needed for this new analysis

## Takeaways:
- Cluster 10 is the only one with Respiratory Epithelium, and Smooth muscle. (NORMAL LUNG)
- Tumor well represented in clusters 1, 11, 4, 5, 6, 9 (TUMOR)
- Lots of NaNs in the annotations
- Only VA is annotated


```R
tissue_slices_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2,
    var = "Tissue_Slice",
    group.by = "Navin_Clusters") + 
    theme(text = element_text(size=20))

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/Tissue_Slices_in_Seurat_Clusters_AB1_A2.png",
    width = 1100,
    height = 700)
print(tissue_slices_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/Tissue_Slices_in_Seurat_Clusters_AB1_A2.pdf",
    width = 11,
    height = 7)
print(tissue_slices_in_seurat_clusters_plot_AB1_A2)
dev.off()

tissue_slices_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_80_2.png)
    



```R
# by navin cluster annotation

tissue_slices_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2,
    var = "Tissue_Slice",
    group.by = "Navin_Clusters2") + 
    theme(text = element_text(size=20))

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/Tissue_Slices_in_Navin_Cluster_Annotations2_AB1_A2.png",
    width = 1100,
    height = 700)
print(tissue_slices_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/Tissue_Slices_in_Navin_Cluster_Annotations2_AB1_A2.pdf",
    width = 11,
    height = 7)
print(tissue_slices_in_seurat_clusters_plot_AB1_A2)
dev.off()

tissue_slices_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_81_2.png)
    



```R
colnames(AB1_A2@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Barcode'</li><li>'Navin_annotations'</li><li>'Tissue_Slice'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li></ol>




```R
Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2,
    var = "Tumor_MHCI_label",
    group.by = "Navin_Clusters") + 
    theme(text = element_text(size=20))

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/Tumor_MHCI_Label_in_Seurat_Clusters_AB1_A2.png",
    width = 1100,
    height = 700)
print(Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/Tumor_MHCI_Label_in_Seurat_Clusters_AB1_A2.pdf",
    width = 11,
    height = 7)
print(Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_83_2.png)
    



```R
vasc_Label_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2,
    var = "vasc_label",
    group.by = "Navin_Clusters") + 
    theme(text = element_text(size=20))

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/vasc_Label_in_Seurat_Clusters_AB1_A2.png",
    width = 1100,
    height = 700)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/vasc_Label_in_Seurat_Clusters_AB1_A2.pdf",
    width = 11,
    height = 7)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

vasc_Label_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_84_2.png)
    


# by Navin Cluster Annotation


```R
Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2,
    var = "Tumor_MHCI_label",
    group.by = "Navin_Clusters2") + 
    theme(text = element_text(size=20))

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/Tumor_MHCI_Label_in_Navin_Cluster_Annotation_AB1_A2.png",
    width = 1100,
    height = 700)
print(Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/Tumor_MHCI_Label_in_Navin_Cluster_Annotation_AB1_A2.pdf",
    width = 11,
    height = 7)
print(Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_86_2.png)
    



```R
vasc_Label_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2,
    var = "vasc_label",
    group.by = "Navin_Clusters2") + 
    theme(text = element_text(size=20))

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/vasc_Label_in_Navin_Cluster_Annotation_AB1_A2.png",
    width = 1100,
    height = 700)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/vasc_Label_in_Navin_Cluster_Annotation_AB1_A2.pdf",
    width = 11,
    height = 7)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

vasc_Label_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_87_2.png)
    



```R
# define cluster colors based on dittoColors...
ditto_colors <- c('0_Lymphocytes_B_Cells' = dittoColors()[1], 
                  '1_SCLC' = dittoColors()[2], 
                  '2_SCLC_TSI' = dittoColors()[3], 
                  '3_SCLC_Stroma' = dittoColors()[4], 
                  '4_SCLC' = dittoColors()[5], 
                  '5_SCLC' = dittoColors()[6], 
                  '6_SCLC_Effector_Immune' = dittoColors()[7], 
                  '7_Normal_Lung_Stroma' = dittoColors()[8],
                  '8_Fibroblasts' = dittoColors()[9], 
                  '9_SCLC_Apoptotic' = dittoColors()[10], 
                  '10_Respiratory_Epithelium_Lung' = dittoColors()[11], 
                  '11_SCLC_Necrotic' = dittoColors()[19])
```


```R
options(repr.plot.width=28, repr.plot.height=8)

Idents(AB1_A2) <- "Navin_Clusters"

plot0 <- SpatialDimPlot(AB1_A2, cols = ditto_colors, pt.size.factor=0)
plot0 

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/HnEs.png",
    width = 2400,
    height = 700)
print(plot0)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/HnEs.pdf",
    width = 24,
    height = 7)
print(plot0)
dev.off()

plot1 <- SpatialDimPlot(AB1_A2, cols = ditto_colors, pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=8)))
plot1 

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/DimPlot.png",
    width = 2400,
    height = 700)
print(plot1)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/DimPlot.pdf",
    width = 24,
    height = 7)
print(plot1)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_89_2.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_89_5.png)
    



```R
# plot in one image
options(repr.plot.width=28, repr.plot.height=16)

Idents(AB1_A2) <- "Navin_Clusters"

plot0 <- SpatialDimPlot(AB1_A2, cols = ditto_colors, pt.size.factor=0)

plot1 <- SpatialDimPlot(AB1_A2, cols = ditto_colors, pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=8)))

plot0 / plot1

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/HnEs_and_DimPlot.png",
    width = 2400,
    height = 1400)
print(plot0 / plot1)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/HnEs_and_DimPlot.pdf",
    width = 24,
    height = 14)
print(plot0 / plot1)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_90_2.png)
    


# Look at collapsed cluster labels


```R
# define cluster colors based on dittoColors...
ditto_colors2 <- c('0_Lymphocytes_B_Cells' = dittoColors()[1], 
                  '1_4_5_11_SCLC' = dittoColors()[2], 
                  '2_SCLC_TSI' = dittoColors()[3], 
                  '3_SCLC_Stroma' = dittoColors()[4], 
                  '6_SCLC_Effector_Immune' = dittoColors()[7], 
                  '7_10_Normal_Lung' = dittoColors()[8],
                  '8_Fibroblasts' = dittoColors()[9], 
                  '9_SCLC_Apoptotic' = dittoColors()[10])
```


```R
options(repr.plot.width=28, repr.plot.height=8)

Idents(AB1_A2) <- "Navin_Clusters2"

plot1 <- SpatialDimPlot(AB1_A2, cols = ditto_colors2, pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=8)))
plot1 

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/DimPlot2.png",
    width = 2400,
    height = 700)
print(plot1)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/DimPlot2.pdf",
    width = 24,
    height = 7)
print(plot1)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_93_2.png)
    



```R
# plot in one image
options(repr.plot.width=28, repr.plot.height=16)

Idents(AB1_A2) <- "Navin_Clusters2"

plot0 <- SpatialDimPlot(AB1_A2, cols = ditto_colors2, pt.size.factor=0)

plot1 <- SpatialDimPlot(AB1_A2, cols = ditto_colors2, pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=8)))

plot0 / plot1

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/HnEs_and_DimPlot2.png",
    width = 2400,
    height = 1400)
print(plot0 / plot1)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/HnEs_and_DimPlot2.pdf",
    width = 24,
    height = 14)
print(plot0 / plot1)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_94_2.png)
    



```R
sort(unique(data_VA@meta.data$Navin_annotations))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>''</li><li>'Alveoli'</li><li>'Endothelial cells'</li><li>'Fibrosis'</li><li>'Immune cells_lymphocytes'</li><li>'Immune cells_myeloid cells'</li><li>'Immune cells+vasculature'</li><li>'Respiratory epithelium'</li><li>'SCLC'</li><li>'SCLC_TSI'</li><li>'SCLC+fibrosis'</li><li>'SCLC+immune cells'</li><li>'SCLC+vasculature'</li><li>'Smooth muscle'</li><li>'Tumor necrosis'</li></ol>




```R
length(unique(data_VA@meta.data$Navin_annotations))
```


15



```R
ditto_colors_navin <- setNames(as.list(dittoColors()[1:15]),
                               as.list(sort(unique(data_VA@meta.data$Navin_annotations))))
```


```R
ditto_colors_navin
```


<dl>
	<dt>[[1]]</dt>
		<dd>'#E69F00'</dd>
	<dt>$Alveoli</dt>
		<dd>'#56B4E9'</dd>
	<dt>$`Endothelial cells`</dt>
		<dd>'#009E73'</dd>
	<dt>$Fibrosis</dt>
		<dd>'#F0E442'</dd>
	<dt>$`Immune cells_lymphocytes`</dt>
		<dd>'#0072B2'</dd>
	<dt>$`Immune cells_myeloid cells`</dt>
		<dd>'#D55E00'</dd>
	<dt>$`Immune cells+vasculature`</dt>
		<dd>'#CC79A7'</dd>
	<dt>$`Respiratory epithelium`</dt>
		<dd>'#666666'</dd>
	<dt>$SCLC</dt>
		<dd>'#AD7700'</dd>
	<dt>$SCLC_TSI</dt>
		<dd>'#1C91D4'</dd>
	<dt>$`SCLC+fibrosis`</dt>
		<dd>'#007756'</dd>
	<dt>$`SCLC+immune cells`</dt>
		<dd>'#D5C711'</dd>
	<dt>$`SCLC+vasculature`</dt>
		<dd>'#005685'</dd>
	<dt>$`Smooth muscle`</dt>
		<dd>'#A04700'</dd>
	<dt>$`Tumor necrosis`</dt>
		<dd>'#B14380'</dd>
</dl>




```R
# spatial dimplot colored by navin annotations
Idents(data_VA) <- "Navin_annotations"

options(repr.plot.width=15, repr.plot.height=8)

plot2 <- SpatialDimPlot(data_VA, cols = ditto_colors_navin, pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=8)))
plot3 <- plot2[1]
plot2

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/DimPlot_Navin_Annot.png",
    width = 1500,
    height = 700)
print(plot2)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/DimPlot_Navin_Annot.pdf",
    width = 15,
    height = 7)
print(plot2)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_99_2.png)
    


# plot simplified annotations


```R
length(unique(data_VA@meta.data$Navin_annotations_simplified))
```


10



```R
ditto_colors_navin_simplified <- setNames(as.list(dittoColors()[1:10]),
                               as.list(sort(unique(data_VA@meta.data$Navin_annotations_simplified))))
```


```R
ditto_colors_navin_simplified
```


<dl>
	<dt>[[1]]</dt>
		<dd>'#E69F00'</dd>
	<dt>$Alveoli</dt>
		<dd>'#56B4E9'</dd>
	<dt>$`Endothelial cells`</dt>
		<dd>'#009E73'</dd>
	<dt>$Fibrosis</dt>
		<dd>'#F0E442'</dd>
	<dt>$Immune</dt>
		<dd>'#0072B2'</dd>
	<dt>$`Respiratory epithelium`</dt>
		<dd>'#D55E00'</dd>
	<dt>$SCLC</dt>
		<dd>'#CC79A7'</dd>
	<dt>$`SCLC+Immune`</dt>
		<dd>'#666666'</dd>
	<dt>$`Smooth muscle`</dt>
		<dd>'#AD7700'</dd>
	<dt>$`Tumor necrosis`</dt>
		<dd>'#1C91D4'</dd>
</dl>




```R
# spatial dimplot colored by navin annotations simplified
Idents(data_VA) <- "Navin_annotations_simplified"

options(repr.plot.width=15, repr.plot.height=8)

plot2 <- SpatialDimPlot(data_VA, cols = ditto_colors_navin_simplified, pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=8)))
plot3 <- plot2[1]
plot2

png(file = "Output/03_Incorporate_Navin_Annotations_VA1/DimPlot_Navin_Annot_simplified.png",
    width = 1500,
    height = 700)
print(plot2)
dev.off()

pdf(file = "Output/03_Incorporate_Navin_Annotations_VA1/DimPlot_Navin_Annot_simplified.pdf",
    width = 15,
    height = 7)
print(plot2)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_104_2.png)
    


# Future work: repeat analyses from previous submission
    - Measure cell type enrichment (GSVA, or STdeconvolve, or something similar)
    - Measure enrichment of various pathways
        - Does enrichment of pathways correspond to Tumor High MHCI High spots?
        - Does enrichment of pathways correspond to vasc_high spots?


```R

```


```R

```


```R

```
