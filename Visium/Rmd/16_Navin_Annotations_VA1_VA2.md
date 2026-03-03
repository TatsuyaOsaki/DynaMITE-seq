## July 2025
## Round 1 REVISIONS for SCLC 10X Visium
## I annotated spots for vasc, Tumor, and MHCI. 
## Add Navin's Annotations for VA1 and see which clusters and groups they represent. 
## Update July 2025, Navin also annotated VA2. Let's add those in and make some plots for the final paper figures.


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
library(RColorBrewer)
library(GSVA)
library(ggplot2)
library(scales)
#library('scCustomize')

options(repr.matrix.max.cols=1000, repr.matrix.max.rows=100)
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
AB1_A2 <- readRDS('Processing/AB1_A2_Annotated_GSVA_by_Tissue_Slice_06_Label_MN4_20250425.rds')
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
AB1_A2@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 61</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td></tr>
</tbody>
</table>




```R
table(AB1_A2@meta.data$Tissue_Slice, AB1_A2@meta.data$Navin_annotations, useNA = "always")
```


          
                Alveoli Endothelial cells Fibrosis Immune cells_lymphocytes
      VA1   413      12                52      155                      192
      VA2     0       0                 0        0                        0
      VB1     0       0                 0        0                        0
      <NA>    0       0                 0        0                        0
          
           Immune cells_myeloid cells Immune cells+vasculature
      VA1                          29                       68
      VA2                           0                        0
      VB1                           0                        0
      <NA>                          0                        0
          
           Respiratory epithelium SCLC SCLC_TSI SCLC+fibrosis SCLC+immune cells
      VA1                      30 1233      112            76                69
      VA2                       0    0        0             0                 0
      VB1                       0    0        0             0                 0
      <NA>                      0    0        0             0                 0
          
           SCLC+vasculature Smooth muscle Tumor necrosis <NA>
      VA1               309            20             36    0
      VA2                 0             0              0 4316
      VB1                 0             0              0 2589
      <NA>                0             0              0    0


# Load Navin's manual annotations


```R
navin_annot_VA <- read.csv('/immuno/ian/Projects/2023_02_SCLC_10X_Visium/Annotated_cloupe/Navin_annotations_VA.csv')
navin_annot_VA2 <- read.csv('/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/Navin_annotations_VA2.csv')
```


```R
# Add "VA1_" to Barcode to match Seurat's barcode naming
navin_annot_VA$Barcode <- lapply(navin_annot_VA$Barcode, function(x) paste0('VA1_', trimws(x)))
navin_annot_VA2$Barcode <- lapply(navin_annot_VA2$Barcode, function(x) paste0('VA2_', trimws(x)))
```


```R
# convert to character
navin_annot_VA$Barcode <- as.character(navin_annot_VA$Barcode)
navin_annot_VA2$Barcode <- as.character(navin_annot_VA2$Barcode)
```


```R
# rename VA2 colname
navin_annot_VA2$Navin_annotations <- navin_annot_VA2$Navin.s.annotations_VA2
navin_annot_VA2 <- navin_annot_VA2 %>% select(!c(`Navin.s.annotations_VA2`))
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




```R
head(navin_annot_VA2,3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 2</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Navin_annotations</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>VA2_AAACAAGTATCTCCCA-1</td><td>SCLC                    </td></tr>
	<tr><th scope=row>2</th><td>VA2_AAACACCAATAACTGC-1</td><td>Immune cells_lymphocytes</td></tr>
	<tr><th scope=row>3</th><td>VA2_AAACAGAGCGACTCCT-1</td><td>SCLC                    </td></tr>
</tbody>
</table>




```R
# concatenate the annotation dataframes
dim(navin_annot_VA)
dim(navin_annot_VA2)
temp <- rbind(navin_annot_VA, navin_annot_VA2)
navin_annot_VA <- temp
dim(navin_annot_VA)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>2806</li><li>2</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>4316</li><li>2</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>7122</li><li>2</li></ol>




```R
# rename to manual annotation
navin_annot_VA$Manual_Annotation <- navin_annot_VA$Navin_annotations
navin_annot_VA <- navin_annot_VA %>% select(!c(Navin_annotations))
```

# Try to merge navin's annotations in with the metadata...


```R
dim(navin_annot_VA)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>7122</li><li>2</li></ol>




```R
dim(AB1_A2[[]])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>61</li></ol>




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
<caption>A data.frame: 5 × 61</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.15846527</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.11080007</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.53207656</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.20627210</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.51946615</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.11225417</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1</td><td>VA1_AAACACCAATAACTGC-1</td><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td>4 </td><td>4 </td><td>VA1</td><td>0.907884</td><td>vasc_low </td><td> 9</td><td>Tumor_low </td><td>Tumor_low </td><td>1</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC+vasculature      </td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td>-0.29499671</td><td>-0.1592542</td><td>-0.6246499</td><td>-0.70918797</td><td>-0.3126876</td><td>-0.128532387</td><td> 0.54796712</td><td>-0.28738253</td><td> 0.05564459</td><td>-0.4488276</td><td>-0.3094928</td><td>-0.48102913</td><td>-0.23406975</td><td>-0.03504189</td><td>-0.1584708</td><td>-0.2611306</td><td>-0.3051610</td><td> 0.03743911</td><td> 0.3381445</td><td>-0.4111955</td><td>-0.3300223</td><td>-0.14181196</td><td>-0.19025652</td><td>-0.7206696</td><td>-0.21902693</td><td>-0.45491178</td><td>-0.3728931</td><td>-0.4603460</td><td>-0.23406975</td><td> 0.07776597</td><td>Low </td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td></tr>
	<tr><th scope=row>VA1_AAACAGCTTTCAGAAG-1</th><td>VA1</td><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>5132</td><td>3219</td><td>round1</td><td>4797</td><td>3218</td><td>4 </td><td>4 </td><td>VA1</td><td>3.751279</td><td>vasc_high</td><td>14</td><td>Tumor_low </td><td>Tumor_low </td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low          </td><td>vasc_high</td><td>SCLC+vasculature      </td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td> 0.14865678</td><td> 0.6701079</td><td> 0.1810290</td><td>-0.68264193</td><td> 0.5112075</td><td> 0.351200349</td><td> 0.29886669</td><td> 0.53696284</td><td>-0.22259553</td><td> 0.2016089</td><td>-0.2723653</td><td> 0.26234173</td><td>-0.24403536</td><td>-0.02474443</td><td>-0.0993218</td><td> 0.5413331</td><td> 0.7253772</td><td>-0.13408735</td><td>-0.4164351</td><td> 0.1162332</td><td> 0.1290473</td><td> 0.15984342</td><td> 0.43378600</td><td> 0.6304333</td><td> 0.59860816</td><td>-0.09760237</td><td> 0.8719730</td><td> 0.9183918</td><td>-0.24403536</td><td>-0.15730040</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
</tbody>
</table>




```R
# how many barcodes actually overlap?
overlap_list <- Reduce(intersect,list(AB1_A2[[]]$Barcode,navin_annot_VA$Barcode))
print(length(overlap_list))
```

    [1] 7122



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
<ol class=list-inline><li>9711</li><li>62</li></ol>




```R
AB1_A2[[]][0:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 62</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Manual_Annotation</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.15846527</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.11080007</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Respiratory epithelium</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.53207656</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.20627210</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.51946615</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.11225417</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1</td><td>VA1_AAACACCAATAACTGC-1</td><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td>4 </td><td>4 </td><td>VA1</td><td>0.907884</td><td>vasc_low </td><td> 9</td><td>Tumor_low </td><td>Tumor_low </td><td>1</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC+vasculature      </td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td>-0.29499671</td><td>-0.1592542</td><td>-0.6246499</td><td>-0.70918797</td><td>-0.3126876</td><td>-0.128532387</td><td> 0.54796712</td><td>-0.28738253</td><td> 0.05564459</td><td>-0.4488276</td><td>-0.3094928</td><td>-0.48102913</td><td>-0.23406975</td><td>-0.03504189</td><td>-0.1584708</td><td>-0.2611306</td><td>-0.3051610</td><td> 0.03743911</td><td> 0.3381445</td><td>-0.4111955</td><td>-0.3300223</td><td>-0.14181196</td><td>-0.19025652</td><td>-0.7206696</td><td>-0.21902693</td><td>-0.45491178</td><td>-0.3728931</td><td>-0.4603460</td><td>-0.23406975</td><td> 0.07776597</td><td>Low </td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC+vasculature      </td></tr>
	<tr><th scope=row>VA1_AAACAGCTTTCAGAAG-1</th><td>VA1</td><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>5132</td><td>3219</td><td>round1</td><td>4797</td><td>3218</td><td>4 </td><td>4 </td><td>VA1</td><td>3.751279</td><td>vasc_high</td><td>14</td><td>Tumor_low </td><td>Tumor_low </td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low          </td><td>vasc_high</td><td>SCLC+vasculature      </td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td> 0.14865678</td><td> 0.6701079</td><td> 0.1810290</td><td>-0.68264193</td><td> 0.5112075</td><td> 0.351200349</td><td> 0.29886669</td><td> 0.53696284</td><td>-0.22259553</td><td> 0.2016089</td><td>-0.2723653</td><td> 0.26234173</td><td>-0.24403536</td><td>-0.02474443</td><td>-0.0993218</td><td> 0.5413331</td><td> 0.7253772</td><td>-0.13408735</td><td>-0.4164351</td><td> 0.1162332</td><td> 0.1290473</td><td> 0.15984342</td><td> 0.43378600</td><td> 0.6304333</td><td> 0.59860816</td><td>-0.09760237</td><td> 0.8719730</td><td> 0.9183918</td><td>-0.24403536</td><td>-0.15730040</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>SCLC+vasculature      </td></tr>
</tbody>
</table>




```R
table(AB1_A2@meta.data$Manual_Annotation, useNA = "always")
```


    
                                                  Alveoli 
                          1033                        609 
             Endothelial cells                   Fibrosis 
                            70                        505 
      Immune cells_lymphocytes Immune cells_myeloid cells 
                           587                         41 
      Immune cells+vasculature     Respiratory epithelium 
                           107                         47 
                          SCLC                   SCLC_TSI 
                          2504                        211 
                 SCLC+fibrosis          SCLC+immune cells 
                           399                        173 
              SCLC+vasculature              Smooth muscle 
                           507                         53 
                Tumor necrosis                       <NA> 
                           276                       2589 



```R
colSums(!is.na(AB1_A2@meta.data))
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>Tissue_Slice</dt><dd>9711</dd><dt>Barcode</dt><dd>9711</dd><dt>orig.ident</dt><dd>9711</dd><dt>nCount_Spatial</dt><dd>9711</dd><dt>nFeature_Spatial</dt><dd>9711</dd><dt>visium_round</dt><dd>9711</dd><dt>nCount_SCT</dt><dd>9711</dd><dt>nFeature_SCT</dt><dd>9711</dd><dt>SCT_snn_res.0.8</dt><dd>9711</dd><dt>seurat_clusters</dt><dd>9711</dd><dt>slice_ident</dt><dd>9711</dd><dt>vasc_summed_expr</dt><dd>9711</dd><dt>vasc_label</dt><dd>9711</dd><dt>Tumor_score</dt><dd>9711</dd><dt>Tumor_label1</dt><dd>9711</dd><dt>Tumor_label</dt><dd>9711</dd><dt>MHCI_score</dt><dd>9711</dd><dt>MHCI_label1</dt><dd>9711</dd><dt>MHCI_label</dt><dd>9711</dd><dt>Tumor_MHCI_label</dt><dd>9711</dd><dt>Tumor_MHCI_label2</dt><dd>9711</dd><dt>vasc_label1</dt><dd>9711</dd><dt>Navin_annotations</dt><dd>2806</dd><dt>Navin_annotations_simplified</dt><dd>2806</dd><dt>Navin_Clusters1</dt><dd>9711</dd><dt>Navin_Clusters</dt><dd>9711</dd><dt>Navin_Clusters2</dt><dd>9711</dd><dt>All_Immune</dt><dd>9711</dd><dt>Angiogenesis</dt><dd>9711</dd><dt>AyersIFNG</dt><dd>9711</dd><dt>B_Cell</dt><dd>9711</dd><dt>CD8_T_Cell</dt><dd>9711</dd><dt>DC</dt><dd>9711</dd><dt>Endothelial_Activation</dt><dd>9711</dd><dt>Endothelial_Cell</dt><dd>9711</dd><dt>Endothelial_Chemokines</dt><dd>9711</dd><dt>Epithelial_cell</dt><dd>9711</dd><dt>Exhaustion</dt><dd>9711</dd><dt>Fibroblast</dt><dd>9711</dd><dt>GOBP_LEUKOCYTE_ADHESION_TO_VASC</dt><dd>9711</dd><dt>GOBP_LEUKOCYTE_MIGRATION</dt><dd>9711</dd><dt>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</dt><dd>9711</dd><dt>ICB_Targets</dt><dd>9711</dd><dt>LEC</dt><dd>9711</dd><dt>M1_Macrophage</dt><dd>9711</dd><dt>M1_Macrophage_old</dt><dd>9711</dd><dt>M2_Macrophage</dt><dd>9711</dd><dt>M2_Macrophage_old</dt><dd>9711</dd><dt>MN4_EC_Phenotype</dt><dd>9711</dd><dt>MN4_EC_Phenotype_Top30</dt><dd>9711</dd><dt>Macrophage</dt><dd>9711</dd><dt>NK_Cell</dt><dd>9711</dd><dt>Proliferation</dt><dd>9711</dd><dt>T_Cell</dt><dd>9711</dd><dt>T_Reg</dt><dd>9711</dd><dt>Upregulated_by_2_3_CGAMP</dt><dd>9711</dd><dt>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</dt><dd>9711</dd><dt>MN4_EC_Phenotype_Label</dt><dd>9711</dd><dt>MN4_EC_Percentile_Label</dt><dd>9711</dd><dt>vasc_MN4_label</dt><dd>9711</dd><dt>vasc_MN4_percentile_label</dt><dd>9711</dd><dt>Manual_Annotation</dt><dd>7122</dd></dl>




```R
AB1_A2@meta.data[0:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 62</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Manual_Annotation</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.15846527</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.11080007</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Respiratory epithelium</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.53207656</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.20627210</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.51946615</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.11225417</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1</td><td>VA1_AAACACCAATAACTGC-1</td><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td>4 </td><td>4 </td><td>VA1</td><td>0.907884</td><td>vasc_low </td><td> 9</td><td>Tumor_low </td><td>Tumor_low </td><td>1</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC+vasculature      </td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td>-0.29499671</td><td>-0.1592542</td><td>-0.6246499</td><td>-0.70918797</td><td>-0.3126876</td><td>-0.128532387</td><td> 0.54796712</td><td>-0.28738253</td><td> 0.05564459</td><td>-0.4488276</td><td>-0.3094928</td><td>-0.48102913</td><td>-0.23406975</td><td>-0.03504189</td><td>-0.1584708</td><td>-0.2611306</td><td>-0.3051610</td><td> 0.03743911</td><td> 0.3381445</td><td>-0.4111955</td><td>-0.3300223</td><td>-0.14181196</td><td>-0.19025652</td><td>-0.7206696</td><td>-0.21902693</td><td>-0.45491178</td><td>-0.3728931</td><td>-0.4603460</td><td>-0.23406975</td><td> 0.07776597</td><td>Low </td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC+vasculature      </td></tr>
	<tr><th scope=row>VA1_AAACAGCTTTCAGAAG-1</th><td>VA1</td><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>5132</td><td>3219</td><td>round1</td><td>4797</td><td>3218</td><td>4 </td><td>4 </td><td>VA1</td><td>3.751279</td><td>vasc_high</td><td>14</td><td>Tumor_low </td><td>Tumor_low </td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low          </td><td>vasc_high</td><td>SCLC+vasculature      </td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td> 0.14865678</td><td> 0.6701079</td><td> 0.1810290</td><td>-0.68264193</td><td> 0.5112075</td><td> 0.351200349</td><td> 0.29886669</td><td> 0.53696284</td><td>-0.22259553</td><td> 0.2016089</td><td>-0.2723653</td><td> 0.26234173</td><td>-0.24403536</td><td>-0.02474443</td><td>-0.0993218</td><td> 0.5413331</td><td> 0.7253772</td><td>-0.13408735</td><td>-0.4164351</td><td> 0.1162332</td><td> 0.1290473</td><td> 0.15984342</td><td> 0.43378600</td><td> 0.6304333</td><td> 0.59860816</td><td>-0.09760237</td><td> 0.8719730</td><td> 0.9183918</td><td>-0.24403536</td><td>-0.15730040</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>SCLC+vasculature      </td></tr>
</tbody>
</table>




```R
print(unique(AB1_A2@meta.data$Tissue_Slice))
```

    [1] "VA1" "VA2" "VB1"


# Looks like navin's annotations have successfully been integrated with the metadata. Now we can try to plot some things...

# Since we only have navin annotations for VA1 and VA2, let's subset data to just look at VA1 and VA2 for now...


```R
# remove troublesome graph 
AB1_A2@graphs <- list() 

data_VA <- subset(AB1_A2, subset = Tissue_Slice %in% c('VA1', 'VA2'))
print(dim(AB1_A2))
print(dim(data_VA))
```

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
    [1] 17943  7122



```R
print(dim(AB1_A2[[]]))
print(dim(data_VA[[]]))
```

    [1] 9711   62
    [1] 7122   62



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
    35796 features across 7122 samples within 2 assays 
    Active assay: Spatial (17943 features, 2000 variable features)
     3 layers present: data, counts, scale.data
     1 other assay present: SCT
     3 dimensional reductions calculated: pca, integrated.cca, umap
     2 spatial fields of view present: slice1 slice1.2



```R
colnames(AB1_A2[[]])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li><li>'Manual_Annotation'</li></ol>




```R
colnames(data_VA[[]])
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li><li>'Manual_Annotation'</li></ol>




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


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'slice1'</li><li>'slice1.2'</li></ol>




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
<ol class=list-inline><li>17943</li><li>7122</li></ol>



# Let's simplify and group some of Navin's annotations together


```R
print(unique(data_VA[[]]$Manual_Annotation))
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
# Relabel empty string so I can recode it
data_VA$Manual_Annotation[data_VA$Manual_Annotation == ""] <- "Unannotated"
AB1_A2$Manual_Annotation[AB1_A2$Manual_Annotation == ""] <- "Unannotated"
```
# old
data_VA[[]]$Manual_Annotation_simplified <- recode(data_VA[[]]$Manual_Annotation, 
                                                   `Unannotated` = "Unannotated",
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

```R
data_VA[[]]$Manual_Annotation_simplified <- recode(data_VA[[]]$Manual_Annotation, 
                                                   `Unannotated` = "Unannotated",
                                                   `Respiratory epithelium` = "Respiratory epithelium", 
                                                   `SCLC+vasculature` = "SCLC", 
                                                   `SCLC` = "SCLC", 
                                                   `SCLC_TSI` = "SCLC", 
                                                   `SCLC+fibrosis` = "SCLC", 
                                                   `SCLC+immune cells` = "SCLC", 
                                                   `Immune cells_lymphocytes` = "Immune", 
                                                   `Immune cells+vasculature` = "Immune", 
                                                   `Immune cells_myeloid cells` = "Immune", 
                                                   `Smooth muscle` = "Respiratory epithelium", 
                                                   `Endothelial cells` = "Endothelial cells", 
                                                   `Tumor necrosis` = "SCLC", 
                                                   `Alveoli` = "Respiratory epithelium", 
                                                   `Fibrosis` = "Fibrosis")
```


```R
print(table(data_VA[[]]$Manual_Annotation, useNA = "always"))
print(table(data_VA[[]]$Manual_Annotation_simplified, useNA = "always"))
```

    
                       Alveoli          Endothelial cells 
                           609                         70 
                      Fibrosis   Immune cells_lymphocytes 
                           505                        587 
    Immune cells_myeloid cells   Immune cells+vasculature 
                            41                        107 
        Respiratory epithelium                       SCLC 
                            47                       2504 
                      SCLC_TSI              SCLC+fibrosis 
                           211                        399 
             SCLC+immune cells           SCLC+vasculature 
                           173                        507 
                 Smooth muscle             Tumor necrosis 
                            53                        276 
                   Unannotated                       <NA> 
                          1033                          0 
    
         Endothelial cells               Fibrosis                 Immune 
                        70                    505                    735 
    Respiratory epithelium                   SCLC            Unannotated 
                       709                   4070                   1033 
                      <NA> 
                         0 



```R
print(dim(data_VA[[]]))
```

    [1] 7122   63


# Create as minimal annotation as possible


```R
data_VA[[]]$Manual_Annotation_minimal <- recode(data_VA[[]]$Manual_Annotation, 
                                                   `Unannotated` = "Unannotated",
                                                   `Respiratory epithelium` = "Stroma/Other", 
                                                   `SCLC+vasculature` = "SCLC", 
                                                   `SCLC` = "SCLC", 
                                                   `SCLC_TSI` = "SCLC", 
                                                   `SCLC+fibrosis` = "SCLC", 
                                                   `SCLC+immune cells` = "SCLC", 
                                                   `Immune cells_lymphocytes` = "Immune", 
                                                   `Immune cells+vasculature` = "Immune", 
                                                   `Immune cells_myeloid cells` = "Immune", 
                                                   `Smooth muscle` = "Stroma/Other", 
                                                   `Endothelial cells` = "Stroma/Other", 
                                                   `Tumor necrosis` = "Stroma/Other", 
                                                   `Alveoli` = "Stroma/Other", 
                                                   `Fibrosis` = "Stroma/Other")
```


```R
print(table(data_VA[[]]$Manual_Annotation_minimal, useNA = "always"))
```

    
          Immune         SCLC Stroma/Other  Unannotated         <NA> 
             735         3794         1560         1033            0 



```R
print(dim(data_VA[[]]))
```

    [1] 7122   64



```R
unique(AB1_A2[[]]$Manual_Annotation)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Respiratory epithelium'</li><li>'SCLC'</li><li>'SCLC+vasculature'</li><li>'Fibrosis'</li><li>'Unannotated'</li><li>'Immune cells+vasculature'</li><li>'Immune cells_lymphocytes'</li><li>'SCLC_TSI'</li><li>'Immune cells_myeloid cells'</li><li>'SCLC+immune cells'</li><li>'Smooth muscle'</li><li>'Endothelial cells'</li><li>'Tumor necrosis'</li><li>'SCLC+fibrosis'</li><li>'Alveoli'</li><li>NA</li></ol>



# Also annotate AB1_A2


```R
AB1_A2[[]]$Manual_Annotation_simplified <- recode(AB1_A2[[]]$Manual_Annotation, 
                                                   `Unannotated` = "Unannotated",
                                                   `Respiratory epithelium` = "Respiratory epithelium", 
                                                   `SCLC+vasculature` = "SCLC", 
                                                   `SCLC` = "SCLC", 
                                                   `SCLC_TSI` = "SCLC", 
                                                   `SCLC+fibrosis` = "SCLC", 
                                                   `SCLC+immune cells` = "SCLC", 
                                                   `Immune cells_lymphocytes` = "Immune", 
                                                   `Immune cells+vasculature` = "Immune", 
                                                   `Immune cells_myeloid cells` = "Immune", 
                                                   `Smooth muscle` = "Respiratory epithelium", 
                                                   `Endothelial cells` = "Endothelial cells", 
                                                   `Tumor necrosis` = "SCLC", 
                                                   `Alveoli` = "Respiratory epithelium", 
                                                   `Fibrosis` = "Fibrosis")

AB1_A2[[]]$Manual_Annotation_minimal <- recode(AB1_A2[[]]$Manual_Annotation, 
                                                   `Unannotated` = "Unannotated",
                                                   `Respiratory epithelium` = "Stroma/Other", 
                                                   `SCLC+vasculature` = "SCLC", 
                                                   `SCLC` = "SCLC", 
                                                   `SCLC_TSI` = "SCLC", 
                                                   `SCLC+fibrosis` = "SCLC", 
                                                   `SCLC+immune cells` = "SCLC", 
                                                   `Immune cells_lymphocytes` = "Immune", 
                                                   `Immune cells+vasculature` = "Immune", 
                                                   `Immune cells_myeloid cells` = "Immune", 
                                                   `Smooth muscle` = "Stroma/Other", 
                                                   `Endothelial cells` = "Stroma/Other", 
                                                   `Tumor necrosis` = "Stroma/Other", 
                                                   `Alveoli` = "Stroma/Other", 
                                                   `Fibrosis` = "Stroma/Other")
```

# Navin also annotated seurat clusters using marker genes and looking at the spots overlapped on the H&E images using the loupe browser. Let's load those annotations in too. 


```R
navin_cluster_annots <- read.csv('Processing/Navin_Annotations_20250331.csv')
```


```R
navin_cluster_annots %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 2</caption>
<thead>
	<tr><th></th><th scope=col>Cluster</th><th scope=col>Navin_Annotation</th></tr>
	<tr><th></th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>0</td><td>Lymphocytes B Cells (and Plasma Cells)</td></tr>
	<tr><th scope=row>2</th><td>1</td><td>SCLC                                  </td></tr>
</tbody>
</table>




```R
navin_cluster_annots <- navin_cluster_annots %>% rename('seurat_clusters' = 'Cluster', 'Manual_Cluster_Annotation'='Navin_Annotation')
```


```R
navin_cluster_annots
```


<table class="dataframe">
<caption>A data.frame: 12 × 2</caption>
<thead>
	<tr><th scope=col>seurat_clusters</th><th scope=col>Manual_Cluster_Annotation</th></tr>
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
<caption>A data.frame: 3 × 64</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Manual_Annotation</th><th scope=col>Manual_Annotation_simplified</th><th scope=col>Manual_Annotation_minimal</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Stroma/Other</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC        </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC        </td></tr>
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
# minimal
AB1_A2[[]]$Labeled_Clusters <- recode(AB1_A2[[]]$seurat_clusters, 
                                                   "0" = "Lymphocytes_B_Cells",
                                                   "1" = "SCLC", 
                                                   "2" = "SCLC_TSI", 
                                                   "3" = "SCLC_Stroma", 
                                                   "4" = "SCLC", 
                                                   "5" = "SCLC", 
                                                   "6" = "SCLC_Effector_Immune", 
                                                   "7" = "Normal_Lung", 
                                                   "8" = "Fibroblasts", 
                                                   "9" = "SCLC_Apoptotic", 
                                                   "10" = "Normal_Lung", 
                                                   "11" = "SCLC")
```


```R
# most minimal
AB1_A2[[]]$Labeled_Clusters_minimal <- recode(AB1_A2[[]]$seurat_clusters, 
                                                   "0" = "Immune",
                                                   "1" = "SCLC", 
                                                   "2" = "SCLC", 
                                                   "3" = "SCLC", 
                                                   "4" = "SCLC", 
                                                   "5" = "SCLC", 
                                                   "6" = "SCLC", 
                                                   "7" = "Stroma/Other", 
                                                   "8" = "Stroma/Other", 
                                                   "9" = "SCLC", 
                                                   "10" = "Stroma/Other", 
                                                   "11" = "SCLC")
```


```R
AB1_A2@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 66</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Manual_Annotation</th><th scope=col>Manual_Annotation_simplified</th><th scope=col>Manual_Annotation_minimal</th><th scope=col>Labeled_Clusters</th><th scope=col>Labeled_Clusters_minimal</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Stroma/Other</td><td>Normal_Lung         </td><td>Stroma/Other</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC        </td><td>SCLC_Effector_Immune</td><td>SCLC        </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC        </td><td>SCLC_Effector_Immune</td><td>SCLC        </td></tr>
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

# minimal
data_VA[[]]$Labeled_Clusters <- recode(data_VA[[]]$seurat_clusters, 
                                                   "0" = "Lymphocytes_B_Cells",
                                                   "1" = "SCLC", 
                                                   "2" = "SCLC_TSI", 
                                                   "3" = "SCLC_Stroma", 
                                                   "4" = "SCLC", 
                                                   "5" = "SCLC", 
                                                   "6" = "SCLC_Effector_Immune", 
                                                   "7" = "Normal_Lung", 
                                                   "8" = "Fibroblasts", 
                                                   "9" = "SCLC_Apoptotic", 
                                                   "10" = "Normal_Lung", 
                                                   "11" = "SCLC")

# most minimal
data_VA[[]]$Labeled_Clusters_minimal <- recode(data_VA[[]]$seurat_clusters, 
                                                   "0" = "Immune",
                                                   "1" = "SCLC", 
                                                   "2" = "SCLC", 
                                                   "3" = "SCLC", 
                                                   "4" = "SCLC", 
                                                   "5" = "SCLC", 
                                                   "6" = "SCLC", 
                                                   "7" = "Stroma/Other", 
                                                   "8" = "Stroma/Other", 
                                                   "9" = "SCLC", 
                                                   "10" = "Stroma/Other", 
                                                   "11" = "SCLC")
```

# Save data now annotated with Navin's manual annotations


```R
saveRDS(AB1_A2, 'Processing/AB1_A2_Annotated_step2Navin_20250714.rds')
saveRDS(data_VA, 'Processing/A1A2_Annotated_step2Navin_20250714.rds')
```

# Dimplots


```R
#options(repr.plot.width=17, repr.plot.height=14)
options(repr.plot.width=20, repr.plot.height=27)

navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Navin_Clusters", "Manual_Annotation", 
                                                                "Labeled_Clusters", "Manual_Annotation_simplified", 
                                                                "Labeled_Clusters_minimal", "Manual_Annotation_minimal",
                                                                "Tumor_MHCI_label", "vasc_label",
                                                                "Tissue_Slice"), 
                      ncol = 2,
                      cols = dittoColors()[1:15],
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap

# save
# 1300, 1100 old
png('Output/16_Navin_Annotations_VA1_VA2/Mulitple_UMAPs.png',
    width = 1400,
    height = 2100)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_75_1.png)
    


# Let's plot a bunch individually for the paper


```R
ditto_colors_Labeled_Clusters <- c(`Normal_Lung` = '#56B4E9',
                                  `SCLC` = '#D55E00',
                                  `Fibroblasts` = '#005685',
                                  `SCLC_Apoptotic` = '#666666',
                                  `Lymphocytes_B_Cells` = '#F0E442',
                                  `SCLC_TSI` = '#009E73',
                                  `SCLC_Stroma` = '#005685',
                                  `SCLC_Effector_Immune` = '#E69F00')
```


```R
#group.by = c("Navin_Clusters", "Manual_Annotation", "Labeled_Clusters", "Manual_Annotation_simplified", 
#            "Tissue_Slice", "Manual_Annotation_minimal","Tumor_MHCI_label", "vasc_label"), 
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Labeled_Clusters"), 
                      ncol = 1,
                      #cols = dittoColors()[1:15],
                      cols = ditto_colors_Labeled_Clusters,
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Labeled_Clusters_UMAP.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()
pdf('Output/16_Navin_Annotations_VA1_VA2/Labeled_Clusters_UMAP.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_78_2.png)
    



```R
#group.by = c("Navin_Clusters", "Manual_Annotation", "Labeled_Clusters", "Manual_Annotation_simplified", 
#            "Tissue_Slice", "Manual_Annotation_minimal","Tumor_MHCI_label", "vasc_label"), 
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Labeled_Clusters_minimal"), 
                      ncol = 1,
                      cols = dittoColors()[1:15],
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Labeled_Clusters_minimal_UMAP.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()
pdf('Output/16_Navin_Annotations_VA1_VA2/Labeled_Clusters_minimal_UMAP.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_79_2.png)
    



```R
ditto_colors_navin_simplified <- c(`Respiratory epithelium` = '#56B4E9',
                                  `SCLC` = '#D55E00',
                                  `Fibrosis` = '#005685',
                                  `Unannotated` = '#666666',
                                  `Immune` = '#F0E442',
                                  `Endothelial cells` = '#009E73')
```


```R
#group.by = c("Navin_Clusters", "Manual_Annotation", "Labeled_Clusters", "Manual_Annotation_simplified", 
#            "Tissue_Slice", "Manual_Annotation_minimal","Tumor_MHCI_label", "vasc_label"), 
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Manual_Annotation_simplified"), 
                      ncol = 1,
                      #cols = dittoColors()[1:15],
                      cols=ditto_colors_navin_simplified,
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()

pdf('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_81_2.png)
    


# Plot each group on its own


```R
Idents(data_VA) <- "Manual_Annotation_simplified"
# get lists of cells in each group
SCLC_cells <- WhichCells(data_VA, idents = c("SCLC"))
Fibrosis_cells <- WhichCells(data_VA, idents = c("Fibrosis"))
Immune_cells <- WhichCells(data_VA, idents = c("Immune"))
Endothelial_cells <- WhichCells(data_VA, idents = c("Endothelial cells"))
Respiratory_epithelium_cells <- WhichCells(data_VA, idents = c("Respiratory epithelium"))
Unannotated_cells <- WhichCells(data_VA, idents = c("Unannotated"))
```


```R
ditto_colors_navin_simplified
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>Respiratory epithelium</dt><dd>'#56B4E9'</dd><dt>SCLC</dt><dd>'#D55E00'</dd><dt>Fibrosis</dt><dd>'#005685'</dd><dt>Unannotated</dt><dd>'#666666'</dd><dt>Immune</dt><dd>'#F0E442'</dd><dt>Endothelial cells</dt><dd>'#009E73'</dd></dl>


# example of how to do it
# DimPlot(integrated, label=T, group.by="Treat", cells.highlight= list(g1_treat, g1_untreat), cols.highlight = c("darkblue", "darkred"), cols= "grey")

```R
# SCLC
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Manual_Annotation_simplified"), 
                      ncol = 1,
                      #cols = dittoColors()[1:15],
                      #cols=ditto_colors_navin_simplified,
                      cells.highlight= list("Endothelial cells" = Endothelial_cells), 
                      cols.highlight = c('#009E73'), 
                      cols= "#d5d5d5",
                      sizes.highlight = 1.25,
                      pt.size=1.25) +
    theme(text = element_text(size=20)) +
    labs(title = "Endothelial cells")
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_Endothelial_focus.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()

pdf('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_Endothelial_focus.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_86_2.png)
    



```R
# SCLC
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Manual_Annotation_simplified"), 
                      ncol = 1,
                      #cols = dittoColors()[1:15],
                      #cols=ditto_colors_navin_simplified,
                      cells.highlight= list("Immune" = Immune_cells), 
                      cols.highlight = c('#d5d205'), 
                      cols= "#d5d5d5",
                      sizes.highlight = 1.25,
                      pt.size=1.25) +
    theme(text = element_text(size=20)) +
    labs(title = "Immune")
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_Immune_focus.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()

pdf('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_Immune_focus.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_87_2.png)
    



```R
# SCLC
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Manual_Annotation_simplified"), 
                      ncol = 1,
                      #cols = dittoColors()[1:15],
                      #cols=ditto_colors_navin_simplified,
                      cells.highlight= list("Unannotated" = Unannotated_cells), 
                      cols.highlight = c('#666666'), 
                      cols= "#d5d5d5",
                      sizes.highlight = 1.25,
                      pt.size=1.25) +
    theme(text = element_text(size=20)) +
    labs(title = "Unannotated")
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_Unannotated_focus.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()

pdf('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_Unannotated_focus.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_88_2.png)
    



```R
# SCLC
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Manual_Annotation_simplified"), 
                      ncol = 1,
                      #cols = dittoColors()[1:15],
                      #cols=ditto_colors_navin_simplified,
                      cells.highlight= list("Respiratory epithelium" = Respiratory_epithelium_cells), 
                      cols.highlight = c('#56B4E9'), 
                      cols= "#d5d5d5",
                      sizes.highlight = 1.25,
                      pt.size=1.25) +
    theme(text = element_text(size=20)) +
    labs(title = "Respiratory epithelium")
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_Resp_epithelium_focus.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()

pdf('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_Resp_epithelium_focus.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_89_2.png)
    



```R
# SCLC
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Manual_Annotation_simplified"), 
                      ncol = 1,
                      #cols = dittoColors()[1:15],
                      #cols=ditto_colors_navin_simplified,
                      cells.highlight= list("Fibrosis" = Fibrosis_cells), 
                      cols.highlight = c('#005685'), 
                      cols= "#d5d5d5",
                      sizes.highlight = 1.25,
                      pt.size=1.25) +
    theme(text = element_text(size=20)) +
    labs(title = "Fibrosis")
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_Fibrosis_focus.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()

pdf('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_Fibrosis_focus.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_90_2.png)
    



```R
# SCLC
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Manual_Annotation_simplified"), 
                      ncol = 1,
                      #cols = dittoColors()[1:15],
                      #cols=ditto_colors_navin_simplified,
                      cells.highlight= list("SCLC" = SCLC_cells), 
                      cols.highlight = c('#D55E00'), 
                      cols= "#d5d5d5",
                      sizes.highlight = 1.25,
                      pt.size=1.25) +
    theme(text = element_text(size=20)) +
    labs(title = "SCLC")
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_SCLC_focus.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()

pdf('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_simplified_UMAP_SCLC_focus.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_91_2.png)
    


# Continue with other plots


```R
ditto_colors_navin_minimal <- c(`Stroma/Other` = '#56B4E9',
                                  `SCLC` = '#D55E00',
                                  `Unannotated` = '#666666',
                                  `Immune` = '#F0E442')
```


```R
#group.by = c("Navin_Clusters", "Manual_Annotation", "Labeled_Clusters", "Manual_Annotation_simplified", 
#            "Tissue_Slice", "Manual_Annotation_minimal","Tumor_MHCI_label", "vasc_label"), 
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Manual_Annotation_minimal"), 
                      ncol = 1,
                      #cols = dittoColors()[1:15],
                      cols=ditto_colors_navin_minimal,
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_minimal_UMAP.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()

pdf('Output/16_Navin_Annotations_VA1_VA2/Manual_Annotation_minimal_UMAP.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_94_2.png)
    



```R
#group.by = c("Navin_Clusters", "Manual_Annotation", "Labeled_Clusters", "Manual_Annotation_simplified", 
#            "Tissue_Slice", "Manual_Annotation_minimal","Tumor_MHCI_label", "vasc_label"), 
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Tissue_Slice"), 
                      ncol = 1,
                      cols = dittoColors()[1:15],
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Tissue_Slice_UMAP.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()

pdf('Output/16_Navin_Annotations_VA1_VA2/Tissue_Slice_UMAP.pdf',
    width = 9,
    height = 7.5)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_95_2.png)
    



```R
#group.by = c("Navin_Clusters", "Manual_Annotation", "Labeled_Clusters", "Manual_Annotation_simplified", 
#            "Tissue_Slice", "Manual_Annotation_minimal","Tumor_MHCI_label", "vasc_label"), 
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Tumor_MHCI_label"), 
                      ncol = 1,
                      cols = dittoColors()[1:15],
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/Tumor_MHCI_label_UMAP.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_96_1.png)
    



```R
#group.by = c("Navin_Clusters", "Manual_Annotation", "Labeled_Clusters", "Manual_Annotation_simplified", 
#            "Tissue_Slice", "Manual_Annotation_minimal","Tumor_MHCI_label", "vasc_label"), 
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("vasc_label"), 
                      ncol = 1,
                      cols = dittoColors()[1:15],
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/vasc_label_UMAP.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_97_1.png)
    



```R
colnames(data_VA@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M2_Macrophage'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li><li>'Manual_Annotation'</li><li>'Manual_Annotation_simplified'</li><li>'Manual_Annotation_minimal'</li><li>'Labeled_Clusters'</li><li>'Labeled_Clusters_minimal'</li></ol>




```R
#group.by = c("Navin_Clusters", "Manual_Annotation", "Labeled_Clusters", "Manual_Annotation_simplified", 
#            "Tissue_Slice", "Manual_Annotation_minimal","Tumor_MHCI_label", "vasc_label"), 
options(repr.plot.width=10, repr.plot.height=8)
navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("vasc_MN4_label"), 
                      ncol = 1,
                      cols = dittoColors()[1:15],
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap
# save
png('Output/16_Navin_Annotations_VA1_VA2/vasc_MN4_label_UMAP.png',
    width = 900,
    height = 750)
print(navin_umap)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_99_1.png)
    


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


    
![png](output_102_0.png)
    



```R
options(repr.plot.width=17, repr.plot.height=8)
```


```R
navin_annot_percent_in_seurat_clusters_plot_VA <- dittoBarPlot(
    object = data_VA,
    var = "Manual_Annotation",
    group.by = "Navin_Clusters") + 
    theme(text = element_text(size=20))

png(file = "Output/16_Navin_Annotations_VA1_VA2/Navin_Annot_Percent_in_Seurat_Clusters_VA1.png",
    width = 1100,
    height = 700)
print(navin_annot_percent_in_seurat_clusters_plot_VA)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/Navin_Annot_Percent_in_Seurat_Clusters_VA1.pdf",
    width = 11,
    height = 7)
print(navin_annot_percent_in_seurat_clusters_plot_VA)
dev.off()

navin_annot_percent_in_seurat_clusters_plot_VA
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_104_2.png)
    



```R
navin_annot_percent_in_seurat_clusters_plot_VA_2 <- dittoBarPlot(
    object = data_VA,
    var = "Manual_Annotation_simplified",
    group.by = "Labeled_Clusters") + 
    theme(text = element_text(size=20))

png(file = "Output/16_Navin_Annotations_VA1_VA2/Navin_Annot_simplified_Percent_in_Navin_Annotated_Clusters_VA1.png",
    width = 1100,
    height = 700)
print(navin_annot_percent_in_seurat_clusters_plot_VA_2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/Navin_Annot_simplified_Percent_in_Navin_Annotated_Clusters_VA1.pdf",
    width = 11,
    height = 7)
print(navin_annot_percent_in_seurat_clusters_plot_VA_2)
dev.off()

navin_annot_percent_in_seurat_clusters_plot_VA_2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_105_2.png)
    


# Let's map very simplified and perform stats test (chi squared or fisher exact)


```R
ditto_colors_navin_minimal <- c(`Stroma/Other` = '#56B4E9',
                                  `SCLC` = '#D55E00',
                                  `Unannotated` = '#666666',
                                  `Immune` = '#F0E442')
```


```R
navin_annot_percent_in_seurat_clusters_plot_VA_2 <- dittoBarPlot(
    object = data_VA,
    var = "Manual_Annotation_minimal",
    color.panel = ditto_colors_navin_minimal,
    group.by = "Labeled_Clusters_minimal") + 
    theme(text = element_text(size=20))

png(file = "Output/16_Navin_Annotations_VA1_VA2/Navin_Annot_minimal_Percent_in_Labeled_Clusters_minimal_VA1.png",
    width = 700,
    height = 700)
print(navin_annot_percent_in_seurat_clusters_plot_VA_2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/Navin_Annot_minimal_Percent_in_Labeled_Clusters_minimal_VA1.pdf",
    width = 7,
    height = 7)
print(navin_annot_percent_in_seurat_clusters_plot_VA_2)
dev.off()

navin_annot_percent_in_seurat_clusters_plot_VA_2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_108_2.png)
    



```R
unique(data_VA@meta.data$Manual_Annotation_simplified)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Respiratory epithelium'</li><li>'SCLC'</li><li>'Fibrosis'</li><li>'Unannotated'</li><li>'Immune'</li><li>'Endothelial cells'</li></ol>




```R
colnames(data_VA@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_Cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M2_Macrophage'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li><li>'Manual_Annotation'</li><li>'Manual_Annotation_simplified'</li><li>'Manual_Annotation_minimal'</li><li>'Labeled_Clusters'</li><li>'Labeled_Clusters_minimal'</li></ol>




```R
ditto_colors_navin_simplified <- c(`Respiratory epithelium` = '#56B4E9',
                                  `SCLC` = '#D55E00',
                                  `Fibrosis` = '#005685',
                                  `Unannotated` = '#666666',
                                  `Immune` = '#F0E442',
                                  `Endothelial cells` = '#009E73')
```


```R
navin_annot_percent_in_seurat_clusters_plot_VA_2 <- dittoBarPlot(
    object = data_VA,
    var = "Manual_Annotation_simplified",
    color.panel = ditto_colors_navin_simplified,
    group.by = "Labeled_Clusters_minimal") + 
    theme(text = element_text(size=20))

png(file = "Output/16_Navin_Annotations_VA1_VA2/Navin_Annot_minimal_Percent_in_Labeled_Clusters_simplified_VA1.png",
    width = 700,
    height = 700)
print(navin_annot_percent_in_seurat_clusters_plot_VA_2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/Navin_Annot_minimal_Percent_in_Labeled_Clusters_simplified_VA1.pdf",
    width = 7,
    height = 7)
print(navin_annot_percent_in_seurat_clusters_plot_VA_2)
dev.off()

navin_annot_percent_in_seurat_clusters_plot_VA_2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_112_2.png)
    


# test if groups are significantly more represented in different annotations


```R
# Create a contingency table
table_data <- table(data_VA$Labeled_Clusters_minimal, data_VA$Manual_Annotation_minimal)
print('actual counts per group: ')
table_data

# Chi-squared test of independence
chi_result <- chisq.test(table_data)
print(chi_result)
print('standardized residuals: ')
chi_result$stdres
```

    [1] "actual counts per group: "



                  
                   Immune SCLC Stroma/Other Unannotated
      Immune          341  400          211         179
      SCLC            156 3232          647         630
      Stroma/Other    238  162          702         224


    
    	Pearson's Chi-squared test
    
    data:  table_data
    X-squared = 2173.7, df = 6, p-value < 2.2e-16
    
    [1] "standardized residuals: "



                  
                       Immune       SCLC Stroma/Other Unannotated
      Immune        23.901162 -13.158592    -2.879452    1.376917
      SCLC         -26.665254  37.315241   -22.590154   -3.300644
      Stroma/Other  10.121906 -33.214634    30.291132    2.737936



```R
names(chi_result)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'statistic'</li><li>'parameter'</li><li>'p.value'</li><li>'method'</li><li>'data.name'</li><li>'observed'</li><li>'expected'</li><li>'residuals'</li><li>'stdres'</li></ol>




```R
chi_result$expected
chi_result$residuals
```


<table class="dataframe">
<caption>A matrix: 3 × 4 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>Immune</th><th scope=col>SCLC</th><th scope=col>Stroma/Other</th><th scope=col>Unannotated</th></tr>
</thead>
<tbody>
	<tr><th scope=row>Immune</th><td>116.7207</td><td> 602.5013</td><td> 247.7338</td><td>164.0442</td></tr>
	<tr><th scope=row>SCLC</th><td>481.4343</td><td>2485.1179</td><td>1021.8197</td><td>676.6281</td></tr>
	<tr><th scope=row>Stroma/Other</th><td>136.8450</td><td> 706.3808</td><td> 290.4465</td><td>192.3277</td></tr>
</tbody>
</table>




                  
                       Immune       SCLC Stroma/Other Unannotated
      Immune        20.759417  -8.249901    -2.333851    1.167692
      SCLC         -14.831832  14.982301   -11.725606   -1.792554
      Stroma/Other   8.647149 -20.482518    24.148670    2.283802


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
    group.by = "Labeled_Clusters") + 
    theme(text = element_text(size=20))

png(file = "Output/16_Navin_Annotations_VA1_VA2/Tissue_Slices_in_Navin_Annotated_Clusters_AB1_A2.png",
    width = 1100,
    height = 700)
print(tissue_slices_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/Tissue_Slices_in_Navin_Annotated_Clusters_AB1_A2.pdf",
    width = 11,
    height = 7)
print(tissue_slices_in_seurat_clusters_plot_AB1_A2)
dev.off()

tissue_slices_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_120_2.png)
    



```R
colnames(AB1_A2@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M2_Macrophage'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li><li>'Manual_Annotation'</li><li>'Manual_Annotation_simplified'</li><li>'Manual_Annotation_minimal'</li><li>'Labeled_Clusters'</li><li>'Labeled_Clusters_minimal'</li></ol>



# by Navin Cluster Annotation


```R
Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2,
    var = "Tumor_MHCI_label",
    group.by = "Labeled_Clusters") + 
    theme(text = element_text(size=20))

png(file = "Output/16_Navin_Annotations_VA1_VA2/Tumor_MHCI_Label_in_Navin_Cluster_Annotation_AB1_A2.png",
    width = 1100,
    height = 700)
print(Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/Tumor_MHCI_Label_in_Navin_Cluster_Annotation_AB1_A2.pdf",
    width = 11,
    height = 7)
print(Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

Tumor_MHCI_Label_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_123_2.png)
    



```R
vasc_Label_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2,
    var = "vasc_label",
    group.by = "Labeled_Clusters") + 
    theme(text = element_text(size=20))

png(file = "Output/16_Navin_Annotations_VA1_VA2/vasc_Label_in_Navin_Cluster_Annotation_AB1_A2.png",
    width = 1100,
    height = 700)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/vasc_Label_in_Navin_Cluster_Annotation_AB1_A2.pdf",
    width = 11,
    height = 7)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

vasc_Label_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_124_2.png)
    


# Vasc MN4 label in Tumor MHCI


```R
vasc_Label_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2,
    var = "vasc_MN4_label",
    group.by = "Tumor_MHCI_label") + 
    theme(text = element_text(size=20))

png(file = "Output/16_Navin_Annotations_VA1_VA2/vasc_MN4_Label_in_Tumor_MHCI.png",
    width = 1100,
    height = 700)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/vasc_MN4_Label_in_Tumor_MHCI.pdf",
    width = 11,
    height = 7)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

vasc_Label_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_126_2.png)
    


# Vasc label in Tumor MHCI


```R
vasc_Label_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2,
    var = "vasc_label",
    group.by = "Tumor_MHCI_label") + 
    theme(text = element_text(size=20))

png(file = "Output/16_Navin_Annotations_VA1_VA2/vasc_Label_in_Tumor_MHCI.png",
    width = 1100,
    height = 700)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/vasc_Label_in_Tumor_MHCI.pdf",
    width = 11,
    height = 7)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

vasc_Label_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_128_2.png)
    



```R
AB1_A2@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 64</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_Cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Manual_Annotation</th><th scope=col>Manual_Annotation_simplified</th><th scope=col>Manual_Annotation_minimal</th><th scope=col>Labeled_Clusters</th><th scope=col>Labeled_Clusters_minimal</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Stroma/Other</td><td>Normal_Lung         </td><td>Stroma/Other</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC        </td><td>SCLC_Effector_Immune</td><td>SCLC        </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC        </td><td>SCLC_Effector_Immune</td><td>SCLC        </td></tr>
</tbody>
</table>




```R
unique(AB1_A2@meta.data$Tumor_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high'</li></ol>




```R
# subset
AB1_A2_tumor_high_subset = subset(AB1_A2,subset = Tumor_label == 'Tumor_high')
```

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



```R
vasc_Label_in_seurat_clusters_plot_AB1_A2 <- dittoBarPlot(
    object = AB1_A2_tumor_high_subset,
    var = "vasc_label",
    group.by = "Tumor_MHCI_label") + 
    theme(text = element_text(size=20))

png(file = "Output/16_Navin_Annotations_VA1_VA2/vasc_Label_in_Tumor_High_MHCI.png",
    width = 600,
    height = 700)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/vasc_Label_in_Tumor_High_MHCI.pdf",
    width = 6,
    height = 7)
print(vasc_Label_in_seurat_clusters_plot_AB1_A2)
dev.off()

vasc_Label_in_seurat_clusters_plot_AB1_A2
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_132_2.png)
    


# test if groups are significantly more represented in different annotations


```R
# Create a contingency table
table_data <- table(AB1_A2_tumor_high_subset$vasc_label, AB1_A2_tumor_high_subset$Tumor_MHCI_label)
print('actual counts per group: ')
table_data

# Chi-squared test of independence
chi_result <- chisq.test(table_data)
print(chi_result)
print('standardized residuals: ')
chi_result$stdres
```

    [1] "actual counts per group: "



               
                Tumor_high_MHCI_high Tumor_high_MHCI_low
      vasc_high                  590                 619
      vasc_low                  1313                2961


    
    	Pearson's Chi-squared test with Yates' continuity correction
    
    data:  table_data
    X-squared = 135.15, df = 1, p-value < 2.2e-16
    
    [1] "standardized residuals: "



               
                Tumor_high_MHCI_high Tumor_high_MHCI_low
      vasc_high             11.65944           -11.65944
      vasc_low             -11.65944            11.65944



```R
names(chi_result)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'statistic'</li><li>'parameter'</li><li>'p.value'</li><li>'method'</li><li>'data.name'</li><li>'observed'</li><li>'expected'</li><li>'residuals'</li><li>'stdres'</li></ol>




```R
chi_result$expected
chi_result$residuals
```


<table class="dataframe">
<caption>A matrix: 2 × 2 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>Tumor_high_MHCI_high</th><th scope=col>Tumor_high_MHCI_low</th></tr>
</thead>
<tbody>
	<tr><th scope=row>vasc_high</th><td> 419.611</td><td> 789.389</td></tr>
	<tr><th scope=row>vasc_low</th><td>1483.389</td><td>2790.611</td></tr>
</tbody>
</table>




               
                Tumor_high_MHCI_high Tumor_high_MHCI_low
      vasc_high             8.317986           -6.064515
      vasc_low             -4.423989            3.225462


# Check if vasc score is higher in Tumor high MHCI high spots


```R
colnames(AB1_A2_tumor_high_subset@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_Cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M2_Macrophage'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li><li>'Manual_Annotation'</li><li>'Manual_Annotation_simplified'</li><li>'Manual_Annotation_minimal'</li><li>'Labeled_Clusters'</li><li>'Labeled_Clusters_minimal'</li></ol>




```R
Idents(AB1_A2_tumor_high_subset) <- 'Tumor_MHCI_label'
p1 <- VlnPlot(AB1_A2_tumor_high_subset, features=c('vasc_summed_expr'), split.by='Tumor_MHCI_label')
p1
```


    
![png](output_139_0.png)
    


# SpatialDimPlots
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
                  '11_SCLC_Necrotic' = dittoColors()[19])options(repr.plot.width=28, repr.plot.height=8)

Idents(AB1_A2) <- "Navin_Clusters"

plot0 <- SpatialDimPlot(AB1_A2, cols = ditto_colors, pt.size.factor=0)
plot0 

png(file = "Output/16_Navin_Annotations_VA1_VA2/HnEs.png",
    width = 2400,
    height = 700)
print(plot0)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/HnEs.pdf",
    width = 24,
    height = 7)
print(plot0)
dev.off()

plot1 <- SpatialDimPlot(AB1_A2, cols = ditto_colors, pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=8)))
plot1 

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot.png",
    width = 2400,
    height = 700)
print(plot1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot.pdf",
    width = 24,
    height = 7)
print(plot1)
dev.off()# plot in one image
options(repr.plot.width=28, repr.plot.height=16)

Idents(AB1_A2) <- "Navin_Clusters"

plot0 <- SpatialDimPlot(AB1_A2, cols = ditto_colors, pt.size.factor=0)

plot1 <- SpatialDimPlot(AB1_A2, cols = ditto_colors, pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=8)))

plot0 / plot1

png(file = "Output/16_Navin_Annotations_VA1_VA2/HnEs_and_DimPlot.png",
    width = 2400,
    height = 1400)
print(plot0 / plot1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/HnEs_and_DimPlot.pdf",
    width = 24,
    height = 14)
print(plot0 / plot1)
dev.off()
# Look at collapsed cluster labels


```R
show_col(dittoColors()[1:24])
```


    
![png](output_145_0.png)
    



```R
length(unique(data_VA@meta.data$Labeled_Clusters))
```


8



```R
unique(data_VA@meta.data$Labeled_Clusters)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>Normal_Lung</li><li>SCLC_Effector_Immune</li><li>SCLC</li><li>SCLC_TSI</li><li>SCLC_Stroma</li><li>Lymphocytes_B_Cells</li><li>Fibroblasts</li><li>SCLC_Apoptotic</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'Lymphocytes_B_Cells'</li><li>'SCLC'</li><li>'SCLC_TSI'</li><li>'SCLC_Stroma'</li><li>'SCLC_Effector_Immune'</li><li>'Normal_Lung'</li><li>'Fibroblasts'</li><li>'SCLC_Apoptotic'</li></ol>
</details>



```R
ditto_colors_Labeled_Clusters <- c(`Normal_Lung` = '#56B4E9',
                                  `SCLC` = '#D55E00',
                                  `Fibroblasts` = '#005685',
                                  `SCLC_Apoptotic` = '#666666',
                                  `Lymphocytes_B_Cells` = '#F0E442',
                                  `SCLC_TSI` = '#009E73',
                                  `SCLC_Stroma` = '#005685',
                                  `SCLC_Effector_Immune` = '#E69F00')
```


```R
ditto_colors_Labeled_Clusters
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>Normal_Lung</dt><dd>'#56B4E9'</dd><dt>SCLC</dt><dd>'#D55E00'</dd><dt>Fibroblasts</dt><dd>'#005685'</dd><dt>SCLC_Apoptotic</dt><dd>'#666666'</dd><dt>Lymphocytes_B_Cells</dt><dd>'#F0E442'</dd><dt>SCLC_TSI</dt><dd>'#009E73'</dd><dt>SCLC_Stroma</dt><dd>'#005685'</dd><dt>SCLC_Effector_Immune</dt><dd>'#E69F00'</dd></dl>


# define cluster colors based on dittoColors...
ditto_colors2 <- c('0_Lymphocytes_B_Cells' = dittoColors()[1], 
                  '1_4_5_11_SCLC' = dittoColors()[2], 
                  '2_SCLC_TSI' = dittoColors()[3], 
                  '3_SCLC_Stroma' = dittoColors()[4], 
                  '6_SCLC_Effector_Immune' = dittoColors()[7], 
                  '7_10_Normal_Lung' = dittoColors()[8],
                  '8_Fibroblasts' = dittoColors()[9], 
                  '9_SCLC_Apoptotic' = dittoColors()[10])

ditto_colors3 <- c('Lymphocytes_B_Cells' = dittoColors()[1], 
                  'SCLC' = dittoColors()[2], 
                  'SCLC_TSI' = dittoColors()[3], 
                  'SCLC_Stroma' = dittoColors()[4], 
                  'SCLC_Effector_Immune' = dittoColors()[7], 
                  'Normal_Lung' = dittoColors()[8],
                  'Fibroblasts' = dittoColors()[9], 
                  'SCLC_Apoptotic' = dittoColors()[10])

```R
show_col(ditto_colors_Labeled_Clusters)
```


    
![png](output_151_0.png)
    



```R
options(repr.plot.width=28, repr.plot.height=8)

Idents(AB1_A2) <- "Labeled_Clusters"

plot1 <- SpatialDimPlot(AB1_A2, cols = ditto_colors_Labeled_Clusters, pt.size.factor=2.2, image.alpha = 0.9, stroke = 0.3) & guides(fill = guide_legend(override.aes = list(size=8)))
plot1 

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Labeled_Clusters.png",
    width = 2400,
    height = 700)
print(plot1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Labeled_Clusters.pdf",
    width = 24,
    height = 7)
print(plot1)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_152_2.png)
    



```R
# plot in one image
options(repr.plot.width=28, repr.plot.height=16)

Idents(AB1_A2) <- "Labeled_Clusters"

plot0 <- SpatialDimPlot(AB1_A2, cols = ditto_colors_Labeled_Clusters, pt.size.factor=0)

plot1 <- SpatialDimPlot(AB1_A2, cols = ditto_colors_Labeled_Clusters, pt.size.factor=2.2, image.alpha = 0.9, stroke = 0.3) & guides(fill = guide_legend(override.aes = list(size=8)))

plot0 / plot1

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Labeled_Clusters_HnE.pdf",
    width = 2400,
    height = 1400)
print(plot0 / plot1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Labeled_Clusters_HnE.pdf",
    width = 24,
    height = 14)
print(plot0 / plot1)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_153_2.png)
    


# Individually


```R
options(repr.plot.width=8.5, repr.plot.height=8)

Idents(AB1_A2) <- "Labeled_Clusters"

plot1 <- SpatialDimPlot(AB1_A2, cols = ditto_colors_Labeled_Clusters, pt.size.factor=2.5, image.alpha = 0.9, stroke = 0.3) & guides(fill = guide_legend(override.aes = list(size=8)))
plot_A1 <- plot1[[2]]
plot_B1 <- plot1[[3]]
plot_A2 <- plot1[[1]]

plot_A1
plot_B1
plot_A2

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Labeled_Clusters_VA1.png",width = 800,height = 700)
print(plot_A1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Labeled_Clusters_VA1.pdf",width = 8,height = 7)
print(plot_A1)
dev.off()

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Labeled_Clusters_VB1.png",width = 800,height = 700)
print(plot_B1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Labeled_Clusters_VB1.pdf",width = 8,height = 7)
print(plot_B1)
dev.off()

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Labeled_Clusters_VA2.png",width = 800,height = 700)
print(plot_A2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Labeled_Clusters_VA2.pdf",width = 8,height = 7)
print(plot_A2)
dev.off()
```


    
![png](output_155_0.png)
    



    
![png](output_155_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_155_8.png)
    



```R

```


```R
show_col(dittoColors()[1:36])
```


    
![png](output_157_0.png)
    



```R
length(unique(data_VA@meta.data$Manual_Annotation))
```


15



```R
unique(data_VA@meta.data$Manual_Annotation)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Respiratory epithelium'</li><li>'SCLC'</li><li>'SCLC+vasculature'</li><li>'Fibrosis'</li><li>'Unannotated'</li><li>'Immune cells+vasculature'</li><li>'Immune cells_lymphocytes'</li><li>'SCLC_TSI'</li><li>'Immune cells_myeloid cells'</li><li>'SCLC+immune cells'</li><li>'Smooth muscle'</li><li>'Endothelial cells'</li><li>'Tumor necrosis'</li><li>'SCLC+fibrosis'</li><li>'Alveoli'</li></ol>




```R
ditto_colors_navin <- setNames(as.list(dittoColors()[1:15]),
                               as.list(sort(unique(data_VA@meta.data$Manual_Annotation))))
```


```R
ditto_colors_navin
```


<dl>
	<dt>$Alveoli</dt>
		<dd>'#E69F00'</dd>
	<dt>$`Endothelial cells`</dt>
		<dd>'#56B4E9'</dd>
	<dt>$Fibrosis</dt>
		<dd>'#009E73'</dd>
	<dt>$`Immune cells_lymphocytes`</dt>
		<dd>'#F0E442'</dd>
	<dt>$`Immune cells_myeloid cells`</dt>
		<dd>'#0072B2'</dd>
	<dt>$`Immune cells+vasculature`</dt>
		<dd>'#D55E00'</dd>
	<dt>$`Respiratory epithelium`</dt>
		<dd>'#CC79A7'</dd>
	<dt>$SCLC</dt>
		<dd>'#666666'</dd>
	<dt>$SCLC_TSI</dt>
		<dd>'#AD7700'</dd>
	<dt>$`SCLC+fibrosis`</dt>
		<dd>'#1C91D4'</dd>
	<dt>$`SCLC+immune cells`</dt>
		<dd>'#007756'</dd>
	<dt>$`SCLC+vasculature`</dt>
		<dd>'#D5C711'</dd>
	<dt>$`Smooth muscle`</dt>
		<dd>'#005685'</dd>
	<dt>$`Tumor necrosis`</dt>
		<dd>'#A04700'</dd>
	<dt>$Unannotated</dt>
		<dd>'#B14380'</dd>
</dl>




```R
# spatial dimplot colored by navin annotations
Idents(data_VA) <- "Manual_Annotation"

options(repr.plot.width=17, repr.plot.height=8)

plot2 <- SpatialDimPlot(data_VA, cols = ditto_colors_navin, pt.size.factor=2.5, image.alpha = 0.9, stroke = 0.3) & guides(fill = guide_legend(override.aes = list(size=8)))
plot_VA1 <- plot2[[2]]
plot_VA2 <- plot2[[1]]
plot2
plot_VA1
plot_VA2

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_VA1_VA2.png",width = 1700,height = 700)
print(plot2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_VA1_VA2.pdf",width = 17,height = 7)
print(plot2)
dev.off()

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_VA1.png",width = 850,height = 700)
print(plot_VA1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_VA1.pdf",width = 8.5,height = 7)
print(plot_VA1)
dev.off()

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_VA2.png",width = 850,height = 700)
print(plot_VA2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_VA2.pdf",width = 8.5,height = 7)
print(plot_VA2)
dev.off()
```


    
![png](output_162_0.png)
    



    
![png](output_162_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_162_8.png)
    


# plot simplified annotations


```R
show_col(dittoColors()[1:16])
```


    
![png](output_164_0.png)
    



```R
length(unique(data_VA@meta.data$Manual_Annotation_simplified))
```


6



```R
unique(data_VA@meta.data$Manual_Annotation_simplified)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Respiratory epithelium'</li><li>'SCLC'</li><li>'Fibrosis'</li><li>'Unannotated'</li><li>'Immune'</li><li>'Endothelial cells'</li></ol>




```R
ditto_colors_navin_simplified <- c(`Respiratory epithelium` = '#56B4E9',
                                  `SCLC` = '#D55E00',
                                  `Fibrosis` = '#005685',
                                  `Unannotated` = '#666666',
                                  `Immune` = '#F0E442',
                                  `Endothelial cells` = '#009E73')
```


```R
ditto_colors_navin_simplified
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>Respiratory epithelium</dt><dd>'#56B4E9'</dd><dt>SCLC</dt><dd>'#D55E00'</dd><dt>Fibrosis</dt><dd>'#005685'</dd><dt>Unannotated</dt><dd>'#666666'</dd><dt>Immune</dt><dd>'#F0E442'</dd><dt>Endothelial cells</dt><dd>'#009E73'</dd></dl>




```R
# spatial dimplot colored by navin annotations simplified
Idents(data_VA) <- "Manual_Annotation_simplified"

options(repr.plot.width=17, repr.plot.height=8)

plot2 <- SpatialDimPlot(data_VA, cols = ditto_colors_navin_simplified, pt.size.factor=2.5, image.alpha = 0.9, stroke = 0.3) & guides(fill = guide_legend(override.aes = list(size=8)))
plot_VA1 <- plot2[[2]]
plot_VA2 <- plot2[[1]]
plot2
plot_VA1
plot_VA2

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_simplified_VA1_VA2.png",width = 1700,height = 700)
print(plot2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_simplified_VA1_VA2.pdf",width = 17,height = 7)
print(plot2)
dev.off()

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_simplified_VA1.png",width = 850,height = 700)
print(plot_VA1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_simplified_VA1.pdf",width = 8.5,height = 7)
print(plot_VA1)
dev.off()

png(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_simplified_VA2.png",width = 850,height = 700)
print(plot_VA2)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/DimPlot_Navin_Annot_simplified_VA2.pdf",width = 8.5,height = 7)
print(plot_VA2)
dev.off()
```


    
![png](output_169_0.png)
    



    
![png](output_169_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_169_8.png)
    


# Let's try to make the GSVA heatmap vs labeled clusters (or Navin annotated clusters)


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
    35796 features across 7122 samples within 2 assays 
    Active assay: Spatial (17943 features, 2000 variable features)
     3 layers present: data, counts, scale.data
     1 other assay present: SCT
     3 dimensional reductions calculated: pca, integrated.cca, umap
     2 spatial fields of view present: slice1 slice1.2



```R
AB1_A2@meta.data %>% head(3)
data_VA@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 64</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Manual_Annotation</th><th scope=col>Manual_Annotation_simplified</th><th scope=col>Manual_Annotation_minimal</th><th scope=col>Labeled_Clusters</th><th scope=col>Labeled_Clusters_minimal</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Stroma/Other</td><td>Normal_Lung         </td><td>Stroma/Other</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC        </td><td>SCLC_Effector_Immune</td><td>SCLC        </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC        </td><td>SCLC_Effector_Immune</td><td>SCLC        </td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 3 × 64</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Manual_Annotation</th><th scope=col>Manual_Annotation_simplified</th><th scope=col>Manual_Annotation_minimal</th><th scope=col>Labeled_Clusters</th><th scope=col>Labeled_Clusters_minimal</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Stroma/Other</td><td>Normal_Lung         </td><td>Stroma/Other</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC        </td><td>SCLC_Effector_Immune</td><td>SCLC        </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC        </td><td>SCLC_Effector_Immune</td><td>SCLC        </td></tr>
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
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M2_Macrophage'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li><li>'Manual_Annotation'</li><li>'Manual_Annotation_simplified'</li><li>'Manual_Annotation_minimal'</li><li>'Labeled_Clusters'</li><li>'Labeled_Clusters_minimal'</li></ol>




```R
# Rename Epithelial_cell to Epithelial_Cell
AB1_A2@meta.data <- AB1_A2@meta.data %>% rename('Epithelial_Cell'='Epithelial_cell')
data_VA@meta.data <- data_VA@meta.data %>% rename('Epithelial_Cell'='Epithelial_cell')
```


```R
colnames(AB1_A2@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_Cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M2_Macrophage'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li><li>'Manual_Annotation'</li><li>'Manual_Annotation_simplified'</li><li>'Manual_Annotation_minimal'</li><li>'Labeled_Clusters'</li><li>'Labeled_Clusters_minimal'</li></ol>




```R
cols_pivot1 <- c('Barcode','Labeled_Clusters','DC','NK_Cell','T_Cell','B_Cell','Macrophage','CD8_T_Cell','Endothelial_Cell','Epithelial_Cell','Fibroblast')
cols_pivot2 <- c('Barcode','Manual_Annotation','DC','NK_Cell','T_Cell','B_Cell','Macrophage','CD8_T_Cell','Endothelial_Cell','Epithelial_Cell','Fibroblast')
cols_pivot2_simp <- c('Barcode','Manual_Annotation_simplified','DC','NK_Cell','T_Cell','B_Cell','Macrophage','CD8_T_Cell','Endothelial_Cell','Epithelial_Cell','Fibroblast')
cols_pivot2_min <- c('Barcode','Manual_Annotation_minimal','DC','NK_Cell','T_Cell','B_Cell','Macrophage','CD8_T_Cell','Endothelial_Cell','Epithelial_Cell','Fibroblast')
cols_keep1 <- c('DC','NK_Cell','T_Cell','B_Cell','Macrophage','CD8_T_Cell','Endothelial_Cell','Epithelial_Cell','Fibroblast')
```


```R
AB1_A2_Labeled_Clusters <- AB1_A2@meta.data %>% select(all_of(cols_pivot1)) %>% pivot_longer(!c(Barcode, Labeled_Clusters), names_to = "cell_type", values_to = "GSVA_Enrichment")
data_VA_Manual_Annotation <- data_VA@meta.data %>% select(all_of(cols_pivot2)) %>% pivot_longer(!c(Barcode, Manual_Annotation), names_to = "cell_type", values_to = "GSVA_Enrichment")
data_VA_Manual_Annotation_simp <- data_VA@meta.data %>% select(all_of(cols_pivot2_simp)) %>% pivot_longer(!c(Barcode, Manual_Annotation_simplified), names_to = "cell_type", values_to = "GSVA_Enrichment")
data_VA_Manual_Annotation_min <- data_VA@meta.data %>% select(all_of(cols_pivot2_min)) %>% pivot_longer(!c(Barcode, Manual_Annotation_minimal), names_to = "cell_type", values_to = "GSVA_Enrichment")
```


```R
AB1_A2_Labeled_Clusters %>% head(5)
```


<table class="dataframe">
<caption>A tibble: 5 × 4</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>Labeled_Clusters</th><th scope=col>cell_type</th><th scope=col>GSVA_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>Normal_Lung</td><td>DC        </td><td>-0.007136273</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>Normal_Lung</td><td>NK_Cell   </td><td> 0.015558701</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>Normal_Lung</td><td>T_Cell    </td><td> 0.633032413</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>Normal_Lung</td><td>B_Cell    </td><td>-0.038430718</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>Normal_Lung</td><td>Macrophage</td><td>-0.165127095</td></tr>
</tbody>
</table>




```R
AB1_A2_Labeled_Clusters_mat <- AB1_A2_Labeled_Clusters %>% 
                group_by(cell_type,Labeled_Clusters) %>% 
                summarise(mean_GSVA_pergroup = mean(GSVA_Enrichment)) %>% 
                spread(Labeled_Clusters,mean_GSVA_pergroup) %>% 
                column_to_rownames('cell_type')
#AB1_A2_Labeled_Clusters_mat <- AB1_A2_Labeled_Clusters_mat %>% rename(Unlabeled = V1)

data_VA_Manual_Annotation_mat <- data_VA_Manual_Annotation %>% 
                group_by(cell_type,Manual_Annotation) %>% 
                summarise(mean_GSVA_pergroup = mean(GSVA_Enrichment)) %>% 
                spread(Manual_Annotation,mean_GSVA_pergroup) %>% 
                column_to_rownames('cell_type')
                #rename(Unlabeled = V1)
                #select(!V1)
#data_VA_Manual_Annotation_mat <- data_VA_Manual_Annotation_mat %>% rename(Unlabeled = V1)

data_VA_Manual_Annotation_simp_mat <- data_VA_Manual_Annotation_simp %>% 
                group_by(cell_type,Manual_Annotation_simplified) %>% 
                summarise(mean_GSVA_pergroup = mean(GSVA_Enrichment)) %>% 
                spread(Manual_Annotation_simplified,mean_GSVA_pergroup) %>% 
                column_to_rownames('cell_type')
                #rename(Unlabeled = V1)
                #select(!V1)
#data_VA_Manual_Annotation_simp_mat <- data_VA_Manual_Annotation_simp_mat %>% rename(Unlabeled = V1)

data_VA_Manual_Annotation_min_mat <- data_VA_Manual_Annotation_min %>% 
                group_by(cell_type,Manual_Annotation_minimal) %>% 
                summarise(mean_GSVA_pergroup = mean(GSVA_Enrichment)) %>% 
                spread(Manual_Annotation_minimal,mean_GSVA_pergroup) %>% 
                column_to_rownames('cell_type')
                #rename(Unlabeled = V1)
                #select(!V1)
#data_VA_Manual_Annotation_min_mat <- data_VA_Manual_Annotation_min_mat %>% rename(Unlabeled = V1)
```

    [1m[22m`summarise()` has grouped output by 'cell_type'. You can override using the
    `.groups` argument.
    [1m[22m`summarise()` has grouped output by 'cell_type'. You can override using the
    `.groups` argument.
    [1m[22m`summarise()` has grouped output by 'cell_type'. You can override using the
    `.groups` argument.
    [1m[22m`summarise()` has grouped output by 'cell_type'. You can override using the
    `.groups` argument.



```R
colnames(AB1_A2_Labeled_Clusters_mat)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Lymphocytes_B_Cells'</li><li>'SCLC'</li><li>'SCLC_TSI'</li><li>'SCLC_Stroma'</li><li>'SCLC_Effector_Immune'</li><li>'Normal_Lung'</li><li>'Fibroblasts'</li><li>'SCLC_Apoptotic'</li></ol>




```R
AB1_A2_Labeled_Clusters_mat
data_VA_Manual_Annotation_mat
data_VA_Manual_Annotation_simp_mat
data_VA_Manual_Annotation_min_mat
```


<table class="dataframe">
<caption>A data.frame: 9 × 8</caption>
<thead>
	<tr><th></th><th scope=col>Lymphocytes_B_Cells</th><th scope=col>SCLC</th><th scope=col>SCLC_TSI</th><th scope=col>SCLC_Stroma</th><th scope=col>SCLC_Effector_Immune</th><th scope=col>Normal_Lung</th><th scope=col>Fibroblasts</th><th scope=col>SCLC_Apoptotic</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B_Cell</th><td> 0.46915645</td><td>-0.40129361</td><td> 0.012871339</td><td>-0.12645085</td><td>-0.176410302</td><td> 0.28198183</td><td> 0.08304008</td><td>-0.324447367</td></tr>
	<tr><th scope=row>CD8_T_Cell</th><td> 0.14926360</td><td>-0.24414396</td><td> 0.023927991</td><td> 0.06850047</td><td>-0.042777073</td><td> 0.27749976</td><td> 0.13113426</td><td>-0.147719230</td></tr>
	<tr><th scope=row>DC</th><td> 0.28466514</td><td>-0.27949359</td><td>-0.049873615</td><td> 0.03916977</td><td>-0.155406625</td><td> 0.06427179</td><td> 0.07730561</td><td>-0.288146373</td></tr>
	<tr><th scope=row>Endothelial_Cell</th><td> 0.09109064</td><td>-0.26230273</td><td>-0.071015097</td><td>-0.08262259</td><td>-0.129507477</td><td> 0.33734812</td><td> 0.05974505</td><td>-0.139393670</td></tr>
	<tr><th scope=row>Epithelial_Cell</th><td>-0.16567499</td><td> 0.02333796</td><td>-0.213228685</td><td>-0.09440080</td><td>-0.081650635</td><td>-0.11180690</td><td>-0.23401081</td><td> 0.007676266</td></tr>
	<tr><th scope=row>Fibroblast</th><td> 0.35651810</td><td>-0.42574978</td><td>-0.074324351</td><td> 0.12138585</td><td>-0.249806157</td><td> 0.31837832</td><td> 0.59493137</td><td>-0.382076232</td></tr>
	<tr><th scope=row>Macrophage</th><td> 0.06524891</td><td>-0.31398641</td><td>-0.240676274</td><td> 0.10324018</td><td>-0.203921708</td><td> 0.01960697</td><td> 0.15294566</td><td>-0.322508157</td></tr>
	<tr><th scope=row>NK_Cell</th><td> 0.10847699</td><td>-0.10498383</td><td> 0.034025903</td><td> 0.09963871</td><td> 0.005818302</td><td> 0.10076710</td><td> 0.10824206</td><td>-0.096980105</td></tr>
	<tr><th scope=row>T_Cell</th><td> 0.16856428</td><td>-0.28385005</td><td>-0.007903851</td><td> 0.02358209</td><td>-0.111682947</td><td> 0.19945756</td><td> 0.18999155</td><td>-0.225536122</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 9 × 15</caption>
<thead>
	<tr><th></th><th scope=col>Alveoli</th><th scope=col>Endothelial cells</th><th scope=col>Fibrosis</th><th scope=col>Immune cells_lymphocytes</th><th scope=col>Immune cells_myeloid cells</th><th scope=col>Immune cells+vasculature</th><th scope=col>Respiratory epithelium</th><th scope=col>SCLC</th><th scope=col>SCLC_TSI</th><th scope=col>SCLC+fibrosis</th><th scope=col>SCLC+immune cells</th><th scope=col>SCLC+vasculature</th><th scope=col>Smooth muscle</th><th scope=col>Tumor necrosis</th><th scope=col>Unannotated</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B_Cell</th><td> 0.30979455</td><td> 0.086804129</td><td> 0.07992840</td><td> 0.39207599</td><td> 0.089069636</td><td> 0.2949419</td><td> 2.618262e-01</td><td>-0.327329978</td><td>-0.17390821</td><td>-0.08771357</td><td> 0.0767659087</td><td>-0.20773076</td><td> 0.26090008</td><td>-0.072803391</td><td> 0.046646311</td></tr>
	<tr><th scope=row>CD8_T_Cell</th><td> 0.30271168</td><td> 0.038881182</td><td> 0.01832951</td><td> 0.29720881</td><td> 0.200984240</td><td> 0.1923273</td><td> 8.004192e-02</td><td>-0.251707066</td><td>-0.12605817</td><td>-0.18891002</td><td>-0.0400619888</td><td>-0.10509266</td><td> 0.14408044</td><td>-0.008325300</td><td> 0.020076269</td></tr>
	<tr><th scope=row>DC</th><td>-0.04880404</td><td> 0.083306169</td><td> 0.01558226</td><td> 0.28692848</td><td> 0.199222405</td><td> 0.2995857</td><td> 5.942674e-02</td><td>-0.285076066</td><td>-0.14670477</td><td>-0.18260870</td><td>-0.0001557487</td><td>-0.09224352</td><td>-0.09752169</td><td>-0.188632104</td><td>-0.012589191</td></tr>
	<tr><th scope=row>Endothelial_Cell</th><td> 0.28102364</td><td> 0.165710529</td><td>-0.04806745</td><td> 0.17287354</td><td> 0.078189858</td><td> 0.1270248</td><td> 9.423644e-02</td><td>-0.265252538</td><td>-0.17044989</td><td>-0.19868387</td><td>-0.0530467914</td><td>-0.08606615</td><td> 0.18755837</td><td>-0.054899254</td><td>-0.001813333</td></tr>
	<tr><th scope=row>Epithelial_Cell</th><td>-0.14570260</td><td>-0.158975910</td><td>-0.16776065</td><td>-0.20087356</td><td>-0.209292908</td><td>-0.1578715</td><td> 1.114332e-02</td><td> 0.004411367</td><td>-0.03684068</td><td>-0.10344673</td><td>-0.0916099336</td><td>-0.03556302</td><td>-0.31254515</td><td>-0.093054755</td><td>-0.129368982</td></tr>
	<tr><th scope=row>Fibroblast</th><td> 0.21655789</td><td> 0.276767003</td><td> 0.43941113</td><td> 0.46573851</td><td> 0.414811234</td><td> 0.3018786</td><td> 1.203406e-01</td><td>-0.421354910</td><td>-0.14869220</td><td> 0.08931514</td><td> 0.1172035761</td><td>-0.12808118</td><td> 0.35882612</td><td>-0.115077796</td><td> 0.033637049</td></tr>
	<tr><th scope=row>Macrophage</th><td>-0.03236368</td><td> 0.007488114</td><td> 0.10237376</td><td> 0.09321139</td><td>-0.004107695</td><td> 0.1217390</td><td>-5.767658e-02</td><td>-0.337663746</td><td>-0.21820860</td><td>-0.07410176</td><td>-0.0299543644</td><td>-0.17132880</td><td>-0.25571139</td><td>-0.114365617</td><td>-0.105458514</td></tr>
	<tr><th scope=row>NK_Cell</th><td> 0.01872419</td><td> 0.192857999</td><td> 0.11207807</td><td> 0.18016854</td><td> 0.209916531</td><td> 0.1757864</td><td> 8.100107e-02</td><td>-0.074555249</td><td> 0.02109971</td><td>-0.04277419</td><td> 0.0762688773</td><td> 0.03766665</td><td> 0.03137024</td><td>-0.006218108</td><td> 0.077933986</td></tr>
	<tr><th scope=row>T_Cell</th><td> 0.19224392</td><td> 0.009909457</td><td> 0.08105488</td><td> 0.30903705</td><td> 0.140945163</td><td> 0.2262118</td><td> 7.979539e-06</td><td>-0.260795864</td><td>-0.12860793</td><td>-0.14833769</td><td> 0.0290415527</td><td>-0.11547943</td><td> 0.12110411</td><td>-0.045533466</td><td> 0.026337183</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 9 × 6</caption>
<thead>
	<tr><th></th><th scope=col>Endothelial cells</th><th scope=col>Fibrosis</th><th scope=col>Immune</th><th scope=col>Respiratory epithelium</th><th scope=col>SCLC</th><th scope=col>Unannotated</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B_Cell</th><td> 0.086804129</td><td> 0.07992840</td><td> 0.3610330</td><td> 0.30295968</td><td>-0.24655021</td><td> 0.046646311</td></tr>
	<tr><th scope=row>CD8_T_Cell</th><td> 0.038881182</td><td> 0.01832951</td><td> 0.2765727</td><td> 0.27609259</td><td>-0.19527232</td><td> 0.020076269</td></tr>
	<tr><th scope=row>DC</th><td> 0.083306169</td><td> 0.01558226</td><td> 0.2838786</td><td>-0.04527117</td><td>-0.22518499</td><td>-0.012589191</td></tr>
	<tr><th scope=row>Endothelial_Cell</th><td> 0.165710529</td><td>-0.04806745</td><td> 0.1609173</td><td> 0.26165459</td><td>-0.20820564</td><td>-0.001813333</td></tr>
	<tr><th scope=row>Epithelial_Cell</th><td>-0.158975910</td><td>-0.16776065</td><td>-0.1950830</td><td>-0.14777721</td><td>-0.02397166</td><td>-0.129368982</td></tr>
	<tr><th scope=row>Fibroblast</th><td> 0.276767003</td><td> 0.43941113</td><td> 0.4390432</td><td> 0.22081459</td><td>-0.27696128</td><td> 0.033637049</td></tr>
	<tr><th scope=row>Macrophage</th><td> 0.007488114</td><td> 0.10237376</td><td> 0.0919357</td><td>-0.05073764</td><td>-0.25669026</td><td>-0.105458514</td></tr>
	<tr><th scope=row>NK_Cell</th><td> 0.192857999</td><td> 0.11207807</td><td> 0.1811900</td><td> 0.02379789</td><td>-0.04145599</td><td> 0.077933986</td></tr>
	<tr><th scope=row>T_Cell</th><td> 0.009909457</td><td> 0.08105488</td><td> 0.2876029</td><td> 0.17418257</td><td>-0.19789852</td><td> 0.026337183</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 9 × 4</caption>
<thead>
	<tr><th></th><th scope=col>Immune</th><th scope=col>SCLC</th><th scope=col>Stroma/Other</th><th scope=col>Unannotated</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B_Cell</th><td> 0.3610330</td><td>-0.25918967</td><td> 0.154580007</td><td> 0.046646311</td></tr>
	<tr><th scope=row>CD8_T_Cell</th><td> 0.2765727</td><td>-0.20887205</td><td> 0.131685863</td><td> 0.020076269</td></tr>
	<tr><th scope=row>DC</th><td> 0.2838786</td><td>-0.22784408</td><td>-0.045166181</td><td>-0.012589191</td></tr>
	<tr><th scope=row>Endothelial_Cell</th><td> 0.1609173</td><td>-0.21935814</td><td> 0.101081144</td><td>-0.001813333</td></tr>
	<tr><th scope=row>Epithelial_Cell</th><td>-0.1950830</td><td>-0.01894611</td><td>-0.145067048</td><td>-0.129368982</td></tr>
	<tr><th scope=row>Fibroblast</th><td> 0.4390432</td><td>-0.28873773</td><td> 0.234661787</td><td> 0.033637049</td></tr>
	<tr><th scope=row>Macrophage</th><td> 0.0919357</td><td>-0.26704387</td><td>-0.009817295</td><td>-0.105458514</td></tr>
	<tr><th scope=row>NK_Cell</th><td> 0.1811900</td><td>-0.04401943</td><td> 0.054651278</td><td> 0.077933986</td></tr>
	<tr><th scope=row>T_Cell</th><td> 0.2876029</td><td>-0.20898253</td><td> 0.097791399</td><td> 0.026337183</td></tr>
</tbody>
</table>



# Make the heatmaps


```R
options(repr.plot.width=10, repr.plot.height=7)

Labeled_Clusters_hm1 <- pheatmap(AB1_A2_Labeled_Clusters_mat,
         col = rev(brewer.pal(10, 'RdBu')),
         cluster_rows = T, cluster_cols = T,
         fontsize = 16,
         angle_col = 315,
         main='GSVA cell type enrichments vs Labeled Clusters')

Labeled_Clusters_hm1

png(file = "Output/16_Navin_Annotations_VA1_VA2/GSVA_vs_Labeled_Clusters.png",
    width = 1000,
    height = 600)
print(Labeled_Clusters_hm1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/GSVA_vs_Labeled_Clusters.pdf",
    width = 10,
    height = 6)
print(Labeled_Clusters_hm1)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_184_2.png)
    



```R
options(repr.plot.width=11, repr.plot.height=7)

Manual_Annotation_hm1 <- pheatmap(data_VA_Manual_Annotation_mat,
         col = rev(brewer.pal(10, 'RdBu')),
         cluster_rows = T, cluster_cols = T,
         fontsize = 16,
         angle_col = 315,
         main='GSVA cell type enrichments vs Manual Annotation')

Manual_Annotation_hm1

png(file = "Output/16_Navin_Annotations_VA1_VA2/GSVA_vs_Manual_Annotation.png",
    width = 1100,
    height = 600)
print(Manual_Annotation_hm1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/GSVA_vs_Manual_Annotation.pdf",
    width = 11,
    height = 6)
print(Manual_Annotation_hm1)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_185_2.png)
    



```R
options(repr.plot.width=11, repr.plot.height=7)

Manual_Annotation_simp_hm1 <- pheatmap(data_VA_Manual_Annotation_simp_mat,
         col = rev(brewer.pal(10, 'RdBu')),
         cluster_rows = T, cluster_cols = T,
         fontsize = 16,
         angle_col = 315,
         main='GSVA cell type enrichments vs Manual Annotation')

Manual_Annotation_simp_hm1

png(file = "Output/16_Navin_Annotations_VA1_VA2/GSVA_vs_Manual_Annotation_Simplified.png",
    width = 1100,
    height = 600)
print(Manual_Annotation_simp_hm1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/GSVA_vs_Manual_Annotation_Simlified.pdf",
    width = 11,
    height = 6)
print(Manual_Annotation_simp_hm1)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_186_2.png)
    



```R
options(repr.plot.width=9, repr.plot.height=7)

Manual_Annotation_Minimal_hm1 <- pheatmap(data_VA_Manual_Annotation_min_mat,
         col = rev(brewer.pal(10, 'RdBu')),
         cluster_rows = T, cluster_cols = T,
         fontsize = 16,
         angle_col = 315,
         main='GSVA cell type enrichments vs Manual Annotation')

Manual_Annotation_Minimal_hm1

png(file = "Output/16_Navin_Annotations_VA1_VA2/GSVA_vs_Manual_Annotation_Minimal.png",
    width = 900,
    height = 600)
print(Manual_Annotation_Minimal_hm1)
dev.off()

pdf(file = "Output/16_Navin_Annotations_VA1_VA2/GSVA_vs_Manual_Annotation_Minimal.pdf",
    width = 9,
    height = 6)
print(Manual_Annotation_Minimal_hm1)
dev.off()
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_187_2.png)
    



```R
sessionInfo()
```


    R version 4.5.0 (2025-04-11)
    Platform: x86_64-pc-linux-gnu
    Running under: Ubuntu 24.04.2 LTS
    
    Matrix products: default
    BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
    LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
    
    locale:
     [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
     [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
     [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
     [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
     [9] LC_ADDRESS=C               LC_TELEPHONE=C            
    [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
    
    time zone: Etc/UTC
    tzcode source: system (glibc)
    
    attached base packages:
    [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
    [8] methods   base     
    
    other attached packages:
     [1] scales_1.4.0         GSVA_2.2.0           RColorBrewer_1.1-3  
     [4] pheatmap_1.0.12      dittoSeq_1.20.0      lubridate_1.9.4     
     [7] forcats_1.0.0        stringr_1.5.1        dplyr_1.1.4         
    [10] purrr_1.0.4          readr_2.1.5          tidyr_1.3.1         
    [13] tibble_3.2.1         tidyverse_2.0.0      org.Hs.eg.db_3.21.0 
    [16] AnnotationDbi_1.70.0 IRanges_2.42.0       S4Vectors_0.46.0    
    [19] Biobase_2.68.0       BiocGenerics_0.54.0  generics_0.1.4      
    [22] ggplot2_3.5.2        doParallel_1.0.17    iterators_1.0.14    
    [25] foreach_1.5.2        Matrix_1.7-3         Seurat_5.3.0        
    [28] SeuratObject_5.1.0   sp_2.2-0            
    
    loaded via a namespace (and not attached):
      [1] RcppAnnoy_0.0.22            splines_4.5.0              
      [3] later_1.4.2                 pbdZMQ_0.3-14              
      [5] polyclip_1.10-7             graph_1.86.0               
      [7] XML_3.99-0.18               fastDummies_1.7.5          
      [9] lifecycle_1.0.4             globals_0.18.0             
     [11] lattice_0.22-7              MASS_7.3-65                
     [13] magrittr_2.0.3              plotly_4.10.4              
     [15] httpuv_1.6.16               sctransform_0.4.2          
     [17] spam_2.11-1                 spatstat.sparse_3.1-0      
     [19] reticulate_1.42.0           cowplot_1.1.3              
     [21] pbapply_1.7-2               DBI_1.2.3                  
     [23] abind_1.4-8                 Rtsne_0.17                 
     [25] GenomicRanges_1.60.0        GenomeInfoDbData_1.2.14    
     [27] ggrepel_0.9.6               irlba_2.3.5.1              
     [29] listenv_0.9.1               spatstat.utils_3.1-4       
     [31] goftest_1.2-3               RSpectra_0.16-2            
     [33] annotate_1.86.0             spatstat.random_3.4-1      
     [35] fitdistrplus_1.2-2          parallelly_1.44.0          
     [37] codetools_0.2-20            DelayedArray_0.34.1        
     [39] tidyselect_1.2.1            UCSC.utils_1.4.0           
     [41] farver_2.1.2                ScaledMatrix_1.16.0        
     [43] matrixStats_1.5.0           base64enc_0.1-3            
     [45] spatstat.explore_3.4-3      jsonlite_2.0.0             
     [47] progressr_0.15.1            ggridges_0.5.6             
     [49] survival_3.8-3              tools_4.5.0                
     [51] ica_1.0-3                   Rcpp_1.0.14                
     [53] glue_1.8.0                  gridExtra_2.3              
     [55] SparseArray_1.8.0           MatrixGenerics_1.20.0      
     [57] GenomeInfoDb_1.44.0         IRdisplay_1.1              
     [59] HDF5Array_1.36.0            withr_3.0.2                
     [61] fastmap_1.2.0               rhdf5filters_1.20.0        
     [63] rsvd_1.0.5                  digest_0.6.37              
     [65] timechange_0.3.0            R6_2.6.1                   
     [67] mime_0.13                   colorspace_2.1-1           
     [69] Cairo_1.6-2                 scattermore_1.2            
     [71] tensor_1.5                  spatstat.data_3.1-6        
     [73] RSQLite_2.3.11              h5mread_1.0.1              
     [75] data.table_1.17.4           httr_1.4.7                 
     [77] htmlwidgets_1.6.4           S4Arrays_1.8.0             
     [79] uwot_0.2.3                  pkgconfig_2.0.3            
     [81] gtable_0.3.6                blob_1.2.4                 
     [83] lmtest_0.9-40               SingleCellExperiment_1.30.1
     [85] XVector_0.48.0              htmltools_0.5.8.1          
     [87] dotCall64_1.2               GSEABase_1.70.0            
     [89] png_0.1-8                   SpatialExperiment_1.18.1   
     [91] spatstat.univar_3.1-3       rjson_0.2.23               
     [93] tzdb_0.5.0                  reshape2_1.4.4             
     [95] uuid_1.2-1                  nlme_3.1-168               
     [97] rhdf5_2.52.0                repr_1.1.7                 
     [99] zoo_1.8-14                  cachem_1.1.0               
    [101] KernSmooth_2.23-26          miniUI_0.1.2               
    [103] pillar_1.10.2               grid_4.5.0                 
    [105] vctrs_0.6.5                 RANN_2.6.2                 
    [107] promises_1.3.3              BiocSingular_1.24.0        
    [109] beachmat_2.24.0             xtable_1.8-4               
    [111] cluster_2.1.8.1             evaluate_1.0.3             
    [113] magick_2.8.6                cli_3.6.5                  
    [115] compiler_4.5.0              rlang_1.1.6                
    [117] crayon_1.5.3                future.apply_1.11.3        
    [119] labeling_0.4.3              plyr_1.8.9                 
    [121] stringi_1.8.7               viridisLite_0.4.2          
    [123] deldir_2.0-4                BiocParallel_1.42.0        
    [125] Biostrings_2.76.0           lazyeval_0.2.2             
    [127] spatstat.geom_3.4-1         IRkernel_1.3.2             
    [129] RcppHNSW_0.6.0              hms_1.1.3                  
    [131] patchwork_1.3.0             sparseMatrixStats_1.20.0   
    [133] bit64_4.6.0-1               future_1.49.0              
    [135] Rhdf5lib_1.30.0             KEGGREST_1.48.0            
    [137] shiny_1.10.0                SummarizedExperiment_1.38.1
    [139] ROCR_1.0-11                 igraph_2.1.4               
    [141] memoise_2.0.1               bit_4.6.0                  



```R

```


```R

```


```R

```
