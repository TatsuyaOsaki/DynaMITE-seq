# Plot enrichments vs our various labels/groups
# Use GSVA by Tissue Slice


```R
R.home()
```


'/usr/local/lib/R'



```R
library(Seurat, quietly = T)
library(Matrix, quietly = T)
library(doParallel, quietly = T)
library(ggplot2, quietly = T)
library(org.Hs.eg.db, quietly = T)
library(msigdbr, quietly = T)
library(tidyverse, quietly = T)
library(dittoSeq, quietly = T)
library(pheatmap, quietly = T)
library(RColorBrewer, quietly = T)
library(scales, quietly = T)
library(GGally, quietly = T)
library(ggbeeswarm)
library(paletteer)
library(scales)

options(repr.matrix.max.cols=1000, repr.matrix.max.rows=100)

```

    
    Attaching package: ‘SeuratObject’
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, t
    
    
    
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
    
    
    Welcome to Bioconductor
    
        Vignettes contain introductory material; view with
        'browseVignettes()'. To cite Bioconductor, see
        'citation("Biobase")', and for packages 'citation("pkgname")'.
    
    
    
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
    
    
    Registered S3 method overwritten by 'GGally':
      method from   
      +.gg   ggplot2
    


# load most up to date RDS object, and separate metadata csv too


```R
annots <- read.csv('Processing/Annotations_05_do_GSVA_by_Tissue_Slice_20250425.csv')
annots_melt <- read.csv('Processing/Annotations_melted_05_do_GSVA_by_Tissue_Slice_20250425.csv')
data <- readRDS('Processing/AB1_A2_Annotated_05_do_GSVA_by_Tissue_Slice_20250425.rds')
```


```R
# remove X column
annots <- annots %>% select(-c(X))
annots_melt <- annots_melt %>% select(-c(X))
```


```R
# set Barcode to rownames for annots
rownames(annots) <- annots$Barcode
```

# check data


```R
dim(annots)
annots %>% head(2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>57</li></ol>




<table class="dataframe">
<caption>A data.frame: 2 × 57</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.3752852</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.6186324</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td></tr>
</tbody>
</table>




```R
dim(annots_melt)
annots_melt %>% head(2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>291330</li><li>59</li></ol>




<table class="dataframe">
<caption>A data.frame: 2 × 59</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>GeneSet</th><th scope=col>GSVA_enrichment</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>-0.03843072</td><td>0.5135266</td><td>-0.007136273</td><td>0.4121959</td><td>0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.282326</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.1651271</td><td>0.0155587</td><td>-0.7552172</td><td>0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>All_Immune  </td><td>-0.05680737</td></tr>
	<tr><th scope=row>2</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>-0.03843072</td><td>0.5135266</td><td>-0.007136273</td><td>0.4121959</td><td>0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.282326</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.1651271</td><td>0.0155587</td><td>-0.7552172</td><td>0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>Angiogenesis</td><td>-0.25984750</td></tr>
</tbody>
</table>




```R
data
```


    An object of class Seurat 
    35796 features across 9711 samples within 2 assays 
    Active assay: Spatial (17943 features, 2000 variable features)
     3 layers present: data, counts, scale.data
     1 other assay present: SCT
     3 dimensional reductions calculated: pca, integrated.cca, umap
     3 spatial fields of view present: slice1 slice1.2 slice1.3



```R
data@meta.data[1:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 57</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.15846527</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.11080007</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.53207656</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.20627210</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.51946615</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.11225417</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1_AAACACCAATAACTGC-1</td><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td>4 </td><td>4 </td><td>VA1</td><td>0.907884</td><td>vasc_low </td><td> 9</td><td>Tumor_low </td><td>Tumor_low </td><td>1</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC+vasculature      </td><td>VA1</td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td>-0.29499671</td><td>-0.1592542</td><td>-0.6246499</td><td>-0.70918797</td><td>-0.3126876</td><td>-0.128532387</td><td> 0.54796712</td><td>-0.28738253</td><td> 0.05564459</td><td>-0.4488276</td><td>-0.3094928</td><td>-0.48102913</td><td>-0.23406975</td><td>-0.03504189</td><td>-0.1584708</td><td>-0.2611306</td><td>-0.3051610</td><td> 0.03743911</td><td> 0.3381445</td><td>-0.4111955</td><td>-0.3300223</td><td>-0.14181196</td><td>-0.19025652</td><td>-0.7206696</td><td>-0.21902693</td><td>-0.45491178</td><td>-0.3728931</td><td>-0.4603460</td><td>-0.23406975</td><td> 0.07776597</td></tr>
	<tr><th scope=row>VA1_AAACAGCTTTCAGAAG-1</th><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>5132</td><td>3219</td><td>round1</td><td>4797</td><td>3218</td><td>4 </td><td>4 </td><td>VA1</td><td>3.751279</td><td>vasc_high</td><td>14</td><td>Tumor_low </td><td>Tumor_low </td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low          </td><td>vasc_high</td><td>SCLC+vasculature      </td><td>VA1</td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td> 0.14865678</td><td> 0.6701079</td><td> 0.1810290</td><td>-0.68264193</td><td> 0.5112075</td><td> 0.351200349</td><td> 0.29886669</td><td> 0.53696284</td><td>-0.22259553</td><td> 0.2016089</td><td>-0.2723653</td><td> 0.26234173</td><td>-0.24403536</td><td>-0.02474443</td><td>-0.0993218</td><td> 0.5413331</td><td> 0.7253772</td><td>-0.13408735</td><td>-0.4164351</td><td> 0.1162332</td><td> 0.1290473</td><td> 0.15984342</td><td> 0.43378600</td><td> 0.6304333</td><td> 0.59860816</td><td>-0.09760237</td><td> 0.8719730</td><td> 0.9183918</td><td>-0.24403536</td><td>-0.15730040</td></tr>
</tbody>
</table>




```R
colnames(data@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Tissue_Slice'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li></ol>




```R
Assays(data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Spatial'</li><li>'SCT'</li></ol>



# unique pathways


```R
GeneSetList <- unique(annots_melt$GeneSet)
GeneSetList
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li></ol>




```R
celltypes_list <- c('All_Immune','T_Cell','CD8_T_Cell','T_Reg','NK_Cell','B_Cell','DC',
                    'Endothelial_Cell','LEC','Macrophage','M2_Macrophage','Macrophage','M1_Macrophage')
pathways_list <- c('Endothelial_Activation','Endothelial_Chemokines','Angiogenesis','GOBP_LEUKOCYTE_ADHESION_TO_VASC',
                   'MN4_EC_Phenotype','Upregulated_by_2_3_CGAMP','Upregulated_by_2_3_CGAMP_IFNb_OVERLAP','MN4_EC_Phenotype_Top30',
                   'AyersIFNG','Exhaustion','Proliferation','ICB_Targets','GOBP_LEUKOCYTE_MIGRATION','GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION')
```

# Recreating paper figures:
- Plot Immune Cell enrichments for Tumor MHCI Low vs High
    - NK cells, Pan T cells, CD8 T cells, B cells, Macrophages, Dendritic Cells
- Plot GOBP Leukocyte adhesion to vasculature, split by Tumor MHCI Low vs High
- Plot GOBP Positive regulation of Leukocyte migration, split by Tumor MHCI Low vs High
- Neighborhoods: Immune Cell enrichments for Tumor MHCI Low vs High
- Combine Endothelial Cell score (vasc_label) with MN4_EC label (called MITE vascular activation signature in the paper)
- Plot MITE (MN4_EC) enrichment, for Tumor MHCI Low vs High
- Plot Immune Cell enrichments, for MITE (MN4_EC) Low vs High
- Neighborhoods: Immune Cell enrichments for MITE (MN4_EC) Low vs High
- Supp Figure 3...
    - plots of Navin's manual annotations
    - cross-referencing barplot (Navin annotations vs clusters)
    - Heatmap of Immune Cell enrichments vs Navin Annotations
    - 
threshold_enrichment_at_median <- function(df, col) {
    curr_median <- median(df[[col]])
    new_col_name <- paste0(col,'_Label')
    #df[[new_col_name]] <- paste0(col,'_Low')
    #df[[new_col_name]][df[[col]] > curr_median] <- paste0(col,'_High')
    df[[new_col_name]] <- 'Low'
    df[[new_col_name]][df[[col]] > curr_median] <- 'High'
    return(df)
}# make a function that splits at the median per group (e.g. Tissue Slice)
threshold_enrichment_at_median <- function(df, col, group_col) {
  col <- enquo(col)
  group_col <- enquo(group_col)
  new_col_name <- paste0(rlang::as_name(col), "_Label")
  
  df %>%
    group_by(!!group_col) %>%
    mutate(
      median_val = median(!!col, na.rm = TRUE),
      !!new_col_name := ifelse(!!col > median_val, "High", "Low")
    ) %>%
    select(-median_val) %>%
    ungroup()
}# make a function that splits at the median per group (e.g. Tissue Slice)
threshold_enrichment_at_median <- function(df, col, group_col) {
  new_col_name <- paste0(col, "_Label")

    # get median per group and put in a new temporary column 
    df <- df %>%
    group_by(!!group_col) %>%
    mutate(median_val = median(.data[[col]], na.rm = TRUE)) %>%
    ungroup()

    print(df %>% select(c(col, group_col, median_val)) %>% head(2))
    print(table(df$median_val))
    
    # Threshold and label
    df[[new_col_name]] <- 'Low'
    df[[new_col_name]][df[[col]] > df$median_val] <- 'High'

    # Remove temporary median column
    group_col_name_remove <- paste0('\"', group_col, '\"')
    df <- df %>% select(-c(median_val, group_col_name_remove))

    # return df
    return(df)
}threshold_enrichment_at_median <- function(df, col, group_col) {
  new_col_name <- paste0(deparse(substitute(col)), "_Label")
  
  df %>%
    group_by({{ group_col }}) %>%
    mutate(
      median_val = median({{ col }}, na.rm = TRUE),
      !!new_col_name := ifelse({{ col }} > median_val, "High", "Low")
    ) %>%
    select(-median_val) %>%
    ungroup()
}# make a function that splits at the median per group (e.g. Tissue Slice)
threshold_enrichment_at_median <- function(df, col, group_col) {
  new_col_name <- paste0(col, "_Label")

    # get median per group and put in a new temporary column 
    df <- df %>%
    group_by({{ group_col }}) %>%
    mutate(median_val = median({{ col }}, na.rm = TRUE),
           !!new_col_name := ifelse({{ col }} > median_val, "High", "Low")) %>%
    ungroup()

    print(df %>% select(c(col, group_col, median_val)) %>% head(2))
    print(table(df$median_val))

    # Remove temporary median column
    group_col_name_remove <- paste0('\"', group_col, '\"')
    df <- df %>% select(-c(median_val, group_col_name_remove))

    # return df
    return(df)
}# make a function that splits at the median per group (e.g. Tissue Slice)
threshold_enrichment_at_median <- function(df, col, group_col) {
  new_col_name <- paste0(col, "_Label")

    # get median per group and put in a new temporary column 
    df <- df %>%
    group_by(!!group_col) %>%
    mutate(median_val = median(.data[[col]], na.rm = TRUE)) %>%
    ungroup()

    median_per_group <- df %>%
      group_by(!!group_col) %>%
      summarise(median_value = median(!!col, na.rm = TRUE))
    median_per_group <- as.data.frame(median_per_group)
    print(median_per_group)

    # merge with df
    df <- merge(df,median_per_group, by=group_col)
    
    print(df %>% select(c(col, group_col, median_value)) %>% head(2))
    
    # Threshold and label
    df[[new_col_name]] <- 'Low'
    df[[new_col_name]][df[[col]] > df$median_value] <- 'High'

    # Remove temporary median column
    group_col_name_remove <- paste0('\"', group_col, '\"')
    df <- df %>% select(-c(median_value, group_col_name_remove))

    # return df
    return(df)
}# make a function that splits at the median per group (e.g. Tissue Slice)
threshold_enrichment_at_median <- function(df, col, group_col) {
  new_col_name <- paste0(col, "_Label")

    median_per_group <- df %>%
      group_by({{group_col}}) %>%
      summarise(median_value = median({{col}}, na.rm = TRUE))
    median_per_group <- as.data.frame(median_per_group)
    print(median_per_group)

    # merge with df
    df <- merge(df,median_per_group, by=group_col)
    
    print(df %>% select(c(col, group_col, median_value)) %>% head(2))
    
    # Threshold and label
    df[[new_col_name]] <- 'Low'
    df[[new_col_name]][df[[col]] > df$median_value] <- 'High'

    # Remove temporary median column
    group_col_name_remove <- paste0('\"', group_col, '\"')
    df <- df %>% select(-c(median_value, group_col_name_remove))

    # return df
    return(df)
}
# Just find the medians per group manually outside of a function...


```R
# find medians
temp_median_testing <- annots
median_testing_result <- temp_median_testing %>%
  group_by(Tissue_Slice) %>%
  summarise(median_value = median(MN4_EC_Phenotype, na.rm = TRUE))
median_testing_result
```


<table class="dataframe">
<caption>A tibble: 3 × 2</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>median_value</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>-0.11315038</td></tr>
	<tr><td>VA2</td><td>-0.07112625</td></tr>
	<tr><td>VB1</td><td>-0.09340693</td></tr>
</tbody>
</table>




```R
# Merge medians with datasets
annots <- merge(annots, median_testing_result, by='Tissue_Slice')

temp_metdata <- merge(data@meta.data, median_testing_result, by='Tissue_Slice')
```


```R
# Threshold and apply labels
annots$MN4_EC_Phenotype_Label <- 'Low'
annots$MN4_EC_Phenotype_Label[annots$MN4_EC_Phenotype > annots$median_value] <- 'High'

temp_metdata$MN4_EC_Phenotype_Label <- 'Low'
temp_metdata$MN4_EC_Phenotype_Label[temp_metdata$MN4_EC_Phenotype > temp_metdata$median_value] <- 'High'
```


```R
# drop extra columns
annots <- annots %>% select(-c(median_value))
temp_metdata <- temp_metdata %>% select(-c(median_value))
```

# Let's annotate top 20% and bot 20% enrichment for MN4 EC Phenotype


```R
percentile_testing_result <- annots %>%
  group_by(Tissue_Slice) %>%
  summarise(
    p20_MN4_EC = quantile(MN4_EC_Phenotype, probs = 0.2, na.rm = TRUE),
    p80_MN4_EC = quantile(MN4_EC_Phenotype, probs = 0.8, na.rm = TRUE)
  )
percentile_testing_result
```


<table class="dataframe">
<caption>A tibble: 3 × 3</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>p20_MN4_EC</th><th scope=col>p80_MN4_EC</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>-0.2438415</td><td>0.02563792</td></tr>
	<tr><td>VA2</td><td>-0.3246353</td><td>0.17477807</td></tr>
	<tr><td>VB1</td><td>-0.2000513</td><td>0.03130858</td></tr>
</tbody>
</table>




```R
# Merge medians with datasets
annots <- merge(annots, percentile_testing_result, by='Tissue_Slice')

temp_metdata <- merge(temp_metdata, percentile_testing_result, by='Tissue_Slice')
```


```R
# Threshold and apply labels
annots$MN4_EC_Percentile_Label <- 'Medium'
annots$MN4_EC_Percentile_Label[annots$MN4_EC_Phenotype > annots$p80_MN4_EC] <- 'Above80Percentile'
annots$MN4_EC_Percentile_Label[annots$MN4_EC_Phenotype < annots$p20_MN4_EC] <- 'Below20Percentile'

temp_metdata$MN4_EC_Percentile_Label <- 'Low'
temp_metdata$MN4_EC_Percentile_Label[temp_metdata$MN4_EC_Phenotype > temp_metdata$p80_MN4_EC] <- 'Above80Percentile'
temp_metdata$MN4_EC_Percentile_Label[temp_metdata$MN4_EC_Phenotype < temp_metdata$p20_MN4_EC] <- 'Below20Percentile'
```


```R
# drop extra columns
annots <- annots %>% select(-c(p20_MN4_EC))
temp_metdata <- temp_metdata %>% select(-c(p20_MN4_EC))

annots <- annots %>% select(-c(p80_MN4_EC))
temp_metdata <- temp_metdata %>% select(-c(p80_MN4_EC))
```

# Need to set rownames


```R
annots <- annots %>% as.data.frame()
rownames(annots) <- annots$Barcode
```


```R
temp_metadata <- temp_metdata %>% as.data.frame()
rownames(temp_metadata) <- temp_metadata$Barcode
data@meta.data <- temp_metadata
```


```R
annots %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 59</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.3752852</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.6186324</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td></tr>
</tbody>
</table>




```R
data@meta.data %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 59</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.3752852</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.6186324</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td></tr>
</tbody>
</table>


# threshold annots df
annots <- threshold_enrichment_at_median(annots, 'MN4_EC_Phenotype', 'Tissue_Slice')
table(annots$MN4_EC_Phenotype_Label)
annots %>% head(2)

# threshold data seurat object
data@meta.data <- threshold_enrichment_at_median(data@meta.data, 'MN4_EC_Phenotype', 'Tissue_Slice')
table(data@meta.data$MN4_EC_Phenotype_Label)
data@meta.data %>% head(2)

```R
unique(annots$vasc_label)
unique(annots$vasc_label1)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_high'</li><li>'vasc_low'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_high'</li><li>'vasc_low'</li><li>'vasc_neg'</li></ol>



# Combine vasc_label and MN4_EC_Phenotype_Label


```R
annots$vasc_MN4_label <- 'vasc_low'
annots$vasc_MN4_label[annots$vasc_label == 'vasc_high' & annots$MN4_EC_Phenotype_Label == 'Low'] <- 'vasc_high_MN4_low'
annots$vasc_MN4_label[annots$vasc_label == 'vasc_high' & annots$MN4_EC_Phenotype_Label == 'High'] <- 'vasc_high_MN4_high'

data@meta.data$vasc_MN4_label <- 'vasc_low'
data@meta.data$vasc_MN4_label[data@meta.data$vasc_label == 'vasc_high' & data@meta.data$MN4_EC_Phenotype_Label == 'Low'] <- 'vasc_high_MN4_low'
data@meta.data$vasc_MN4_label[data@meta.data$vasc_label == 'vasc_high' & data@meta.data$MN4_EC_Phenotype_Label == 'High'] <- 'vasc_high_MN4_high'
```

# check counts...


```R
table(annots$vasc_MN4_label)
```


    
    vasc_high_MN4_high  vasc_high_MN4_low           vasc_low 
                  2410                759               6542 



```R
dplyr::count(data@meta.data, Tissue_Slice, vasc_label)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>vasc_label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>vasc_high</td><td> 771</td></tr>
	<tr><td>VA1</td><td>vasc_low </td><td>2035</td></tr>
	<tr><td>VA2</td><td>vasc_high</td><td>1733</td></tr>
	<tr><td>VA2</td><td>vasc_low </td><td>2583</td></tr>
	<tr><td>VB1</td><td>vasc_high</td><td> 665</td></tr>
	<tr><td>VB1</td><td>vasc_low </td><td>1924</td></tr>
</tbody>
</table>




```R
dplyr::count(data@meta.data, Tissue_Slice, MN4_EC_Phenotype_Label)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>High</td><td>1403</td></tr>
	<tr><td>VA1</td><td>Low </td><td>1403</td></tr>
	<tr><td>VA2</td><td>High</td><td>2158</td></tr>
	<tr><td>VA2</td><td>Low </td><td>2158</td></tr>
	<tr><td>VB1</td><td>High</td><td>1294</td></tr>
	<tr><td>VB1</td><td>Low </td><td>1295</td></tr>
</tbody>
</table>




```R
dplyr::count(data@meta.data, Tissue_Slice, vasc_MN4_label)
```


<table class="dataframe">
<caption>A data.frame: 9 × 3</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>vasc_MN4_label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>vasc_high_MN4_high</td><td> 511</td></tr>
	<tr><td>VA1</td><td>vasc_high_MN4_low </td><td> 260</td></tr>
	<tr><td>VA1</td><td>vasc_low          </td><td>2035</td></tr>
	<tr><td>VA2</td><td>vasc_high_MN4_high</td><td>1473</td></tr>
	<tr><td>VA2</td><td>vasc_high_MN4_low </td><td> 260</td></tr>
	<tr><td>VA2</td><td>vasc_low          </td><td>2583</td></tr>
	<tr><td>VB1</td><td>vasc_high_MN4_high</td><td> 426</td></tr>
	<tr><td>VB1</td><td>vasc_high_MN4_low </td><td> 239</td></tr>
	<tr><td>VB1</td><td>vasc_low          </td><td>1924</td></tr>
</tbody>
</table>



# Combine vasc_label and MN4_EC_Percentile_Label


```R
annots$vasc_MN4_percentile_label <- 'vasc_low'
annots$vasc_MN4_percentile_label[annots$vasc_label == 'vasc_high' & annots$MN4_EC_Percentile_Label == 'Below20Percentile'] <- 'vasc_high_MN4_below20'
annots$vasc_MN4_percentile_label[annots$vasc_label == 'vasc_high' & annots$MN4_EC_Percentile_Label == 'Above80Percentile'] <- 'vasc_high_MN4_upper80'

data@meta.data$vasc_MN4_percentile_label <- 'vasc_low'
data@meta.data$vasc_MN4_percentile_label[data@meta.data$vasc_label == 'vasc_high' & data@meta.data$MN4_EC_Percentile_Label == 'Below20Percentile'] <- 'vasc_high_MN4_below20'
data@meta.data$vasc_MN4_percentile_label[data@meta.data$vasc_label == 'vasc_high' & data@meta.data$MN4_EC_Percentile_Label == 'Above80Percentile'] <- 'vasc_high_MN4_upper80'
```

# check counts...


```R
table(annots$vasc_MN4_percentile_label)
```


    
    vasc_high_MN4_below20 vasc_high_MN4_upper80              vasc_low 
                      170                  1179                  8362 



```R
dplyr::count(data@meta.data, Tissue_Slice, vasc_label)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>vasc_label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>vasc_high</td><td> 771</td></tr>
	<tr><td>VA1</td><td>vasc_low </td><td>2035</td></tr>
	<tr><td>VA2</td><td>vasc_high</td><td>1733</td></tr>
	<tr><td>VA2</td><td>vasc_low </td><td>2583</td></tr>
	<tr><td>VB1</td><td>vasc_high</td><td> 665</td></tr>
	<tr><td>VB1</td><td>vasc_low </td><td>1924</td></tr>
</tbody>
</table>




```R
dplyr::count(data@meta.data, Tissue_Slice, MN4_EC_Phenotype_Label)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>High</td><td>1403</td></tr>
	<tr><td>VA1</td><td>Low </td><td>1403</td></tr>
	<tr><td>VA2</td><td>High</td><td>2158</td></tr>
	<tr><td>VA2</td><td>Low </td><td>2158</td></tr>
	<tr><td>VB1</td><td>High</td><td>1294</td></tr>
	<tr><td>VB1</td><td>Low </td><td>1295</td></tr>
</tbody>
</table>




```R
dplyr::count(data@meta.data, Tissue_Slice, vasc_MN4_percentile_label)
```


<table class="dataframe">
<caption>A data.frame: 9 × 3</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>vasc_high_MN4_below20</td><td>  86</td></tr>
	<tr><td>VA1</td><td>vasc_high_MN4_upper80</td><td> 267</td></tr>
	<tr><td>VA1</td><td>vasc_low             </td><td>2453</td></tr>
	<tr><td>VA2</td><td>vasc_high_MN4_below20</td><td>  11</td></tr>
	<tr><td>VA2</td><td>vasc_high_MN4_upper80</td><td> 707</td></tr>
	<tr><td>VA2</td><td>vasc_low             </td><td>3598</td></tr>
	<tr><td>VB1</td><td>vasc_high_MN4_below20</td><td>  73</td></tr>
	<tr><td>VB1</td><td>vasc_high_MN4_upper80</td><td> 205</td></tr>
	<tr><td>VB1</td><td>vasc_low             </td><td>2311</td></tr>
</tbody>
</table>



# Save updated enrichment df and seurat object


```R
write.csv(annots, 'Processing/Annotations_GSVA_by_Tissue_Slice_06_Label_MN4_20250425.csv')
saveRDS(data, 'Processing/AB1_A2_Annotated_GSVA_by_Tissue_Slice_06_Label_MN4_20250425.rds')
```

# Plot cell type enrichments vs vasc_MN4_label


```R
colnames(annots)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li></ol>




```R
annots_melt_immune <- annots %>% pivot_longer(c('Macrophage','M1_Macrophage','M2_Macrophage','T_Cell','CD8_T_Cell','NK_Cell','B_Cell','DC','Endothelial_Cell'), names_to = "Cell_Type_Label", values_to = "Cell_Type_Enrichment")

```


```R
annots_melt_immune[1:3,]
```


<table class="dataframe">
<caption>A tibble: 3 × 54</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Proliferation</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Macrophage   </td><td>-0.16512709</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M1_Macrophage</td><td>-0.06733032</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M2_Macrophage</td><td>-0.28232601</td></tr>
</tbody>
</table>




```R
dplyr::count(annots_melt_immune, orig.ident)
```


<table class="dataframe">
<caption>A tibble: 3 × 2</caption>
<thead>
	<tr><th scope=col>orig.ident</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>25254</td></tr>
	<tr><td>VA2</td><td>38844</td></tr>
	<tr><td>VB1</td><td>23301</td></tr>
</tbody>
</table>



# Plot PCA of enrichments
BiocManager::install("PCAtools")

```R
#GeneSetList
library(PCAtools)
```

    Loading required package: ggrepel
    
    
    Attaching package: ‘PCAtools’
    
    
    The following objects are masked from ‘package:stats’:
    
        biplot, screeplot
    
    



```R
annots %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 61</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.3752852</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.6186324</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td></tr>
</tbody>
</table>




```R
# set up pca...
annots_pre_pca <- annots %>% select(all_of(GeneSetList)) %>% t()
dim(annots_pre_pca)
annots_pre_pca %>% head(2)

annots_meta_pre_pca <- annots %>% select(-any_of(GeneSetList))
dim(annots_meta_pre_pca)
annots_meta_pre_pca %>% head(2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>30</li><li>9711</li></ol>




<table class="dataframe">
<caption>A matrix: 2 × 9711 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>VA1_AAACAACGAATAGTTC-1</th><th scope=col>VA1_AAACAAGTATCTCCCA-1</th><th scope=col>VA1_AAACAATCTACTAGCA-1</th><th scope=col>VA1_AAACACCAATAACTGC-1</th><th scope=col>VA1_AAACAGCTTTCAGAAG-1</th><th scope=col>VA1_AAACAGGGTCTATATT-1</th><th scope=col>VA1_AAACCCGAACGAAATC-1</th><th scope=col>VA1_AAACCGGAAATGTTAA-1</th><th scope=col>VA1_AAACCGGGTAGGTACC-1</th><th scope=col>VA1_AAACCGTTCGTCCAGG-1</th><th scope=col>VA1_AAACCTCATGAAGTTG-1</th><th scope=col>VA1_AAACGAGACGGTTGAT-1</th><th scope=col>VA1_AAACGGGTTGGTATCC-1</th><th scope=col>VA1_AAACTTGCAAACGTAT-1</th><th scope=col>VA1_AAAGACTGGGCGCTTT-1</th><th scope=col>VA1_AAAGCTTGCCTACATA-1</th><th scope=col>VA1_AAAGGCTACGGACCAT-1</th><th scope=col>VA1_AAAGGCTCTCGCGCCG-1</th><th scope=col>VA1_AAAGGGATGTAGCAAG-1</th><th scope=col>VA1_AAAGGGCAGCTTGAAT-1</th><th scope=col>VA1_AAAGTAGCATTGCTCA-1</th><th scope=col>VA1_AAAGTCGACCCTCAGT-1</th><th scope=col>VA1_AAAGTGTGATTTATCT-1</th><th scope=col>VA1_AAAGTTGACTCCCGTA-1</th><th scope=col>VA1_AAATAAGGTAGTGCCC-1</th><th scope=col>VA1_AAATAGCTTAGACTTT-1</th><th scope=col>VA1_AAATAGGGTGCTATTG-1</th><th scope=col>VA1_AAATCGCGGAAGGAGT-1</th><th scope=col>VA1_AAATCGTGTACCACAA-1</th><th scope=col>VA1_AAATCTAGCCCTGCTA-1</th><th scope=col>VA1_AAATGCTCGTTACGTT-1</th><th scope=col>VA1_AAATGGCCCGTGCCCT-1</th><th scope=col>VA1_AAATGGTCAATGTGCC-1</th><th scope=col>VA1_AAATGTATCTTATCCC-1</th><th scope=col>VA1_AAATTAACGGGTAGCT-1</th><th scope=col>VA1_AAATTAATAAGCGCGA-1</th><th scope=col>VA1_AAATTACCTATCGATG-1</th><th scope=col>VA1_AAATTCCAGGTCCAAA-1</th><th scope=col>VA1_AAATTGATAGTCCTTT-1</th><th scope=col>VA1_AAATTTGCGGGTGTGG-1</th><th scope=col>VA1_AACAACTGGTAGTTGC-1</th><th scope=col>VA1_AACAATACATTGTCGA-1</th><th scope=col>VA1_AACAATTACTCTACGC-1</th><th scope=col>VA1_AACACACGCTCGCCGC-1</th><th scope=col>VA1_AACACGAGACGCGGCC-1</th><th scope=col>VA1_AACACGCGGCCGCGAA-1</th><th scope=col>VA1_AACAGCTGTGTGGCAA-1</th><th scope=col>VA1_AACAGGATGGGCCGCG-1</th><th scope=col>VA1_AACAGGTAGTATGGAT-1</th><th scope=col>VA1_AACATATCAACTGGTG-1</th><th scope=col>VA1_AACATCGATACGTCTA-1</th><th scope=col>VA1_AACCAAGACTTCTCTG-1</th><th scope=col>VA1_AACCATGGGATCGCTA-1</th><th scope=col>VA1_AACCCAGAGACGGAGA-1</th><th scope=col>VA1_AACCCATCCCATGATC-1</th><th scope=col>VA1_AACCCGAGCAGAATCG-1</th><th scope=col>VA1_AACCCTACTGTCAATA-1</th><th scope=col>VA1_AACCGCTAAGGGATGC-1</th><th scope=col>VA1_AACCGTTGTGTTTGCT-1</th><th scope=col>VA1_AACCTGTCACGGAATT-1</th><th scope=col>VA1_AACGAAAGTCGTCCCA-1</th><th scope=col>VA1_AACGCATGATCTGGGT-1</th><th scope=col>VA1_AACGCGAACGGCAACA-1</th><th scope=col>VA1_AACGCTGTTGCTGAAA-1</th><th scope=col>VA1_AACGGACGTACGTATA-1</th><th scope=col>VA1_AACGGCCATCTCCGGT-1</th><th scope=col>VA1_AACGTACTGTGGGTAC-1</th><th scope=col>VA1_AACGTAGTCTACCCAT-1</th><th scope=col>VA1_AACGTCAGACTAGTGG-1</th><th scope=col>VA1_AACGTTATCAGCACCT-1</th><th scope=col>VA1_AACTACCCGTTTGTCA-1</th><th scope=col>VA1_AACTAGGCTTGGGTGT-1</th><th scope=col>VA1_AACTCAAGTTAATTGC-1</th><th scope=col>VA1_AACTCCAGAGCGTGTT-1</th><th scope=col>VA1_AACTCGATAAACACGT-1</th><th scope=col>VA1_AACTCGATGGCGCAGT-1</th><th scope=col>VA1_AACTCTCAATAGAGCG-1</th><th scope=col>VA1_AACTCTCAGTGTGCTC-1</th><th scope=col>VA1_AACTGAGTTATACTGA-1</th><th scope=col>VA1_AACTTCTGCGTCTATC-1</th><th scope=col>VA1_AACTTGCCCGTATGCA-1</th><th scope=col>VA1_AACTTGCGTTCTCGCG-1</th><th scope=col>VA1_AACTTTACGGGAGCTT-1</th><th scope=col>VA1_AACTTTCTCGATCATG-1</th><th scope=col>VA1_AAGAAGGATCAGTTAG-1</th><th scope=col>VA1_AAGACCCAACTGAACA-1</th><th scope=col>VA1_AAGACTAACCCGTTGT-1</th><th scope=col>VA1_AAGAGATGAATCGGTA-1</th><th scope=col>VA1_AAGAGCTCTTTATCGG-1</th><th scope=col>VA1_AAGAGGATGTACGCGA-1</th><th scope=col>VA1_AAGAGGCATGGATCGC-1</th><th scope=col>VA1_AAGAGGCCCTTTGGAA-1</th><th scope=col>VA1_AAGATGGCACCGGACC-1</th><th scope=col>VA1_AAGCACCCTGCGTATC-1</th><th scope=col>VA1_AAGCCGAAGCGGTTTA-1</th><th scope=col>VA1_AAGCGTCCCTCATCGA-1</th><th scope=col>VA1_AAGCTCTTTCATGGTG-1</th><th scope=col>VA1_AAGGAGCGGTTGGTGC-1</th><th scope=col>VA1_AAGGATCGATCGCTTG-1</th><th scope=col>VA1_AAGGCGCGTAAAGCTT-1</th><th scope=col>VA1_AAGGCTGTGCTCATCG-1</th><th scope=col>VA1_AAGGGACTATGCATTC-1</th><th scope=col>VA1_AAGGTGATAAACCAGC-1</th><th scope=col>VA1_AAGTAAGCTTCCAAAC-1</th><th scope=col>VA1_AAGTAGAAGACCGGGT-1</th><th scope=col>VA1_AAGTAGTGACGCGAGG-1</th><th scope=col>VA1_AAGTCAATTGTCGTCA-1</th><th scope=col>VA1_AAGTCTAGTAGCTGCC-1</th><th scope=col>VA1_AAGTCTTCTGTGGCCT-1</th><th scope=col>VA1_AAGTGACGACCGAATT-1</th><th scope=col>VA1_AAGTGAGTCGGGTTTA-1</th><th scope=col>VA1_AAGTGCAAAGGTAGAC-1</th><th scope=col>VA1_AAGTGCCTTGACTGTA-1</th><th scope=col>VA1_AAGTGCGTTAGAATCT-1</th><th scope=col>VA1_AAGTGTTTGGAGACGG-1</th><th scope=col>VA1_AAGTTCACTCCAAGCT-1</th><th scope=col>VA1_AAGTTCAGTCTGCGTA-1</th><th scope=col>VA1_AAGTTCGGCCAACAGG-1</th><th scope=col>VA1_AATAACAACGCTCGGC-1</th><th scope=col>VA1_AATAACACTAGAACAA-1</th><th scope=col>VA1_AATAATCTTCGTATCG-1</th><th scope=col>VA1_AATACAATGTTTCAGG-1</th><th scope=col>VA1_AATAGAACAGAGTGGC-1</th><th scope=col>VA1_AATAGAATCTGTTTCA-1</th><th scope=col>VA1_AATAGCTACCGCGTGC-1</th><th scope=col>VA1_AATAGTCCGTCCCGAC-1</th><th scope=col>VA1_AATAGTCGCGAGTCGG-1</th><th scope=col>VA1_AATATCAAGGTCGGAT-1</th><th scope=col>VA1_AATATCCTAGCAAACT-1</th><th scope=col>VA1_AATATCGAATCAATGC-1</th><th scope=col>VA1_AATATCGAGGGTTCTC-1</th><th scope=col>VA1_AATATTGGAGTATTGA-1</th><th scope=col>VA1_AATCAGACTGCAGGAC-1</th><th scope=col>VA1_AATCATGTAAAGACTC-1</th><th scope=col>VA1_AATCCATGCAAGGGTG-1</th><th scope=col>VA1_AATCGCCTCAGCGCCA-1</th><th scope=col>VA1_AATCGCGCAGAGGACT-1</th><th scope=col>VA1_AATCTAGGTTTACTTG-1</th><th scope=col>VA1_AATCTATGCCGGAGCC-1</th><th scope=col>VA1_AATCTCTACTGTGGTT-1</th><th scope=col>VA1_AATCTGCGTTGGGACG-1</th><th scope=col>VA1_AATCTGGCTTTCTAGT-1</th><th scope=col>VA1_AATGACAGCAATGTCT-1</th><th scope=col>VA1_AATGACGTAGGATGTC-1</th><th scope=col>VA1_AATGACTGTCAGCCGG-1</th><th scope=col>VA1_AATGAGTTCGCATATG-1</th><th scope=col>VA1_AATGATGCGACTCCTG-1</th><th scope=col>VA1_AATGCAACCGGGTACC-1</th><th scope=col>VA1_AATGGTCCACCGTTCA-1</th><th scope=col>VA1_AATGGTTCTCACAAGC-1</th><th scope=col>VA1_AATGTGCCCGAGGTGT-1</th><th scope=col>VA1_AATTAAAGGTCGGCGT-1</th><th scope=col>VA1_AATTAACGGATTTCCA-1</th><th scope=col>VA1_AATTACTCGTACGCTC-1</th><th scope=col>VA1_AATTAGCGCTGCAGCG-1</th><th scope=col>VA1_AATTCCAACTTGGTGA-1</th><th scope=col>VA1_AATTCTAGAGTTAGGC-1</th><th scope=col>VA1_AATTGAACGCTCTGGT-1</th><th scope=col>VA1_ACAAAGAAGGTAGGCC-1</th><th scope=col>VA1_ACAAATGGTAGTGTTT-1</th><th scope=col>VA1_ACAACGGTCCCTGCGA-1</th><th scope=col>VA1_ACAAGCAGTGCCTAGC-1</th><th scope=col>VA1_ACAAGCTATATGGAAG-1</th><th scope=col>VA1_ACAAGGACAAGAGGTT-1</th><th scope=col>VA1_ACAATAGTCGTACGTT-1</th><th scope=col>VA1_ACAATCCATTTAAACC-1</th><th scope=col>VA1_ACAATGAATACGGAGA-1</th><th scope=col>VA1_ACAATTGTGTCTCTTT-1</th><th scope=col>VA1_ACAATTTAGGAGGCTC-1</th><th scope=col>VA1_ACACAAAGACGGGTGG-1</th><th scope=col>VA1_ACACACCAGGACCAGT-1</th><th scope=col>VA1_ACACATGATCAAATCT-1</th><th scope=col>VA1_ACACATTGACGCAACA-1</th><th scope=col>VA1_ACACCACATAATTAGC-1</th><th scope=col>VA1_ACACCCAGCATGCAGC-1</th><th scope=col>VA1_ACACCCGAGAAATCCG-1</th><th scope=col>VA1_ACACCGAGCGCTCTTT-1</th><th scope=col>VA1_ACACCTTAAGTAGGGC-1</th><th scope=col>VA1_ACACCTTACTACTTGC-1</th><th scope=col>VA1_ACACGTAGGCCACAAG-1</th><th scope=col>VA1_ACACTCCAATGTCACT-1</th><th scope=col>VA1_ACACTGATCAAGGTGT-1</th><th scope=col>VA1_ACACTGGGACAGTCGT-1</th><th scope=col>VA1_ACAGAACTGAGAACAA-1</th><th scope=col>VA1_ACAGCATAGAGCCAGT-1</th><th scope=col>VA1_ACAGCGACATTCTCAT-1</th><th scope=col>VA1_ACAGGAGGCGCAGCCG-1</th><th scope=col>VA1_ACAGGCTTGCCCGACT-1</th><th scope=col>VA1_ACAGGTGGAGGTGAGG-1</th><th scope=col>VA1_ACAGGTGTGTTGTTGC-1</th><th scope=col>VA1_ACAGTAATACAACTTG-1</th><th scope=col>VA1_ACATAAGTCGTGGTGA-1</th><th scope=col>VA1_ACATCAGCTGGGACGC-1</th><th scope=col>VA1_ACATCCCGGCCATACG-1</th><th scope=col>VA1_ACATCCTGGTAACTGT-1</th><th scope=col>VA1_ACATCGATCGTTTACC-1</th><th scope=col>VA1_ACATCGGTCAGCCGCG-1</th><th scope=col>VA1_ACATCGTATGCAATGG-1</th><th scope=col>VA1_ACATGGCGCCAAAGTA-1</th><th scope=col>VA1_ACATGGCTCAATTTAG-1</th><th scope=col>VA1_ACATTAGTTTATATCC-1</th><th scope=col>VA1_ACCAAACACCCAGCGA-1</th><th scope=col>VA1_ACCAACACCACACACT-1</th><th scope=col>VA1_ACCAACGCTTATTTAT-1</th><th scope=col>VA1_ACCAATATGCAAGTTA-1</th><th scope=col>VA1_ACCACAAGTTTCTATC-1</th><th scope=col>VA1_ACCACCAATGTAACAA-1</th><th scope=col>VA1_ACCACCCTCTCTTCTA-1</th><th scope=col>VA1_ACCACTGTTCAAGAAG-1</th><th scope=col>VA1_ACCAGACCATAACAAC-1</th><th scope=col>VA1_ACCAGTGCGGGAGACG-1</th><th scope=col>VA1_ACCATATCCGCAATAA-1</th><th scope=col>VA1_ACCATCCGCCAACTAG-1</th><th scope=col>VA1_ACCATCGTATATGGTA-1</th><th scope=col>VA1_ACCATTAAGGGTGTCA-1</th><th scope=col>VA1_ACCCACCTACATGCTC-1</th><th scope=col>VA1_ACCCATCTTGAGGGTA-1</th><th scope=col>VA1_ACCCATTTGTCCCTCT-1</th><th scope=col>VA1_ACCCGAGCGAAATTAC-1</th><th scope=col>VA1_ACCCGGAAACTCCCAG-1</th><th scope=col>VA1_ACCCGGAGCACCACAA-1</th><th scope=col>VA1_ACCCGGATGACGCATC-1</th><th scope=col>VA1_ACCCGGTTACACTTCC-1</th><th scope=col>VA1_ACCCGTGTCATCAGTA-1</th><th scope=col>VA1_ACCCTATAGGACTGAG-1</th><th scope=col>VA1_ACCCTATGCCATATCG-1</th><th scope=col>VA1_ACCCTCCCTTGCTATT-1</th><th scope=col>VA1_ACCCTTCATCTGCGAA-1</th><th scope=col>VA1_ACCGAAAGGGCCCTGC-1</th><th scope=col>VA1_ACCGAAGAGTCTGGTT-1</th><th scope=col>VA1_ACCGACACATCTCCCA-1</th><th scope=col>VA1_ACCGATAGGCATAACC-1</th><th scope=col>VA1_ACCGATATTTAATCAT-1</th><th scope=col>VA1_ACCGATGGTAGCATCG-1</th><th scope=col>VA1_ACCGCGGTGGAAGTCG-1</th><th scope=col>VA1_ACCGGGCCTTTGTTGA-1</th><th scope=col>VA1_ACCGGTCAGGTACACC-1</th><th scope=col>VA1_ACCGTGACCACGTGGG-1</th><th scope=col>VA1_ACCTAAGTACCTTTCA-1</th><th scope=col>VA1_ACCTAATCGACTTCCT-1</th><th scope=col>VA1_ACCTACAGTATGTGGT-1</th><th scope=col>VA1_ACCTCCGCCCTCGCTG-1</th><th scope=col>VA1_ACCTCCGTTATTCACC-1</th><th scope=col>VA1_ACCTGCGTGTCATGTT-1</th><th scope=col>VA1_ACGAAATGGGCGGCAC-1</th><th scope=col>VA1_ACGACCCATGAGTTGC-1</th><th scope=col>VA1_ACGAGAACCCATCACG-1</th><th scope=col>VA1_ACGAGGATACCACTCT-1</th><th scope=col>VA1_ACGAGGTTTACAACGT-1</th><th scope=col>VA1_ACGAGTACGGATGCCC-1</th><th scope=col>VA1_ACGAGTCGCCGGCGTT-1</th><th scope=col>VA1_ACGATACCTATCCTGA-1</th><th scope=col>VA1_ACGATCATACATAGAG-1</th><th scope=col>VA1_ACGATCATCTTGTAAA-1</th><th scope=col>VA1_ACGATGCATATGTTAT-1</th><th scope=col>VA1_ACGCAAACTAATAGAT-1</th><th scope=col>VA1_ACGCAATCACTACAGC-1</th><th scope=col>VA1_ACGCATACGTTTACTA-1</th><th scope=col>VA1_ACGCATTCGTGAGTAC-1</th><th scope=col>VA1_ACGCCTGACACGCGCT-1</th><th scope=col>VA1_ACGCGAAGTCAGACGA-1</th><th scope=col>VA1_ACGCGCTACACAGGGT-1</th><th scope=col>VA1_ACGCGGGCCAAGGACA-1</th><th scope=col>VA1_ACGCGTTTCTTAAGAG-1</th><th scope=col>VA1_ACGCTAGTGATACACT-1</th><th scope=col>VA1_ACGCTGTGAGGCGTAG-1</th><th scope=col>VA1_ACGGAATTTAGCAAAT-1</th><th scope=col>VA1_ACGGACTCTCAAAGCG-1</th><th scope=col>VA1_ACGGAGCGCAAATTAC-1</th><th scope=col>VA1_ACGGATGGTGCGGATA-1</th><th scope=col>VA1_ACGGCACTTGCTTGGG-1</th><th scope=col>VA1_ACGGCCAACATGGACT-1</th><th scope=col>VA1_ACGGCGGGTTGCCCTG-1</th><th scope=col>VA1_ACGGGAGTGTCGGCCC-1</th><th scope=col>VA1_ACGGGTCATGTGACTT-1</th><th scope=col>VA1_ACGTATTACTCCGATC-1</th><th scope=col>VA1_ACGTCTCGTTCCGGGA-1</th><th scope=col>VA1_ACGTGACAAAGTAAGT-1</th><th scope=col>VA1_ACGTTAGATTTGCCCG-1</th><th scope=col>VA1_ACGTTATTGGTCACTC-1</th><th scope=col>VA1_ACGTTCGCAATCAATT-1</th><th scope=col>VA1_ACGTTCGTTCAGGAAA-1</th><th scope=col>VA1_ACTACCAGCTCTCTGG-1</th><th scope=col>VA1_ACTACGCGTTAGAATT-1</th><th scope=col>VA1_ACTAGTTGCGATCGTC-1</th><th scope=col>VA1_ACTATCCAGGGCATGG-1</th><th scope=col>VA1_ACTATCGCCGGCTAAA-1</th><th scope=col>VA1_ACTATCTGCCCGCGTA-1</th><th scope=col>VA1_ACTATTTCCGGGCCCA-1</th><th scope=col>VA1_ACTCAACGAATGTATT-1</th><th scope=col>VA1_ACTCAAGTGCAAGGCT-1</th><th scope=col>VA1_ACTCAGACCTGCTTCT-1</th><th scope=col>VA1_ACTCCCATTCCTAAAG-1</th><th scope=col>VA1_ACTCCCGAATTCGTTT-1</th><th scope=col>VA1_ACTCCCGTAGACTAGG-1</th><th scope=col>VA1_ACTCCCTAATGCTAAA-1</th><th scope=col>VA1_ACTCGATGTATTTCAT-1</th><th scope=col>VA1_ACTCGTAACCCGTCCT-1</th><th scope=col>VA1_ACTCGTCAGTAATCCC-1</th><th scope=col>VA1_ACTCTAAACCTGGGAT-1</th><th scope=col>VA1_ACTCTCTGACTTAGGT-1</th><th scope=col>VA1_ACTCTTCAGCTCCCGC-1</th><th scope=col>VA1_ACTGAAACGCCGTTAG-1</th><th scope=col>VA1_ACTGAATGGCGAAAGT-1</th><th scope=col>VA1_ACTGCGGACACACCGT-1</th><th scope=col>VA1_ACTGCTCGGAAGGATG-1</th><th scope=col>VA1_ACTGTAGCACTTTGGA-1</th><th scope=col>VA1_ACTGTATACGCGAGCA-1</th><th scope=col>VA1_ACTGTCCAGGATTATA-1</th><th scope=col>VA1_ACTGTGCTAGTAGATC-1</th><th scope=col>VA1_ACTTACGCATCCACGC-1</th><th scope=col>VA1_ACTTAGTACGACAAGA-1</th><th scope=col>VA1_ACTTATACTTACCCGG-1</th><th scope=col>VA1_ACTTATCTGATCTATA-1</th><th scope=col>VA1_ACTTATTAGGATCGGT-1</th><th scope=col>VA1_ACTTATTTATGTGCCA-1</th><th scope=col>VA1_ACTTCCAGTGGAAGCT-1</th><th scope=col>VA1_ACTTCCATGCGGGACA-1</th><th scope=col>VA1_ACTTCGCTAGCGAGTG-1</th><th scope=col>VA1_ACTTGGGACCCGGTGG-1</th><th scope=col>VA1_ACTTGTAGTCCCTTCA-1</th><th scope=col>VA1_ACTTGTGGATGGAACG-1</th><th scope=col>VA1_ACTTTACCCTCATGAA-1</th><th scope=col>VA1_ACTTTCCTATAGCTTC-1</th><th scope=col>VA1_ACTTTGACTGCATCCT-1</th><th scope=col>VA1_ACTTTGTCGACGCACT-1</th><th scope=col>VA1_AGAAATTATGACTCGC-1</th><th scope=col>VA1_AGAAGAGCGCCGTTCC-1</th><th scope=col>VA1_AGAAGGTACACTTCAC-1</th><th scope=col>VA1_AGAAGGTTGCCGAATT-1</th><th scope=col>VA1_AGAAGTGATTCGTGAT-1</th><th scope=col>VA1_AGAATAAATCTTCAGG-1</th><th scope=col>VA1_AGAATTATGGATTCGA-1</th><th scope=col>VA1_AGAATTGTTTGACATA-1</th><th scope=col>VA1_AGACAGCTCAGAATCC-1</th><th scope=col>VA1_AGACAGGCATCTCAGC-1</th><th scope=col>VA1_AGACCAAACCACACCT-1</th><th scope=col>VA1_AGACCATGGGATACAA-1</th><th scope=col>VA1_AGACCCGCCCTCCTCG-1</th><th scope=col>VA1_AGACCGCTCCGCGGTT-1</th><th scope=col>VA1_AGACCGGGAAACCCTG-1</th><th scope=col>VA1_AGACGAAGTGCCGGTC-1</th><th scope=col>VA1_AGACGACGATGCCGCT-1</th><th scope=col>VA1_AGACGCCCACTTCGCC-1</th><th scope=col>VA1_AGACTAGCCTTCCAGA-1</th><th scope=col>VA1_AGACTGTTACCGGGTC-1</th><th scope=col>VA1_AGAGAAGGAGTACAAT-1</th><th scope=col>VA1_AGAGATCTCTAAAGCG-1</th><th scope=col>VA1_AGAGCAGTTATGAGAC-1</th><th scope=col>VA1_AGAGCCGCCGAGATTT-1</th><th scope=col>VA1_AGAGCGTACAAGCTCG-1</th><th scope=col>VA1_AGAGCTACGAAAGCAT-1</th><th scope=col>VA1_AGAGGCTTCGGAAACC-1</th><th scope=col>VA1_AGAGTAAACTTCACTA-1</th><th scope=col>VA1_AGAGTATAGTGTTACG-1</th><th scope=col>VA1_AGATAACTTCAGGGCC-1</th><th scope=col>VA1_AGATACCGGTGTTCAC-1</th><th scope=col>VA1_AGATACGACTTCATAT-1</th><th scope=col>VA1_AGATACTCAAGATCGA-1</th><th scope=col>VA1_AGATGATGGAGTCTGG-1</th><th scope=col>VA1_AGATGCAAGACGTGCA-1</th><th scope=col>VA1_AGATGCATCCTGTGTC-1</th><th scope=col>VA1_AGATGCTATAACGAGC-1</th><th scope=col>VA1_AGATTATAGGACGTTT-1</th><th scope=col>VA1_AGATTCACAACCGATA-1</th><th scope=col>VA1_AGCAAAGGCCGCTAGT-1</th><th scope=col>VA1_AGCAACATATCTTATT-1</th><th scope=col>VA1_AGCAACCGAAAGTAAT-1</th><th scope=col>VA1_AGCACTACCGGCCTGT-1</th><th scope=col>VA1_AGCACTACCTCACCAG-1</th><th scope=col>VA1_AGCACTTAAGGACGCC-1</th><th scope=col>VA1_AGCAGCCAGATGAATA-1</th><th scope=col>VA1_AGCATATCAATATGCT-1</th><th scope=col>VA1_AGCATCATTTCGAAAG-1</th><th scope=col>VA1_AGCATCGTCGATAATT-1</th><th scope=col>VA1_AGCCCATACATGTAAG-1</th><th scope=col>VA1_AGCCCGGTAGCCTGTA-1</th><th scope=col>VA1_AGCCCTTCTAATCCGA-1</th><th scope=col>VA1_AGCCCTTGGACATCCC-1</th><th scope=col>VA1_AGCCGCAAATTCAAAT-1</th><th scope=col>VA1_AGCCGCTTGATTAGCG-1</th><th scope=col>VA1_AGCCGTGGCTAAATGT-1</th><th scope=col>VA1_AGCCTAATACCCACGT-1</th><th scope=col>VA1_AGCCTTAAAGCGGAAG-1</th><th scope=col>VA1_AGCGACAGGAACGGTC-1</th><th scope=col>VA1_AGCGACCAACGATATT-1</th><th scope=col>VA1_AGCGCATAATGAATCG-1</th><th scope=col>VA1_AGCGGACACTTCGTAG-1</th><th scope=col>VA1_AGCGGCGGTTAGCGGT-1</th><th scope=col>VA1_AGCGGGAAGGGTCCAT-1</th><th scope=col>VA1_AGCGTACGAGAGCTAG-1</th><th scope=col>VA1_AGCGTAGCGCTAGACC-1</th><th scope=col>VA1_AGCGTCTGAACCCGCA-1</th><th scope=col>VA1_AGCTAACAAGCAATGT-1</th><th scope=col>VA1_AGCTAGAAGCAGAAGT-1</th><th scope=col>VA1_AGCTCCATATATGTTC-1</th><th scope=col>VA1_AGCTCCTTCGCACATC-1</th><th scope=col>VA1_AGCTCTAGACGTTCCA-1</th><th scope=col>VA1_AGCTCTTCGTAACCTT-1</th><th scope=col>VA1_AGCTCTTTACTCAGTT-1</th><th scope=col>VA1_AGCTGTAACCTCAATC-1</th><th scope=col>VA1_AGCTTGATCTTAACTT-1</th><th scope=col>VA1_AGGAAAGCCTCTGATG-1</th><th scope=col>VA1_AGGAAGCTGTCCGCCG-1</th><th scope=col>VA1_AGGACAGTCGAATCCC-1</th><th scope=col>VA1_AGGACATCGCACGTCG-1</th><th scope=col>VA1_AGGACGACCCATTAGA-1</th><th scope=col>VA1_AGGACGCTCGATGTTG-1</th><th scope=col>VA1_AGGACTTATAGGAGAA-1</th><th scope=col>VA1_AGGAGAGTCTGGCTAC-1</th><th scope=col>VA1_AGGAGGCCTTCGCGCG-1</th><th scope=col>VA1_AGGATAAAGTCGGGAT-1</th><th scope=col>VA1_AGGATATAGGGATTTA-1</th><th scope=col>VA1_AGGATCACGCGATCTG-1</th><th scope=col>VA1_AGGATTGCTTACGACA-1</th><th scope=col>VA1_AGGCAAAGAGGAATCA-1</th><th scope=col>VA1_AGGCAATACGGAGGAC-1</th><th scope=col>VA1_AGGCACGTGACTGTCC-1</th><th scope=col>VA1_AGGCAGGGAGCGTACT-1</th><th scope=col>VA1_AGGCATTGTCGTAGGG-1</th><th scope=col>VA1_AGGCCACATTGGTTAC-1</th><th scope=col>VA1_AGGCCACCCGTTATGA-1</th><th scope=col>VA1_AGGCCCAGTGACTGGT-1</th><th scope=col>VA1_AGGCCCTAGAACGCCA-1</th><th scope=col>VA1_AGGCGATAACTGGCGT-1</th><th scope=col>VA1_AGGCGGTTTGTCCCGC-1</th><th scope=col>VA1_AGGCGTCTATGGACGG-1</th><th scope=col>VA1_AGGCTATGGTTAGCTT-1</th><th scope=col>VA1_AGGCTTCCCGAAGAAG-1</th><th scope=col>VA1_AGGCTTGCTAGACACC-1</th><th scope=col>VA1_AGGGACCGGCTGCGTT-1</th><th scope=col>VA1_AGGGACTCTACGCGAC-1</th><th scope=col>VA1_AGGGAGACATACTTCG-1</th><th scope=col>VA1_AGGGCGAGCAGCTGAT-1</th><th scope=col>VA1_AGGGCGTGATCGGCTA-1</th><th scope=col>VA1_AGGGCTGCAGTTACAG-1</th><th scope=col>VA1_AGGGTCAGTAACCCTA-1</th><th scope=col>VA1_AGGGTCGATGCGAACT-1</th><th scope=col>VA1_AGGGTGCTCTCGAGGG-1</th><th scope=col>VA1_AGGGTGGATAGTGCAT-1</th><th scope=col>VA1_AGGGTTTAGTTCGGGA-1</th><th scope=col>VA1_AGGTAACCTCCTATTC-1</th><th scope=col>VA1_AGGTAGGTACAAAGCT-1</th><th scope=col>VA1_AGGTATAATTGATAGT-1</th><th scope=col>VA1_AGGTCAGGTGAGAGTG-1</th><th scope=col>VA1_AGGTCGCCACTTCGGT-1</th><th scope=col>VA1_AGGTCGCGGAGTTACT-1</th><th scope=col>VA1_AGGTGTATCGCCATGA-1</th><th scope=col>VA1_AGGTTACACCATGCCG-1</th><th scope=col>VA1_AGTACATCATTTATCA-1</th><th scope=col>VA1_AGTACGGCCCGTATCG-1</th><th scope=col>VA1_AGTACTCTTATGCCCA-1</th><th scope=col>VA1_AGTATAATACTAGGCA-1</th><th scope=col>VA1_AGTATCCATAATAACG-1</th><th scope=col>VA1_AGTATGCTGGAGACCA-1</th><th scope=col>VA1_AGTATTTGGCACGACC-1</th><th scope=col>VA1_AGTCAACACCACCATC-1</th><th scope=col>VA1_AGTCAAGATGACACTT-1</th><th scope=col>VA1_AGTCACTAGCTCTCGA-1</th><th scope=col>VA1_AGTCACTCCGCCTCAT-1</th><th scope=col>VA1_AGTCAGCCACCGCCTG-1</th><th scope=col>VA1_AGTCCAGCGGGTACGT-1</th><th scope=col>VA1_AGTCCATTGGCTGATG-1</th><th scope=col>VA1_AGTCCCTCGCAGAAAG-1</th><th scope=col>VA1_AGTCGACGGTCTCAAG-1</th><th scope=col>VA1_AGTCGGCCCAAACGAC-1</th><th scope=col>VA1_AGTCGGCTCAACTTTA-1</th><th scope=col>VA1_AGTCGGTTGCGTGAGA-1</th><th scope=col>VA1_AGTCGTCGACCACCAA-1</th><th scope=col>VA1_AGTCGTGGGCATTACG-1</th><th scope=col>VA1_AGTCTAAAGTATACTC-1</th><th scope=col>VA1_AGTCTCACAAGACTAC-1</th><th scope=col>VA1_AGTCTGGACATCCTTG-1</th><th scope=col>VA1_AGTGAACAAACTTCTC-1</th><th scope=col>VA1_AGTGAAGATGGTGTCC-1</th><th scope=col>VA1_AGTGACCTACTTTACG-1</th><th scope=col>VA1_AGTGAGACTTCCAGTA-1</th><th scope=col>VA1_AGTGAGCCTCGCCGCC-1</th><th scope=col>VA1_AGTGATAACCTGCGCG-1</th><th scope=col>VA1_AGTGATATGAGTAGTT-1</th><th scope=col>VA1_AGTGATTCAAGCAGGA-1</th><th scope=col>VA1_AGTGCACGCTTAAGAA-1</th><th scope=col>VA1_AGTGCGTAGCTCGTAA-1</th><th scope=col>VA1_AGTGCTTGCACGAATA-1</th><th scope=col>VA1_AGTGGCTCCGTCGGCC-1</th><th scope=col>VA1_AGTGGGAGTATACACG-1</th><th scope=col>VA1_AGTGGTGTTACCCGTG-1</th><th scope=col>VA1_AGTGGTTGCGTATAGG-1</th><th scope=col>VA1_AGTGTGCTAAGATCGC-1</th><th scope=col>VA1_AGTGTGGTCTATTGTG-1</th><th scope=col>VA1_AGTTAAACACTTGCGA-1</th><th scope=col>VA1_AGTTAAGTCAACCGCT-1</th><th scope=col>VA1_AGTTACCCTTAAGACT-1</th><th scope=col>VA1_AGTTACCGCACATGGT-1</th><th scope=col>VA1_AGTTCACCGGTTGGAC-1</th><th scope=col>VA1_AGTTCCTACAGAATTA-1</th><th scope=col>VA1_AGTTCTGCGTTGTATC-1</th><th scope=col>VA1_AGTTGACATCGGCTGG-1</th><th scope=col>VA1_AGTTGCTGACTGATAT-1</th><th scope=col>VA1_AGTTTATGTAAAGACA-1</th><th scope=col>⋯</th><th scope=col>VB1_TCAAAGAGCTATCTGT-1</th><th scope=col>VB1_TCAACAAAGATAATTC-1</th><th scope=col>VB1_TCAACACATTGGGTAA-1</th><th scope=col>VB1_TCAACATAGCGCCCTA-1</th><th scope=col>VB1_TCAACATCGACCGAGA-1</th><th scope=col>VB1_TCAACGCAGGAAATAA-1</th><th scope=col>VB1_TCAACGCGACCGGCAG-1</th><th scope=col>VB1_TCAACTAACGTATAAC-1</th><th scope=col>VB1_TCAAGCTGCCTTGAAA-1</th><th scope=col>VB1_TCAAGGTTACTACACC-1</th><th scope=col>VB1_TCAATACAATTGCTGC-1</th><th scope=col>VB1_TCAATACGCCGTCATG-1</th><th scope=col>VB1_TCACAGATCCTCAAAC-1</th><th scope=col>VB1_TCACAGCAAACTCGAA-1</th><th scope=col>VB1_TCACAGGAGAATAAGA-1</th><th scope=col>VB1_TCACAGGGAATCGCAA-1</th><th scope=col>VB1_TCACATCTTATCTGAT-1</th><th scope=col>VB1_TCACCGCTCGGCACTC-1</th><th scope=col>VB1_TCACGATGTCCGTGGA-1</th><th scope=col>VB1_TCACGCATTGTAGATC-1</th><th scope=col>VB1_TCACGTGCCCGATTCA-1</th><th scope=col>VB1_TCACTCGTGCAACGGC-1</th><th scope=col>VB1_TCACTCTTCGTCTGTC-1</th><th scope=col>VB1_TCAGAACCTCCACAGG-1</th><th scope=col>VB1_TCAGCAAATGCATCTC-1</th><th scope=col>VB1_TCAGCAGTAGGCCCTG-1</th><th scope=col>VB1_TCAGCTTGAGCTTTCG-1</th><th scope=col>VB1_TCAGGGCGACTTCCTT-1</th><th scope=col>VB1_TCAGGGCGCAAACTCG-1</th><th scope=col>VB1_TCAGGTTCTTTGAGAA-1</th><th scope=col>VB1_TCAGTACTGACCCGCG-1</th><th scope=col>VB1_TCAGTAGGGACTATAA-1</th><th scope=col>VB1_TCAGTGTATACGTCAT-1</th><th scope=col>VB1_TCATACTTACAGATCC-1</th><th scope=col>VB1_TCATATGAGCTTTGTT-1</th><th scope=col>VB1_TCATCGATGGTCCCAA-1</th><th scope=col>VB1_TCATGCAGGTTCTCAT-1</th><th scope=col>VB1_TCATTTAAGTCTCCGA-1</th><th scope=col>VB1_TCCAACTTTAAATTCT-1</th><th scope=col>VB1_TCCAATAAAGGCTACC-1</th><th scope=col>VB1_TCCAATGCGTCGCCGC-1</th><th scope=col>VB1_TCCACAATGGTTTACG-1</th><th scope=col>VB1_TCCACCTCTAGCCTTT-1</th><th scope=col>VB1_TCCACTTTATCTAGGT-1</th><th scope=col>VB1_TCCAGAGCACCGGTTC-1</th><th scope=col>VB1_TCCAGATGTACGCCAA-1</th><th scope=col>VB1_TCCAGCGCTATAAGCG-1</th><th scope=col>VB1_TCCAGGGTATATACGA-1</th><th scope=col>VB1_TCCATCAATACTAATC-1</th><th scope=col>VB1_TCCATTCCCACTAGAG-1</th><th scope=col>VB1_TCCCAAACAGACAACG-1</th><th scope=col>VB1_TCCCAAAGCCCTAAAT-1</th><th scope=col>VB1_TCCCACTCTCTTCCGG-1</th><th scope=col>VB1_TCCCGCCTATGTGCGT-1</th><th scope=col>VB1_TCCCGGGTGTGCTGCT-1</th><th scope=col>VB1_TCCCGTGTGCAATTTG-1</th><th scope=col>VB1_TCCCTGGCGTATTAAC-1</th><th scope=col>VB1_TCCCTGGCTCGCTGGA-1</th><th scope=col>VB1_TCCCTTAGATTACTCG-1</th><th scope=col>VB1_TCCGAACGTTGCCGCT-1</th><th scope=col>VB1_TCCGAACTTGGCTTAC-1</th><th scope=col>VB1_TCCGATGACTGAGCTC-1</th><th scope=col>VB1_TCCGATTACATTGCCG-1</th><th scope=col>VB1_TCCGCAGCCACCTAGC-1</th><th scope=col>VB1_TCCGCGGCAGCATCTG-1</th><th scope=col>VB1_TCCGCGGCCCAATGAA-1</th><th scope=col>VB1_TCCGCTGTCATCCCGG-1</th><th scope=col>VB1_TCCGGCCTAGCGTACA-1</th><th scope=col>VB1_TCCGTAACCACAATCC-1</th><th scope=col>VB1_TCCGTTAAGCTAATAT-1</th><th scope=col>VB1_TCCGTTTAGCCTTGAA-1</th><th scope=col>VB1_TCCTAACCGTCGGGCA-1</th><th scope=col>VB1_TCCTACATCCACGGCC-1</th><th scope=col>VB1_TCCTAGCAAAGAAGCT-1</th><th scope=col>VB1_TCCTATCATAGGTAAC-1</th><th scope=col>VB1_TCCTCCTAAGACATTC-1</th><th scope=col>VB1_TCCTCGGGCTGGGCTT-1</th><th scope=col>VB1_TCCTCTACGAGATGGC-1</th><th scope=col>VB1_TCCTCTCCAGTTGTCC-1</th><th scope=col>VB1_TCCTCTGGCCCATTAG-1</th><th scope=col>VB1_TCCTGCCAACTGGAGA-1</th><th scope=col>VB1_TCCTGCGTTGATACTC-1</th><th scope=col>VB1_TCCTTACGACGGTCCG-1</th><th scope=col>VB1_TCCTTCAATCCCTACG-1</th><th scope=col>VB1_TCCTTTAAATCCGCTT-1</th><th scope=col>VB1_TCGAAATTTAGGACCA-1</th><th scope=col>VB1_TCGAAGAACCGAGCAC-1</th><th scope=col>VB1_TCGAATATCCCGCAGG-1</th><th scope=col>VB1_TCGACAACTGAACCCG-1</th><th scope=col>VB1_TCGAGCCAGGCAGGCC-1</th><th scope=col>VB1_TCGAGTCTACGATTCG-1</th><th scope=col>VB1_TCGCAAAGATGCATTT-1</th><th scope=col>VB1_TCGCACCAGGAGGCAG-1</th><th scope=col>VB1_TCGCACTAACGTTTGT-1</th><th scope=col>VB1_TCGCATCCCTAAGTGT-1</th><th scope=col>VB1_TCGCCCACTGCGAGAG-1</th><th scope=col>VB1_TCGCCGAAGTTGCGTC-1</th><th scope=col>VB1_TCGCCGCACCGCGTGA-1</th><th scope=col>VB1_TCGCCGGAGAGTCTTA-1</th><th scope=col>VB1_TCGCCGGATGGGCAAG-1</th><th scope=col>VB1_TCGCCTCGACCTGTTG-1</th><th scope=col>VB1_TCGCGTCCAGAAGGTC-1</th><th scope=col>VB1_TCGCTCGGCACCAGCG-1</th><th scope=col>VB1_TCGCTGGGCGGATTGT-1</th><th scope=col>VB1_TCGCTTAATTACGAAG-1</th><th scope=col>VB1_TCGCTTTAAACGTTTG-1</th><th scope=col>VB1_TCGGAATGCGCTCTGA-1</th><th scope=col>VB1_TCGGACGCCCAGCCCA-1</th><th scope=col>VB1_TCGGAGAGTATCGGGA-1</th><th scope=col>VB1_TCGGAGTACATGAGTA-1</th><th scope=col>VB1_TCGGCGAACCCAAACC-1</th><th scope=col>VB1_TCGGCGTACTGCACAA-1</th><th scope=col>VB1_TCGGGAACGTGCCTAG-1</th><th scope=col>VB1_TCGGGATTCAAACATA-1</th><th scope=col>VB1_TCGGGCCGTCGTGGTA-1</th><th scope=col>VB1_TCGGTCCCGACAATAG-1</th><th scope=col>VB1_TCGGTCCCTGACTCCA-1</th><th scope=col>VB1_TCGTAAGACGACATTG-1</th><th scope=col>VB1_TCGTAAGCTCCGAGGA-1</th><th scope=col>VB1_TCGTATTACCCATTGC-1</th><th scope=col>VB1_TCGTCAAGTACGCGCA-1</th><th scope=col>VB1_TCGTCACACTGTTAGC-1</th><th scope=col>VB1_TCGTCTTAGGCGTTAA-1</th><th scope=col>VB1_TCGTGTACTATGGATG-1</th><th scope=col>VB1_TCGTGTCACGCTGACA-1</th><th scope=col>VB1_TCGTTAGGAGTCCCTA-1</th><th scope=col>VB1_TCGTTGACAGGGTCCC-1</th><th scope=col>VB1_TCGTTGCTATCCGGTC-1</th><th scope=col>VB1_TCTAAAGAACAGTCTC-1</th><th scope=col>VB1_TCTAACTGTATGTAAA-1</th><th scope=col>VB1_TCTAATACTGCCTCAG-1</th><th scope=col>VB1_TCTACGGGCTCAGTTG-1</th><th scope=col>VB1_TCTAGCAATCTCCGCC-1</th><th scope=col>VB1_TCTAGGTGGCGACGCT-1</th><th scope=col>VB1_TCTAGTGATATCGTGG-1</th><th scope=col>VB1_TCTATTACGCTGGCGA-1</th><th scope=col>VB1_TCTATTACTAGAGGAT-1</th><th scope=col>VB1_TCTCATGAGATAGGGT-1</th><th scope=col>VB1_TCTCCAACGTAGGTTA-1</th><th scope=col>VB1_TCTCCCTGGGCAGCGT-1</th><th scope=col>VB1_TCTCGACGTATCGCCG-1</th><th scope=col>VB1_TCTCGAGGAGGTTCGC-1</th><th scope=col>VB1_TCTCGTGTTACGAGGA-1</th><th scope=col>VB1_TCTCTTACCGCGAACC-1</th><th scope=col>VB1_TCTGAACTCGTACCCG-1</th><th scope=col>VB1_TCTGAAGCACGTGGTC-1</th><th scope=col>VB1_TCTGAATTCCGTACAA-1</th><th scope=col>VB1_TCTGACGGGCTAACCC-1</th><th scope=col>VB1_TCTGACTGTAATGGTT-1</th><th scope=col>VB1_TCTGATCGGGTGCTAG-1</th><th scope=col>VB1_TCTGATGTATTCTGTC-1</th><th scope=col>VB1_TCTGATTGGAAATGGA-1</th><th scope=col>VB1_TCTGCACCATTAGTAA-1</th><th scope=col>VB1_TCTGCAGATTCGAGTC-1</th><th scope=col>VB1_TCTGCCAGAAACTGCA-1</th><th scope=col>VB1_TCTGCGAATCGTTCGC-1</th><th scope=col>VB1_TCTGCGTCCGGTTTCT-1</th><th scope=col>VB1_TCTGCTTAGAACAAGC-1</th><th scope=col>VB1_TCTGGAGCGTAAGAGT-1</th><th scope=col>VB1_TCTGGCCGTTCAAGTT-1</th><th scope=col>VB1_TCTGGCGCAAGCCGGG-1</th><th scope=col>VB1_TCTGGGAACCTTTGAA-1</th><th scope=col>VB1_TCTGGGTAGCGCTCAT-1</th><th scope=col>VB1_TCTGTGACTGACCGTT-1</th><th scope=col>VB1_TCTGTGCCATCATAGT-1</th><th scope=col>VB1_TCTGTTACCCAGCATA-1</th><th scope=col>VB1_TCTTAACTCGGATGTA-1</th><th scope=col>VB1_TCTTACCGGAACTCGT-1</th><th scope=col>VB1_TCTTACGGCATCCGAC-1</th><th scope=col>VB1_TCTTACTTATGCCTCT-1</th><th scope=col>VB1_TCTTAGAGCTCCAATT-1</th><th scope=col>VB1_TCTTAGAGTGAACTCT-1</th><th scope=col>VB1_TCTTCCCATGGGCACA-1</th><th scope=col>VB1_TCTTCGAATAGACGTT-1</th><th scope=col>VB1_TCTTCGATACCAATAA-1</th><th scope=col>VB1_TCTTCTATAACCCGCC-1</th><th scope=col>VB1_TCTTGATGCGTAGCGA-1</th><th scope=col>VB1_TCTTTAAGACTATGAA-1</th><th scope=col>VB1_TCTTTAGAGTCTAACA-1</th><th scope=col>VB1_TCTTTAGCAGGCGAAC-1</th><th scope=col>VB1_TCTTTCATCCGTCCTT-1</th><th scope=col>VB1_TCTTTCGGCGGGACAC-1</th><th scope=col>VB1_TGAAAGGACCTGACTC-1</th><th scope=col>VB1_TGAACAAGCAGGGACT-1</th><th scope=col>VB1_TGAACACCCGAAGCAG-1</th><th scope=col>VB1_TGAACTGCTATGACTT-1</th><th scope=col>VB1_TGAAGAGCGGTCCTAG-1</th><th scope=col>VB1_TGAAGTAGCTTACGGA-1</th><th scope=col>VB1_TGAATGAGATACAGCA-1</th><th scope=col>VB1_TGAATGTCAGCCGGCC-1</th><th scope=col>VB1_TGAATTTCACTTGCCT-1</th><th scope=col>VB1_TGACAACGCATGTCGC-1</th><th scope=col>VB1_TGACACTTCTCTTTGC-1</th><th scope=col>VB1_TGACAGGACAAGTCCA-1</th><th scope=col>VB1_TGACATATATGACGAT-1</th><th scope=col>VB1_TGACATCGAGCGGACC-1</th><th scope=col>VB1_TGACCCACGTTAGACA-1</th><th scope=col>VB1_TGACGAATATTTCCCT-1</th><th scope=col>VB1_TGACGATGCACTAGAA-1</th><th scope=col>VB1_TGACTCCGAATCATAC-1</th><th scope=col>VB1_TGAGACGTACCTCTCA-1</th><th scope=col>VB1_TGAGAGATTTACCACG-1</th><th scope=col>VB1_TGAGCCATACAGTCTC-1</th><th scope=col>VB1_TGAGCTTTAATGACGC-1</th><th scope=col>VB1_TGAGGCATGTACTGTG-1</th><th scope=col>VB1_TGAGTGCCTCTTAAAT-1</th><th scope=col>VB1_TGAGTGGTCCGTGACG-1</th><th scope=col>VB1_TGAGTTAAAGACATTC-1</th><th scope=col>VB1_TGATACATTTAGCCGT-1</th><th scope=col>VB1_TGATAGCGGGATTCTA-1</th><th scope=col>VB1_TGATCAGGGAACTGCT-1</th><th scope=col>VB1_TGATCCCAGCATTAGT-1</th><th scope=col>VB1_TGATGGGACTAAGTCA-1</th><th scope=col>VB1_TGATTCAGGTCCCGCG-1</th><th scope=col>VB1_TGATTCCCGGTTACCT-1</th><th scope=col>VB1_TGATTCGTCTATCACT-1</th><th scope=col>VB1_TGATTCTGTCGCCGGT-1</th><th scope=col>VB1_TGATTTATTAGCTGTG-1</th><th scope=col>VB1_TGATTTCCTCCTGACG-1</th><th scope=col>VB1_TGCAATTTGGGCACGG-1</th><th scope=col>VB1_TGCACAGTGAAGTTAT-1</th><th scope=col>VB1_TGCACTATGTGAGTGC-1</th><th scope=col>VB1_TGCAGAACTATATCGT-1</th><th scope=col>VB1_TGCAGAGTACCGAGCA-1</th><th scope=col>VB1_TGCAGATCGTCCTAGG-1</th><th scope=col>VB1_TGCAGCTACGTACTTC-1</th><th scope=col>VB1_TGCAGTGAGGCTCGGG-1</th><th scope=col>VB1_TGCAGTTTCCTCCCAT-1</th><th scope=col>VB1_TGCATATGTCTGTCAC-1</th><th scope=col>VB1_TGCATGAGTAGATTCG-1</th><th scope=col>VB1_TGCATGTGGTAATCTA-1</th><th scope=col>VB1_TGCCAAAGTCAGACTT-1</th><th scope=col>VB1_TGCCAATGGGTACTCT-1</th><th scope=col>VB1_TGCCACACTAGAGGAA-1</th><th scope=col>VB1_TGCCATTACTAAAGAA-1</th><th scope=col>VB1_TGCCCGTACCGTTAAA-1</th><th scope=col>VB1_TGCCGAAAGCGTATTC-1</th><th scope=col>VB1_TGCCGATGTCATCAAT-1</th><th scope=col>VB1_TGCCGGATGTACGAGC-1</th><th scope=col>VB1_TGCCGTTCTTAATCGG-1</th><th scope=col>VB1_TGCCTAAATTTAATAG-1</th><th scope=col>VB1_TGCCTGACATCGGTCA-1</th><th scope=col>VB1_TGCCTTGGCCAGGCAA-1</th><th scope=col>VB1_TGCGAATATGGGATTT-1</th><th scope=col>VB1_TGCGACACCCTAGTGC-1</th><th scope=col>VB1_TGCGACGGCCGAACGT-1</th><th scope=col>VB1_TGCGAGAAACGTTACG-1</th><th scope=col>VB1_TGCGAGAATATTACCC-1</th><th scope=col>VB1_TGCGATGCTAATGGCT-1</th><th scope=col>VB1_TGCGCAAAGCATTTGG-1</th><th scope=col>VB1_TGCGCCGTTAATAACG-1</th><th scope=col>VB1_TGCGCGATTAACGGAG-1</th><th scope=col>VB1_TGCGTAAGAACCTGAT-1</th><th scope=col>VB1_TGCGTCATGACTGAGC-1</th><th scope=col>VB1_TGCGTGATTGGGTGTC-1</th><th scope=col>VB1_TGCGTTTGTTGACACT-1</th><th scope=col>VB1_TGCTAAGTGTCTATTT-1</th><th scope=col>VB1_TGCTCCACAGTTCTTA-1</th><th scope=col>VB1_TGCTCTGCCGGTTCAC-1</th><th scope=col>VB1_TGCTGTTGAAGAACTC-1</th><th scope=col>VB1_TGCTTAGAGAGAATGC-1</th><th scope=col>VB1_TGCTTGAAACCATGCA-1</th><th scope=col>VB1_TGGAAACGGAGTGAAC-1</th><th scope=col>VB1_TGGAACCACTGACACA-1</th><th scope=col>VB1_TGGAAGAAGGGAACGT-1</th><th scope=col>VB1_TGGAATATCCTTGACC-1</th><th scope=col>VB1_TGGAATTAGACGCTTT-1</th><th scope=col>VB1_TGGACACCGTTGCTTG-1</th><th scope=col>VB1_TGGACCACGGCGTTGA-1</th><th scope=col>VB1_TGGACGCAATCCAGCC-1</th><th scope=col>VB1_TGGACGTAGGCGAATC-1</th><th scope=col>VB1_TGGACTGTTCGCTCAA-1</th><th scope=col>VB1_TGGAGGGAAACACCTC-1</th><th scope=col>VB1_TGGAGTGATGCGATGA-1</th><th scope=col>VB1_TGGATAGAGTAACAGA-1</th><th scope=col>VB1_TGGCAAACTAAATTAC-1</th><th scope=col>VB1_TGGCAACTCGCGCGCC-1</th><th scope=col>VB1_TGGCAATGGGACGGCG-1</th><th scope=col>VB1_TGGCAGCAGTAATAGT-1</th><th scope=col>VB1_TGGCATGAAGTTTGGG-1</th><th scope=col>VB1_TGGCCAAACTGAAGTA-1</th><th scope=col>VB1_TGGCCAATTTGGTACT-1</th><th scope=col>VB1_TGGCCGTATATTGACC-1</th><th scope=col>VB1_TGGCGACTGCTCCAAA-1</th><th scope=col>VB1_TGGCTCTTGTCGCGTA-1</th><th scope=col>VB1_TGGCTTATGTATAATG-1</th><th scope=col>VB1_TGGCTTGTACAAGCTT-1</th><th scope=col>VB1_TGGGAAATGCCTTTCC-1</th><th scope=col>VB1_TGGGACCATTGGGAGT-1</th><th scope=col>VB1_TGGGCAATAGTTGGGT-1</th><th scope=col>VB1_TGGGCCACAAGAGCGC-1</th><th scope=col>VB1_TGGGCCCATACTAATT-1</th><th scope=col>VB1_TGGGCCTTGCCTGCAT-1</th><th scope=col>VB1_TGGGTAAGGTTCCCGC-1</th><th scope=col>VB1_TGGGTTCCCGGACGGA-1</th><th scope=col>VB1_TGGTAAGCAGGATTGA-1</th><th scope=col>VB1_TGGTAGAATATATGGG-1</th><th scope=col>VB1_TGGTCCCACGCTACGG-1</th><th scope=col>VB1_TGGTCGTGCAAGGCAA-1</th><th scope=col>VB1_TGGTCTAGCTTACATG-1</th><th scope=col>VB1_TGGTCTGTTGGGCGTA-1</th><th scope=col>VB1_TGGTGATCGTATTTGT-1</th><th scope=col>VB1_TGGTTAACTTACATTT-1</th><th scope=col>VB1_TGGTTATGCTTGCGGT-1</th><th scope=col>VB1_TGGTTGGAGGATCCTG-1</th><th scope=col>VB1_TGGTTTAAACGTGGGT-1</th><th scope=col>VB1_TGTAAGACTGATAAGA-1</th><th scope=col>VB1_TGTAATGACCACAATA-1</th><th scope=col>VB1_TGTACCCGACCCTAAT-1</th><th scope=col>VB1_TGTACGAACAAATCCG-1</th><th scope=col>VB1_TGTACTACTCTCACGG-1</th><th scope=col>VB1_TGTAGGAGAAATTTCC-1</th><th scope=col>VB1_TGTAGTGATCTATAAT-1</th><th scope=col>VB1_TGTATAACAGATCCTG-1</th><th scope=col>VB1_TGTATACGGATGATGA-1</th><th scope=col>VB1_TGTATGGCGCAGACAG-1</th><th scope=col>VB1_TGTCCGTGGCGCCTTT-1</th><th scope=col>VB1_TGTCCTAAGTCACCGC-1</th><th scope=col>VB1_TGTCGTTATCACATAT-1</th><th scope=col>VB1_TGTGAGACTAGCCCAA-1</th><th scope=col>VB1_TGTGCCAGAGGCAAAG-1</th><th scope=col>VB1_TGTGCCGGTGCCGGAA-1</th><th scope=col>VB1_TGTGCTTTACGTAAGA-1</th><th scope=col>VB1_TGTGGCAAAGCGTATG-1</th><th scope=col>VB1_TGTGGCTCCCACCAAC-1</th><th scope=col>VB1_TGTGGTTGCTAAAGCT-1</th><th scope=col>VB1_TGTGTCGAAGTCGAGG-1</th><th scope=col>VB1_TGTGTGACCATGAATC-1</th><th scope=col>VB1_TGTGTTCGTATCCAAG-1</th><th scope=col>VB1_TGTTATTGTATGTGGC-1</th><th scope=col>VB1_TGTTCATAAATGTGCT-1</th><th scope=col>VB1_TGTTCCGCTTCCATGA-1</th><th scope=col>VB1_TGTTCGTATTGCGGTG-1</th><th scope=col>VB1_TGTTCTTCCATTGACT-1</th><th scope=col>VB1_TGTTGCCGTCGTCCCA-1</th><th scope=col>VB1_TGTTGTCAAGAAGTCT-1</th><th scope=col>VB1_TGTTTCTGAAGCGTGC-1</th><th scope=col>VB1_TGTTTGAGATCGTCAG-1</th><th scope=col>VB1_TTAAACCTGGTTCCTT-1</th><th scope=col>VB1_TTAAACTCGAATTCAT-1</th><th scope=col>VB1_TTAACACCTCGAACAT-1</th><th scope=col>VB1_TTAACGAACAAGCAGT-1</th><th scope=col>VB1_TTAACGTTAAAGCCTG-1</th><th scope=col>VB1_TTAACTCACGCGTGGA-1</th><th scope=col>VB1_TTAACTTCAGGTAGGA-1</th><th scope=col>VB1_TTAAGATAGGATTGAC-1</th><th scope=col>VB1_TTAAGCGCCTGACCCA-1</th><th scope=col>VB1_TTAAGGCCCGTACTTT-1</th><th scope=col>VB1_TTAATCAGTACGTCAG-1</th><th scope=col>VB1_TTAATGTAGACCAGGT-1</th><th scope=col>VB1_TTAATGTGTTTGCAGG-1</th><th scope=col>VB1_TTACAGACCTAAATGA-1</th><th scope=col>VB1_TTACATGCCACAACTA-1</th><th scope=col>VB1_TTACCATTGATTACCC-1</th><th scope=col>VB1_TTACCCATTGCCGGGT-1</th><th scope=col>VB1_TTACCCTAACAGTCCT-1</th><th scope=col>VB1_TTACCCTAGGGATTGG-1</th><th scope=col>VB1_TTACCGCCTTAGGGAA-1</th><th scope=col>VB1_TTACGGATGGTTCGAG-1</th><th scope=col>VB1_TTACTAAAGGACTTTA-1</th><th scope=col>VB1_TTACTATCGGCTTCTC-1</th><th scope=col>VB1_TTACTCCGGCCGGGAA-1</th><th scope=col>VB1_TTACTCTGGTACGTAC-1</th><th scope=col>VB1_TTACTGTCTAGAGCTC-1</th><th scope=col>VB1_TTAGAATAAGGGTCGG-1</th><th scope=col>VB1_TTAGACACGATCGTTG-1</th><th scope=col>VB1_TTAGACGAGTCACCTC-1</th><th scope=col>VB1_TTAGAGGGATATACAG-1</th><th scope=col>VB1_TTAGAGTTTAGAAGGA-1</th><th scope=col>VB1_TTAGATAGGTCGATAC-1</th><th scope=col>VB1_TTAGCTCTGTAATCCG-1</th><th scope=col>VB1_TTAGCTGATTTGCCGT-1</th><th scope=col>VB1_TTAGGATGGGAGGGTA-1</th><th scope=col>VB1_TTAGGTCATAACCGAC-1</th><th scope=col>VB1_TTAGGTGTGACTGGTC-1</th><th scope=col>VB1_TTAGTAGGGCGGCGGG-1</th><th scope=col>VB1_TTAGTTATTCGTGGCA-1</th><th scope=col>VB1_TTAGTTCAAGTGTTCG-1</th><th scope=col>VB1_TTATATACGCTGTCAC-1</th><th scope=col>VB1_TTATCCAATCGAACTC-1</th><th scope=col>VB1_TTATCCGGGATCTATA-1</th><th scope=col>VB1_TTATCCTCAAGGAATA-1</th><th scope=col>VB1_TTATCGCCTGCGAAGC-1</th><th scope=col>VB1_TTATCTGACATTAGGA-1</th><th scope=col>VB1_TTATGAATGAAAGGGA-1</th><th scope=col>VB1_TTATGACAAACTGGAT-1</th><th scope=col>VB1_TTATGATCTTAACGAA-1</th><th scope=col>VB1_TTATGTTTGCGATAGA-1</th><th scope=col>VB1_TTATTAGAGCGTGTTC-1</th><th scope=col>VB1_TTATTTCATCCCAAAC-1</th><th scope=col>VB1_TTCAAAGTCTCTAGCC-1</th><th scope=col>VB1_TTCAATACTCTGAATC-1</th><th scope=col>VB1_TTCACGAAAGGATCAC-1</th><th scope=col>VB1_TTCAGAGTAACCTGAC-1</th><th scope=col>VB1_TTCAGCCCTGGTCCAC-1</th><th scope=col>VB1_TTCAGCTGGCGTGCCC-1</th><th scope=col>VB1_TTCAGGCGTCAAAGCC-1</th><th scope=col>VB1_TTCAGTTTGTGGCAGC-1</th><th scope=col>VB1_TTCATAGCCTTGTAAC-1</th><th scope=col>VB1_TTCATAGGGTGTCCAT-1</th><th scope=col>VB1_TTCCAATCAGAGCTAG-1</th><th scope=col>VB1_TTCCAATCTGGCTATC-1</th><th scope=col>VB1_TTCCACACAGATTTGA-1</th><th scope=col>VB1_TTCCAGACGAGATTTA-1</th><th scope=col>VB1_TTCCATCATGCGGTGA-1</th><th scope=col>VB1_TTCCCGGCGCCAATAG-1</th><th scope=col>VB1_TTCCGCAGAGAAATAT-1</th><th scope=col>VB1_TTCCGGCCTTGAGGCT-1</th><th scope=col>VB1_TTCCGGTATCTGTGTC-1</th><th scope=col>VB1_TTCCGGTTACCCACTT-1</th><th scope=col>VB1_TTCCTCGGACTAACCA-1</th><th scope=col>VB1_TTCCTCTGCCCGAATA-1</th><th scope=col>VB1_TTCGCACTGTACGACA-1</th><th scope=col>VB1_TTCGCATCCGGAAGCA-1</th><th scope=col>VB1_TTCGCCGCTCGCGCTA-1</th><th scope=col>VB1_TTCGCGCGCCATACGA-1</th><th scope=col>VB1_TTCGCTAGGAAGTTGT-1</th><th scope=col>VB1_TTCGGCAACCCGCTGA-1</th><th scope=col>VB1_TTCGGCTAGAGATGGT-1</th><th scope=col>VB1_TTCGGGACTAATCGCG-1</th><th scope=col>VB1_TTCGGGCGCTAGTCTT-1</th><th scope=col>VB1_TTCGGTGGAGACGCCC-1</th><th scope=col>VB1_TTCGTAATCCCAGCGG-1</th><th scope=col>VB1_TTCGTACTCCAGAACG-1</th><th scope=col>VB1_TTCGTTCAACGAAGTT-1</th><th scope=col>VB1_TTCTAACCGAAGCTTA-1</th><th scope=col>VB1_TTCTACTTGCGAGGGC-1</th><th scope=col>VB1_TTCTAGAAAGTCTTAT-1</th><th scope=col>VB1_TTCTATGCCTTTCGCA-1</th><th scope=col>VB1_TTCTCTTACAGGTGAT-1</th><th scope=col>VB1_TTCTGACCGGGCTCAA-1</th><th scope=col>VB1_TTCTGCGAGCGCCCTT-1</th><th scope=col>VB1_TTCTGCGGGTTAGCGG-1</th><th scope=col>VB1_TTCTGCTAGACTCCAA-1</th><th scope=col>VB1_TTCTGTTTCCTGTCGC-1</th><th scope=col>VB1_TTCTTATCCGCTGGGT-1</th><th scope=col>VB1_TTCTTGAGCCGCGCTA-1</th><th scope=col>VB1_TTCTTGCTAGCATCTC-1</th><th scope=col>VB1_TTCTTGGACGATCTGC-1</th><th scope=col>VB1_TTCTTGGAGTAATGAG-1</th><th scope=col>VB1_TTCTTGTGTCCATCAG-1</th><th scope=col>VB1_TTCTTTGGTCGCGACG-1</th><th scope=col>VB1_TTGAAAGGTGTAAAGG-1</th><th scope=col>VB1_TTGAACGAATCCTTTG-1</th><th scope=col>VB1_TTGAAGGATGGGCGCC-1</th><th scope=col>VB1_TTGAATATGGACTTTC-1</th><th scope=col>VB1_TTGAATCGTTGTATAA-1</th><th scope=col>VB1_TTGAATTCACGTGAGG-1</th><th scope=col>VB1_TTGACAGGAGCTCCCG-1</th><th scope=col>VB1_TTGACATGAACGTGGA-1</th><th scope=col>VB1_TTGACCATGTTCTCCG-1</th><th scope=col>VB1_TTGACCGTGTTAATGA-1</th><th scope=col>VB1_TTGACGCTCCATGAGC-1</th><th scope=col>VB1_TTGACTACCATATGGT-1</th><th scope=col>VB1_TTGACTATTGTCCGGC-1</th><th scope=col>VB1_TTGAGAAGTTTAGCAT-1</th><th scope=col>VB1_TTGAGAGTACTGCTAA-1</th><th scope=col>VB1_TTGAGCAGCCCACGGT-1</th><th scope=col>VB1_TTGAGCGCCACGTGAT-1</th><th scope=col>VB1_TTGAGTCCCGCTGCTG-1</th><th scope=col>VB1_TTGATCTAACTTTGTC-1</th><th scope=col>VB1_TTGATGTGTAGTCCCG-1</th><th scope=col>VB1_TTGATTAGCTGTTTCT-1</th><th scope=col>VB1_TTGATTATGCAGATGA-1</th><th scope=col>VB1_TTGCATGCTGATCACG-1</th><th scope=col>VB1_TTGCCAAGCAGAACCC-1</th><th scope=col>VB1_TTGCCATAGCCCGCTC-1</th><th scope=col>VB1_TTGCCCTGATCACGGG-1</th><th scope=col>VB1_TTGCCGCAGACCTACA-1</th><th scope=col>VB1_TTGCGCTCTCTCGCTT-1</th><th scope=col>VB1_TTGCGGCATCAGAAAG-1</th><th scope=col>VB1_TTGCGTAGTTTGAGGA-1</th><th scope=col>VB1_TTGCGTCGGCCAACCG-1</th><th scope=col>VB1_TTGCTGCACCTATCCA-1</th><th scope=col>VB1_TTGGAAGAATACAGTC-1</th><th scope=col>VB1_TTGGACATGTGGCTTA-1</th><th scope=col>VB1_TTGGACCATCTGGCAA-1</th><th scope=col>VB1_TTGGACCTATAACAGT-1</th><th scope=col>VB1_TTGGCCATCTTGCGCT-1</th><th scope=col>VB1_TTGGCCTAGAATTTCG-1</th><th scope=col>VB1_TTGGCGATCCGAATAT-1</th><th scope=col>VB1_TTGGGACGTAAGAGTT-1</th><th scope=col>VB1_TTGGGCGGCGGTTGCC-1</th><th scope=col>VB1_TTGGTCACACTCGTAA-1</th><th scope=col>VB1_TTGGTTCGCTCAAAGG-1</th><th scope=col>VB1_TTGGTTGCGGTGCGCG-1</th><th scope=col>VB1_TTGTAAGGACCTAAGT-1</th><th scope=col>VB1_TTGTAAGGCCAGTTGG-1</th><th scope=col>VB1_TTGTACACCTCGAACA-1</th><th scope=col>VB1_TTGTCACCGCGGTATC-1</th><th scope=col>VB1_TTGTGAACCTAATCCG-1</th><th scope=col>VB1_TTGTGAGGCATGACGC-1</th><th scope=col>VB1_TTGTGCGGAAGCGGAT-1</th><th scope=col>VB1_TTGTGGTAGGAGGGAT-1</th><th scope=col>VB1_TTGTGGTATAGGTATG-1</th><th scope=col>VB1_TTGTTCTAGATACGCT-1</th><th scope=col>VB1_TTGTTGTGTGTCAAGA-1</th><th scope=col>VB1_TTGTTTCACATCCAGG-1</th><th scope=col>VB1_TTGTTTCATTAGTCTA-1</th><th scope=col>VB1_TTGTTTCCATACAACT-1</th></tr>
</thead>
<tbody>
	<tr><th scope=row>All_Immune</th><td>-0.05680737</td><td>-0.4320402</td><td>-0.3122538</td><td>-0.2949967</td><td>0.1486568</td><td>-0.3237087</td><td>-0.09053114</td><td>-0.2726288931</td><td>-0.02927195</td><td>0.2272304</td><td>-0.30533013</td><td>-0.09327672</td><td> 0.23709902</td><td>-0.2045095</td><td>-0.4812392</td><td>-0.32047356</td><td>-0.35327980</td><td>0.4838026</td><td>-0.2073011</td><td>-0.03469408</td><td>-0.2743828</td><td>-0.3542925</td><td>-0.09232902</td><td>0.4385228</td><td>0.05718737</td><td>-0.01091304</td><td> 0.238856534</td><td> 0.29074356</td><td> 0.1345291</td><td>-0.21412090</td><td>-0.1148071</td><td>-0.14994788</td><td> 0.04422197</td><td>-0.2702715</td><td>-0.1906549</td><td>-0.48872274</td><td>-0.005195924</td><td> 0.07597530</td><td>0.4637180</td><td>0.006989429</td><td>-0.13753233</td><td>0.1763031</td><td>-0.1672269</td><td>0.4987293</td><td>-0.04657804</td><td>-0.1632526</td><td>-0.06586516</td><td>-0.00495179</td><td>-0.33884358</td><td>-0.4092121</td><td>0.1621810</td><td>-0.2234214</td><td>-0.11026609</td><td> 0.2153078</td><td>-0.133797391</td><td>-0.2135907</td><td>-0.3517899</td><td>0.1234798</td><td>-0.3243628</td><td>0.2002890</td><td>-0.04985154</td><td>-0.008122756</td><td> 0.2318517</td><td>0.02700073</td><td>0.5013315</td><td>0.2932873</td><td>0.5036274</td><td>-0.006819626</td><td>0.07229338</td><td>-0.1626247</td><td>-0.1727530</td><td>0.20111033</td><td>0.20023388</td><td>0.1647004</td><td>-0.1418727</td><td>-0.2083585</td><td>-0.2562545</td><td>-0.07275779</td><td>-0.3083723</td><td>0.5178613</td><td>-0.29783665</td><td>0.3413110</td><td>-0.2283483</td><td>-0.3099277</td><td> 0.2804614</td><td>0.09717966</td><td>-0.2605779</td><td>-0.1602114</td><td>-0.1265614</td><td>-0.3650418</td><td>-0.34978905</td><td>-0.03342306</td><td>-0.1495886</td><td> 0.2111308</td><td> 0.3691197</td><td>-0.36311415</td><td>0.2906142</td><td> 0.3043971</td><td>-0.2216314</td><td>0.09697745</td><td>0.26740324</td><td>-0.4198759</td><td>0.1913189</td><td> 0.1751910</td><td>0.03892444</td><td>-0.33401566</td><td>-0.1514095</td><td>-0.3094709</td><td>-0.2773055</td><td>0.01613779</td><td>-0.3357305</td><td>-0.4375212</td><td>-0.3993235</td><td>-0.10989747</td><td>0.3719291</td><td>0.14575928</td><td>-0.1511215</td><td>-0.05479382</td><td>-0.3760440</td><td>0.1822646</td><td>-0.09470229</td><td>-0.2519540</td><td>0.39750486</td><td>0.1740756365</td><td>-0.03440722</td><td> 0.38131700</td><td>-0.3725317</td><td>-0.1138695</td><td>-0.3014305</td><td> 0.266644379</td><td>-0.1413543</td><td>-0.09065147</td><td>0.2145880</td><td>-0.08463717</td><td>-0.10433462</td><td> 0.2987798</td><td>-0.3685592</td><td>0.1588956</td><td>-0.22060337</td><td>-0.3843525</td><td>-0.2447416</td><td>-0.3612612</td><td>-0.3515896</td><td> 0.01404673</td><td>0.1465583</td><td>0.2040514</td><td>-0.11202497</td><td>-0.1634989</td><td>-0.2980107</td><td>-0.06094555</td><td>0.1870253</td><td> 0.1963632</td><td>0.1323246</td><td>-0.2477115</td><td>0.22866705</td><td> 0.1579011</td><td>0.44849005</td><td>0.1531529</td><td>-0.4167674</td><td>-0.2525756</td><td>-0.125925773</td><td>0.1921504</td><td>-0.3907248</td><td>-0.08905283</td><td>-0.1999226</td><td>-0.15195161</td><td>-0.3366642</td><td>0.3090893</td><td>-0.2599865</td><td>-0.04987366</td><td>-0.10468954</td><td>-0.3218797</td><td>-0.3751124</td><td>-0.3066825</td><td>-0.2880849</td><td> 0.1106334</td><td>0.225502974</td><td>0.06321371</td><td>-0.40169368</td><td>-0.1541610</td><td>0.3181454</td><td>-0.312967</td><td>-0.1041995</td><td>-0.2309157</td><td> 0.02846773</td><td>0.01689617</td><td>0.09466947</td><td>0.34109461</td><td>-0.02670114</td><td>-0.4510007</td><td>-0.1484833</td><td>-0.3446877</td><td>-0.3924737</td><td> 0.06633051</td><td>0.09357155</td><td> 0.1496268</td><td>-0.03428584</td><td>-0.1374192</td><td>-0.2697457</td><td> 0.04532911</td><td>0.2551975</td><td>0.03371713</td><td>-0.11389825</td><td>-0.26719668</td><td>0.1881671</td><td> 0.04751356</td><td>-0.4429910</td><td>-0.1663062</td><td>-0.14316841</td><td>0.04707504</td><td>0.08857542</td><td>-0.1121169</td><td> 0.06640465</td><td>0.1297974</td><td>-0.3859977</td><td>-0.2088221</td><td>-0.2299124</td><td>-0.2121040</td><td>-0.2083052</td><td>0.1361060</td><td>0.1042946</td><td>-0.07619976</td><td>-0.2955252</td><td> 0.3427874</td><td>0.02147337</td><td>-0.02333834</td><td>0.2099141</td><td> 0.1628724</td><td>-0.1048633</td><td>-0.0574130</td><td> 0.08005199</td><td>-0.008591061</td><td>0.0359427</td><td>-0.3095405</td><td>0.004611903</td><td>-0.2315231</td><td>-0.3677786</td><td>-0.08670818</td><td>-0.4006289</td><td>-0.3013862</td><td>-0.3480400</td><td>-0.05320688</td><td>-0.2530860</td><td>-0.003050881</td><td>-0.09547891</td><td>-0.1322741</td><td> 0.0599449156</td><td>0.2960318</td><td>-0.05729111</td><td> 0.03415177</td><td> 0.35973398</td><td>-0.08032796</td><td>-0.25941455</td><td>-0.1402787</td><td> 0.308501258</td><td>-0.1530283</td><td> 0.2631800</td><td> 0.05489949</td><td>-0.1884276</td><td>0.03410911</td><td>-0.01978387</td><td>-0.30752147</td><td> 0.02157453</td><td>-0.2931447</td><td>0.12261848</td><td>-0.1969595</td><td>0.0366359</td><td>0.2512319</td><td>0.2127055</td><td>0.1654879</td><td>0.3324959</td><td> 0.12794105</td><td>-0.07336339</td><td>-0.01620169</td><td>-0.2931645</td><td>-0.07830054</td><td>-0.11741478</td><td>-0.4276686</td><td> 0.09273912</td><td>0.2058952</td><td>0.1336438</td><td>-0.01090172</td><td>-0.21244503</td><td>0.5352924</td><td> 0.1538791</td><td>-0.2934441</td><td> 0.41886770</td><td>-0.07718126</td><td>-0.06498525</td><td>0.2483888</td><td>-0.02361509</td><td> 0.08057546</td><td>-0.2515057</td><td>-0.339675881</td><td>-0.3242714</td><td> 0.10113235</td><td>-0.1865536</td><td>-0.2633231</td><td>-0.2889749</td><td>-0.3571617</td><td>0.05678134</td><td>-0.3095882</td><td>0.37605923</td><td>0.003838381</td><td>0.09577287</td><td> 0.3103853</td><td>0.2555810</td><td>-0.1524959</td><td>-0.1266704</td><td>0.2552103</td><td>0.5333913</td><td>0.4348059</td><td>-0.03721258</td><td> 0.02933871</td><td>0.2124162</td><td>0.12310057</td><td>-0.1342795</td><td> 0.03830715</td><td>-0.2487820</td><td>-0.005472197</td><td>-0.1567211</td><td> 0.1599959</td><td>-0.2051320</td><td>0.03464732</td><td>-0.2837288</td><td>0.11383795</td><td>-0.04501895</td><td>-0.2108367</td><td> 0.29823152</td><td>0.01334924</td><td>-0.44301757</td><td>-0.4559371</td><td>0.1923826</td><td>-0.1898745</td><td>-0.3292752</td><td>-0.1913707</td><td>-0.11669427</td><td>-0.3511826</td><td>-0.3672774</td><td>-0.043306726</td><td>-0.1853263</td><td>0.1438177</td><td>-0.1537074</td><td>0.2550648</td><td>-0.03524294</td><td>-0.3187641</td><td>-0.06438519</td><td>-0.2717979</td><td>-0.4751055</td><td>-0.02929902</td><td>0.2685364</td><td>-0.03586882</td><td>-0.2183029</td><td>0.2739386</td><td>-0.2886498</td><td>-0.09301603</td><td>0.1204538</td><td>0.4315001</td><td>0.3214277</td><td>0.1483025</td><td>-0.08633239</td><td>-0.3805011</td><td>-0.4267998</td><td>-0.2201739</td><td>-0.29056816</td><td>0.3843505</td><td>0.09264109</td><td>-0.06767316</td><td>-0.3407686</td><td>0.1030633</td><td>-0.2388037</td><td> 0.2308503</td><td>0.4010376</td><td>-0.4192040</td><td>-0.2471875</td><td> 0.2033636</td><td>-0.3755958</td><td> 0.09311522</td><td>-0.002573741</td><td>-0.3767622</td><td> 0.2209938</td><td>-0.4132854</td><td>-0.3872059</td><td>0.15557849</td><td>-0.1404971</td><td>0.1894501</td><td> 0.05837784</td><td>-0.1641069</td><td>-0.4304419</td><td> 0.3296439</td><td>0.3226914</td><td> 0.05398832</td><td>-0.08717873</td><td>-0.3040831</td><td>-0.2208123</td><td>-0.19020838</td><td>0.09456765</td><td>-0.07534609</td><td>-0.4226652</td><td>-0.11185539</td><td>-0.08965617</td><td>-0.1952422</td><td>0.05758277</td><td>0.004762314</td><td>-0.1964032</td><td>0.1895437</td><td>0.24933858</td><td>-0.03527213</td><td>-0.4760470</td><td>-0.3099167</td><td>0.2464180</td><td>0.2400421</td><td>-0.1644267</td><td>-0.3641763</td><td>-0.2166708</td><td>0.1315210</td><td>-0.02280512</td><td>-0.009367156</td><td>-0.1945582</td><td>-0.136449580</td><td>-0.3098445</td><td> 0.04896163</td><td>-0.3075252</td><td>0.04972026</td><td>0.2480508</td><td>-0.2863842</td><td>0.1155203</td><td>0.2648070</td><td>-0.0521151</td><td>-0.3093807</td><td>-0.30064091</td><td>-0.1276637</td><td>-0.2037646</td><td>-0.09284332</td><td>0.08045442</td><td>-0.1490994</td><td>-0.11645977</td><td> 0.07215251</td><td>-0.05720342</td><td>-0.3740565</td><td>-0.11731536</td><td>0.09530779</td><td>-0.36706489</td><td>-0.3433258</td><td>-0.17452875</td><td>0.1680964</td><td>-0.0600713</td><td>-0.2060481</td><td> 0.1025282</td><td>-0.05558927</td><td>0.07665226</td><td> 0.1833969</td><td>0.1940869</td><td>0.2489087</td><td>-0.2747439</td><td> 0.1270747</td><td>0.17577214</td><td>0.1807008</td><td> 0.04279143</td><td> 0.39460038</td><td>-0.1741859</td><td> 0.13042175</td><td>0.0001290622</td><td>-0.3099167</td><td> 0.1203429</td><td>-0.03140071</td><td>0.05294388</td><td>-0.3052213</td><td> 0.02076767</td><td>-0.3543945</td><td>-0.09103041</td><td>-0.2309566</td><td> 0.2224596</td><td>0.07138237</td><td>-0.2954934</td><td>0.2712612</td><td>-0.04273489</td><td>0.2819382</td><td>0.03289402</td><td>0.002305699</td><td>-0.06404339</td><td>-0.003608402</td><td>0.1568993</td><td> 0.01233974</td><td>-0.3340874</td><td>-0.2117121</td><td>-0.4735961</td><td>-0.4434453</td><td>-0.4282333</td><td>-0.3209603</td><td>0.3765564</td><td>-0.3719606</td><td>-0.1485771</td><td>-0.39659213</td><td>-0.2805517</td><td> 0.06797252</td><td>0.1568422</td><td>-0.2954867</td><td> 0.3575270</td><td>-0.3088493</td><td>⋯</td><td>-0.1867634</td><td> 0.4788550</td><td> 0.1530436</td><td> 0.4767225</td><td>-0.3758204</td><td>-0.2709445</td><td>0.24696047</td><td> 0.1815945</td><td>-0.1029377</td><td>-0.04035266</td><td>-0.4379857</td><td> 0.3272288</td><td>-0.05055972</td><td>0.2431261</td><td>0.2605170</td><td>0.2795922</td><td> 0.05849899</td><td>-0.05374773</td><td>0.2353706</td><td>-0.08403543</td><td>-0.22998002</td><td>0.077529167</td><td>-0.1130794</td><td> 0.1850373</td><td>-0.2161740</td><td> 0.04999117</td><td>-0.1976540</td><td>-0.2676224</td><td>0.3459654</td><td>-0.2538753</td><td>-0.1868426</td><td>0.201937813</td><td>-0.12193045</td><td>-0.2580797</td><td>-0.2704434</td><td>-0.3120529</td><td> 0.04369227</td><td>-0.4190505</td><td> 0.1482068</td><td>-0.3625550</td><td>-0.3660676</td><td>-0.3748491</td><td>-0.01291808</td><td> 0.01459989</td><td>-0.4513064</td><td>-0.1429345</td><td> 0.1661141</td><td>0.2458993</td><td>-0.06533643</td><td>-0.1295380</td><td>-0.1822065</td><td>-0.1998893</td><td>-0.3578749</td><td>-0.0386328</td><td>-0.314490</td><td>-0.08806295</td><td>-0.24626976</td><td>-0.05542134</td><td>-0.3414229</td><td>0.2219668</td><td>-0.2922476</td><td>-0.00179332</td><td>-0.1135965</td><td>-0.08218692</td><td>0.1542454</td><td>-0.3979012</td><td>0.2678977</td><td>-0.2675547</td><td>0.3359187</td><td>-0.2974123</td><td>-0.08740255</td><td> 0.5112075</td><td> 0.21804445</td><td> 0.1326490</td><td>0.07168898</td><td>-0.42377995</td><td>0.01961548</td><td>-0.09446445</td><td>-0.04006233</td><td> 0.134997</td><td> 0.1152315</td><td>-0.05249716</td><td>-0.06943997</td><td>-0.1114678</td><td> 0.09342948</td><td>-0.2700468</td><td> 0.1077615</td><td>-0.008133184</td><td>-0.1011371</td><td>-0.1208799</td><td>-0.4493862</td><td>-0.4037542</td><td>0.1249917</td><td> 0.3761847</td><td>0.02830934</td><td>-0.4168332</td><td> 0.08497293</td><td>-0.3881808</td><td>0.2128561</td><td>-0.02981664</td><td>0.4830124</td><td>-0.0416154</td><td>-0.08952562</td><td>-0.04717037</td><td>-0.1712284</td><td>-0.003889946</td><td> 0.3725608</td><td>-0.02220208</td><td>0.2603706</td><td>-0.1038253</td><td>-0.3221950</td><td>-0.21501092</td><td> 0.02488805</td><td>-0.04693871</td><td>-0.1712358</td><td>0.2788991</td><td>-0.04173324</td><td>-0.1296062</td><td>-0.1463775</td><td> 0.2431421</td><td>-0.3331820</td><td>-0.3058837</td><td>-0.2679635</td><td>-0.1307840</td><td>-0.4110108</td><td>0.3404154</td><td>0.0478274</td><td>0.17154349</td><td>-0.02292978</td><td>-0.3555711</td><td>0.2248661</td><td> 0.07872112</td><td>-0.16932970</td><td>-0.24150183</td><td> 0.189078731</td><td> 0.05250383</td><td>-0.40903544</td><td>-0.331409302</td><td>-0.2207622</td><td>-0.02102189</td><td>-0.01434016</td><td>-0.15264498</td><td>-0.2823520</td><td>-0.2799390</td><td>0.2570737</td><td>-0.3250829</td><td>-0.3595741</td><td>-0.0769582052</td><td>0.2586078</td><td>-0.2311522</td><td>0.3928050</td><td>-0.1285662</td><td>-0.1753022</td><td>-0.4504785</td><td>-0.31590700</td><td>-0.1547715</td><td> 0.1560738</td><td>0.1495561</td><td>0.1442271</td><td>-0.01742539</td><td>-0.2279759</td><td>-0.1331947</td><td>-0.1100318</td><td>-0.1585935</td><td>-0.2438267</td><td>-0.3958475</td><td>-0.2910635</td><td>-0.30136470</td><td> 0.01173644</td><td>0.01026568</td><td>0.0716369</td><td>0.1124081</td><td>-0.4041440</td><td>-0.1350245</td><td>-0.06057534</td><td>-0.1609034</td><td>0.03561598</td><td>-0.41089336</td><td>-0.1750870</td><td> 0.180930648</td><td>-0.4592652</td><td>0.07835832</td><td> 0.06508965</td><td>-0.07552363</td><td>-0.3225593</td><td>-0.1677291</td><td> 0.3512619</td><td>-0.3204145</td><td>-0.1709429</td><td>-0.09960947</td><td>-0.4409674</td><td> 0.07894281</td><td>-0.3116489</td><td>-0.3639286</td><td>-0.21087035</td><td>-0.06955116</td><td> 0.2742827</td><td>-0.05444851</td><td>0.06749174</td><td>-0.2651614</td><td>-0.3958076</td><td> 0.003335668</td><td> 0.3998135</td><td>-0.16050048</td><td>-0.07827426</td><td>0.2212427</td><td>-0.06808805</td><td>-0.1329919</td><td>0.4419936</td><td>-0.19777949</td><td>-0.07029857</td><td> 0.10664154</td><td> 0.07933052</td><td>0.15578896</td><td>0.33827494</td><td> 0.12175110</td><td>-0.1171267</td><td>0.4384580</td><td> 0.1148695</td><td>-0.1835288</td><td>-0.361643308</td><td>-0.1553061</td><td>-0.3186185</td><td> 0.04109448</td><td>-0.1059865</td><td>-0.2427387</td><td>-0.2760731</td><td>-0.2385194</td><td>-0.12021035</td><td>-0.4143417</td><td>-0.1518727</td><td>-0.26182718</td><td>0.3892487</td><td>-0.02973562</td><td>-0.09875486</td><td>0.1025138</td><td>-0.1643691</td><td>-0.3836949</td><td> 0.1192650</td><td>-0.4060676</td><td>0.2537932</td><td>0.3471116</td><td> 0.08060762</td><td>-0.07647283</td><td>0.2929339</td><td>-0.11235892</td><td>0.1695442</td><td> 0.1835923</td><td>-0.2356944</td><td> 0.009843151</td><td>0.1269156</td><td>-0.4201098</td><td> 0.31678774</td><td>-0.3978131</td><td>-0.2157045</td><td>-0.03490663</td><td>0.02447295</td><td>-0.03929059</td><td>-0.1675013</td><td>0.38264060</td><td>-0.2147420</td><td>0.02955894</td><td>-0.3619161</td><td>-0.05623514</td><td>-0.01676065</td><td> 0.1497797</td><td> 0.2772710</td><td>0.16189149</td><td> 0.04800534</td><td>-0.2745364</td><td> 0.1693959</td><td>-0.3039321</td><td>-0.2983708</td><td>-0.1323475</td><td> 0.30964315</td><td>0.2314451</td><td>-0.0139266</td><td>-0.09291987</td><td>-0.2587566</td><td>-0.3382832</td><td>-0.01589071</td><td>0.001714197</td><td>0.03909666</td><td>0.4131416</td><td>-0.36692512</td><td>-0.33691858</td><td>-0.2807307</td><td>-0.01658442</td><td>-0.0234535</td><td>-0.3324207</td><td>-0.08247293</td><td>0.2839719</td><td>0.4382558</td><td>-0.1035817</td><td>0.3359403</td><td>-0.3726248</td><td>-0.43697356</td><td>0.52612948</td><td>-0.4460085</td><td>-0.3153409</td><td>-0.3452522</td><td>0.293858</td><td>-0.3708351</td><td>0.04039368</td><td> 0.01320762</td><td>-0.05379969</td><td>-0.4126428</td><td>-0.05972381</td><td>-0.3980399</td><td>0.26909186</td><td> 0.1131846</td><td>-0.189109085</td><td>-0.1158904</td><td>-0.3377623</td><td>-0.1059766</td><td>0.04681667</td><td>0.22766870</td><td>0.1302709</td><td>-0.5595581</td><td>-0.07152642</td><td>-0.3951038</td><td>-0.3373137</td><td>-0.2678920</td><td>-0.37903305</td><td>-0.004550087</td><td> 0.179811780</td><td>-0.3099022</td><td>-0.04854573</td><td> 0.2251317</td><td>0.08552927</td><td>-0.009245791</td><td>-0.1798787</td><td>-0.26881563</td><td>-0.02711064</td><td> 0.1012830</td><td>-0.2192352</td><td>-0.4082294</td><td>0.06047531</td><td>-0.011167889</td><td>-0.37845537</td><td>-0.05018487</td><td>-0.13281380</td><td>-0.4435681</td><td>0.2985744</td><td>-0.3764525</td><td>0.05362786</td><td>-0.09927756</td><td>-0.4627661</td><td> 0.02207465</td><td> 0.07294996</td><td>-0.2320804</td><td>-0.1782898</td><td>-0.05359784</td><td>-0.19232055</td><td> 0.1475623</td><td>-0.1886234</td><td> 0.2301737</td><td>0.2342421</td><td>-0.2839519</td><td>-0.1695806</td><td>-0.1871436</td><td>-0.2120458</td><td>0.45880546</td><td> 0.0530361</td><td> 0.20436942</td><td>-0.3973664</td><td>0.03575068</td><td> 0.1967032</td><td> 0.01355138</td><td>-0.002604579</td><td>0.4341153</td><td>-0.1668098</td><td>-0.06336949</td><td>-0.4347483</td><td>0.1927845</td><td>-0.1160828</td><td>0.2912478</td><td>-0.09396141</td><td>-0.12771137</td><td>-0.1654102</td><td> 0.01926317</td><td>-0.2649523</td><td>-0.13232947</td><td>-0.3124342</td><td>-0.08780874</td><td> 0.02274332</td><td> 0.20458957</td><td> 0.02938472</td><td>0.1510418</td><td> 0.2491153</td><td>-0.34277889</td><td>-0.08575543</td><td>0.07957073</td><td> 0.0934102</td><td>0.2659089</td><td> 0.08195431</td><td>-0.4361261</td><td>-0.1277941</td><td> 0.06608009</td><td>-0.08299027</td><td>-0.1246469</td><td>0.119697989</td><td>0.5373233</td><td>-0.2499181</td><td>-0.3461719</td><td>-0.1029538</td><td>-0.2999278</td><td>0.4198123</td><td>-0.2126493</td><td>0.04753214</td><td>-0.04205889</td><td>-0.08409886</td><td>0.3042187</td><td>0.1606759</td><td>-0.1431611</td><td>-0.08137913</td><td>0.02835725</td><td> 0.37608373</td><td> 0.1054727</td><td>-0.04061472</td><td>0.2693739</td><td> 0.03287404</td><td>-0.1057462</td><td>-0.1722411</td><td>-0.220961332</td><td> 0.1279675</td><td>-0.1954983</td><td>-0.02159357</td><td>-0.06668337</td><td>-0.4126998</td><td> 0.05307225</td><td> 0.08825069</td><td>-0.4187567</td><td> 0.14149411</td><td> 0.4672388</td><td>-0.3400920</td><td> 0.2914611</td><td>-0.10230253</td><td>0.4705167</td><td>-0.112970142</td><td>-0.09549635</td><td>-0.1889675</td><td>-0.1983383</td><td> 0.09670054</td><td>-0.4104028</td><td>-0.2012322</td><td>-0.4296237</td><td>0.02753232</td><td>-0.01893576</td><td>-0.3270123</td><td>-0.2221356</td><td>-0.06358689</td><td>-0.19617519</td><td> 0.1700201</td><td>-0.06619192</td><td>-0.1446155</td><td>-0.3209004</td><td>-0.2057163</td><td>-0.1182971</td><td>0.08331947</td><td>-0.2036227</td><td>-0.3809979</td><td> 0.32525563</td><td> 0.05790675</td><td>-0.03960102</td><td>-0.01721894</td><td>0.06612352</td><td>-0.278505</td><td>-0.4159777</td><td>-0.05612973</td><td>0.3624167</td><td> 0.3665413</td><td>0.2060093</td><td>0.12799406</td><td>-0.36692701</td><td>0.3230871</td><td>0.007297871</td><td>-0.1040083</td><td>-0.04166428</td><td>-0.4272843</td><td> 0.1134232</td><td> 0.1339973</td><td>-0.1298473</td><td>0.3075708</td><td>0.1394638</td><td>-0.3207854</td><td>-0.1742502</td><td>0.1140986</td><td>0.1189916</td><td>-0.4418161</td><td>0.04587000</td><td>-0.2617494</td><td>-0.2568759</td><td>-0.3723166</td><td>-0.2127830</td><td>-0.02391508</td><td>-0.2477704</td><td>-0.3752086</td><td>0.2921108</td><td>-0.09576479</td></tr>
	<tr><th scope=row>Angiogenesis</th><td>-0.25984750</td><td>-0.3636956</td><td>-0.1422460</td><td>-0.1592542</td><td>0.6701079</td><td>-0.4711628</td><td>-0.17020046</td><td>-0.0003359958</td><td>-0.04696236</td><td>0.4667481</td><td>-0.04242067</td><td> 0.11057705</td><td>-0.07343843</td><td>-0.1535245</td><td>-0.2085669</td><td>-0.06393342</td><td>-0.04797819</td><td>0.2628455</td><td>-0.4874150</td><td>-0.01397154</td><td>-0.1809975</td><td>-0.3770225</td><td> 0.01444596</td><td>0.4560055</td><td>0.11433245</td><td> 0.61165446</td><td>-0.001063347</td><td>-0.03768425</td><td>-0.2768315</td><td> 0.01258583</td><td> 0.2140868</td><td>-0.07065291</td><td>-0.22365934</td><td>-0.4863134</td><td>-0.2287146</td><td> 0.02370807</td><td>-0.240516559</td><td>-0.04283355</td><td>0.6907327</td><td>0.396864672</td><td> 0.01099594</td><td>0.3814791</td><td>-0.3625070</td><td>0.7135791</td><td> 0.27613940</td><td>-0.3672031</td><td> 0.14704333</td><td> 0.23477637</td><td> 0.09344184</td><td> 0.3074853</td><td>0.1826995</td><td>-0.2058980</td><td>-0.06749687</td><td>-0.6499562</td><td>-0.006384953</td><td>-0.6206675</td><td> 0.1279143</td><td>0.5372368</td><td>-0.4168516</td><td>0.1076207</td><td> 0.58565151</td><td>-0.367754347</td><td>-0.2085319</td><td>0.02740239</td><td>0.1686087</td><td>0.4487901</td><td>0.1020619</td><td> 0.094354333</td><td>0.21763470</td><td>-0.3082402</td><td> 0.2555036</td><td>0.01507229</td><td>0.05536619</td><td>0.3260233</td><td> 0.3826466</td><td>-0.6637662</td><td>-0.3141366</td><td> 0.16673535</td><td>-0.5719394</td><td>0.6015446</td><td> 0.05026585</td><td>0.5492422</td><td>-0.4585718</td><td>-0.5732724</td><td>-0.1157847</td><td>0.28321425</td><td>-0.5225046</td><td>-0.1523671</td><td>-0.1509056</td><td>-0.3757087</td><td> 0.04724007</td><td>-0.14152041</td><td>-0.6278948</td><td>-0.4243095</td><td>-0.2376023</td><td> 0.08075442</td><td>0.4282780</td><td>-0.2423978</td><td>-0.3745077</td><td>0.10814100</td><td>0.00230075</td><td>-0.1653922</td><td>0.5162211</td><td>-0.5073745</td><td>0.48993776</td><td>-0.02407317</td><td> 0.2581412</td><td>-0.5728008</td><td>-0.3960099</td><td>0.02009229</td><td>-0.1175501</td><td>-0.3902330</td><td>-0.3670949</td><td> 0.02849531</td><td>0.4618849</td><td>0.06591504</td><td> 0.2997496</td><td> 0.55482390</td><td>-0.5620277</td><td>0.1282902</td><td> 0.52177948</td><td>-0.3624435</td><td>0.07253998</td><td>0.0006543838</td><td>-0.25638925</td><td>-0.01351337</td><td> 0.1853688</td><td> 0.1053426</td><td> 0.2635640</td><td>-0.001490455</td><td>-0.1183433</td><td>-0.13246527</td><td>0.1095951</td><td>-0.37234996</td><td>-0.08675361</td><td>-0.3257736</td><td>-0.1960821</td><td>0.1526852</td><td> 0.05432829</td><td>-0.3253852</td><td>-0.1346982</td><td>-0.1195909</td><td>-0.4204694</td><td>-0.21360517</td><td>0.1637281</td><td>0.1868878</td><td>-0.01860376</td><td> 0.0957722</td><td>-0.4551103</td><td> 0.05011947</td><td>0.1719030</td><td>-0.1243636</td><td>0.3431719</td><td>-0.2632384</td><td>0.03841197</td><td>-0.1836877</td><td>0.09101893</td><td>0.0354472</td><td>-0.4107937</td><td>-0.4403284</td><td>-0.007476868</td><td>0.2138696</td><td>-0.4330532</td><td>-0.45805870</td><td>-0.4656138</td><td> 0.04163774</td><td>-0.3292435</td><td>0.4142333</td><td>-0.4230567</td><td> 0.25687756</td><td>-0.06838872</td><td>-0.4048449</td><td>-0.6527865</td><td>-0.5192462</td><td>-0.3269120</td><td>-0.1827541</td><td>0.002104178</td><td>0.31587959</td><td> 0.05083448</td><td>-0.3940874</td><td>0.2231452</td><td> 0.261876</td><td>-0.4073699</td><td>-0.4787080</td><td>-0.41212495</td><td>0.03065182</td><td>0.17517157</td><td>0.07535777</td><td> 0.44671168</td><td>-0.3176447</td><td> 0.4127027</td><td>-0.5546477</td><td>-0.4446784</td><td>-0.16263060</td><td>0.49236853</td><td>-0.1453150</td><td>-0.03648807</td><td> 0.6561875</td><td>-0.2251917</td><td>-0.56131852</td><td>0.2993442</td><td>0.61747439</td><td> 0.01430326</td><td> 0.01296601</td><td>0.0600290</td><td>-0.51300439</td><td>-0.1570327</td><td>-0.1190901</td><td> 0.01839513</td><td>0.14179510</td><td>0.29429791</td><td>-0.6706975</td><td>-0.22967887</td><td>0.3210804</td><td>-0.4638339</td><td>-0.3161120</td><td>-0.1694196</td><td>-0.1972168</td><td>-0.3133599</td><td>0.2510939</td><td>0.4257403</td><td> 0.30092100</td><td>-0.2695913</td><td>-0.1415836</td><td>0.20648382</td><td>-0.07598554</td><td>0.7727161</td><td>-0.2332948</td><td> 0.4802928</td><td> 0.3654595</td><td>-0.12088754</td><td> 0.414434083</td><td>0.1537230</td><td> 0.1481233</td><td>0.337811379</td><td>-0.2380994</td><td>-0.4289141</td><td> 0.14240438</td><td>-0.2457093</td><td> 0.2802964</td><td> 0.4871435</td><td>-0.44100913</td><td> 0.2866335</td><td> 0.523455236</td><td> 0.31051372</td><td>-0.5766261</td><td>-0.0005736743</td><td>0.5233992</td><td> 0.33514176</td><td>-0.48942444</td><td>-0.01317547</td><td> 0.13304379</td><td> 0.06251385</td><td> 0.1004945</td><td>-0.001697443</td><td>-0.3978925</td><td>-0.1879647</td><td>-0.14385531</td><td>-0.2543577</td><td>0.21621210</td><td> 0.28560014</td><td>-0.09525923</td><td>-0.20645404</td><td>-0.1596103</td><td>0.06196806</td><td>-0.3783880</td><td>0.3961926</td><td>0.2178788</td><td>0.5343948</td><td>0.1286802</td><td>0.4970720</td><td>-0.05419484</td><td> 0.11537147</td><td> 0.17842386</td><td>-0.4569987</td><td> 0.17961701</td><td>-0.08888919</td><td>-0.3823266</td><td>-0.53334887</td><td>0.4159262</td><td>0.3136098</td><td>-0.42383488</td><td>-0.06621759</td><td>0.1972199</td><td>-0.1274751</td><td>-0.3042103</td><td>-0.09028521</td><td>-0.19527835</td><td>-0.01164516</td><td>0.2560857</td><td>-0.19949592</td><td>-0.17877276</td><td>-0.2930005</td><td>-0.005822685</td><td>-0.5875706</td><td>-0.01002931</td><td>-0.5390341</td><td>-0.4966469</td><td> 0.1264815</td><td>-0.3151843</td><td>0.02726519</td><td>-0.5729324</td><td>0.08995127</td><td>0.202007648</td><td>0.34766154</td><td>-0.3481231</td><td>0.2015457</td><td>-0.3742436</td><td>-0.3429816</td><td>0.3412271</td><td>0.3526454</td><td>0.1126705</td><td>-0.02751105</td><td>-0.12337573</td><td>0.2070949</td><td>0.07664133</td><td> 0.2237549</td><td>-0.17447329</td><td>-0.5737808</td><td>-0.094377752</td><td>-0.4010789</td><td>-0.1882925</td><td> 0.1479069</td><td>0.12957561</td><td>-0.3115511</td><td>0.08251391</td><td> 0.25970381</td><td>-0.5852084</td><td>-0.01629256</td><td>0.33415612</td><td>-0.09765444</td><td>-0.2707552</td><td>0.6014395</td><td>-0.1479389</td><td>-0.1618588</td><td>-0.5761008</td><td> 0.09938865</td><td>-0.4077140</td><td>-0.6777638</td><td>-0.002964503</td><td>-0.1035770</td><td>0.6083108</td><td>-0.2726706</td><td>0.3169254</td><td>-0.07108634</td><td>-0.3860530</td><td>-0.38500765</td><td>-0.4911162</td><td>-0.3567562</td><td>-0.07969108</td><td>0.3316270</td><td>-0.34388221</td><td>-0.2786365</td><td>0.1590021</td><td>-0.3794368</td><td> 0.33393767</td><td>0.3614352</td><td>0.6014407</td><td>0.5687580</td><td>0.2087130</td><td> 0.31722013</td><td>-0.4305445</td><td>-0.3911505</td><td>-0.2132914</td><td>-0.03583908</td><td>0.4940780</td><td>0.64278284</td><td> 0.13021716</td><td>-0.5502058</td><td>0.3106238</td><td> 0.1146198</td><td>-0.1772358</td><td>0.3948534</td><td>-0.5619103</td><td>-0.5075470</td><td>-0.2381862</td><td>-0.4042711</td><td>-0.41869658</td><td> 0.153941079</td><td>-0.4026765</td><td>-0.5711471</td><td>-0.2168016</td><td>-0.4653136</td><td>0.04389784</td><td> 0.1714884</td><td>0.4008581</td><td>-0.33049223</td><td> 0.0141049</td><td>-0.1548210</td><td>-0.1929613</td><td>0.2509835</td><td>-0.25792093</td><td> 0.06632242</td><td> 0.1037984</td><td> 0.0789143</td><td> 0.05733791</td><td>0.07088974</td><td> 0.34632939</td><td>-0.1325203</td><td>-0.07412579</td><td> 0.33189791</td><td>-0.5691118</td><td>0.18418883</td><td>0.091821948</td><td> 0.2218705</td><td>0.4733507</td><td>0.09622021</td><td>-0.39473312</td><td>-0.2970306</td><td>-0.5737476</td><td>0.3364499</td><td>0.6017285</td><td>-0.3018409</td><td>-0.5891914</td><td>-0.2076279</td><td>0.1578641</td><td>-0.38012324</td><td>-0.638585752</td><td>-0.5717998</td><td>-0.008010502</td><td>-0.5732485</td><td>-0.29324451</td><td>-0.4317382</td><td>0.09221203</td><td>0.4131254</td><td>-0.3501003</td><td>0.3647816</td><td>0.2526828</td><td>-0.3228010</td><td>-0.5538268</td><td> 0.09832949</td><td>-0.1601513</td><td> 0.3469274</td><td> 0.32979024</td><td>0.09953012</td><td> 0.6117070</td><td>-0.09630877</td><td>-0.26549581</td><td>-0.19407524</td><td>-0.3041007</td><td>-0.09032109</td><td>0.24971848</td><td>-0.06902876</td><td>-0.2101776</td><td> 0.02375782</td><td>0.1909901</td><td>-0.6401896</td><td>-0.2352884</td><td>-0.1083845</td><td> 0.32587672</td><td>0.20420617</td><td>-0.2461510</td><td>0.5240396</td><td>0.2880494</td><td> 0.1422201</td><td>-0.1230288</td><td>0.01410649</td><td>0.4819849</td><td>-0.10710857</td><td>-0.01018997</td><td>-0.1458161</td><td>-0.01352768</td><td>0.4918288615</td><td>-0.5730847</td><td>-0.2159404</td><td> 0.14613546</td><td>0.62646184</td><td> 0.3317283</td><td>-0.09949803</td><td>-0.2937368</td><td> 0.08840018</td><td>-0.1440285</td><td>-0.0715974</td><td>0.28913461</td><td>-0.3961021</td><td>0.3497378</td><td> 0.07229490</td><td>0.3236912</td><td>0.07903127</td><td>0.081012806</td><td> 0.06780445</td><td> 0.249056137</td><td>0.2546789</td><td>-0.23828801</td><td>-0.2697891</td><td>-0.1395339</td><td>-0.2107730</td><td>-0.2703367</td><td>-0.4367136</td><td>-0.3755254</td><td>0.3853345</td><td>-0.2335643</td><td>-0.2508135</td><td> 0.09773947</td><td>-0.4968358</td><td>-0.55452051</td><td>0.4432255</td><td>-0.2632692</td><td>-0.3156793</td><td>-0.3149645</td><td>⋯</td><td>-0.3390147</td><td>-0.1227305</td><td>-0.5294721</td><td>-0.4063255</td><td>-0.2696010</td><td>-0.1611574</td><td>0.04812393</td><td>-0.1905383</td><td>-0.2107007</td><td> 0.03931479</td><td>-0.4710200</td><td>-0.5025380</td><td>-0.35550901</td><td>0.2399132</td><td>0.6833896</td><td>0.4964123</td><td>-0.36809234</td><td> 0.04354493</td><td>0.5029911</td><td>-0.31397810</td><td>-0.02598467</td><td>0.005668143</td><td>-0.1602844</td><td>-0.0493400</td><td> 0.1655819</td><td>-0.08618036</td><td>-0.1589374</td><td>-0.4543803</td><td>0.3557728</td><td>-0.3799799</td><td>-0.3216259</td><td>0.005984613</td><td> 0.09573245</td><td>-0.2861227</td><td>-0.3690004</td><td>-0.1501640</td><td>-0.34937186</td><td>-0.3478297</td><td>-0.1781148</td><td>-0.3869951</td><td> 0.3140536</td><td>-0.1915203</td><td> 0.29305694</td><td>-0.06231968</td><td>-0.1658483</td><td>-0.1736317</td><td>-0.1366788</td><td>0.1768220</td><td>-0.52369924</td><td>-0.4784521</td><td> 0.2751477</td><td> 0.0984491</td><td>-0.5137177</td><td>-0.2491935</td><td>-0.395493</td><td> 0.02660855</td><td> 0.07327944</td><td>-0.46547397</td><td>-0.2313559</td><td>0.1515850</td><td>-0.3066094</td><td>-0.38505009</td><td>-0.4074883</td><td>-0.35909834</td><td>0.1426220</td><td>-0.5938501</td><td>0.1700843</td><td> 0.4067966</td><td>0.2916024</td><td>-0.4903922</td><td>-0.01700463</td><td>-0.0774881</td><td>-0.04141737</td><td>-0.1433562</td><td>0.26993163</td><td> 0.01586061</td><td>0.12109550</td><td>-0.52329804</td><td>-0.52501635</td><td>-0.240485</td><td>-0.4919003</td><td> 0.23829851</td><td> 0.63166921</td><td>-0.4915326</td><td>-0.05899697</td><td>-0.2874448</td><td>-0.5509113</td><td>-0.155682972</td><td>-0.1230804</td><td> 0.1988217</td><td>-0.4755008</td><td>-0.3043791</td><td>0.3019636</td><td>-0.7033401</td><td>0.56243350</td><td>-0.2883838</td><td>-0.10450374</td><td>-0.4045382</td><td>0.1376480</td><td> 0.17729134</td><td>0.4717843</td><td>-0.5098666</td><td> 0.02086567</td><td> 0.22817417</td><td>-0.2947724</td><td>-0.387473921</td><td>-0.3818490</td><td>-0.31405947</td><td>0.5187994</td><td>-0.3048920</td><td>-0.4291686</td><td>-0.02235258</td><td>-0.51631671</td><td> 0.40910824</td><td>-0.1515091</td><td>0.5446293</td><td> 0.57734101</td><td> 0.3935537</td><td>-0.4479107</td><td>-0.2352606</td><td>-0.3405165</td><td>-0.3888115</td><td>-0.4577648</td><td>-0.2432291</td><td>-0.3290535</td><td>0.5397457</td><td>0.3075494</td><td>0.09749611</td><td>-0.02562862</td><td>-0.1339351</td><td>0.3957596</td><td>-0.54413284</td><td> 0.07439447</td><td>-0.06263923</td><td>-0.003168958</td><td>-0.35072935</td><td>-0.04586119</td><td> 0.001802348</td><td>-0.1820981</td><td> 0.40473870</td><td>-0.55875771</td><td> 0.01446457</td><td>-0.1554227</td><td>-0.4439322</td><td>0.5823722</td><td>-0.3544284</td><td>-0.2333885</td><td> 0.0008232575</td><td>0.5377679</td><td>-0.4772320</td><td>0.4037898</td><td> 0.4528646</td><td> 0.2122021</td><td>-0.6533889</td><td>-0.07677771</td><td>-0.2578457</td><td>-0.1907226</td><td>0.2940598</td><td>0.4928764</td><td> 0.27531493</td><td>-0.5166856</td><td>-0.3163683</td><td>-0.5373713</td><td>-0.2346069</td><td> 0.2362494</td><td>-0.3389715</td><td>-0.2185356</td><td> 0.03334188</td><td>-0.12125087</td><td>0.21871671</td><td>0.4582868</td><td>0.3353924</td><td>-0.3913554</td><td>-0.1051620</td><td>-0.09611414</td><td> 0.3401315</td><td>0.03856625</td><td>-0.05725103</td><td>-0.1845809</td><td>-0.007964884</td><td>-0.6371507</td><td>0.01365092</td><td>-0.57579949</td><td> 0.07446113</td><td>-0.6060727</td><td>-0.2041980</td><td>-0.2112953</td><td>-0.3928505</td><td>-0.2429969</td><td> 0.31923595</td><td>-0.6110997</td><td>-0.19855913</td><td>-0.5111058</td><td> 0.0575339</td><td>-0.01143839</td><td>-0.20932478</td><td>-0.5979168</td><td> 0.32923586</td><td>0.49756728</td><td>-0.4483616</td><td>-0.4763712</td><td>-0.015358799</td><td>-0.1314418</td><td>-0.07469618</td><td>-0.25302903</td><td>0.1809130</td><td>-0.07205173</td><td>-0.2993059</td><td>0.4839097</td><td> 0.07552712</td><td>-0.14908292</td><td>-0.02179004</td><td>-0.37654779</td><td>0.04643748</td><td>0.03614418</td><td>-0.05925172</td><td> 0.3391374</td><td>0.4100915</td><td>-0.2698218</td><td>-0.4776318</td><td> 0.009418023</td><td>-0.4618295</td><td>-0.6046176</td><td>-0.26654483</td><td>-0.1181522</td><td>-0.0883438</td><td>-0.3002060</td><td>-0.2114746</td><td> 0.04593501</td><td>-0.2102116</td><td> 0.1645404</td><td> 0.04796873</td><td>0.1158586</td><td> 0.47794737</td><td> 0.13332868</td><td>0.2872786</td><td>-0.1807839</td><td>-0.2622786</td><td>-0.6952196</td><td>-0.4085504</td><td>0.6443345</td><td>0.6828261</td><td>-0.11713623</td><td>-0.42243410</td><td>0.4951474</td><td>-0.01301023</td><td>0.3024781</td><td>-0.3186378</td><td>-0.5222819</td><td>-0.202233871</td><td>0.2292280</td><td> 0.3017152</td><td>-0.07032425</td><td>-0.3792439</td><td>-0.1641742</td><td> 0.40131311</td><td>0.26300527</td><td> 0.06070517</td><td>-0.3076398</td><td>0.04169458</td><td> 0.2537745</td><td>0.12632188</td><td>-0.4970589</td><td> 0.13691841</td><td>-0.35866232</td><td>-0.1678586</td><td>-0.1439681</td><td>0.09679793</td><td>-0.14299002</td><td>-0.4559607</td><td>-0.2612019</td><td> 0.1232419</td><td>-0.0228524</td><td> 0.1154897</td><td>-0.06999098</td><td>0.4223545</td><td> 0.3228061</td><td>-0.47932174</td><td> 0.2158796</td><td> 0.2898621</td><td>-0.12054943</td><td>0.418125539</td><td>0.25280302</td><td>0.6329719</td><td>-0.08492894</td><td> 0.05850902</td><td>-0.4985332</td><td> 0.28053016</td><td> 0.2350526</td><td>-0.1181108</td><td>-0.66290842</td><td>0.3887361</td><td>0.1014662</td><td>-0.4884211</td><td>0.4103845</td><td>-0.5104947</td><td>-0.01717798</td><td>0.02774256</td><td>-0.7036878</td><td>-0.3823652</td><td>-0.1940602</td><td>0.131812</td><td>-0.6030550</td><td>0.31236602</td><td>-0.31421212</td><td> 0.11841910</td><td>-0.6001392</td><td>-0.04590991</td><td>-0.3415829</td><td>0.07254807</td><td>-0.2834003</td><td>-0.009397029</td><td>-0.4078298</td><td>-0.4285412</td><td>-0.3477204</td><td>0.22567992</td><td>0.06965403</td><td>0.0357814</td><td>-0.5321642</td><td>-0.30042120</td><td>-0.3827237</td><td> 0.0108267</td><td>-0.3335698</td><td> 0.05038181</td><td> 0.306120037</td><td>-0.007431344</td><td>-0.3101018</td><td>-0.63386973</td><td>-0.2043446</td><td>0.17743601</td><td> 0.081655931</td><td> 0.1378592</td><td> 0.03356901</td><td>-0.04835605</td><td>-0.5065171</td><td>-0.3575868</td><td>-0.1333571</td><td>0.29191053</td><td> 0.003450711</td><td>-0.05998573</td><td>-0.37628065</td><td>-0.08841685</td><td>-0.6459033</td><td>0.2927352</td><td>-0.4792835</td><td>0.05409607</td><td> 0.03294957</td><td>-0.2274207</td><td>-0.26294682</td><td>-0.14356862</td><td>-0.4683692</td><td>-0.1721464</td><td> 0.13543621</td><td> 0.06174825</td><td>-0.5097365</td><td>-0.1226915</td><td>-0.1620284</td><td>0.1692557</td><td>-0.4626556</td><td>-0.3073997</td><td> 0.6053792</td><td> 0.1156187</td><td>0.02749207</td><td>-0.4120545</td><td>-0.09829443</td><td>-0.1468669</td><td>0.19731708</td><td>-0.3244289</td><td>-0.11371405</td><td> 0.222947339</td><td>0.1951857</td><td> 0.4681358</td><td> 0.11261977</td><td>-0.3523676</td><td>0.2030592</td><td>-0.3392617</td><td>0.4701438</td><td>-0.03030032</td><td>-0.04157864</td><td>-0.5327607</td><td>-0.06862742</td><td> 0.3808016</td><td> 0.09435969</td><td> 0.1373493</td><td> 0.20309419</td><td>-0.17009612</td><td>-0.03987281</td><td>-0.20223592</td><td>0.2216734</td><td>-0.1537745</td><td>-0.07614608</td><td>-0.06268976</td><td>0.06516105</td><td>-0.3297110</td><td>0.1917183</td><td>-0.22054819</td><td>-0.4715682</td><td>-0.1490165</td><td>-0.03982557</td><td> 0.08688297</td><td> 0.4771506</td><td>0.005699211</td><td>0.1473911</td><td> 0.3705624</td><td>-0.5102403</td><td>-0.1600281</td><td>-0.3527140</td><td>0.3577708</td><td>-0.1933785</td><td>0.28106506</td><td>-0.31284489</td><td>-0.48950307</td><td>0.2009962</td><td>0.2343899</td><td> 0.0480065</td><td> 0.05053740</td><td>0.04916981</td><td>-0.02753923</td><td>-0.2503862</td><td>-0.39039776</td><td>0.4895703</td><td>-0.32853346</td><td> 0.1092684</td><td>-0.2514773</td><td>-0.004770143</td><td>-0.1519982</td><td>-0.2223781</td><td>-0.30871940</td><td>-0.25963421</td><td>-0.2367030</td><td>-0.01579496</td><td>-0.29405391</td><td>-0.2885283</td><td>-0.05358234</td><td>-0.0240433</td><td> 0.2492149</td><td>-0.5475317</td><td>-0.06684206</td><td>0.2085477</td><td>-0.008451096</td><td> 0.37755713</td><td> 0.3569229</td><td>-0.5086309</td><td>-0.39919391</td><td>-0.3818828</td><td>-0.2864245</td><td>-0.2075318</td><td>0.03080272</td><td>-0.28294419</td><td>-0.3956065</td><td>-0.3891204</td><td> 0.17052181</td><td>-0.09495074</td><td>-0.1744794</td><td> 0.58765200</td><td>-0.3751693</td><td> 0.1842594</td><td>-0.4827688</td><td>-0.1960242</td><td>0.23315471</td><td>-0.2259431</td><td>-0.3010001</td><td>-0.03315458</td><td>-0.34345912</td><td>-0.28109768</td><td>-0.33129724</td><td>0.23152291</td><td>-0.142882</td><td> 0.1018331</td><td>-0.34126208</td><td>0.5326904</td><td>-0.6005726</td><td>0.3175016</td><td>0.04804074</td><td> 0.08302104</td><td>0.3083732</td><td>0.380427701</td><td>-0.5336194</td><td> 0.66093932</td><td>-0.2005885</td><td>-0.2213533</td><td>-0.3282462</td><td> 0.4557656</td><td>0.3209752</td><td>0.2601042</td><td>-0.3384646</td><td> 0.2394631</td><td>0.5359387</td><td>0.3672545</td><td>-0.5035941</td><td>0.07933361</td><td>-0.3412286</td><td>-0.5237940</td><td>-0.3509303</td><td>-0.2196897</td><td>-0.53662302</td><td>-0.5375393</td><td>-0.3140421</td><td>0.3118755</td><td>-0.01556957</td></tr>
</tbody>
</table>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>31</li></ol>




<table class="dataframe">
<caption>A data.frame: 2 × 31</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td></tr>
</tbody>
</table>




```R
p <- pca(annots %>% select(all_of(GeneSetList)) %>% t(), metadata = annots %>% select(-any_of(GeneSetList)), removeVar = 0.1)
```

    -- removing the lower 10% of variables based on variance
    



```R
screeplot(p, axisLabSize = 18, titleLabSize = 22)
```


    
![png](output_69_0.png)
    



```R
options(repr.plot.width=9, repr.plot.height=9)
biplot(p, showLoadings = TRUE, labSize = 5, pointSize = 5, sizeLoadingsNames = 5, alpha = 0.25)
```

    Warning message:
    “[1m[22mRemoved 2 rows containing missing values or values outside the scale range
    (`geom_segment()`).”
    Warning message:
    “[1m[22mRemoved 2 rows containing missing values or values outside the scale range
    (`geom_label_repel()`).”
    Warning message:
    “ggrepel: 9711 unlabeled data points (too many overlaps). Consider increasing max.overlaps”



    
![png](output_70_1.png)
    



```R
options(repr.plot.width=9, repr.plot.height=9)
# color by Tissue Slice
biplot(p,
    colby = 'Tissue_Slice',
    hline = 0, vline = 0,
    legendPosition = 'right',
    alpha = 0.5) 
```

    Warning message:
    “ggrepel: 9704 unlabeled data points (too many overlaps). Consider increasing max.overlaps”



    
![png](output_71_1.png)
    



```R
options(repr.plot.width=11, repr.plot.height=9)
# color by Navin Clusters
biplot(p,
    colby = 'Navin_Clusters',
    hline = 0, vline = 0,
    legendPosition = 'right',
    alpha = 0.5)
```

    Warning message:
    “ggrepel: 9704 unlabeled data points (too many overlaps). Consider increasing max.overlaps”



    
![png](output_72_1.png)
    



```R
options(repr.plot.width=11, repr.plot.height=9)
# color by vasc_MN4_label
biplot(p,
    colby = 'vasc_MN4_label',
    hline = 0, vline = 0,
    legendPosition = 'right',
    alpha = 0.5)
```

    Warning message:
    “ggrepel: 9704 unlabeled data points (too many overlaps). Consider increasing max.overlaps”



    
![png](output_73_1.png)
    



```R
options(repr.plot.width=11, repr.plot.height=9)
# color by Tumor_MHCI_label2
biplot(p,
    colby = 'Tumor_MHCI_label2',
    hline = 0, vline = 0,
    legendPosition = 'right',
    alpha = 0.5)
```

    Warning message:
    “ggrepel: 9704 unlabeled data points (too many overlaps). Consider increasing max.overlaps”



    
![png](output_74_1.png)
    


# Do UMAP of GSVA values
install.packages('umap')

```R
library(umap)
```


```R
# prepare data... just need the GSVA values
gsva_values <- annots %>% select(all_of(GeneSetList))
gsva_values[1:3,1:3]
```


<table class="dataframe">
<caption>A data.frame: 3 × 3</caption>
<thead>
	<tr><th></th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td></tr>
</tbody>
</table>




```R
gsva_umap <- umap(gsva_values)
```


```R
gsva_umap
```

    umap embedding of 9711 items in 2 dimensions
    object components: layout, data, knn, config
    



```R
gsva_umap$layout %>% head(3)
```


<table class="dataframe">
<caption>A matrix: 3 × 2 of type dbl</caption>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>-1.609071</td><td> 0.8862022</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td> 3.396564</td><td>-2.2385961</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td> 2.815883</td><td> 2.0775624</td></tr>
</tbody>
</table>




```R
# get some labels
gsva_umap_labels_Tissue_Slice <- annots$Tissue_Slice
gsva_umap_labels_Navin_Clusters2 <- annots$Navin_Clusters2
```


```R
# get embeddings
embedding <- gsva_umap$layout
#gsva_umap_plot_df <- data.frame(X = embedding[,1], Y = embedding[,2], Label_Tissue_Slice = gsva_umap_labels_Tissue_Slice, Label_Navin_Clusters2 = gsva_umap_labels_Navin_Clusters2)
gsva_umap_plot_df <- data.frame(X = embedding[,1], Y = embedding[,2])

# add barcode
gsva_umap_plot_df$Barcode <- rownames(gsva_umap_plot_df)
# merge with annots
gsva_umap_plot_df <- merge(gsva_umap_plot_df, annots, by='Barcode')
# add rownames again
rownames(gsva_umap_plot_df) <- gsva_umap_plot_df$Barcode
```


```R
gsva_umap_plot_df %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 63</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>X</th><th scope=col>Y</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>-1.609071</td><td> 0.8862022</td><td>VA1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td> 3.396564</td><td>-2.2385961</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1_AAACAATCTACTAGCA-1</td><td> 2.815883</td><td> 2.0775624</td><td>VA1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td> 6</td><td> 6</td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Medium           </td><td>vasc_low          </td><td>vasc_low             </td></tr>
</tbody>
</table>




```R
# plot umap
gsva_umap_plot_Tissue_Slice <- ggplot(gsva_umap_plot_df, aes(x = X, y = Y, color = Tissue_Slice)) + geom_point(alpha=0.2)
gsva_umap_plot_Tissue_Slice
```


    
![png](output_85_0.png)
    



```R
# plot umap
gsva_umap_plot_Tissue_Slice <- ggplot(gsva_umap_plot_df, aes(x = X, y = Y, color = Navin_Clusters2)) + geom_point(alpha=0.35)
gsva_umap_plot_Tissue_Slice
```


    
![png](output_86_0.png)
    


# Plot with zoomed in axes


```R
# subset to zoom in
#gsva_umap_plot_df_sub <- gsva_umap_plot_df[(gsva_umap_plot_df$X > -10) & (gsva_umap_plot_df$X < 10),]
gsva_umap_plot_df_sub <- gsva_umap_plot_df[(gsva_umap_plot_df$X < 12) & (gsva_umap_plot_df$Y > -5),]
```


```R
# plot umap
gsva_umap_plot_Tissue_Slice <- ggplot(gsva_umap_plot_df_sub, aes(x = X, y = Y, color = Tissue_Slice)) + geom_point(alpha=0.35)
gsva_umap_plot_Tissue_Slice
```


    
![png](output_89_0.png)
    



```R
# plot umap
gsva_umap_plot_Tissue_Slice <- ggplot(gsva_umap_plot_df_sub, aes(x = X, y = Y, color = Navin_Clusters2)) + geom_point(alpha=0.35)
gsva_umap_plot_Tissue_Slice
```


    
![png](output_90_0.png)
    



```R
# plot umap
gsva_umap_plot_Tissue_Slice <- ggplot(gsva_umap_plot_df_sub, aes(x = X, y = Y, color = vasc_MN4_percentile_label)) + geom_point(alpha=0.35)
gsva_umap_plot_Tissue_Slice
```


    
![png](output_91_0.png)
    



```R
colnames(gsva_umap_plot_df_sub)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Barcode'</li><li>'X'</li><li>'Y'</li><li>'Tissue_Slice'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li></ol>




```R
# plot umap
gsva_umap_plot_Tissue_Slice <- ggplot(gsva_umap_plot_df_sub, aes(x = X, y = Y, color = Tumor_MHCI_label2)) + geom_point(alpha=0.35)
gsva_umap_plot_Tissue_Slice
```


    
![png](output_93_0.png)
    



```R

```


```R

```


```R

```


```R
annots_melt_immune %>% head(4)
```


<table class="dataframe">
<caption>A tibble: 4 × 54</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Proliferation</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Macrophage   </td><td>-0.16512709</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M1_Macrophage</td><td>-0.06733032</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M2_Macrophage</td><td>-0.28232601</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>T_Cell       </td><td> 0.63303241</td></tr>
</tbody>
</table>




```R
options(repr.plot.width=20, repr.plot.height=16)

# make sure to subset for vasc_label = vasc_high
celltype_vs_MN4_EC_vaschigh_plot1 <- ggplot(annots_melt_immune[annots_melt_immune$vasc_label == 'vasc_high',], 
                                          aes(x=MN4_EC_Phenotype_Label, y=Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey', scale = "width", width = 0.9) + 
  geom_beeswarm(cex = 0.5, size=1) +
  #geom_point() + 
  facet_wrap(vars(Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Cell Type GSVA Enrichment vs. MN4_EC_Label, in vasc_high spots")
             
celltype_vs_MN4_EC_vaschigh_plot1

# save plot
png(file = "Output/06_Label_MN4/cell_type_GSVA_per_MN4_EC_Label_for_vasc_high.png",
    width = 1800,
    height = 1400)
print(celltype_vs_MN4_EC_vaschigh_plot1)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_98_1.png)
    



```R
options(repr.plot.width=20, repr.plot.height=16)

# make sure to subset for vasc_label = vasc_high
celltype_vs_MN4_EC_vaschigh_plot1 <- ggplot(annots_melt_immune[annots_melt_immune$vasc_label == 'vasc_high',], 
                                          aes(x=MN4_EC_Phenotype_Label, y=Cell_Type_Enrichment, color=Tissue_Slice)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey', scale = "width", width = 0.9) + 
  geom_beeswarm(cex = 0.5, size=1) +
  #geom_point() + 
  facet_wrap(vars(Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Cell Type GSVA Enrichment vs. MN4_EC_Label, in vasc_high spots")
             
celltype_vs_MN4_EC_vaschigh_plot1

# save plot
png(file = "Output/06_Label_MN4/cell_type_GSVA_per_MN4_EC_Label_for_vasc_high_perTissueSlice.png",
    width = 1800,
    height = 1400)
print(celltype_vs_MN4_EC_vaschigh_plot1)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_99_1.png)
    



```R
options(repr.plot.width=20, repr.plot.height=16)

# subset to only VA1 and VB1
temp_df = annots_melt_immune %>% subset(Tissue_Slice %in% c('VA1', 'VB1'))

# make sure to subset for vasc_label = vasc_high
celltype_vs_MN4_EC_vaschigh_plot1 <- ggplot(temp_df[temp_df$vasc_label == 'vasc_high',], 
                                          aes(x=MN4_EC_Phenotype_Label, y=Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey', scale = "width", width = 0.9) + 
  geom_beeswarm(cex = 0.5, size=1, aes(color=Tissue_Slice)) +
  #geom_point() + 
  facet_wrap(vars(Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Cell Type GSVA Enrichment vs. MN4_EC_Label, in vasc_high spots")
             
celltype_vs_MN4_EC_vaschigh_plot1

# save plot
png(file = "Output/06_Label_MN4/cell_type_GSVA_per_MN4_EC_Label_for_vasc_high_VA1VB1.png",
    width = 1800,
    height = 1400)
print(celltype_vs_MN4_EC_vaschigh_plot1)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_100_1.png)
    



```R
options(repr.plot.width=20, repr.plot.height=16)

# subset to only VA1 and VB1
temp_df = annots_melt_immune %>% subset(Tissue_Slice %in% c('VA2'))

# make sure to subset for vasc_label = vasc_high
celltype_vs_MN4_EC_vaschigh_plot1 <- ggplot(temp_df[temp_df$vasc_label == 'vasc_high',], 
                                          aes(x=MN4_EC_Phenotype_Label, y=Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey', scale = "width", width = 0.9) + 
  geom_beeswarm(cex = 0.5, size=1, aes(color=Tissue_Slice)) +
  #geom_point() + 
  facet_wrap(vars(Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Cell Type GSVA Enrichment vs. MN4_EC_Label, in vasc_high spots")
             
celltype_vs_MN4_EC_vaschigh_plot1

# save plot
png(file = "Output/06_Label_MN4/cell_type_GSVA_per_MN4_EC_Label_for_vasc_high_VA2.png",
    width = 1800,
    height = 1400)
print(celltype_vs_MN4_EC_vaschigh_plot1)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_101_1.png)
    


### Is this because VA2 has mostly vasc_high spots?


```R
dplyr::count(data@meta.data, Tissue_Slice, vasc_label)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>vasc_label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>vasc_high</td><td> 771</td></tr>
	<tr><td>VA1</td><td>vasc_low </td><td>2035</td></tr>
	<tr><td>VA2</td><td>vasc_high</td><td>1733</td></tr>
	<tr><td>VA2</td><td>vasc_low </td><td>2583</td></tr>
	<tr><td>VB1</td><td>vasc_high</td><td> 665</td></tr>
	<tr><td>VB1</td><td>vasc_low </td><td>1924</td></tr>
</tbody>
</table>




```R
colnames(data@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li></ol>




```R
dplyr::count(data@meta.data, Tissue_Slice, MN4_EC_Phenotype_Label)
```


<table class="dataframe">
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>High</td><td>1403</td></tr>
	<tr><td>VA1</td><td>Low </td><td>1403</td></tr>
	<tr><td>VA2</td><td>High</td><td>2158</td></tr>
	<tr><td>VA2</td><td>Low </td><td>2158</td></tr>
	<tr><td>VB1</td><td>High</td><td>1294</td></tr>
	<tr><td>VB1</td><td>Low </td><td>1295</td></tr>
</tbody>
</table>



### I think my median splitting by tissue slice didn't work properly...


```R
dplyr::count(data@meta.data, Tissue_Slice, vasc_MN4_label)
```


<table class="dataframe">
<caption>A data.frame: 9 × 3</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>vasc_MN4_label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>vasc_high_MN4_high</td><td> 511</td></tr>
	<tr><td>VA1</td><td>vasc_high_MN4_low </td><td> 260</td></tr>
	<tr><td>VA1</td><td>vasc_low          </td><td>2035</td></tr>
	<tr><td>VA2</td><td>vasc_high_MN4_high</td><td>1473</td></tr>
	<tr><td>VA2</td><td>vasc_high_MN4_low </td><td> 260</td></tr>
	<tr><td>VA2</td><td>vasc_low          </td><td>2583</td></tr>
	<tr><td>VB1</td><td>vasc_high_MN4_high</td><td> 426</td></tr>
	<tr><td>VB1</td><td>vasc_high_MN4_low </td><td> 239</td></tr>
	<tr><td>VB1</td><td>vasc_low          </td><td>1924</td></tr>
</tbody>
</table>



# Let's subset data to includce only VA2, as a separate analysis. and also VA1 and VB1. (old and new analyses)


```R
annots_melt_immune_VA2 <- annots_melt_immune %>% subset(Tissue_Slice == "VA2")
annots_melt_immune_VAB1 <- annots_melt_immune %>% subset(Tissue_Slice %in% c("VA1", "VB1"))
```


```R
dim(annots_melt_immune)
dim(annots_melt_immune_VA2)
dim(annots_melt_immune_VAB1)
table(annots_melt_immune$Tissue_Slice)
table(annots_melt_immune_VA2$Tissue_Slice)
table(annots_melt_immune_VAB1$Tissue_Slice)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>87399</li><li>54</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>38844</li><li>54</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>48555</li><li>54</li></ol>




    
      VA1   VA2   VB1 
    25254 38844 23301 



    
      VA2 
    38844 



    
      VA1   VB1 
    25254 23301 



```R
annots_melt_immune %>% head(2)
```


<table class="dataframe">
<caption>A tibble: 2 × 54</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Proliferation</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Macrophage   </td><td>-0.16512709</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M1_Macrophage</td><td>-0.06733032</td></tr>
</tbody>
</table>




```R

```

# do stats, use t-test


```R
# do t-test on groups of interest
celltype_vs_MN4_EC_vaschigh_stats1 <- annots_melt_immune[annots_melt_immune$vasc_label %in% c('vasc_high'),] %>% group_by(Cell_Type_Label) %>% summarise(p_VascHigh_MN4_EC_HiVsLo=t.test(Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='High'], Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='Low'], paired=FALSE)$p.value)
celltype_vs_MN4_EC_vaschigh_stats1_VA2 <- annots_melt_immune_VA2[annots_melt_immune_VA2$vasc_label %in% c('vasc_high'),] %>% group_by(Cell_Type_Label) %>% summarise(p_VascHigh_MN4_EC_HiVsLo_VA2=t.test(Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='High'], Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='Low'], paired=FALSE)$p.value)
celltype_vs_MN4_EC_vaschigh_stats1_VAB1 <- annots_melt_immune_VAB1[annots_melt_immune_VAB1$vasc_label %in% c('vasc_high'),] %>% group_by(Cell_Type_Label) %>% summarise(p_VascHigh_MN4_EC_HiVsLo_VAB1=t.test(Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='High'], Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='Low'], paired=FALSE)$p.value)
# add multiple testing corrections
celltype_vs_MN4_EC_vaschigh_stats1 <- celltype_vs_MN4_EC_vaschigh_stats1 %>% mutate(padj_VascHigh_MN4_EC_HiVsLo = p.adjust(p_VascHigh_MN4_EC_HiVsLo, method = "BH"))
celltype_vs_MN4_EC_vaschigh_stats1_VA2 <- celltype_vs_MN4_EC_vaschigh_stats1_VA2 %>% mutate(padj_VascHigh_MN4_EC_HiVsLo_VA2 = p.adjust(p_VascHigh_MN4_EC_HiVsLo_VA2, method = "BH"))
celltype_vs_MN4_EC_vaschigh_stats1_VAB1 <- celltype_vs_MN4_EC_vaschigh_stats1_VAB1 %>% mutate(padj_VascHigh_MN4_EC_HiVsLo_VAB1 = p.adjust(p_VascHigh_MN4_EC_HiVsLo_VAB1, method = "BH"))
# combine results
stats_celltypesigs3_all <- merge(celltype_vs_MN4_EC_vaschigh_stats1, 
                                 celltype_vs_MN4_EC_vaschigh_stats1_VAB1, by.x='Cell_Type_Label', by.y='Cell_Type_Label')
stats_celltypesigs3_all <- merge(stats_celltypesigs3_all, 
                                 celltype_vs_MN4_EC_vaschigh_stats1_VA2, by.x='Cell_Type_Label', by.y='Cell_Type_Label')

print('t-test: MN4_EC high vs MN4_EC low, for vasc_high spots')
stats_celltypesigs3_all
```

    [1] "t-test: MN4_EC high vs MN4_EC low, for vasc_high spots"



<table class="dataframe">
<caption>A data.frame: 9 × 7</caption>
<thead>
	<tr><th scope=col>Cell_Type_Label</th><th scope=col>p_VascHigh_MN4_EC_HiVsLo</th><th scope=col>padj_VascHigh_MN4_EC_HiVsLo</th><th scope=col>p_VascHigh_MN4_EC_HiVsLo_VAB1</th><th scope=col>padj_VascHigh_MN4_EC_HiVsLo_VAB1</th><th scope=col>p_VascHigh_MN4_EC_HiVsLo_VA2</th><th scope=col>padj_VascHigh_MN4_EC_HiVsLo_VA2</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>B_Cell          </td><td> 7.439974e-93</td><td> 3.347988e-92</td><td>8.189999e-30</td><td>3.685499e-29</td><td>5.721564e-37</td><td>1.029882e-36</td></tr>
	<tr><td>CD8_T_Cell      </td><td> 1.987921e-76</td><td> 3.578258e-76</td><td>5.323797e-16</td><td>6.844882e-16</td><td>1.834088e-43</td><td>5.502265e-43</td></tr>
	<tr><td>DC              </td><td> 3.194757e-27</td><td> 3.594102e-27</td><td>9.899605e-25</td><td>2.743982e-24</td><td>3.615764e-12</td><td>3.615764e-12</td></tr>
	<tr><td>Endothelial_Cell</td><td>1.247763e-103</td><td>1.122987e-102</td><td>4.933405e-36</td><td>4.440064e-35</td><td>1.165064e-58</td><td>5.242788e-58</td></tr>
	<tr><td>M1_Macrophage   </td><td> 7.370377e-78</td><td> 1.658335e-77</td><td>3.090578e-22</td><td>5.563040e-22</td><td>1.586618e-37</td><td>3.569891e-37</td></tr>
	<tr><td>M2_Macrophage   </td><td> 5.810371e-72</td><td> 8.715556e-72</td><td>1.219548e-24</td><td>2.743982e-24</td><td>9.136811e-27</td><td>1.370522e-26</td></tr>
	<tr><td>Macrophage      </td><td> 4.077770e-44</td><td> 5.242847e-44</td><td>6.181970e-17</td><td>9.272956e-17</td><td>1.351394e-21</td><td>1.520318e-21</td></tr>
	<tr><td>NK_Cell         </td><td> 6.488350e-22</td><td> 6.488350e-22</td><td>1.167224e-07</td><td>1.167224e-07</td><td>2.615148e-25</td><td>3.362334e-25</td></tr>
	<tr><td>T_Cell          </td><td> 1.122053e-87</td><td> 3.366158e-87</td><td>6.434171e-16</td><td>7.238443e-16</td><td>9.404048e-59</td><td>5.242788e-58</td></tr>
</tbody>
</table>




```R
# Let's get the effect sizes
# Modified t-test code with effect sizes (mean differences)
celltype_vs_MN4_EC_vaschigh_effects1 <- annots_melt_immune[annots_melt_immune$vasc_label %in% c('vasc_high'),] %>% 
  group_by(Cell_Type_Label) %>% 
  summarise(
    effect_size_VascHigh_MN4_EC_HiVsLo = t.test(
      Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='High'],
      Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='Low'],
      paired=FALSE
    )$estimate
  )

# Repeat same pattern for VA2 and VAB1 datasets:
celltype_vs_MN4_EC_vaschigh_effects1_VA2 <- annots_melt_immune_VA2[annots_melt_immune_VA2$vasc_label %in% c('vasc_high'),] %>% 
  group_by(Cell_Type_Label) %>% 
  summarise(
    effect_size_VascHigh_MN4_EC_HiVsLo_VA2 = t.test(
      Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='High'], 
      Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='Low'], 
      paired=FALSE
    )$estimate
  )

celltype_vs_MN4_EC_vaschigh_effects1_VAB1 <- annots_melt_immune_VAB1[annots_melt_immune_VAB1$vasc_label %in% c('vasc_high'),] %>% 
  group_by(Cell_Type_Label) %>% 
  summarise(
    effect_size_VascHigh_MN4_EC_HiVsLo_VAB1 = t.test(
      Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='High'], 
      Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='Low'], 
      paired=FALSE
    )$estimate
  )

# ...repeat for VA2 and VAB1...
# Combine results (will now include effect sizes)
effects_celltypesigs3_all <- merge(celltype_vs_MN4_EC_vaschigh_effects1, 
                                 celltype_vs_MN4_EC_vaschigh_effects1_VAB1, by='Cell_Type_Label')
effects_celltypesigs3_all <- merge(stats_celltypesigs3_all, 
                                 celltype_vs_MN4_EC_vaschigh_effects1_VA2, by='Cell_Type_Label')

print('Effect Sizes: MN4_EC high vs MN4_EC low, for vasc_high spots')
effects_celltypesigs3_all

```

    Warning message:
    “[1m[22mReturning more (or less) than 1 row per `summarise()` group was deprecated in
    dplyr 1.1.0.
    [36mℹ[39m Please use `reframe()` instead.
    [36mℹ[39m When switching from `summarise()` to `reframe()`, remember that `reframe()`
      always returns an ungrouped data frame and adjust accordingly.”
    [1m[22m`summarise()` has grouped output by 'Cell_Type_Label'. You can override using
    the `.groups` argument.
    Warning message:
    “[1m[22mReturning more (or less) than 1 row per `summarise()` group was deprecated in
    dplyr 1.1.0.
    [36mℹ[39m Please use `reframe()` instead.
    [36mℹ[39m When switching from `summarise()` to `reframe()`, remember that `reframe()`
      always returns an ungrouped data frame and adjust accordingly.”
    [1m[22m`summarise()` has grouped output by 'Cell_Type_Label'. You can override using
    the `.groups` argument.
    Warning message:
    “[1m[22mReturning more (or less) than 1 row per `summarise()` group was deprecated in
    dplyr 1.1.0.
    [36mℹ[39m Please use `reframe()` instead.
    [36mℹ[39m When switching from `summarise()` to `reframe()`, remember that `reframe()`
      always returns an ungrouped data frame and adjust accordingly.”
    [1m[22m`summarise()` has grouped output by 'Cell_Type_Label'. You can override using
    the `.groups` argument.


    [1] "Effect Sizes: MN4_EC high vs MN4_EC low, for vasc_high spots"



<table class="dataframe">
<caption>A data.frame: 18 × 8</caption>
<thead>
	<tr><th scope=col>Cell_Type_Label</th><th scope=col>p_VascHigh_MN4_EC_HiVsLo</th><th scope=col>padj_VascHigh_MN4_EC_HiVsLo</th><th scope=col>p_VascHigh_MN4_EC_HiVsLo_VAB1</th><th scope=col>padj_VascHigh_MN4_EC_HiVsLo_VAB1</th><th scope=col>p_VascHigh_MN4_EC_HiVsLo_VA2</th><th scope=col>padj_VascHigh_MN4_EC_HiVsLo_VA2</th><th scope=col>effect_size_VascHigh_MN4_EC_HiVsLo_VA2</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>B_Cell          </td><td> 7.439974e-93</td><td> 3.347988e-92</td><td>8.189999e-30</td><td>3.685499e-29</td><td>5.721564e-37</td><td>1.029882e-36</td><td> 0.33973300</td></tr>
	<tr><td>B_Cell          </td><td> 7.439974e-93</td><td> 3.347988e-92</td><td>8.189999e-30</td><td>3.685499e-29</td><td>5.721564e-37</td><td>1.029882e-36</td><td>-0.05460438</td></tr>
	<tr><td>CD8_T_Cell      </td><td> 1.987921e-76</td><td> 3.578258e-76</td><td>5.323797e-16</td><td>6.844882e-16</td><td>1.834088e-43</td><td>5.502265e-43</td><td> 0.26526060</td></tr>
	<tr><td>CD8_T_Cell      </td><td> 1.987921e-76</td><td> 3.578258e-76</td><td>5.323797e-16</td><td>6.844882e-16</td><td>1.834088e-43</td><td>5.502265e-43</td><td>-0.14418769</td></tr>
	<tr><td>DC              </td><td> 3.194757e-27</td><td> 3.594102e-27</td><td>9.899605e-25</td><td>2.743982e-24</td><td>3.615764e-12</td><td>3.615764e-12</td><td> 0.08341547</td></tr>
	<tr><td>DC              </td><td> 3.194757e-27</td><td> 3.594102e-27</td><td>9.899605e-25</td><td>2.743982e-24</td><td>3.615764e-12</td><td>3.615764e-12</td><td>-0.20807927</td></tr>
	<tr><td>Endothelial_Cell</td><td>1.247763e-103</td><td>1.122987e-102</td><td>4.933405e-36</td><td>4.440064e-35</td><td>1.165064e-58</td><td>5.242788e-58</td><td> 0.35074855</td></tr>
	<tr><td>Endothelial_Cell</td><td>1.247763e-103</td><td>1.122987e-102</td><td>4.933405e-36</td><td>4.440064e-35</td><td>1.165064e-58</td><td>5.242788e-58</td><td> 0.04036004</td></tr>
	<tr><td>M1_Macrophage   </td><td> 7.370377e-78</td><td> 1.658335e-77</td><td>3.090578e-22</td><td>5.563040e-22</td><td>1.586618e-37</td><td>3.569891e-37</td><td> 0.15746554</td></tr>
	<tr><td>M1_Macrophage   </td><td> 7.370377e-78</td><td> 1.658335e-77</td><td>3.090578e-22</td><td>5.563040e-22</td><td>1.586618e-37</td><td>3.569891e-37</td><td>-0.12450132</td></tr>
	<tr><td>M2_Macrophage   </td><td> 5.810371e-72</td><td> 8.715556e-72</td><td>1.219548e-24</td><td>2.743982e-24</td><td>9.136811e-27</td><td>1.370522e-26</td><td> 0.14689041</td></tr>
	<tr><td>M2_Macrophage   </td><td> 5.810371e-72</td><td> 8.715556e-72</td><td>1.219548e-24</td><td>2.743982e-24</td><td>9.136811e-27</td><td>1.370522e-26</td><td>-0.17356731</td></tr>
	<tr><td>Macrophage      </td><td> 4.077770e-44</td><td> 5.242847e-44</td><td>6.181970e-17</td><td>9.272956e-17</td><td>1.351394e-21</td><td>1.520318e-21</td><td> 0.08377533</td></tr>
	<tr><td>Macrophage      </td><td> 4.077770e-44</td><td> 5.242847e-44</td><td>6.181970e-17</td><td>9.272956e-17</td><td>1.351394e-21</td><td>1.520318e-21</td><td>-0.22776030</td></tr>
	<tr><td>NK_Cell         </td><td> 6.488350e-22</td><td> 6.488350e-22</td><td>1.167224e-07</td><td>1.167224e-07</td><td>2.615148e-25</td><td>3.362334e-25</td><td> 0.09192346</td></tr>
	<tr><td>NK_Cell         </td><td> 6.488350e-22</td><td> 6.488350e-22</td><td>1.167224e-07</td><td>1.167224e-07</td><td>2.615148e-25</td><td>3.362334e-25</td><td>-0.14005255</td></tr>
	<tr><td>T_Cell          </td><td> 1.122053e-87</td><td> 3.366158e-87</td><td>6.434171e-16</td><td>7.238443e-16</td><td>9.404048e-59</td><td>5.242788e-58</td><td> 0.24195674</td></tr>
	<tr><td>T_Cell          </td><td> 1.122053e-87</td><td> 3.366158e-87</td><td>6.434171e-16</td><td>7.238443e-16</td><td>9.404048e-59</td><td>5.242788e-58</td><td>-0.19291127</td></tr>
</tbody>
</table>



# save data from this plot for navin and marco


```R
annots_melt_immune %>% head(3)
```


<table class="dataframe">
<caption>A tibble: 3 × 54</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Proliferation</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Macrophage   </td><td>-0.16512709</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M1_Macrophage</td><td>-0.06733032</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M2_Macrophage</td><td>-0.28232601</td></tr>
</tbody>
</table>




```R
vasc_high_data <- annots_melt_immune[annots_melt_immune$vasc_label == 'vasc_high',]
#vasc_high_data_va2 <- annots_melt_immune_VA2[annots_melt_immune_VA2$vasc_label == 'vasc_high',]

vasc_high_data_tosave = vasc_high_data %>% select(Barcode, orig.ident, vasc_label, MN4_EC_Phenotype_Label, 
                                                    GOBP_LEUKOCYTE_ADHESION_TO_VASC, 
                                                    GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION, MN4_EC_Phenotype,
                                                  Cell_Type_Label, Cell_Type_Enrichment)
#vasc_high_data_tosave_va2 = vasc_high_data_va2 %>% select(Barcode, orig.ident, vasc_label, MN4_EC_Phenotype_Label, Cell_Type_Label, Cell_Type_Enrichment)

# drop duplicates
vasc_high_data_tosave <- vasc_high_data_tosave[!duplicated(vasc_high_data_tosave), ]
#vasc_high_data_tosave_va2 <- vasc_high_data_tosave_va2[!duplicated(vasc_high_data_tosave_va2), ]

vasc_high_data_tosave[1:3,]
```


<table class="dataframe">
<caption>A tibble: 3 × 9</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>vasc_label</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>-0.04649791</td><td>0.102126</td><td>0.1894989</td><td>Macrophage   </td><td>-0.16512709</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>-0.04649791</td><td>0.102126</td><td>0.1894989</td><td>M1_Macrophage</td><td>-0.06733032</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>-0.04649791</td><td>0.102126</td><td>0.1894989</td><td>M2_Macrophage</td><td>-0.28232601</td></tr>
</tbody>
</table>




```R
dplyr::count(vasc_high_data_tosave, orig.ident, MN4_EC_Phenotype_Label)
```


<table class="dataframe">
<caption>A tibble: 6 × 3</caption>
<thead>
	<tr><th scope=col>orig.ident</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>High</td><td> 4599</td></tr>
	<tr><td>VA1</td><td>Low </td><td> 2340</td></tr>
	<tr><td>VA2</td><td>High</td><td>13257</td></tr>
	<tr><td>VA2</td><td>Low </td><td> 2340</td></tr>
	<tr><td>VB1</td><td>High</td><td> 3834</td></tr>
	<tr><td>VB1</td><td>Low </td><td> 2151</td></tr>
</tbody>
</table>




```R
# save
write.csv(vasc_high_data_tosave, 'Output/06_Label_MN4/vasc_high_cell_type_enrichment_by_MN4_EC.csv', row.names = FALSE)
#write.csv(vasc_high_data_tosave_va2, 'Output/06_Label_MN4/vasc_high_cell_type_enrichment_by_MN4_EC_VA2only.csv', row.names = FALSE)

```


```R
# Save stats
# stats_celltypesigs3_all
write.csv(stats_celltypesigs3_all, 'Output/06_Label_MN4/STATS_vasc_high_cell_type_enrichment_by_MN4_EC.csv', row.names = FALSE)

```

# save a version that's unmelted


```R
vasc_high_data_tosave_unmelted <- vasc_high_data_tosave %>% pivot_wider(names_from = Cell_Type_Label, values_from = Cell_Type_Enrichment)
```


```R
vasc_high_data_tosave_unmelted <- vasc_high_data_tosave_unmelted[!duplicated(vasc_high_data_tosave_unmelted), ]
```


```R
vasc_high_data_tosave_unmelted[1:3,]
```


<table class="dataframe">
<caption>A tibble: 3 × 16</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>vasc_label</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>Macrophage</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>T_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>NK_Cell</th><th scope=col>B_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Cell</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>-0.04649791</td><td> 0.1021260</td><td> 0.1894989</td><td>-0.1651271</td><td>-0.06733032</td><td>-0.2823260</td><td> 0.63303241</td><td> 0.5135266</td><td>0.0155587</td><td>-0.03843072</td><td>-0.007136273</td><td>0.4046568</td></tr>
	<tr><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>-0.24403536</td><td>-0.0993218</td><td> 0.1598434</td><td> 0.6304333</td><td>-0.13408735</td><td> 0.1162332</td><td> 0.87197296</td><td> 0.5112075</td><td>0.5986082</td><td>-0.68264193</td><td> 0.351200349</td><td>0.5369628</td></tr>
	<tr><td>VA1_AAACCCGAACGAAATC-1</td><td>VA1</td><td>vasc_high</td><td>Low </td><td> 0.10929027</td><td>-0.2076077</td><td>-0.2645908</td><td>-0.6906275</td><td>-0.26293664</td><td>-0.1793283</td><td>-0.07225781</td><td>-0.1138683</td><td>0.5001156</td><td> 0.18982041</td><td>-0.443574962</td><td>0.2452949</td></tr>
</tbody>
</table>




```R
# save
write.csv(vasc_high_data_tosave_unmelted, 'Output/06_Label_MN4/vasc_high_cell_type_enrichment_by_MN4_EC_wideform.csv', row.names = FALSE)

```


```R
# also save the data that's not only vasc high...

vasc_all_data <- annots_melt_immune
#vasc_all_data_va2 <- annots_melt_immune_VA2

vasc_all_data_tosave = vasc_all_data %>% select(Barcode, orig.ident, vasc_label, MN4_EC_Phenotype_Label, MN4_EC_Phenotype, Cell_Type_Label, Cell_Type_Enrichment)
#vasc_all_data_tosave_va2 = vasc_all_data_va2 %>% select(Barcode, orig.ident, vasc_label, MN4_EC_Phenotype_Label, Cell_Type_Label, Cell_Type_Enrichment)

# drop duplicates
vasc_all_data_tosave <- vasc_all_data_tosave[!duplicated(vasc_all_data_tosave), ]
#vasc_all_data_tosave_va2 <- vasc_all_data_tosave_va2[!duplicated(vasc_all_data_tosave_va2), ]


vasc_all_data_tosave[1:3,]

```


<table class="dataframe">
<caption>A tibble: 3 × 7</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>vasc_label</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>0.1894989</td><td>Macrophage   </td><td>-0.16512709</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>0.1894989</td><td>M1_Macrophage</td><td>-0.06733032</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>0.1894989</td><td>M2_Macrophage</td><td>-0.28232601</td></tr>
</tbody>
</table>




```R
unique(vasc_all_data_tosave$vasc_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_high'</li><li>'vasc_low'</li></ol>




```R
# save
write.csv(vasc_all_data_tosave, 'Output/06_Label_MN4/vasc_all_cell_type_enrichment_by_MN4_EC.csv', row.names = FALSE)
#write.csv(vasc_all_data_tosave_va2, 'Output/06_Label_MN4/vasc_all_cell_type_enrichment_by_MN4_EC_VA2only.csv', row.names = FALSE)
```


```R

```

# Let's redo the cell type enrichment plots from Notebook 16, with the Updated GSVA enrichments from April
# Want to see cell type enrichments for Tumor High, comparing MHCI Low vs MHCI High


```R
# check
annots_melt_immune %>% head(3)
```


<table class="dataframe">
<caption>A tibble: 3 × 54</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Proliferation</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Macrophage   </td><td>-0.16512709</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M1_Macrophage</td><td>-0.06733032</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M2_Macrophage</td><td>-0.28232601</td></tr>
</tbody>
</table>




```R
unique(annots_melt_immune$Tumor_label)
unique(annots_melt_immune$MHCI_label)
unique(annots_melt_immune$Tumor_MHCI_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'MHCI_low'</li><li>'MHCI_high'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low_MHCI_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_low_MHCI_high'</li><li>'Tumor_high_MHCI_high'</li></ol>




```R
unique(annots_melt_immune$Cell_Type_Label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Macrophage'</li><li>'M1_Macrophage'</li><li>'M2_Macrophage'</li><li>'T_Cell'</li><li>'CD8_T_Cell'</li><li>'NK_Cell'</li><li>'B_Cell'</li><li>'DC'</li><li>'Endothelial_Cell'</li></ol>




```R
options(repr.plot.width=20, repr.plot.height=16)

# make sure to subset for vasc_label = vasc_high
celltype_vs_MHCI_tumorhigh_plot1 <- ggplot(annots_melt_immune[(annots_melt_immune$Tumor_label == 'Tumor_high') & (annots_melt_immune$MHCI_label %in% c('MHCI_low','MHCI_high')),], 
                                          aes(x=Tumor_MHCI_label, y=Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey', scale = "width", width = 0.9) + 
  geom_beeswarm(cex = 0.4, size=1) +
  #geom_point() + 
  facet_wrap(vars(Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Cell Type GSVA Enrichment for MHCI High vs. MHCI Low, in Tumor High spots")
             
celltype_vs_MHCI_tumorhigh_plot1

# save plot
png(file = "Output/06_Label_MN4/cell_type_GSVA_for_MHCI_High_vs_Low_within_Tumor_High.png",
    width = 1800,
    height = 1400)
print(celltype_vs_MHCI_tumorhigh_plot1)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_135_1.png)
    



```R
options(repr.plot.width=20, repr.plot.height=16)

# make sure to subset for vasc_label = vasc_high
celltype_vs_MHCI_tumorhigh_plot1 <- ggplot(annots_melt_immune[(annots_melt_immune$Tumor_label == 'Tumor_high') & (annots_melt_immune$MHCI_label %in% c('MHCI_low','MHCI_high')),], 
                                          aes(x=Tumor_MHCI_label, y=Cell_Type_Enrichment, color=Tissue_Slice)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey', scale = "width", width = 0.9) + 
  geom_beeswarm(cex = 0.4, size=1) +
  #geom_point() + 
  facet_wrap(vars(Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Cell Type GSVA Enrichment for MHCI High vs. MHCI Low, in Tumor High spots")
             
celltype_vs_MHCI_tumorhigh_plot1

# save plot
png(file = "Output/06_Label_MN4/cell_type_GSVA_for_MHCI_High_vs_Low_within_Tumor_High_perTissueSlice.png",
    width = 1800,
    height = 1400)
print(celltype_vs_MHCI_tumorhigh_plot1)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_136_1.png)
    



```R
options(repr.plot.width=20, repr.plot.height=16)

# make sure to subset for vasc_label = vasc_high
celltype_vs_MHCI_tumorhigh_plot1 <- ggplot(annots_melt_immune_VA2[(annots_melt_immune_VA2$Tumor_label == 'Tumor_high') & (annots_melt_immune_VA2$MHCI_label %in% c('MHCI_low','MHCI_high')),], 
                                          aes(x=Tumor_MHCI_label, y=Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey', scale = "width", width = 0.9) + 
  geom_beeswarm(cex = 0.4, size=1, aes(color=Tissue_Slice)) +
  #geom_point() + 
  facet_wrap(vars(Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Cell Type GSVA Enrichment for MHCI High vs. MHCI Low, in Tumor High spots")
             
celltype_vs_MHCI_tumorhigh_plot1

# save plot
png(file = "Output/06_Label_MN4/cell_type_GSVA_for_MHCI_High_vs_Low_within_Tumor_High_VA2.png",
    width = 1800,
    height = 1400)
print(celltype_vs_MHCI_tumorhigh_plot1)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_137_1.png)
    



```R
options(repr.plot.width=20, repr.plot.height=16)

# make sure to subset for vasc_label = vasc_high
celltype_vs_MHCI_tumorhigh_plot1 <- ggplot(annots_melt_immune_VAB1[(annots_melt_immune_VAB1$Tumor_label == 'Tumor_high') & (annots_melt_immune_VAB1$MHCI_label %in% c('MHCI_low','MHCI_high')),], 
                                          aes(x=Tumor_MHCI_label, y=Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey', scale = "width", width = 0.9) + 
  geom_beeswarm(cex = 0.4, size=1, aes(color=Tissue_Slice)) +
  #geom_point() + 
  facet_wrap(vars(Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Cell Type GSVA Enrichment for MHCI High vs. MHCI Low, in Tumor High spots")
             
celltype_vs_MHCI_tumorhigh_plot1

# save plot
png(file = "Output/06_Label_MN4/cell_type_GSVA_for_MHCI_High_vs_Low_within_Tumor_High_VAB1.png",
    width = 1800,
    height = 1400)
print(celltype_vs_MHCI_tumorhigh_plot1)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_138_1.png)
    


# Do Stats, use t-test


```R
annots_melt_immune %>% head(2)
```


<table class="dataframe">
<caption>A tibble: 2 × 54</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Proliferation</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Macrophage   </td><td>-0.16512709</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M1_Macrophage</td><td>-0.06733032</td></tr>
</tbody>
</table>




```R
unique(annots_melt_immune$MHCI_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'MHCI_low'</li><li>'MHCI_high'</li></ol>




```R
# do t-test on groups of interest
celltype_vs_MHCI_tumorhigh_stats1 <- annots_melt_immune[annots_melt_immune$Tumor_label %in% c('Tumor_high'),] %>% group_by(Cell_Type_Label) %>% summarise(p_TumorHigh_MHCI_HiVsLo=t.test(Cell_Type_Enrichment[MHCI_label=='MHCI_high'], Cell_Type_Enrichment[MHCI_label=='MHCI_low'], paired=FALSE)$p.value)
celltype_vs_MHCI_tumorhigh_stats1_VA2 <- annots_melt_immune_VA2[annots_melt_immune_VA2$Tumor_label %in% c('Tumor_high'),] %>% group_by(Cell_Type_Label) %>% summarise(p_TumorHigh_MHCI_HiVsLo_VA2=t.test(Cell_Type_Enrichment[MHCI_label=='MHCI_high'], Cell_Type_Enrichment[MHCI_label=='MHCI_low'], paired=FALSE)$p.value)
celltype_vs_MHCI_tumorhigh_stats1_VAB1 <- annots_melt_immune_VAB1[annots_melt_immune_VAB1$Tumor_label %in% c('Tumor_high'),] %>% group_by(Cell_Type_Label) %>% summarise(p_TumorHigh_MHCI_HiVsLo_VAB1=t.test(Cell_Type_Enrichment[MHCI_label=='MHCI_high'], Cell_Type_Enrichment[MHCI_label=='MHCI_low'], paired=FALSE)$p.value)
# add multiple testing corrections
celltype_vs_MHCI_tumorhigh_stats1 <- celltype_vs_MHCI_tumorhigh_stats1 %>% mutate(padj_TumorHigh_MHCI_HiVsLo = p.adjust(p_TumorHigh_MHCI_HiVsLo, method = "BH"))
celltype_vs_MHCI_tumorhigh_stats1_VA2 <- celltype_vs_MHCI_tumorhigh_stats1_VA2 %>% mutate(padj_TumorHigh_MHCI_HiVsLo_VA2 = p.adjust(p_TumorHigh_MHCI_HiVsLo_VA2, method = "BH"))
celltype_vs_MHCI_tumorhigh_stats1_VAB1 <- celltype_vs_MHCI_tumorhigh_stats1_VAB1 %>% mutate(padj_TumorHigh_MHCI_HiVsLo_VAB1 = p.adjust(p_TumorHigh_MHCI_HiVsLo_VAB1, method = "BH"))
# combine results
stats_celltypesigs4_all <- merge(celltype_vs_MHCI_tumorhigh_stats1, 
                                 celltype_vs_MHCI_tumorhigh_stats1_VAB1, by.x='Cell_Type_Label', by.y='Cell_Type_Label')
stats_celltypesigs4_all <- merge(stats_celltypesigs4_all, 
                                 celltype_vs_MHCI_tumorhigh_stats1_VA2, by.x='Cell_Type_Label', by.y='Cell_Type_Label')

print('t-test: MHCI high vs MHCI low, for Tumor High spots')
stats_celltypesigs4_all
```

    [1] "t-test: MHCI high vs MHCI low, for Tumor High spots"



<table class="dataframe">
<caption>A data.frame: 9 × 7</caption>
<thead>
	<tr><th scope=col>Cell_Type_Label</th><th scope=col>p_TumorHigh_MHCI_HiVsLo</th><th scope=col>padj_TumorHigh_MHCI_HiVsLo</th><th scope=col>p_TumorHigh_MHCI_HiVsLo_VAB1</th><th scope=col>padj_TumorHigh_MHCI_HiVsLo_VAB1</th><th scope=col>p_TumorHigh_MHCI_HiVsLo_VA2</th><th scope=col>padj_TumorHigh_MHCI_HiVsLo_VA2</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>B_Cell          </td><td> 6.717752e-57</td><td> 1.007663e-56</td><td>1.881302e-14</td><td>1.693172e-13</td><td>1.313121e-57</td><td>1.688298e-57</td></tr>
	<tr><td>CD8_T_Cell      </td><td>3.265949e-101</td><td>1.469677e-100</td><td>2.313809e-07</td><td>2.974897e-07</td><td>2.514214e-89</td><td>7.542642e-89</td></tr>
	<tr><td>DC              </td><td> 1.689213e-46</td><td> 1.689213e-46</td><td>1.801365e-13</td><td>8.106145e-13</td><td>6.694807e-27</td><td>6.694807e-27</td></tr>
	<tr><td>Endothelial_Cell</td><td> 4.572097e-77</td><td> 1.028722e-76</td><td>1.715562e-07</td><td>2.573343e-07</td><td>2.267198e-65</td><td>4.080956e-65</td></tr>
	<tr><td>M1_Macrophage   </td><td> 5.060872e-79</td><td> 1.518262e-78</td><td>6.278616e-11</td><td>1.130151e-10</td><td>2.492140e-84</td><td>5.607315e-84</td></tr>
	<tr><td>M2_Macrophage   </td><td> 1.719620e-71</td><td> 3.095316e-71</td><td>2.927126e-11</td><td>6.586033e-11</td><td>1.321807e-92</td><td>5.948133e-92</td></tr>
	<tr><td>Macrophage      </td><td> 3.208567e-53</td><td> 4.125300e-53</td><td>1.362108e-05</td><td>1.532372e-05</td><td>4.505411e-64</td><td>6.758116e-64</td></tr>
	<tr><td>NK_Cell         </td><td> 2.756668e-48</td><td> 3.101251e-48</td><td>1.493250e-02</td><td>1.493250e-02</td><td>2.333537e-46</td><td>2.625230e-46</td></tr>
	<tr><td>T_Cell          </td><td>3.869846e-118</td><td>3.482861e-117</td><td>5.102459e-13</td><td>1.530738e-12</td><td>6.268058e-98</td><td>5.641252e-97</td></tr>
</tbody>
</table>



# Measure effect size (difference in means)


```R

```


```R
annots_melt_immune %>% head(2)
```


<table class="dataframe">
<caption>A tibble: 2 × 54</caption>
<thead>
	<tr><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Proliferation</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>Macrophage   </td><td>-0.16512709</td></tr>
	<tr><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>0.4121959</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.7552172</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>M1_Macrophage</td><td>-0.06733032</td></tr>
</tbody>
</table>




```R
unique(annots_melt_immune$vasc_label)
unique(annots_melt_immune$vasc_label1)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_high'</li><li>'vasc_low'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_high'</li><li>'vasc_low'</li><li>'vasc_neg'</li></ol>



# save data from this plot for navin and marco


```R
tumor_high_data <- annots_melt_immune[(annots_melt_immune$Tumor_label == 'Tumor_high'),]

tumor_high_data_tosave = tumor_high_data %>% select(Barcode, orig.ident, vasc_label, Tumor_label, MHCI_label, Tumor_MHCI_label, 
                                                    GOBP_LEUKOCYTE_ADHESION_TO_VASC, 
                                                    GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION, MN4_EC_Phenotype, 
                                                    Cell_Type_Label, Cell_Type_Enrichment)

# drop duplicates
tumor_high_data_tosave <- tumor_high_data_tosave[!duplicated(tumor_high_data_tosave), ]

tumor_high_data_tosave[1:3,]
```


<table class="dataframe">
<caption>A tibble: 3 × 11</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>vasc_label</th><th scope=col>Tumor_label</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>vasc_low</td><td>Tumor_high</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>-0.2969164</td><td>-0.1397081</td><td>-0.03932724</td><td>Macrophage   </td><td>-0.7073804</td></tr>
	<tr><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>vasc_low</td><td>Tumor_high</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>-0.2969164</td><td>-0.1397081</td><td>-0.03932724</td><td>M1_Macrophage</td><td>-0.4601281</td></tr>
	<tr><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>vasc_low</td><td>Tumor_high</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>-0.2969164</td><td>-0.1397081</td><td>-0.03932724</td><td>M2_Macrophage</td><td>-0.6018211</td></tr>
</tbody>
</table>




```R
dplyr::count(tumor_high_data_tosave, orig.ident, Tumor_MHCI_label)
```


<table class="dataframe">
<caption>A tibble: 6 × 3</caption>
<thead>
	<tr><th scope=col>orig.ident</th><th scope=col>Tumor_MHCI_label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>Tumor_high_MHCI_high</td><td> 5094</td></tr>
	<tr><td>VA1</td><td>Tumor_high_MHCI_low </td><td> 7389</td></tr>
	<tr><td>VA2</td><td>Tumor_high_MHCI_high</td><td> 6822</td></tr>
	<tr><td>VA2</td><td>Tumor_high_MHCI_low </td><td>17136</td></tr>
	<tr><td>VB1</td><td>Tumor_high_MHCI_high</td><td> 5211</td></tr>
	<tr><td>VB1</td><td>Tumor_high_MHCI_low </td><td> 7695</td></tr>
</tbody>
</table>




```R
dim(tumor_high_data_tosave)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>49347</li><li>11</li></ol>


#ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/06_Label_MN4/DimPlot_MHCI_Label_VA2.pdf

```R
# save
write.csv(tumor_high_data_tosave, 'Output/06_Label_MN4/tumor_high_cell_type_enrichment_by_MHCI.csv', row.names = FALSE)

```


```R
# Save stats
# stats_celltypesigs3_all
write.csv(stats_celltypesigs4_all, 'Output/06_Label_MN4/STATS_tumor_high_cell_type_enrichment_by_MHCI.csv', row.names = FALSE)

```

# save a version that's unmelted


```R
tumor_high_data_tosave_unmelted <- tumor_high_data_tosave %>% pivot_wider(names_from = Cell_Type_Label, values_from = Cell_Type_Enrichment)
```


```R
tumor_high_data_tosave_unmelted <- tumor_high_data_tosave_unmelted[!duplicated(tumor_high_data_tosave_unmelted), ]
```


```R
tumor_high_data_tosave_unmelted[1:3,]
```


<table class="dataframe">
<caption>A tibble: 3 × 18</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>vasc_label</th><th scope=col>Tumor_label</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>Macrophage</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>T_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>NK_Cell</th><th scope=col>B_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Cell</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>vasc_low</td><td>Tumor_high</td><td>MHCI_low </td><td>Tumor_high_MHCI_low </td><td>-0.29691635</td><td>-0.13970812</td><td>-0.03932724</td><td>-0.7073804</td><td>-0.4601281</td><td>-0.6018211</td><td>-0.3041263</td><td>-0.2544987</td><td>-0.09144675</td><td>-0.5228728</td><td> 0.3441700</td><td> 0.01578354</td></tr>
	<tr><td>VA1_AAACAGGGTCTATATT-1</td><td>VA1</td><td>vasc_low</td><td>Tumor_high</td><td>MHCI_low </td><td>Tumor_high_MHCI_low </td><td>-0.04150094</td><td>-0.08681175</td><td>-0.19693860</td><td>-0.3770171</td><td>-0.5632510</td><td>-0.2699566</td><td>-0.3726310</td><td>-0.3717230</td><td> 0.18796917</td><td>-0.6948350</td><td>-0.2468816</td><td>-0.36918291</td></tr>
	<tr><td>VA1_AAACCGGGTAGGTACC-1</td><td>VA1</td><td>vasc_low</td><td>Tumor_high</td><td>MHCI_high</td><td>Tumor_high_MHCI_high</td><td> 0.26347300</td><td>-0.09219470</td><td> 0.08827364</td><td>-0.7079362</td><td> 0.2037043</td><td>-0.3680926</td><td> 0.1315690</td><td>-0.2706203</td><td> 0.42364524</td><td>-0.1603960</td><td> 0.4476608</td><td>-0.24471808</td></tr>
</tbody>
</table>




```R
# save
write.csv(tumor_high_data_tosave_unmelted, 'Output/06_Label_MN4/tumor_high_cell_type_enrichment_by_MHCI_wideform.csv', row.names = FALSE)

```


```R
# also save the data that's not only tumor high, MHCI low or high...

tumor_high_data <- annots_melt_immune

tumor_high_data_tosave = tumor_high_data %>% select(Barcode, orig.ident, vasc_label, Tumor_label, MHCI_label, Tumor_MHCI_label, Cell_Type_Label, Cell_Type_Enrichment)

# drop duplicates
tumor_high_data_tosave <- tumor_high_data_tosave[!duplicated(tumor_high_data_tosave), ]

tumor_high_data_tosave[1:3,]
```


<table class="dataframe">
<caption>A tibble: 3 × 8</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>vasc_label</th><th scope=col>Tumor_label</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Cell_Type_Label</th><th scope=col>Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>vasc_high</td><td>Tumor_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Macrophage   </td><td>-0.16512709</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>vasc_high</td><td>Tumor_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>M1_Macrophage</td><td>-0.06733032</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>vasc_high</td><td>Tumor_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>M2_Macrophage</td><td>-0.28232601</td></tr>
</tbody>
</table>




```R
unique(tumor_high_data_tosave$Tumor_MHCI_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low_MHCI_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_low_MHCI_high'</li><li>'Tumor_high_MHCI_high'</li></ol>




```R
# save
write.csv(tumor_high_data_tosave, 'Output/06_Label_MN4/tumor_all_cell_type_enrichment_by_MHCI.csv', row.names = FALSE)
```


```R

```

# Plot vasc_label, MN4_EC_Phenotype_label, and combined label


```R
data@meta.data %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 61</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.3752852</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.6186324</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td></tr>
</tbody>
</table>




```R
unique(data@meta.data$vasc_label1)
unique(data@meta.data$vasc_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_high'</li><li>'vasc_low'</li><li>'vasc_neg'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_high'</li><li>'vasc_low'</li></ol>




```R
vasc_colors <- c('vasc_neg' = dittoColors()[8], 'vasc_low' = dittoColors()[2], 'vasc_high' = dittoColors()[4])

# spatial dimplot colored by navin annotations
Idents(data) <- "vasc_label1"

options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialDimPlot(data, cols = vasc_colors, pt.size.factor=2.6, image.alpha = 0.95, stroke = 0.3)
plot_VA2 <- plot2[[1]] + 
            labs(title="Vasc Label VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="Vasc Label VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="Vasc Label VB1")
plot_VA2
plot_VA1
plot_VB1

png(file = "Output/06_Label_MN4/DimPlot_Vasc_Label_VA2.png",
    width = 1000,
    height = 800)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_Vasc_Label_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_Vasc_Label_VA1.png",
    width = 1000,
    height = 800)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_Vasc_Label_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_Vasc_Label_VB1.png",
    width = 1000,
    height = 800)
print(plot_VB1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_Vasc_Label_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_166_0.png)
    



    
![png](output_166_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_166_8.png)
    



```R
unique(data@meta.data$MN4_EC_Phenotype_Label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'High'</li><li>'Low'</li></ol>




```R
MN4_EC_colors <- c('Low' = dittoColors()[2], 'High' = dittoColors()[4])

# spatial dimplot colored by navin annotations
Idents(data) <- "MN4_EC_Phenotype_Label"

options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialDimPlot(data, cols = MN4_EC_colors, pt.size.factor=2.6, image.alpha = 0.95, stroke = 0.3)
plot_VA2 <- plot2[[1]] + 
            labs(title="MN4_EC Label VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="MN4_EC Label VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="MN4_EC Label VB1")
plot_VA2
plot_VA1
plot_VB1

png(file = "Output/06_Label_MN4/DimPlot_MN4_EC_Label_VA2.png",
    width = 1000,
    height = 800)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_MN4_EC_Label_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_MN4_EC_Label_VA1.png",
    width = 1000,
    height = 800)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_MN4_EC_Label_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_MN4_EC_Label_VB1.png",
    width = 1000,
    height = 800)
print(plot_VB1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_MN4_EC_Label_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_168_0.png)
    



    
![png](output_168_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_168_8.png)
    



```R
unique(data@meta.data$vasc_MN4_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_high_MN4_high'</li><li>'vasc_low'</li><li>'vasc_high_MN4_low'</li></ol>




```R
Vasc_MN4_colors <- c('vasc_low' = '#9d9d9d',
                     'vasc_high_MN4_low' = '#FC8D59FF', 
                     'vasc_high_MN4_high' = '#B30000FF')

# spatial dimplot colored by navin annotations
Idents(data) <- "vasc_MN4_label"

options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialDimPlot(data, cols = Vasc_MN4_colors, pt.size.factor=2.6, image.alpha = 0.95, stroke = 0.3)
plot_VA2 <- plot2[[1]] + 
            labs(title="vasc_MN4_EC Label VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="vasc_MN4_EC Label VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="vasc_MN4_EC Label VB1")
plot_VA2
plot_VA1
plot_VB1

png(file = "Output/06_Label_MN4/DimPlot_vasc_MN4_EC_Label_VA2.png",
    width = 1000,
    height = 800)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_vasc_MN4_EC_Label_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_vasc_MN4_EC_Label_VA1.png",
    width = 1000,
    height = 800)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_vasc_MN4_EC_Label_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_vasc_MN4_EC_Label_VB1.png",
    width = 1000,
    height = 800)
print(plot_VB1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_vasc_MN4_EC_Label_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_170_0.png)
    



    
![png](output_170_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_170_8.png)
    


# Now let's plot Tumor_label, MHCI_label, and combined


```R
data@meta.data %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 61</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.3752852</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.6186324</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td></tr>
</tbody>
</table>




```R
unique(data@meta.data$Tumor_label1)
unique(data@meta.data$Tumor_label)
unique(data@meta.data$MHCI_label1)
unique(data@meta.data$MHCI_label)
unique(data@meta.data$Tumor_MHCI_label)
unique(data@meta.data$Tumor_MHCI_label2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high'</li><li>'Tumor_neg'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'MHCI_low'</li><li>'MHCI_high'</li><li>'MHCI_neg'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'MHCI_low'</li><li>'MHCI_high'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low_MHCI_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_low_MHCI_high'</li><li>'Tumor_high_MHCI_high'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_high_MHCI_high'</li></ol>




```R
tumor_colors <- c('Tumor_low' = dittoColors()[2], 'Tumor_high' = dittoColors()[4])

# spatial dimplot colored by navin annotations
Idents(data) <- "Tumor_label"

options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialDimPlot(data, cols = tumor_colors, pt.size.factor=2.6, image.alpha = 0.95, stroke = 0.3)
plot_VA2 <- plot2[[1]] + 
            labs(title="Tumor Label VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="Tumor Label VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="Tumor Label VB1")
plot_VA2
plot_VA1
plot_VB1

png(file = "Output/06_Label_MN4/DimPlot_Tumor_Label_VA2.png",
    width = 1000,
    height = 800)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_Tumor_Label_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_Tumor_Label_VA1.png",
    width = 1000,
    height = 800)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_Tumor_Label_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_Tumor_Label_VB1.png",
    width = 1000,
    height = 800)
print(plot_VB1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_Tumor_Label_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_174_0.png)
    



    
![png](output_174_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_174_8.png)
    



```R
unique(data@meta.data$MHCI_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'MHCI_low'</li><li>'MHCI_high'</li></ol>




```R
MHCI_colors <- c('MHCI_low' = dittoColors()[2], 'MHCI_high' = dittoColors()[4])

# spatial dimplot colored by navin annotations
Idents(data) <- "MHCI_label"

options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialDimPlot(data, cols = MHCI_colors, pt.size.factor=2.6, image.alpha = 0.95, stroke = 0.3)
plot_VA2 <- plot2[[1]] + 
            labs(title="MHCI Label VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="MHCI Label VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="MHCI Label VB1")
plot_VA2
plot_VA1
plot_VB1

png(file = "Output/06_Label_MN4/DimPlot_MHCI_Label_VA2.png",
    width = 1000,
    height = 800)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_MHCI_Label_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_MHCI_Label_VA1.png",
    width = 1000,
    height = 800)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_MHCI_Label_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_MHCI_Label_VB1.png",
    width = 1000,
    height = 800)
print(plot_VB1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_MHCI_Label_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_176_0.png)
    



    
![png](output_176_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_176_8.png)
    



```R
unique(data@meta.data$Tumor_MHCI_label2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_high_MHCI_high'</li></ol>




```R
Vasc_MN4_colors <- c('Tumor_low' = '#9d9d9d',
                     'Tumor_high_MHCI_low' = '#FC8D59FF', 
                     'Tumor_high_MHCI_high' = '#B30000FF')

# spatial dimplot colored by navin annotations
Idents(data) <- "Tumor_MHCI_label2"

options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialDimPlot(data, cols = Vasc_MN4_colors, pt.size.factor=2.6, image.alpha = 0.95, stroke = 0.3)
plot_VA2 <- plot2[[1]] + 
            labs(title="Tumor_MHCI Label VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="Tumor_MHCI Label VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="Tumor_MHCI Label VB1")
plot_VA2
plot_VA1
plot_VB1

png(file = "Output/06_Label_MN4/DimPlot_Tumor_MHCI_Label_VA2.png",
    width = 1000,
    height = 800)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_Tumor_MHCI_Label_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_Tumor_MHCI_Label_VA1.png",
    width = 1000,
    height = 800)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_Tumor_MHCI_Label_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

png(file = "Output/06_Label_MN4/DimPlot_Tumor_MHCI_Label_VB1.png",
    width = 1000,
    height = 800)
print(plot_VB1)
dev.off()

pdf(file = "Output/06_Label_MN4/DimPlot_Tumor_MHCI_Label_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_178_0.png)
    



    
![png](output_178_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_178_8.png)
    


# Plot SCORE of Tumor and MHCI


```R
colnames(data@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li></ol>




```R
# spatial feature plot colored by navin annotations
options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialFeaturePlot(data, features=c('Tumor_score'), pt.size.factor=2)
plot_VA2 <- plot2[[1]] + 
            labs(title="Tumor_Score VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="Tumor_Score VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="Tumor_Score VB1")
plot_VA2
plot_VA1
plot_VB1

pdf(file = "Output/06_Label_MN4/FeaturePlot_Tumor_Score_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/FeaturePlot_Tumor_Score_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/FeaturePlot_Tumor_Score_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_181_0.png)
    



    
![png](output_181_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_181_5.png)
    



```R
# spatial feature plot colored by navin annotations
options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialFeaturePlot(data, features=c('MHCI_score'), pt.size.factor=2)
plot_VA2 <- plot2[[1]] + 
            labs(title="MHCI_score VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="MHCI_score VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="MHCI_score VB1")
plot_VA2
plot_VA1
plot_VB1

pdf(file = "Output/06_Label_MN4/FeaturePlot_MHCI_score_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/FeaturePlot_MHCI_score_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/FeaturePlot_MHCI_score_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_182_0.png)
    



    
![png](output_182_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_182_5.png)
    


# Plot IGF2BP1 (navin request 20250508)


```R
# spatial feature plot colored by navin annotations
options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialFeaturePlot(data, features=c('IGF2BP1'), pt.size.factor=2)
plot_VA2 <- plot2[[1]] + 
            labs(title="IGF2BP1 VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="IGF2BP1 VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="IGF2BP1 VB1")
plot_VA2
plot_VA1
plot_VB1

png(file = "Output/06_Label_MN4/FeaturePlot_IGF2BP1_VA2.png",
    width = 1000,
    height = 800)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/FeaturePlot_IGF2BP1_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

png(file = "Output/06_Label_MN4/FeaturePlot_IGF2BP1_VA1.png",
    width = 1000,
    height = 800)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/FeaturePlot_IGF2BP1_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

png(file = "Output/06_Label_MN4/FeaturePlot_IGF2BP1_VB1.png",
    width = 1000,
    height = 800)
print(plot_VB1)
dev.off()

pdf(file = "Output/06_Label_MN4/FeaturePlot_IGF2BP1_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_184_0.png)
    



    
![png](output_184_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_184_8.png)
    


# Plot ENG, VSIR, and FASL (marco request 20250509)


```R
# spatial feature plot colored by navin annotations
options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialFeaturePlot(data, features=c('ENG'), pt.size.factor=2)
plot_VA2 <- plot2[[1]] + 
            labs(title="ENG VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="ENG VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="ENG VB1")
plot_VA2
plot_VA1
plot_VB1

png(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_ENG_VA2.png",
    width = 1000,
    height = 800)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_ENG_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

png(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_ENG_VA1.png",
    width = 1000,
    height = 800)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_ENG_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

png(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_ENG_VB1.png",
    width = 1000,
    height = 800)
print(plot_VB1)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_ENG_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_186_0.png)
    



    
![png](output_186_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_186_8.png)
    



```R
# spatial feature plot colored by navin annotations
options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialFeaturePlot(data, features=c('VSIR'), pt.size.factor=2)
plot_VA2 <- plot2[[1]] + 
            labs(title="VSIR VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="VSIR VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="VSIR VB1")
plot_VA2
plot_VA1
plot_VB1

png(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_VSIR_VA2.png",
    width = 1000,
    height = 800)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_VSIR_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

png(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_VSIR_VA1.png",
    width = 1000,
    height = 800)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_VSIR_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

png(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_VSIR_VB1.png",
    width = 1000,
    height = 800)
print(plot_VB1)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_VSIR_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_187_0.png)
    



    
![png](output_187_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_187_8.png)
    



```R
# spatial feature plot colored by navin annotations
options(repr.plot.width=11, repr.plot.height=9)

plot2 <- SpatialFeaturePlot(data, features=c('FASLG'), pt.size.factor=2)
plot_VA2 <- plot2[[1]] + 
            labs(title="FASLG VA2")
plot_VA1 <- plot2[[2]] + 
            labs(title="FASLG VA1")
plot_VB1 <- plot2[[3]] + 
            labs(title="FASLG VB1")
plot_VA2
plot_VA1
plot_VB1

png(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_FASLG_VA2.png",
    width = 1000,
    height = 800)
print(plot_VA2)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_FASLG_VA2.pdf",
    width = 10,
    height = 8)
print(plot_VA2)
dev.off()

png(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_FASLG_VA1.png",
    width = 1000,
    height = 800)
print(plot_VA1)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_FASLG_VA1.pdf",
    width = 10,
    height = 8)
print(plot_VA1)
dev.off()

png(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_FASLG_VB1.png",
    width = 1000,
    height = 800)
print(plot_VB1)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/FeaturePlot_FASLG_VB1.pdf",
    width = 10,
    height = 8)
print(plot_VB1)
dev.off()
```


    
![png](output_188_0.png)
    



    
![png](output_188_1.png)
    



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_188_8.png)
    



```R

```

# Plot ENG, VSIR, and FASLG for vasc_high spots, comparing MN4 high vs MN4 low


```R
unique(data@meta.data$vasc_MN4_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_high_MN4_high'</li><li>'vasc_low'</li><li>'vasc_high_MN4_low'</li></ol>




```R
library(cowplot)
```

    
    Attaching package: ‘cowplot’
    
    
    The following object is masked from ‘package:lubridate’:
    
        stamp
    
    



```R
# set plot size
options(repr.plot.width=18, repr.plot.height=12)
# set idents to vasc_MN4_label
Idents(data) = 'vasc_MN4_label'
# copy seurat object
data2 <- data
# remove troublesome graph 
data2@graphs <- list() 
# subet data to exclude vasc_low
data_sub <- subset(data2, subset=vasc_MN4_label != 'vasc_low')
# plot
#VlnPlot(data_sub, features=c('ENG', 'VSIR', 'FASLG'), split.by='Tissue_Slice')
p1 <- VlnPlot(data_sub, features=c('ENG'), split.by='Tissue_Slice')
p2 <- VlnPlot(data_sub, features=c('VSIR'), split.by='Tissue_Slice')
p3 <- VlnPlot(data_sub, features=c('FASLG'), split.by='Tissue_Slice')

marco_plot <- plot_grid(p1, p2, p3)
marco_plot

png(file = "Output/06_Label_MN4/Marco_plots_20250509/ENG_VSIR_FASLG_VlnPlot.png",
    width = 1800,
    height = 1200)
print(marco_plot)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/ENG_VSIR_FASLG_VlnPlot.pdf",
    width = 18,
    height = 12)
print(marco_plot)
dev.off()
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
    The default behaviour of split.by has changed.
    Separate violin plots are now plotted side-by-side.
    To restore the old behaviour of a single split violin,
    set split.plot = TRUE.
          
    This message will be shown once per session.
    
    Warning message:
    “[1m[22mThe `slot` argument of `FetchData()` is deprecated as of SeuratObject 5.0.0.
    [36mℹ[39m Please use the `layer` argument instead.
    [36mℹ[39m The deprecated feature was likely used in the [34mSeurat[39m package.
      Please report the issue at [3m[34m<https://github.com/satijalab/seurat/issues>[39m[23m.”
    Warning message:
    “[1m[22m`PackageCheck()` was deprecated in SeuratObject 5.0.0.
    [36mℹ[39m Please use `rlang::check_installed()` instead.
    [36mℹ[39m The deprecated feature was likely used in the [34mSeurat[39m package.
      Please report the issue at [3m[34m<https://github.com/satijalab/seurat/issues>[39m[23m.”



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_193_3.png)
    



```R
install.packages("scCustomize")
```

    Installing package into ‘/usr/local/lib/R/site-library’
    (as ‘lib’ is unspecified)
    
    Warning message in install.packages("scCustomize"):
    “installation of package ‘scCustomize’ had non-zero exit status”



```R
library(scCustomize)
```

    [1m[22mscCustomize v3.0.1
    If you find the scCustomize useful please cite.
    See 'samuel-marsh.github.io/scCustomize/articles/FAQ.html' for citation info.
    



```R
# set plot size
options(repr.plot.width=24, repr.plot.height=10)
# plot with scCustomize
VlnPlot_scCustom(seurat_object = data_sub, 
             features = c('ENG', 'VSIR', 'FASLG'), 
             split.by='Tissue_Slice')
```


    
![png](output_196_0.png)
    


# do again but only in spots that aren't Normal Lung or Fibroblasts


```R
# set plot size
options(repr.plot.width=18, repr.plot.height=12)
# set idents to vasc_MN4_label
Idents(data) = 'vasc_MN4_label'
# copy seurat object
data2 <- data
# remove troublesome graph 
data2@graphs <- list() 
# subet data to exclude vasc_low
data_sub <- subset(data2, subset=vasc_MN4_label != 'vasc_low') 
# subet data to exclude Normal Lung and Fibroblast spots
data_sub <- subset(data_sub, subset=Navin_Clusters2 != '7_10_Normal_Lung')
data_sub <- subset(data_sub, subset=Navin_Clusters2 != '8_Fibroblasts')
# plot
#VlnPlot(data_sub, features=c('ENG', 'VSIR', 'FASLG'), split.by='Tissue_Slice')
p1 <- VlnPlot(data_sub, features=c('ENG'), split.by='Tissue_Slice')
p2 <- VlnPlot(data_sub, features=c('VSIR'), split.by='Tissue_Slice')
p3 <- VlnPlot(data_sub, features=c('FASLG'), split.by='Tissue_Slice')

marco_plot <- plot_grid(p1, p2, p3)
marco_plot

png(file = "Output/06_Label_MN4/Marco_plots_20250509/ENG_VSIR_FASLG_VlnPlot_NotNormalNotFibroblast.png",
    width = 1800,
    height = 1200)
print(marco_plot)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/ENG_VSIR_FASLG_VlnPlot_NotNormalNotFibroblast.pdf",
    width = 18,
    height = 12)
print(marco_plot)
dev.off()
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



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_198_3.png)
    


# Marco request: Look for ENG, VSIR, and FASLG in CD3 and/or CD8 high cells. 


```R
unique(data@meta.data$Tumor_MHCI_label2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_high_MHCI_high'</li></ol>




```R
library(GGally)
```


```R
# set plot size
options(repr.plot.width=18, repr.plot.height=12)
# set idents to vasc_MN4_label
Idents(data) = 'vasc_MN4_label'
# copy seurat object
data2 <- data
# remove troublesome graph 
data2@graphs <- list() 
# subet data to Tumor MHCI High and Low
data_sub_Tumor_MHCI_High <- subset(data2, subset=Tumor_MHCI_label2 == 'Tumor_high_MHCI_high')
data_sub_Tumor_MHCI_Low <- subset(data2, subset=Tumor_MHCI_label2 == 'Tumor_high_MHCI_low')
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
options(repr.plot.width=8, repr.plot.height=8)
FeatureScatter(object = data_sub_Tumor_MHCI_High, feature1 = 'ENG', feature2 = 'CD3E')
FeatureScatter(object = data_sub_Tumor_MHCI_High, feature1 = 'VSIR', feature2 = 'CD3E')
FeatureScatter(object = data_sub_Tumor_MHCI_High, feature1 = 'FASLG', feature2 = 'CD3E')
FeatureScatter(object = data_sub_Tumor_MHCI_High, feature1 = 'ENG', feature2 = 'CD8A')
FeatureScatter(object = data_sub_Tumor_MHCI_High, feature1 = 'VSIR', feature2 = 'CD8A')
FeatureScatter(object = data_sub_Tumor_MHCI_High, feature1 = 'FASLG', feature2 = 'CD8A')
```


    
![png](output_203_0.png)
    



    
![png](output_203_1.png)
    



    
![png](output_203_2.png)
    



    
![png](output_203_3.png)
    



    
![png](output_203_4.png)
    



    
![png](output_203_5.png)
    


# Extract the genes and plot them against the enrichments


```R
data2 <- data
expr_temp <- FetchData(object = data2, vars = c("ENG", "VSIR", "FASLG"))
expr_temp$Barcode <- rownames(expr_temp)
metadata_temp <- data2@meta.data
#colnames(metadata_temp)
temp_merged <- merge(metadata_temp, expr_temp, by='Barcode')
rownames(temp_merged) <- temp_merged$Barcode
data2@meta.data <- temp_merged
data2@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 64</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>ENG</th><th scope=col>VSIR</th><th scope=col>FASLG</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>1.324052</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>0.000000</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>1.101286</td><td>0</td><td>0</td></tr>
</tbody>
</table>




```R
library(GGally)
```


```R
# use GGally package to do a pair plot of marker values
options(repr.plot.width=24, repr.plot.height=24)
plot_list <- c('ENG','VSIR','FASLG','T_Cell','CD8_T_Cell','NK_Cell','All_Immune')
plot1 <- ggpairs(temp_merged %>% select(all_of(plot_list)),
                upper = list(continuous = wrap("cor", size = 9))) 

plot1 <- plot1 + 
            theme_bw() +
            theme(strip.text.x = element_text(size = 20),
                  strip.text.y = element_text(size = 20))
print(plot1)

```


    
![png](output_207_0.png)
    

# Subset to above median T Cell
above_med_TCell <- data2@meta.data[data2@meta.data$T_Cell > median(data2@meta.data$T_Cell),]
above_med_CD8TCell <- data2@meta.data[data2@meta.data$CD8_T_Cell > median(data2@meta.data$CD8_T_Cell),]
above_80p_TCell <- data2@meta.data[data2@meta.data$T_Cell > quantile(data2@meta.data$T_Cell, c(.80))[[1]],]
above_80p_CD8TCell <- data2@meta.data[data2@meta.data$CD8_T_Cell > quantile(data2@meta.data$CD8_T_Cell, c(.80))[[1]],]

```R
# Subset seurat objects to above median T Cell
# remove troublesome graph 
data2@graphs <- list() 
# remove Tumor low
#data_sub <- subset(data2, subset=Tumor_MHCI_label2 != 'Tumor_low') 
data_sub <- data2
Idents(data_sub) <- 'Tumor_MHCI_label2'
above_med_TCell <- subset(data_sub, subset = T_Cell > median(data_sub@meta.data$T_Cell))
above_med_CD8TCell <- subset(data_sub, subset = CD8_T_Cell > median(data_sub@meta.data$CD8_T_Cell))
above_80p_TCell <- subset(data_sub, subset = T_Cell > quantile(data_sub@meta.data$T_Cell, c(.80))[[1]])
above_80p_CD8TCell <- subset(data_sub, subset = CD8_T_Cell > quantile(data_sub@meta.data$CD8_T_Cell, c(.80))[[1]])
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

```


```R
dim(data2@meta.data)
dim(above_med_TCell)
dim(above_med_CD8TCell)
dim(above_80p_TCell)
dim(above_80p_CD8TCell)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>64</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>17943</li><li>4855</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>17943</li><li>4855</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>17943</li><li>1942</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>17943</li><li>1942</li></ol>




```R
# plot
options(repr.plot.width=18, repr.plot.height=12)
#VlnPlot(data_sub, features=c('ENG', 'VSIR', 'FASLG'), split.by='Tissue_Slice')
p1 <- VlnPlot(above_med_TCell, features=c('ENG'), split.by='Tumor_MHCI_label2')
p2 <- VlnPlot(above_med_TCell, features=c('VSIR'), split.by='Tumor_MHCI_label2')
p3 <- VlnPlot(above_med_TCell, features=c('FASLG'), split.by='Tumor_MHCI_label2')

marco_plot <- plot_grid(p1, p2, p3) + ggtitle('Above median TCell Enrichment')
marco_plot

png(file = "Output/06_Label_MN4/Marco_plots_20250509/ENG_VSIR_FASLG_VlnPlot_TCell_high.png",
    width = 1800,
    height = 1200)
print(marco_plot)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/ENG_VSIR_FASLG_VlnPlot_TCell_high.pdf",
    width = 18,
    height = 12)
print(marco_plot)
dev.off()
```

    Warning message:
    “The following variables were found in both object meta data and the default assay: ENG
    Returning meta data; if you want the feature, please use the assay's key (eg. spatial_ENG)”
    Warning message:
    “The following variables were found in both object meta data and the default assay: VSIR
    Returning meta data; if you want the feature, please use the assay's key (eg. spatial_VSIR)”
    Warning message:
    “The following variables were found in both object meta data and the default assay: FASLG
    Returning meta data; if you want the feature, please use the assay's key (eg. spatial_FASLG)”



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_212_3.png)
    



```R
# plot
options(repr.plot.width=18, repr.plot.height=12)
#VlnPlot(data_sub, features=c('ENG', 'VSIR', 'FASLG'), split.by='Tissue_Slice')
p1 <- VlnPlot(above_med_CD8TCell, features=c('ENG'), split.by='Tumor_MHCI_label2')
p2 <- VlnPlot(above_med_CD8TCell, features=c('VSIR'), split.by='Tumor_MHCI_label2')
p3 <- VlnPlot(above_med_CD8TCell, features=c('FASLG'), split.by='Tumor_MHCI_label2')

marco_plot <- plot_grid(p1, p2, p3) + ggtitle('Above median CD8TCell Enrichment')
marco_plot

png(file = "Output/06_Label_MN4/Marco_plots_20250509/ENG_VSIR_FASLG_VlnPlot_CD8TCell_high.png",
    width = 1800,
    height = 1200)
print(marco_plot)
dev.off()

pdf(file = "Output/06_Label_MN4/Marco_plots_20250509/ENG_VSIR_FASLG_VlnPlot_CD8TCell_high.pdf",
    width = 18,
    height = 12)
print(marco_plot)
dev.off()
```

    Warning message:
    “The following variables were found in both object meta data and the default assay: ENG
    Returning meta data; if you want the feature, please use the assay's key (eg. spatial_ENG)”
    Warning message:
    “The following variables were found in both object meta data and the default assay: VSIR
    Returning meta data; if you want the feature, please use the assay's key (eg. spatial_VSIR)”
    Warning message:
    “The following variables were found in both object meta data and the default assay: FASLG
    Returning meta data; if you want the feature, please use the assay's key (eg. spatial_FASLG)”



<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_213_3.png)
    



```R

```


```R

```


```R
options(repr.plot.width=20, repr.plot.height=16)

# subset to only VA1 and VB1
temp_df = annots_melt_immune %>% subset(Tissue_Slice %in% c('VA2'))

# make sure to subset for vasc_label = vasc_high
celltype_vs_MN4_EC_vaschigh_plot1 <- ggplot(temp_df[temp_df$vasc_label == 'vasc_high',], 
                                          aes(x=MN4_EC_Phenotype_Label, y=Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey', scale = "width", width = 0.9) + 
  geom_beeswarm(cex = 0.9, size=1, aes(color=Tissue_Slice)) +
  #geom_point() + 
  facet_wrap(vars(Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Cell Type GSVA Enrichment vs. MN4_EC_Label, in vasc_high spots")
             
celltype_vs_MN4_EC_vaschigh_plot1

# save plot
png(file = "Output/06_Label_MN4/cell_type_GSVA_per_MN4_EC_Label_for_vasc_high_VA2.png",
    width = 1800,
    height = 1400)
print(celltype_vs_MN4_EC_vaschigh_plot1)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_216_1.png)
    



```R

```


```R

```

# test colors


```R
dittoColors(get.names=FALSE)[1:15]
dittoColors(get.names=TRUE)[1:15]
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'#E69F00'</li><li>'#56B4E9'</li><li>'#009E73'</li><li>'#F0E442'</li><li>'#0072B2'</li><li>'#D55E00'</li><li>'#CC79A7'</li><li>'#666666'</li><li>'#AD7700'</li><li>'#1C91D4'</li><li>'#007756'</li><li>'#D5C711'</li><li>'#005685'</li><li>'#A04700'</li><li>'#B14380'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'orange'</li><li>'skyBlue'</li><li>'bluishGreen'</li><li>'yellow'</li><li>'blue'</li><li>'vermillion'</li><li>'reddishPurple'</li><li>'gray40'</li><li>'orange-25%'</li><li>'skyBlue-25%'</li><li>'bluishGreen-25%'</li><li>'yellow-25%'</li><li>'blue-25%'</li><li>'vermillion-25%'</li><li>'reddishPurple-25%'</li></ol>




```R
show_col(dittoColors()[1:32])
```


    
![png](output_221_0.png)
    



```R
#paletteer::paletteer_d("beyonce::X58")
clrs <- paletteer::paletteer_d("RColorBrewer::OrRd")
show_col(clrs)
```


    
![png](output_222_0.png)
    



```R
clrs
```


    <colors>
    [30m[48;5;231m#FFF7ECFF[49m[39m [30m[48;5;230m#FEE8C8FF[49m[39m [30m[48;5;223m#FDD49EFF[49m[39m [30m[48;5;223m#FDBB84FF[49m[39m [30m[48;5;216m#FC8D59FF[49m[39m [30m[48;5;209m#EF6548FF[49m[39m [37m[48;5;202m#D7301FFF[49m[39m [37m[48;5;160m#B30000FF[49m[39m [37m[48;5;88m#7F0000FF[49m[39m 



```R

```


```R

```


```R

```


```R

```


```R

```


```R

```


```R

```
