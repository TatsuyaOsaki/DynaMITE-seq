# Marco and Navin requested GSEA with hallmark pathways, comparing Tumor MHCI High vs Tumor MHCI Low spots. 


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
library(EnhancedVolcano)
library(ComplexHeatmap)
library(circlize)
library(fgsea)

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
    
    Loading required package: ggrepel
    
    Loading required package: grid
    
    ========================================
    ComplexHeatmap version 2.24.0
    Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
    Github page: https://github.com/jokergoo/ComplexHeatmap
    Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
    
    If you use it in published research, please cite either one:
    - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
    - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
        genomic data. Bioinformatics 2016.
    
    
    The new InteractiveComplexHeatmap package can directly export static 
    complex heatmaps into an interactive Shiny app with zero effort. Have a try!
    
    This message can be suppressed by:
      suppressPackageStartupMessages(library(ComplexHeatmap))
    ========================================
    ! pheatmap() has been masked by ComplexHeatmap::pheatmap(). Most of the arguments
       in the original pheatmap() are identically supported in the new function. You 
       can still use the original function by explicitly calling pheatmap::pheatmap().
    
    
    
    Attaching package: ‘ComplexHeatmap’
    
    
    The following object is masked from ‘package:pheatmap’:
    
        pheatmap
    
    
    ========================================
    circlize version 0.4.16
    CRAN page: https://cran.r-project.org/package=circlize
    Github page: https://github.com/jokergoo/circlize
    Documentation: https://jokergoo.github.io/circlize_book/book/
    
    If you use it in published research, please cite:
    Gu, Z. circlize implements and enhances circular visualization
      in R. Bioinformatics 2014.
    
    This message can be suppressed by:
      suppressPackageStartupMessages(library(circlize))
    ========================================
    
    


# Import Seurat Objects


```R
data <- readRDS('Processing/AB1_A2_Annotated_GSVA_by_Tissue_Slice_06_Label_MN4_20250425.rds')
```


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
# check metadata
data@meta.data[1:3,]
```


<table class="dataframe">
<caption>A data.frame: 3 × 57</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td> 0.6906627</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.2739867</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td></tr>
</tbody>
</table>




```R
unique(data@meta.data$Tumor_MHCI_label2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_high_MHCI_high'</li></ol>



# Use FindMarkers to find DE genes between Tumor MHCI High and Tumor MHCI Low


```R
# set ident to Tumor_MHCI_Label
Idents(data) <- 'Tumor_MHCI_label2'
```


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
# find all markers distinguishing Tumor MHCI High and Tumor MHCI Low
Tumor_MHCI_High_Markers <- FindMarkers(data, ident.1 = 'Tumor_high_MHCI_high', ident.2 = 'Tumor_high_MHCI_low', 
                                       logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1, 
                                       min.pct = 0, fc.slot = "counts")

# make new columns from seurat FindMarkers output
Tumor_MHCI_High_Markers$log2FoldChange <- Tumor_MHCI_High_Markers$avg_log2FC
Tumor_MHCI_High_Markers$padj <- Tumor_MHCI_High_Markers$p_val_adj
Tumor_MHCI_High_Markers$gene_name <- rownames(Tumor_MHCI_High_Markers)

head(Tumor_MHCI_High_Markers, n = 5)
```

    Warning message:
    “[1m[22mThe `slot` argument of `GetAssayData()` is deprecated as of SeuratObject 5.0.0.
    [36mℹ[39m Please use the `layer` argument instead.
    [36mℹ[39m The deprecated feature was likely used in the [34mSeurat[39m package.
      Please report the issue at [3m[34m<https://github.com/satijalab/seurat/issues>[39m[23m.”
    Warning message:
    “[1m[22m`PackageCheck()` was deprecated in SeuratObject 5.0.0.
    [36mℹ[39m Please use `rlang::check_installed()` instead.
    [36mℹ[39m The deprecated feature was likely used in the [34mSeurat[39m package.
      Please report the issue at [3m[34m<https://github.com/satijalab/seurat/issues>[39m[23m.”
    For a (much!) faster implementation of the Wilcoxon Rank Sum Test,
    (default method for FindMarkers) please install the presto package
    --------------------------------------------
    install.packages('devtools')
    devtools::install_github('immunogenomics/presto')
    --------------------------------------------
    After installation of presto, Seurat will automatically use the more 
    efficient implementation (no further action necessary).
    This message will be shown once per session
    



<table class="dataframe">
<caption>A data.frame: 5 × 8</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FoldChange</th><th scope=col>padj</th><th scope=col>gene_name</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B2M</th><td> 0.000000e+00</td><td>0.1686889</td><td>1.000</td><td>1.000</td><td> 0.000000e+00</td><td>0.1686889</td><td> 0.000000e+00</td><td>B2M    </td></tr>
	<tr><th scope=row>TAP1</th><td>3.747407e-285</td><td>0.6774325</td><td>0.907</td><td>0.622</td><td>6.723973e-281</td><td>0.6774325</td><td>6.723973e-281</td><td>TAP1   </td></tr>
	<tr><th scope=row>TAP2</th><td>2.070688e-253</td><td>0.7210781</td><td>0.779</td><td>0.485</td><td>3.715436e-249</td><td>0.7210781</td><td>3.715436e-249</td><td>TAP2   </td></tr>
	<tr><th scope=row>CD74</th><td>1.321012e-155</td><td>0.4820377</td><td>0.995</td><td>0.981</td><td>2.370293e-151</td><td>0.4820377</td><td>2.370293e-151</td><td>CD74   </td></tr>
	<tr><th scope=row>HLA-DRA</th><td>4.071785e-136</td><td>0.5069825</td><td>0.970</td><td>0.948</td><td>7.306004e-132</td><td>0.5069825</td><td>7.306004e-132</td><td>HLA-DRA</td></tr>
</tbody>
</table>



# Repeat DE for use in GSEA rankings


```R
# find all markers distinguishing Tumor MHCI High and Tumor MHCI Low
Tumor_MHCI_High_Markers_for_GSEA <- FindMarkers(data, ident.1 = 'Tumor_high_MHCI_high', ident.2 = 'Tumor_high_MHCI_low', 
                                       logfc.threshold = 0, min.cells.feature = 1, min.cells.group = 1, 
                                       min.pct = 0)

# make new columns from seurat FindMarkers output
Tumor_MHCI_High_Markers_for_GSEA$log2FoldChange <- Tumor_MHCI_High_Markers_for_GSEA$avg_log2FC
Tumor_MHCI_High_Markers_for_GSEA$padj <- Tumor_MHCI_High_Markers_for_GSEA$p_val_adj
Tumor_MHCI_High_Markers_for_GSEA$gene_name <- rownames(Tumor_MHCI_High_Markers_for_GSEA)

head(Tumor_MHCI_High_Markers_for_GSEA, n = 5)
```


<table class="dataframe">
<caption>A data.frame: 5 × 8</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FoldChange</th><th scope=col>padj</th><th scope=col>gene_name</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B2M</th><td> 0.000000e+00</td><td>0.6443753</td><td>1.000</td><td>1.000</td><td> 0.000000e+00</td><td>0.6443753</td><td> 0.000000e+00</td><td>B2M    </td></tr>
	<tr><th scope=row>TAP1</th><td>3.747407e-285</td><td>1.5546946</td><td>0.907</td><td>0.622</td><td>6.723973e-281</td><td>1.5546946</td><td>6.723973e-281</td><td>TAP1   </td></tr>
	<tr><th scope=row>TAP2</th><td>2.070688e-253</td><td>1.8187478</td><td>0.779</td><td>0.485</td><td>3.715436e-249</td><td>1.8187478</td><td>3.715436e-249</td><td>TAP2   </td></tr>
	<tr><th scope=row>CD74</th><td>1.321012e-155</td><td>0.7369142</td><td>0.995</td><td>0.981</td><td>2.370293e-151</td><td>0.7369142</td><td>2.370293e-151</td><td>CD74   </td></tr>
	<tr><th scope=row>HLA-DRA</th><td>4.071785e-136</td><td>0.7700200</td><td>0.970</td><td>0.948</td><td>7.306004e-132</td><td>0.7700200</td><td>7.306004e-132</td><td>HLA-DRA</td></tr>
</tbody>
</table>




```R
dim(data[['Spatial']]$counts)
dim(Tumor_MHCI_High_Markers_for_GSEA)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>17943</li><li>9711</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>17943</li><li>8</li></ol>



# Perform GSEA on the Hallmark Pathways (and others)


```R
# Function to perform Gene Set Enrichment Analysis (GSEA) on the current comparison
do_GSEA <- function(df, save_filename, group_name_1, group_name_2, save_dir) {
    
    print("do GSEA")

    # make new columns from seurat FindMarkers output
    df$log2FoldChange <- df$avg_log2FC
    df$padj <- df$p_val_adj
    df$Gene <- rownames(df)
    
    # dropna from log2FoldChange and padj
    df <- df %>% drop_na(log2FoldChange)
    df <- df %>% drop_na(p_val)
    # add a very small number to padj
    df$p_val <- df$p_val + 0.0000000000000000001
    
    # Create a significance score to rank the gene lists for GSEA
    # https://crazyhottommy.blogspot.com/2016/08/gene-set-enrichment-analysis-gsea.html?showComment=1624983978207#c8639673391160120332
    # use Log2 Fold Change * (-1) * log10(adjusted p value)
    # can first undo the log2 transformation by applying 2^x
    df$sigScore <- df$log2FoldChange * (-1) * log10(df$p_val)
    
    # Sort by sigScore and check...
    df <- df %>% arrange(desc(sigScore))
    
    # for normal GSEA, supply the ENTIRE sequencing results... all of the genes in a ranked list. 
    # Complete ranked gene list for GSEA ------------------------------------------
    # Subset gene name and significance score for ranked lists
    # set gene name column as the row names using deframe()
    #ranked <- df[c('gene_name','sigScore')] %>% deframe()
    ranked <- df[c('Gene','sigScore')] %>% deframe()
    
    hallmark_gmt_path <- "/immuno/jason/Projects/2024_01_Giorgia-RNA-seq/Reference/Pathway/h.all.v2023.2.Hs.symbols.gmt"
    # Load the pathways into a named list
    hallmarkGeneSets <- gmtPathways(hallmark_gmt_path)

    # Load the neuroendocrine pathway
    NE_sig_df <- read.csv('Processing/Neuroendocrine_signature.csv')
    NE_sig <- NE_sig_df$Gene

    # add the neuroendocrine pathway to the hallmark pathways
    hallmarkGeneSets$Neuroendocrine <- NE_sig
    
    # Run FGSEA only hallmark for now
    # sort by NES descending and padj in ascending order
    # ------ HALLMARK GENE SETS -------
    H_fgseaRes <- fgsea(pathways = hallmarkGeneSets, ranked, minSize=10, maxSize=2000) %>% arrange(desc(NES),padj)
    
    print("fgsea complete")
    
    # Save excel sheets with results
    #save_filename <- paste(unlist(contrast_item[2]), "vs", unlist(contrast_item[3]), sep="_")
    print(save_filename)
    GSEA_results <- list("HALLMARK" = H_fgseaRes)
    require(openxlsx)
    #gsea_filename <- paste(save_dir, paste0("GSEA_",paste0(save_filename,".xlsx"), sep="/"))
    gsea_filename <- file.path(save_dir, paste0("GSEA_",paste0(save_filename,".xlsx")))
    write.xlsx(GSEA_results, file = gsea_filename)
    
    # Plot NES scores for each
    # function to plot GSEA enrichment scores
    plot_gsea <- function(df, save_dir, save_filename2) {
        # get only top 25 and bottom 25 NES scores
        df2 = rbind(df %>% slice_max(NES, n = 25),df %>% slice_min(NES, n = 26))

        png(paste(save_dir, paste0("GSEA_",save_filename2,".png"), sep="/"), width=1100, height=700)
        plt <- ggplot(df2, aes(reorder(pathway, NES), NES)) +
            geom_col(aes(fill=padj<0.05)) +
            coord_flip() +
            labs(x="Pathway",
                 title=save_filename2) + 
            theme_minimal()
        print(plt)
        dev.off()
    }
    
    # Plot NES scores for each
    plot_gsea(H_fgseaRes, save_dir, paste(save_filename,"_HALLMARK"))
    }
```


```R
# Run GSEA
do_GSEA(Tumor_MHCI_High_Markers_for_GSEA, 'Tumor_high_MHCI_high_vs_low','Tumor_MHCI_High', 'Tumor_MHCI_Low', '/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/15_DE_and_GSEA/GSEA_20250703_TumorHigh_MHCIHigh_vs_MHCILow')
```

    [1] "do GSEA"
    [1] "fgsea complete"
    [1] "Tumor_high_MHCI_high_vs_low"



<strong>pdf:</strong> 2


# Make some plots


```R
# A function to make the volcano plot
plot_volcano <- function(df, group_name_1, group_name_2, save_dir) {
    # Make and save the volcano plot
    save_filename = paste(group_name_1, "vs", group_name_2, sep="_")
    png(paste0(paste(save_dir, save_filename, save_filename, sep="/"),".png"), width=600, height=600, res=100)
    print(ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val), col=diff_exp, label=gene_label)) +
      geom_point() + 
      theme_minimal() +
      geom_text_repel(max.overlaps = 30) +
      scale_color_manual(values=c("blue", "black", "red")) +
      geom_vline(xintercept=c(-0.8, 0.8), col="grey") +
      geom_hline(yintercept=-log10(0.05), col="grey") +
      labs(title=save_filename)
      #guides(fill = guide_legend(override.aes = aes(label = ""))) +
      #xlim(min(df$log2FoldChange), max(df$avg_log2FC)) + 
      #ylim(0, max(-log10(df$p_val)))
         )
    dev.off()
    print('plot_volcano complete')
}
```


```R
df <- Tumor_MHCI_High_Markers

# add a very small number to padj
df$padj <- df$padj + 0.0000000000000000001

# prepare data for custom volcano plots
# add a column of NAs
df$diff_exp <- "NO"
# if log2Foldchange > 0.8 and pvalue < 0.05, set as "UP" 
df$diff_exp[df$log2FoldChange > 0.8 & df$padj < 0.05] <- "UP"
# log2(1.74) = 0.8, log2(1/1.74) = -0.8
# if log2Foldchange < -0.8 and pvalue < 0.05, set as "DOWN"
df$diff_exp[df$log2FoldChange < -0.8 & df$padj < 0.05] <- "DOWN"
# Annotate gene name labels for the differentially expressed genes
df$gene_label <- NA
df$gene_label[df$diff_exp != "NO"] <- df$gene_name[df$diff_exp != "NO"]
```


```R
plot_volcano(df, 'Tumor_MHCI_High', 'Tumor_MHCI_Low', '/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/15_DE_and_GSEA')

```

    Warning message:
    “[1m[22mRemoved 16458 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).”
    Warning message:
    “ggrepel: 1477 unlabeled data points (too many overlaps). Consider increasing max.overlaps”


    [1] "plot_volcano complete"


# save DE genes


```R
save_filename = 'Tumor_MHCI_High_vs_Tumor_MHCI_Low'
save_dir = '/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/15_DE_and_GSEA'

# save csv files for genes shown in the main heatmap
write.csv(as.matrix(df), file=paste0(paste(save_dir, save_filename, save_filename, sep="/"),".csv"))
# now only include significant genes
write.csv(as.matrix(subset(df, diff_exp != "NO")), file=paste0(paste(save_dir, save_filename, save_filename, sep="/"),"_SigDE",".csv"))
```

# Save DE genes used for GSEA rankings


```R
df2 <- Tumor_MHCI_High_Markers_for_GSEA

# add a very small number to padj
df2$padj <- df2$padj + 0.0000000000000000001

# prepare data for custom volcano plots
# add a column of NAs
df2$diff_exp <- "NO"
# if log2Foldchange > 0.8 and pvalue < 0.05, set as "UP" 
df2$diff_exp[df2$log2FoldChange > 0.8 & df2$padj < 0.05] <- "UP"
# log2(1.74) = 0.8, log2(1/1.74) = -0.8
# if log2Foldchange < -0.8 and pvalue < 0.05, set as "DOWN"
df2$diff_exp[df2$log2FoldChange < -0.8 & df2$padj < 0.05] <- "DOWN"
# Annotate gene name labels for the differentially expressed genes
df2$gene_label <- NA
df2$gene_label[df2$diff_exp != "NO"] <- df2$gene_name[df2$diff_exp != "NO"]
```


```R
plot_volcano(df2, 'Tumor_MHCI_High', 'Tumor_MHCI_Low_for_GSEA', '/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/15_DE_and_GSEA')

```

    Warning message:
    “[1m[22mRemoved 17695 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).”
    Warning message:
    “ggrepel: 232 unlabeled data points (too many overlaps). Consider increasing max.overlaps”


    [1] "plot_volcano complete"



```R
save_filename = 'Tumor_MHCI_High_vs_Tumor_MHCI_Low_for_GSEA'
save_dir = '/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/15_DE_and_GSEA'

# save csv files for genes shown in the main heatmap
write.csv(as.matrix(df2), file=paste0(paste(save_dir, save_filename, save_filename, sep="/"),".csv"))
# now only include significant genes
write.csv(as.matrix(subset(df2, diff_exp != "NO")), file=paste0(paste(save_dir, save_filename, save_filename, sep="/"),"_SigDE",".csv"))
```

# Remaking Figure 3E, the volcano plot, for the paper

# Use the DE genes used for GSEA


```R
# load de genes again
save_filename = 'Tumor_MHCI_High_vs_Tumor_MHCI_Low_for_GSEA'
save_dir = '/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/15_DE_and_GSEA'
load_file = paste0(paste(save_dir, save_filename, save_filename, sep="/"),".csv")

df <- read.csv(load_file)
df <- df %>% column_to_rownames('X')
df %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 10</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FoldChange</th><th scope=col>padj</th><th scope=col>gene_name</th><th scope=col>diff_exp</th><th scope=col>gene_label</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B2M</th><td> 0.000000e+00</td><td>0.6443753</td><td>1.000</td><td>1.000</td><td> 0.000000e+00</td><td>0.6443753</td><td>1e-19</td><td>B2M </td><td>NO</td><td>NA  </td></tr>
	<tr><th scope=row>TAP1</th><td>3.747407e-285</td><td>1.5546950</td><td>0.907</td><td>0.622</td><td>6.723973e-281</td><td>1.5546950</td><td>1e-19</td><td>TAP1</td><td>UP</td><td>TAP1</td></tr>
	<tr><th scope=row>TAP2</th><td>2.070688e-253</td><td>1.8187480</td><td>0.779</td><td>0.485</td><td>3.715436e-249</td><td>1.8187480</td><td>1e-19</td><td>TAP2</td><td>UP</td><td>TAP2</td></tr>
</tbody>
</table>




```R
# check genes of interest
old_goi_list <- c('HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'B2M', 'CXCL9', 'CXCL10', 'AXL', 'PRRX1', 'VIM', 'SPARC', 'ZEB1', 'SNAI1', 'SNAI2', 'SNAI3')
#goi_list <- c('TAP1', 'TAP2', 'B2M', 'CXCL9', 'VIM', 'SPARC')
goi_list <- c('TAP1', 'TAP2', 'B2M', 'CXCL9', 'VIM', 'SPARC', 'HLA-A')

gene_name_list <- unique(df$gene_name)

both_list <- intersect(goi_list, gene_name_list)
only_goi_list <- setdiff(goi_list, gene_name_list)

print(both_list)
print(only_goi_list)
```

    [1] "TAP1"  "TAP2"  "B2M"   "CXCL9" "VIM"   "SPARC" "HLA-A"
    character(0)



```R
# label genes of interest
# Annotate gene name labels for the differentially expressed genes
df$gene_label2 <- NA
df$gene_label2[df$gene_name %in% goi_list] <- df$gene_name[df$gene_name %in% goi_list]
```


```R
df %>% head(30)
print(unique(df$gene_label2))
```


<table class="dataframe">
<caption>A data.frame: 30 × 11</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FoldChange</th><th scope=col>padj</th><th scope=col>gene_name</th><th scope=col>diff_exp</th><th scope=col>gene_label</th><th scope=col>gene_label2</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B2M</th><td> 0.000000e+00</td><td> 0.6443753</td><td>1.000</td><td>1.000</td><td> 0.000000e+00</td><td> 0.6443753</td><td>1e-19</td><td>B2M     </td><td>NO</td><td>NA    </td><td>B2M </td></tr>
	<tr><th scope=row>TAP1</th><td>3.747407e-285</td><td> 1.5546950</td><td>0.907</td><td>0.622</td><td>6.723973e-281</td><td> 1.5546950</td><td>1e-19</td><td>TAP1    </td><td>UP</td><td>TAP1  </td><td>TAP1</td></tr>
	<tr><th scope=row>TAP2</th><td>2.070688e-253</td><td> 1.8187480</td><td>0.779</td><td>0.485</td><td>3.715436e-249</td><td> 1.8187480</td><td>1e-19</td><td>TAP2    </td><td>UP</td><td>TAP2  </td><td>TAP2</td></tr>
	<tr><th scope=row>CD74</th><td>1.321012e-155</td><td> 0.7369142</td><td>0.995</td><td>0.981</td><td>2.370293e-151</td><td> 0.7369142</td><td>1e-19</td><td>CD74    </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>HLA-DRA</th><td>4.071785e-136</td><td> 0.7700200</td><td>0.970</td><td>0.948</td><td>7.306004e-132</td><td> 0.7700200</td><td>1e-19</td><td>HLA-DRA </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>IFITM3</th><td>2.012365e-121</td><td> 0.7131448</td><td>0.940</td><td>0.906</td><td>3.610787e-117</td><td> 0.7131448</td><td>1e-19</td><td>IFITM3  </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>VIM</th><td>4.030422e-116</td><td> 0.6374194</td><td>0.971</td><td>0.962</td><td>7.231787e-112</td><td> 0.6374194</td><td>1e-19</td><td>VIM     </td><td>NO</td><td>NA    </td><td>VIM </td></tr>
	<tr><th scope=row>C1QB</th><td>8.725089e-105</td><td> 0.7204158</td><td>0.919</td><td>0.872</td><td>1.565543e-100</td><td> 0.7204158</td><td>1e-19</td><td>C1QB    </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>S100A6</th><td> 1.069319e-97</td><td> 0.4796522</td><td>0.985</td><td>0.977</td><td> 1.918678e-93</td><td> 0.4796522</td><td>1e-19</td><td>S100A6  </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>S100A10</th><td> 1.421101e-91</td><td> 0.7294962</td><td>0.878</td><td>0.845</td><td> 2.549882e-87</td><td> 0.7294962</td><td>1e-19</td><td>S100A10 </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>CTSB</th><td> 5.734065e-91</td><td> 0.5737418</td><td>0.953</td><td>0.925</td><td> 1.028863e-86</td><td> 0.5737418</td><td>1e-19</td><td>CTSB    </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>LYZ</th><td> 4.461857e-86</td><td> 0.8472520</td><td>0.834</td><td>0.809</td><td> 8.005910e-82</td><td> 0.8472520</td><td>1e-19</td><td>LYZ     </td><td>UP</td><td>LYZ   </td><td>NA  </td></tr>
	<tr><th scope=row>LAPTM5</th><td> 1.811341e-81</td><td> 0.6645867</td><td>0.879</td><td>0.851</td><td> 3.250089e-77</td><td> 0.6645867</td><td>1e-19</td><td>LAPTM5  </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>WARS</th><td> 1.436961e-78</td><td> 0.6097824</td><td>0.985</td><td>0.968</td><td> 2.578339e-74</td><td> 0.6097824</td><td>1e-19</td><td>WARS    </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>HLA-DPA1</th><td> 2.486971e-77</td><td> 0.7884652</td><td>0.828</td><td>0.775</td><td> 4.462373e-73</td><td> 0.7884652</td><td>1e-19</td><td>HLA-DPA1</td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>BGN</th><td> 3.774688e-76</td><td> 0.8372316</td><td>0.812</td><td>0.784</td><td> 6.772922e-72</td><td> 0.8372316</td><td>1e-19</td><td>BGN     </td><td>UP</td><td>BGN   </td><td>NA  </td></tr>
	<tr><th scope=row>APOE</th><td> 7.233996e-73</td><td> 0.7209996</td><td>0.863</td><td>0.839</td><td> 1.297996e-68</td><td> 0.7209996</td><td>1e-19</td><td>APOE    </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>MT2A</th><td> 2.330259e-70</td><td> 0.5031570</td><td>0.995</td><td>0.990</td><td> 4.181184e-66</td><td> 0.5031570</td><td>1e-19</td><td>MT2A    </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>MIF</th><td> 2.305805e-63</td><td>-0.2501502</td><td>0.999</td><td>0.997</td><td> 4.137306e-59</td><td>-0.2501502</td><td>1e-19</td><td>MIF     </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>TYROBP</th><td> 2.697552e-60</td><td> 0.7178357</td><td>0.772</td><td>0.714</td><td> 4.840218e-56</td><td> 0.7178357</td><td>1e-19</td><td>TYROBP  </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>HLA-DPB1</th><td> 1.255953e-58</td><td> 0.6789488</td><td>0.800</td><td>0.770</td><td> 2.253557e-54</td><td> 0.6789488</td><td>1e-19</td><td>HLA-DPB1</td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>HIST1H4C</th><td> 1.286125e-58</td><td>-0.2510312</td><td>0.998</td><td>0.997</td><td> 2.307693e-54</td><td>-0.2510312</td><td>1e-19</td><td>HIST1H4C</td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>TAGLN2</th><td> 1.955413e-58</td><td> 0.3812701</td><td>0.977</td><td>0.975</td><td> 3.508598e-54</td><td> 0.3812701</td><td>1e-19</td><td>TAGLN2  </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>LUM</th><td> 2.288667e-58</td><td> 0.9701985</td><td>0.693</td><td>0.630</td><td> 4.106556e-54</td><td> 0.9701985</td><td>1e-19</td><td>LUM     </td><td>UP</td><td>LUM   </td><td>NA  </td></tr>
	<tr><th scope=row>PLAAT4</th><td> 7.329870e-58</td><td> 0.9603808</td><td>0.647</td><td>0.523</td><td> 1.315199e-53</td><td> 0.9603808</td><td>1e-19</td><td>PLAAT4  </td><td>UP</td><td>PLAAT4</td><td>NA  </td></tr>
	<tr><th scope=row>A2M</th><td> 8.142055e-58</td><td> 0.6148800</td><td>0.840</td><td>0.795</td><td> 1.460929e-53</td><td> 0.6148800</td><td>1e-19</td><td>A2M     </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>IFI27</th><td> 3.595947e-57</td><td> 0.6052083</td><td>0.894</td><td>0.782</td><td> 6.452208e-53</td><td> 0.6052083</td><td>1e-19</td><td>IFI27   </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>HLA-F</th><td> 3.697085e-56</td><td> 0.7414117</td><td>0.785</td><td>0.684</td><td> 6.633680e-52</td><td> 0.7414117</td><td>1e-19</td><td>HLA-F   </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>CTSZ</th><td> 6.190173e-55</td><td> 0.4248343</td><td>0.917</td><td>0.909</td><td> 1.110703e-50</td><td> 0.4248343</td><td>1e-19</td><td>CTSZ    </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
	<tr><th scope=row>IFITM2</th><td> 3.306573e-53</td><td> 0.6257996</td><td>0.787</td><td>0.765</td><td> 5.932984e-49</td><td> 0.6257996</td><td>1e-19</td><td>IFITM2  </td><td>NO</td><td>NA    </td><td>NA  </td></tr>
</tbody>
</table>



    [1] "B2M"   "TAP1"  "TAP2"  NA      "VIM"   "HLA-A" "SPARC" "CXCL9"



```R
list2 <- c('ASCL1', 'CHGA', 'INSM1','POU2F3','PROX1', 'ISL1')
test1 <- df %>% subset(subset=(gene_name %in% list2))
test1
```


<table class="dataframe">
<caption>A data.frame: 6 × 11</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FoldChange</th><th scope=col>padj</th><th scope=col>gene_name</th><th scope=col>diff_exp</th><th scope=col>gene_label</th><th scope=col>gene_label2</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ASCL1</th><td>2.581329e-41</td><td>-0.41837860</td><td>0.796</td><td>0.863</td><td>4.631679e-37</td><td>-0.41837860</td><td>1.000000e-19</td><td>ASCL1 </td><td>NO</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>CHGA</th><td>2.542059e-32</td><td>-0.68947040</td><td>0.377</td><td>0.523</td><td>4.561217e-28</td><td>-0.68947040</td><td>1.000000e-19</td><td>CHGA  </td><td>NO</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>ISL1</th><td>1.308378e-25</td><td>-0.75136830</td><td>0.193</td><td>0.335</td><td>2.347623e-21</td><td>-0.75136830</td><td>1.023476e-19</td><td>ISL1  </td><td>NO</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>INSM1</th><td>3.397364e-19</td><td>-0.46457130</td><td>0.379</td><td>0.511</td><td>6.095891e-15</td><td>-0.46457130</td><td>6.095991e-15</td><td>INSM1 </td><td>NO</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>PROX1</th><td>1.510856e-10</td><td>-0.06026615</td><td>0.265</td><td>0.376</td><td>2.710929e-06</td><td>-0.06026615</td><td>2.710929e-06</td><td>PROX1 </td><td>NO</td><td>NA</td><td>NA</td></tr>
	<tr><th scope=row>POU2F3</th><td>6.234851e-01</td><td> 0.18785200</td><td>0.047</td><td>0.044</td><td>1.000000e+00</td><td> 0.18785200</td><td>1.000000e+00</td><td>POU2F3</td><td>NO</td><td>NA</td><td>NA</td></tr>
</tbody>
</table>




```R
old_goi_list <- c('HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'B2M', 'CXCL9', 'CXCL10', 'AXL', 'PRRX1', 'VIM', 'SPARC', 'ZEB1', 'SNAI1', 'SNAI2', 'SNAI3')
test1 <- df %>% subset(subset=(gene_name %in% old_goi_list))
test1
```


<table class="dataframe">
<caption>A data.frame: 15 × 11</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FoldChange</th><th scope=col>padj</th><th scope=col>gene_name</th><th scope=col>diff_exp</th><th scope=col>gene_label</th><th scope=col>gene_label2</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B2M</th><td> 0.000000e+00</td><td> 0.6443753</td><td>1.000</td><td>1.000</td><td> 0.000000e+00</td><td> 0.6443753</td><td>1.000000e-19</td><td>B2M   </td><td>NO</td><td>NA    </td><td>B2M  </td></tr>
	<tr><th scope=row>TAP1</th><td>3.747407e-285</td><td> 1.5546950</td><td>0.907</td><td>0.622</td><td>6.723973e-281</td><td> 1.5546950</td><td>1.000000e-19</td><td>TAP1  </td><td>UP</td><td>TAP1  </td><td>TAP1 </td></tr>
	<tr><th scope=row>TAP2</th><td>2.070688e-253</td><td> 1.8187480</td><td>0.779</td><td>0.485</td><td>3.715436e-249</td><td> 1.8187480</td><td>1.000000e-19</td><td>TAP2  </td><td>UP</td><td>TAP2  </td><td>TAP2 </td></tr>
	<tr><th scope=row>VIM</th><td>4.030422e-116</td><td> 0.6374194</td><td>0.971</td><td>0.962</td><td>7.231787e-112</td><td> 0.6374194</td><td>1.000000e-19</td><td>VIM   </td><td>NO</td><td>NA    </td><td>VIM  </td></tr>
	<tr><th scope=row>HLA-A</th><td> 8.661541e-32</td><td> 0.5908459</td><td>0.597</td><td>0.454</td><td> 1.554140e-27</td><td> 0.5908459</td><td>1.000000e-19</td><td>HLA-A </td><td>NO</td><td>NA    </td><td>HLA-A</td></tr>
	<tr><th scope=row>SPARC</th><td> 2.524304e-31</td><td> 0.6966320</td><td>0.798</td><td>0.784</td><td> 4.529358e-27</td><td> 0.6966320</td><td>1.000000e-19</td><td>SPARC </td><td>NO</td><td>NA    </td><td>SPARC</td></tr>
	<tr><th scope=row>CXCL9</th><td> 4.168324e-28</td><td> 0.8689559</td><td>0.574</td><td>0.522</td><td> 7.479225e-24</td><td> 0.8689559</td><td>1.000075e-19</td><td>CXCL9 </td><td>UP</td><td>CXCL9 </td><td>CXCL9</td></tr>
	<tr><th scope=row>PRRX1</th><td> 3.793255e-15</td><td> 0.8658513</td><td>0.369</td><td>0.301</td><td> 6.806237e-11</td><td> 0.8658513</td><td>6.806237e-11</td><td>PRRX1 </td><td>UP</td><td>PRRX1 </td><td>NA   </td></tr>
	<tr><th scope=row>CXCL10</th><td> 7.037682e-07</td><td> 0.9956671</td><td>0.139</td><td>0.099</td><td> 1.262771e-02</td><td> 0.9956671</td><td>1.262771e-02</td><td>CXCL10</td><td>UP</td><td>CXCL10</td><td>NA   </td></tr>
	<tr><th scope=row>AXL</th><td> 7.267088e-04</td><td> 0.5337296</td><td>0.336</td><td>0.325</td><td> 1.000000e+00</td><td> 0.5337296</td><td>1.000000e+00</td><td>AXL   </td><td>NO</td><td>NA    </td><td>NA   </td></tr>
	<tr><th scope=row>SNAI3</th><td> 4.313011e-03</td><td>-0.2506671</td><td>0.015</td><td>0.028</td><td> 1.000000e+00</td><td>-0.2506671</td><td>1.000000e+00</td><td>SNAI3 </td><td>NO</td><td>NA    </td><td>NA   </td></tr>
	<tr><th scope=row>SNAI2</th><td> 6.340428e-03</td><td> 0.7717588</td><td>0.207</td><td>0.189</td><td> 1.000000e+00</td><td> 0.7717588</td><td>1.000000e+00</td><td>SNAI2 </td><td>NO</td><td>NA    </td><td>NA   </td></tr>
	<tr><th scope=row>SNAI1</th><td> 8.073293e-03</td><td> 0.1270014</td><td>0.077</td><td>0.102</td><td> 1.000000e+00</td><td> 0.1270014</td><td>1.000000e+00</td><td>SNAI1 </td><td>NO</td><td>NA    </td><td>NA   </td></tr>
	<tr><th scope=row>ZEB1</th><td> 2.753795e-01</td><td> 0.4908157</td><td>0.245</td><td>0.252</td><td> 1.000000e+00</td><td> 0.4908157</td><td>1.000000e+00</td><td>ZEB1  </td><td>NO</td><td>NA    </td><td>NA   </td></tr>
	<tr><th scope=row>HLA-C</th><td> 7.444542e-01</td><td> 0.3052066</td><td>0.019</td><td>0.021</td><td> 1.000000e+00</td><td> 0.3052066</td><td>1.000000e+00</td><td>HLA-C </td><td>NO</td><td>NA    </td><td>NA   </td></tr>
</tbody>
</table>




```R

```

# try plotting volcano plot again differently


```R
# A function to make the volcano plot
plot_volcano_2 <- function(df, group_name_1, group_name_2, save_dir) {

    # add a very small number to padj
    df$padj <- df$padj + 0.0000000000000000001
    df$p_val <- df$p_val + 0.0000000000000000001
    
    # Make the volcano plot
    save_filename = paste(group_name_1, "vs", group_name_2, sep="_")

    # genes of interest
    old_goi_list <- c('HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'B2M', 'CXCL9', 'CXCL10', 'AXL', 'PRRX1', 'VIM', 'SPARC', 'ZEB1', 'SNAI1', 'SNAI2', 'SNAI3')
    goi_list <- c('TAP1', 'TAP2', 'B2M', 'CXCL9', 'VIM', 'SPARC', 'HLA-A')
    
    #vplot <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val), col=diff_exp, label=gene_label2)) +
    #vplot <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val), label=gene_label2)) +
    vplot <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(padj), label=gene_label2)) +
      geom_point(data = df[!df$gene_label %in% goi_list,], color = "grey", alpha=0.65) + 
      theme_minimal() +
      geom_vline(xintercept=c(-0.5, 0.5), col="grey") +
      geom_hline(yintercept=-log10(0.05), col="grey") +
      geom_text_repel(max.overlaps = 30, color = "black", box.padding=0.4, size=4, force=2) +
      #geom_label_repel(max.overlaps = 30, color = "black") +
      geom_point(data = df[df$gene_label2 %in% goi_list,], color = "black") +
      labs(title=save_filename) +
      #guides(fill = guide_legend(override.aes = aes(label = ""))) +
      xlim(min(df$log2FoldChange), max(df$avg_log2FC)) + 
      ylim(0, max(-log10(df$padj))+3)
    
    # Save the volcano plot
    png(paste0(paste(save_dir, save_filename, save_filename, sep="/"),".png"), width=1400, height=1400, res=300)
    print(vplot)
    dev.off()
    print('plot_volcano complete')

    # label instead of text
    vplot2 <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(padj), label=gene_label2)) +
      geom_point(data = df[!df$gene_label %in% goi_list,], color = "grey", alpha=0.65) + 
      theme_minimal() +
      geom_vline(xintercept=c(-0.5, 0.5), col="grey") +
      geom_hline(yintercept=-log10(0.05), col="grey") +
      #geom_text_repel(max.overlaps = 30, color = "black", box.padding=0.4, size=4) +
      geom_label_repel(max.overlaps = 30, color = "black", box.padding=0.4, size=4, force=2) +
      geom_point(data = df[df$gene_label2 %in% goi_list,], color = "black") +
      labs(title=save_filename) +
      #guides(fill = guide_legend(override.aes = aes(label = ""))) +
      xlim(min(df$log2FoldChange), max(df$avg_log2FC)) + 
      ylim(0, max(-log10(df$padj))+3)

    # Save the volcano plot
    png(paste0(paste(save_dir, save_filename, save_filename, sep="/"),"_label.png"), width=1400, height=1400, res=300)
    print(vplot2)
    dev.off()

    # Return the 1st volcano plot
    return(vplot)
}

vplot1 <- plot_volcano_2(df, 'Tumor_MHCI_High', 'Tumor_MHCI_Low_for_GSEA', '/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/15_DE_and_GSEA')
# Display the volcano plot
vplot1
```

    Warning message:
    “[1m[22mRemoved 17936 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).”


    [1] "plot_volcano complete"


    Warning message:
    “[1m[22mRemoved 17936 rows containing missing values or values outside the scale range
    (`geom_label_repel()`).”
    Warning message:
    “[1m[22mRemoved 17936 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).”



    
![png](output_37_3.png)
    


# Remaking Figure 3E, the volcano plot, for the paper

# Use the DE genes specifying counts


```R
# load de genes again
save_filename = 'Tumor_MHCI_High_vs_Tumor_MHCI_Low'
save_dir = '/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/15_DE_and_GSEA'
load_file = paste0(paste(save_dir, save_filename, save_filename, sep="/"),".csv")

df <- read.csv(load_file)
df <- df %>% column_to_rownames('X')
df %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 10</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FoldChange</th><th scope=col>padj</th><th scope=col>gene_name</th><th scope=col>diff_exp</th><th scope=col>gene_label</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B2M</th><td> 0.000000e+00</td><td>0.1686889</td><td>1.000</td><td>1.000</td><td> 0.000000e+00</td><td>0.1686889</td><td>1e-19</td><td>B2M </td><td>NO</td><td>NA</td></tr>
	<tr><th scope=row>TAP1</th><td>3.747407e-285</td><td>0.6774325</td><td>0.907</td><td>0.622</td><td>6.723973e-281</td><td>0.6774325</td><td>1e-19</td><td>TAP1</td><td>NO</td><td>NA</td></tr>
	<tr><th scope=row>TAP2</th><td>2.070688e-253</td><td>0.7210781</td><td>0.779</td><td>0.485</td><td>3.715436e-249</td><td>0.7210781</td><td>1e-19</td><td>TAP2</td><td>NO</td><td>NA</td></tr>
</tbody>
</table>




```R
# check genes of interest
old_goi_list <- c('HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'B2M', 'CXCL9', 'CXCL10', 'AXL', 'PRRX1', 'VIM', 'SPARC', 'ZEB1', 'SNAI1', 'SNAI2', 'SNAI3')
#goi_list <- c('TAP1', 'TAP2', 'B2M', 'CXCL9', 'VIM', 'SPARC')
goi_list <- c('TAP1', 'TAP2', 'B2M', 'CXCL9', 'VIM', 'SPARC', 'HLA-A')

gene_name_list <- unique(df$gene_name)

both_list <- intersect(goi_list, gene_name_list)
only_goi_list <- setdiff(goi_list, gene_name_list)

print(both_list)
print(only_goi_list)
```

    [1] "TAP1"  "TAP2"  "B2M"   "CXCL9" "VIM"   "SPARC" "HLA-A"
    character(0)



```R
# label genes of interest
# Annotate gene name labels for the differentially expressed genes
df$gene_label2 <- NA
df$gene_label2[df$gene_name %in% goi_list] <- df$gene_name[df$gene_name %in% goi_list]
```


```R
df %>% head(30)
print(unique(df$gene_label2))
```


<table class="dataframe">
<caption>A data.frame: 30 × 11</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FoldChange</th><th scope=col>padj</th><th scope=col>gene_name</th><th scope=col>diff_exp</th><th scope=col>gene_label</th><th scope=col>gene_label2</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B2M</th><td> 0.000000e+00</td><td> 0.16868886</td><td>1.000</td><td>1.000</td><td> 0.000000e+00</td><td> 0.16868886</td><td>1e-19</td><td>B2M     </td><td>NO  </td><td>NA      </td><td>B2M </td></tr>
	<tr><th scope=row>TAP1</th><td>3.747407e-285</td><td> 0.67743252</td><td>0.907</td><td>0.622</td><td>6.723973e-281</td><td> 0.67743252</td><td>1e-19</td><td>TAP1    </td><td>NO  </td><td>NA      </td><td>TAP1</td></tr>
	<tr><th scope=row>TAP2</th><td>2.070688e-253</td><td> 0.72107809</td><td>0.779</td><td>0.485</td><td>3.715436e-249</td><td> 0.72107809</td><td>1e-19</td><td>TAP2    </td><td>NO  </td><td>NA      </td><td>TAP2</td></tr>
	<tr><th scope=row>CD74</th><td>1.321012e-155</td><td> 0.48203770</td><td>0.995</td><td>0.981</td><td>2.370293e-151</td><td> 0.48203770</td><td>1e-19</td><td>CD74    </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>HLA-DRA</th><td>4.071785e-136</td><td> 0.50698246</td><td>0.970</td><td>0.948</td><td>7.306004e-132</td><td> 0.50698246</td><td>1e-19</td><td>HLA-DRA </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>IFITM3</th><td>2.012365e-121</td><td> 0.47289614</td><td>0.940</td><td>0.906</td><td>3.610787e-117</td><td> 0.47289614</td><td>1e-19</td><td>IFITM3  </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>VIM</th><td>4.030422e-116</td><td> 0.43505663</td><td>0.971</td><td>0.962</td><td>7.231787e-112</td><td> 0.43505663</td><td>1e-19</td><td>VIM     </td><td>NO  </td><td>NA      </td><td>VIM </td></tr>
	<tr><th scope=row>C1QB</th><td>8.725089e-105</td><td> 0.40864262</td><td>0.919</td><td>0.872</td><td>1.565543e-100</td><td> 0.40864262</td><td>1e-19</td><td>C1QB    </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>S100A6</th><td> 1.069319e-97</td><td> 0.15015020</td><td>0.985</td><td>0.977</td><td> 1.918678e-93</td><td> 0.15015020</td><td>1e-19</td><td>S100A6  </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>S100A10</th><td> 1.421101e-91</td><td> 0.58528430</td><td>0.878</td><td>0.845</td><td> 2.549882e-87</td><td> 0.58528430</td><td>1e-19</td><td>S100A10 </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>CTSB</th><td> 5.734065e-91</td><td> 0.14103873</td><td>0.953</td><td>0.925</td><td> 1.028863e-86</td><td> 0.14103873</td><td>1e-19</td><td>CTSB    </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>LYZ</th><td> 4.461857e-86</td><td> 0.57384290</td><td>0.834</td><td>0.809</td><td> 8.005910e-82</td><td> 0.57384290</td><td>1e-19</td><td>LYZ     </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>LAPTM5</th><td> 1.811341e-81</td><td> 0.34508149</td><td>0.879</td><td>0.851</td><td> 3.250089e-77</td><td> 0.34508149</td><td>1e-19</td><td>LAPTM5  </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>WARS</th><td> 1.436961e-78</td><td> 0.01142795</td><td>0.985</td><td>0.968</td><td> 2.578339e-74</td><td> 0.01142795</td><td>1e-19</td><td>WARS    </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>HLA-DPA1</th><td> 2.486971e-77</td><td> 0.40467400</td><td>0.828</td><td>0.775</td><td> 4.462373e-73</td><td> 0.40467400</td><td>1e-19</td><td>HLA-DPA1</td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>BGN</th><td> 3.774688e-76</td><td> 0.57389508</td><td>0.812</td><td>0.784</td><td> 6.772922e-72</td><td> 0.57389508</td><td>1e-19</td><td>BGN     </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>APOE</th><td> 7.233996e-73</td><td> 0.39663131</td><td>0.863</td><td>0.839</td><td> 1.297996e-68</td><td> 0.39663131</td><td>1e-19</td><td>APOE    </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>MT2A</th><td> 2.330259e-70</td><td> 0.31906573</td><td>0.995</td><td>0.990</td><td> 4.181184e-66</td><td> 0.31906573</td><td>1e-19</td><td>MT2A    </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>MIF</th><td> 2.305805e-63</td><td>-0.86214103</td><td>0.999</td><td>0.997</td><td> 4.137306e-59</td><td>-0.86214103</td><td>1e-19</td><td>MIF     </td><td>DOWN</td><td>MIF     </td><td>NA  </td></tr>
	<tr><th scope=row>TYROBP</th><td> 2.697552e-60</td><td> 0.45337595</td><td>0.772</td><td>0.714</td><td> 4.840218e-56</td><td> 0.45337595</td><td>1e-19</td><td>TYROBP  </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>HLA-DPB1</th><td> 1.255953e-58</td><td> 0.26864238</td><td>0.800</td><td>0.770</td><td> 2.253557e-54</td><td> 0.26864238</td><td>1e-19</td><td>HLA-DPB1</td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>HIST1H4C</th><td> 1.286125e-58</td><td>-0.84892703</td><td>0.998</td><td>0.997</td><td> 2.307693e-54</td><td>-0.84892703</td><td>1e-19</td><td>HIST1H4C</td><td>DOWN</td><td>HIST1H4C</td><td>NA  </td></tr>
	<tr><th scope=row>TAGLN2</th><td> 1.955413e-58</td><td> 0.05015211</td><td>0.977</td><td>0.975</td><td> 3.508598e-54</td><td> 0.05015211</td><td>1e-19</td><td>TAGLN2  </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>LUM</th><td> 2.288667e-58</td><td> 0.82351239</td><td>0.693</td><td>0.630</td><td> 4.106556e-54</td><td> 0.82351239</td><td>1e-19</td><td>LUM     </td><td>UP  </td><td>LUM     </td><td>NA  </td></tr>
	<tr><th scope=row>PLAAT4</th><td> 7.329870e-58</td><td> 0.68634823</td><td>0.647</td><td>0.523</td><td> 1.315199e-53</td><td> 0.68634823</td><td>1e-19</td><td>PLAAT4  </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>A2M</th><td> 8.142055e-58</td><td> 0.28418797</td><td>0.840</td><td>0.795</td><td> 1.460929e-53</td><td> 0.28418797</td><td>1e-19</td><td>A2M     </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>IFI27</th><td> 3.595947e-57</td><td> 0.61636481</td><td>0.894</td><td>0.782</td><td> 6.452208e-53</td><td> 0.61636481</td><td>1e-19</td><td>IFI27   </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>HLA-F</th><td> 3.697085e-56</td><td> 0.53293466</td><td>0.785</td><td>0.684</td><td> 6.633680e-52</td><td> 0.53293466</td><td>1e-19</td><td>HLA-F   </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>CTSZ</th><td> 6.190173e-55</td><td> 0.06233561</td><td>0.917</td><td>0.909</td><td> 1.110703e-50</td><td> 0.06233561</td><td>1e-19</td><td>CTSZ    </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
	<tr><th scope=row>IFITM2</th><td> 3.306573e-53</td><td> 0.38775314</td><td>0.787</td><td>0.765</td><td> 5.932984e-49</td><td> 0.38775314</td><td>1e-19</td><td>IFITM2  </td><td>NO  </td><td>NA      </td><td>NA  </td></tr>
</tbody>
</table>



    [1] "B2M"   "TAP1"  "TAP2"  NA      "VIM"   "HLA-A" "SPARC" "CXCL9"



```R
list2 <- c('ASCL1', 'CHGA', 'INSM1','POU2F3','PROX1', 'ISL1')
test1 <- df %>% subset(subset=(gene_name %in% list2))
test1
```


<table class="dataframe">
<caption>A data.frame: 6 × 11</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FoldChange</th><th scope=col>padj</th><th scope=col>gene_name</th><th scope=col>diff_exp</th><th scope=col>gene_label</th><th scope=col>gene_label2</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ASCL1</th><td>2.581329e-41</td><td>-0.8329161</td><td>0.796</td><td>0.863</td><td>4.631679e-37</td><td>-0.8329161</td><td>1.000000e-19</td><td>ASCL1 </td><td>DOWN</td><td>ASCL1</td><td>NA</td></tr>
	<tr><th scope=row>CHGA</th><td>2.542059e-32</td><td>-0.9855884</td><td>0.377</td><td>0.523</td><td>4.561217e-28</td><td>-0.9855884</td><td>1.000000e-19</td><td>CHGA  </td><td>DOWN</td><td>CHGA </td><td>NA</td></tr>
	<tr><th scope=row>ISL1</th><td>1.308378e-25</td><td>-1.0720142</td><td>0.193</td><td>0.335</td><td>2.347623e-21</td><td>-1.0720142</td><td>1.023476e-19</td><td>ISL1  </td><td>DOWN</td><td>ISL1 </td><td>NA</td></tr>
	<tr><th scope=row>INSM1</th><td>3.397364e-19</td><td>-0.8329417</td><td>0.379</td><td>0.511</td><td>6.095891e-15</td><td>-0.8329417</td><td>6.095991e-15</td><td>INSM1 </td><td>DOWN</td><td>INSM1</td><td>NA</td></tr>
	<tr><th scope=row>PROX1</th><td>1.510856e-10</td><td>-0.7364917</td><td>0.265</td><td>0.376</td><td>2.710929e-06</td><td>-0.7364917</td><td>2.710929e-06</td><td>PROX1 </td><td>NO  </td><td>NA   </td><td>NA</td></tr>
	<tr><th scope=row>POU2F3</th><td>6.234851e-01</td><td> 0.1239876</td><td>0.047</td><td>0.044</td><td>1.000000e+00</td><td> 0.1239876</td><td>1.000000e+00</td><td>POU2F3</td><td>NO  </td><td>NA   </td><td>NA</td></tr>
</tbody>
</table>




```R
old_goi_list <- c('HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'B2M', 'CXCL9', 'CXCL10', 'AXL', 'PRRX1', 'VIM', 'SPARC', 'ZEB1', 'SNAI1', 'SNAI2', 'SNAI3')
test1 <- df %>% subset(subset=(gene_name %in% old_goi_list))
test1
```


<table class="dataframe">
<caption>A data.frame: 15 × 11</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>log2FoldChange</th><th scope=col>padj</th><th scope=col>gene_name</th><th scope=col>diff_exp</th><th scope=col>gene_label</th><th scope=col>gene_label2</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B2M</th><td> 0.000000e+00</td><td> 0.16868886</td><td>1.000</td><td>1.000</td><td> 0.000000e+00</td><td> 0.16868886</td><td>1.000000e-19</td><td>B2M   </td><td>NO</td><td>NA</td><td>B2M  </td></tr>
	<tr><th scope=row>TAP1</th><td>3.747407e-285</td><td> 0.67743252</td><td>0.907</td><td>0.622</td><td>6.723973e-281</td><td> 0.67743252</td><td>1.000000e-19</td><td>TAP1  </td><td>NO</td><td>NA</td><td>TAP1 </td></tr>
	<tr><th scope=row>TAP2</th><td>2.070688e-253</td><td> 0.72107809</td><td>0.779</td><td>0.485</td><td>3.715436e-249</td><td> 0.72107809</td><td>1.000000e-19</td><td>TAP2  </td><td>NO</td><td>NA</td><td>TAP2 </td></tr>
	<tr><th scope=row>VIM</th><td>4.030422e-116</td><td> 0.43505663</td><td>0.971</td><td>0.962</td><td>7.231787e-112</td><td> 0.43505663</td><td>1.000000e-19</td><td>VIM   </td><td>NO</td><td>NA</td><td>VIM  </td></tr>
	<tr><th scope=row>HLA-A</th><td> 8.661541e-32</td><td> 0.65483235</td><td>0.597</td><td>0.454</td><td> 1.554140e-27</td><td> 0.65483235</td><td>1.000000e-19</td><td>HLA-A </td><td>NO</td><td>NA</td><td>HLA-A</td></tr>
	<tr><th scope=row>SPARC</th><td> 2.524304e-31</td><td> 0.35980863</td><td>0.798</td><td>0.784</td><td> 4.529358e-27</td><td> 0.35980863</td><td>1.000000e-19</td><td>SPARC </td><td>NO</td><td>NA</td><td>SPARC</td></tr>
	<tr><th scope=row>CXCL9</th><td> 4.168324e-28</td><td> 0.54323738</td><td>0.574</td><td>0.522</td><td> 7.479225e-24</td><td> 0.54323738</td><td>1.000075e-19</td><td>CXCL9 </td><td>NO</td><td>NA</td><td>CXCL9</td></tr>
	<tr><th scope=row>PRRX1</th><td> 3.793255e-15</td><td> 0.58226368</td><td>0.369</td><td>0.301</td><td> 6.806237e-11</td><td> 0.58226368</td><td>6.806237e-11</td><td>PRRX1 </td><td>NO</td><td>NA</td><td>NA   </td></tr>
	<tr><th scope=row>CXCL10</th><td> 7.037682e-07</td><td> 0.56731854</td><td>0.139</td><td>0.099</td><td> 1.262771e-02</td><td> 0.56731854</td><td>1.262771e-02</td><td>CXCL10</td><td>NO</td><td>NA</td><td>NA   </td></tr>
	<tr><th scope=row>AXL</th><td> 7.267088e-04</td><td> 0.25174332</td><td>0.336</td><td>0.325</td><td> 1.000000e+00</td><td> 0.25174332</td><td>1.000000e+00</td><td>AXL   </td><td>NO</td><td>NA</td><td>NA   </td></tr>
	<tr><th scope=row>SNAI3</th><td> 4.313011e-03</td><td>-0.77481650</td><td>0.015</td><td>0.028</td><td> 1.000000e+00</td><td>-0.77481650</td><td>1.000000e+00</td><td>SNAI3 </td><td>NO</td><td>NA</td><td>NA   </td></tr>
	<tr><th scope=row>SNAI2</th><td> 6.340428e-03</td><td> 0.18847159</td><td>0.207</td><td>0.189</td><td> 1.000000e+00</td><td> 0.18847159</td><td>1.000000e+00</td><td>SNAI2 </td><td>NO</td><td>NA</td><td>NA   </td></tr>
	<tr><th scope=row>SNAI1</th><td> 8.073293e-03</td><td>-0.45726424</td><td>0.077</td><td>0.102</td><td> 1.000000e+00</td><td>-0.45726424</td><td>1.000000e+00</td><td>SNAI1 </td><td>NO</td><td>NA</td><td>NA   </td></tr>
	<tr><th scope=row>ZEB1</th><td> 2.753795e-01</td><td> 0.04256594</td><td>0.245</td><td>0.252</td><td> 1.000000e+00</td><td> 0.04256594</td><td>1.000000e+00</td><td>ZEB1  </td><td>NO</td><td>NA</td><td>NA   </td></tr>
	<tr><th scope=row>HLA-C</th><td> 7.444542e-01</td><td>-0.06970030</td><td>0.019</td><td>0.021</td><td> 1.000000e+00</td><td>-0.06970030</td><td>1.000000e+00</td><td>HLA-C </td><td>NO</td><td>NA</td><td>NA   </td></tr>
</tbody>
</table>




```R

```

# try plotting volcano plot again differently


```R
# A function to make the volcano plot
plot_volcano_2 <- function(df, group_name_1, group_name_2, save_dir) {

    # add a very small number to padj
    df$padj <- df$padj + 0.0000000000000000001
    df$p_val <- df$p_val + 0.0000000000000000001
    
    # Make the volcano plot
    save_filename = paste(group_name_1, "vs", group_name_2, sep="_")

    # genes of interest
    old_goi_list <- c('HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'B2M', 'CXCL9', 'CXCL10', 'AXL', 'PRRX1', 'VIM', 'SPARC', 'ZEB1', 'SNAI1', 'SNAI2', 'SNAI3')
    goi_list <- c('TAP1', 'TAP2', 'B2M', 'CXCL9', 'VIM', 'SPARC', 'HLA-A')
    
    #vplot <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val), col=diff_exp, label=gene_label2)) +
    #vplot <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val), label=gene_label2)) +
    vplot <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(padj), label=gene_label2)) +
      geom_point(data = df[!df$gene_label %in% goi_list,], color = "grey", alpha=0.65) + 
      theme_minimal() +
      geom_vline(xintercept=c(-0.5, 0.5), col="grey") +
      geom_hline(yintercept=-log10(0.05), col="grey") +
      geom_text_repel(max.overlaps = 30, color = "black", box.padding=0.4, size=4, force=2) +
      #geom_label_repel(max.overlaps = 30, color = "black") +
      geom_point(data = df[df$gene_label2 %in% goi_list,], color = "black") +
      labs(title=save_filename) +
      #guides(fill = guide_legend(override.aes = aes(label = ""))) +
      xlim(min(df$log2FoldChange), max(df$avg_log2FC)) + 
      ylim(0, max(-log10(df$padj))+3)
    
    # Save the volcano plot
    png(paste0(paste(save_dir, save_filename, save_filename, sep="/"),".png"), width=1400, height=1400, res=300)
    print(vplot)
    dev.off()
    print('plot_volcano complete')

    # label instead of text
    vplot2 <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(padj), label=gene_label2)) +
      geom_point(data = df[!df$gene_label %in% goi_list,], color = "grey", alpha=0.65) + 
      theme_minimal() +
      geom_vline(xintercept=c(-0.5, 0.5), col="grey") +
      geom_hline(yintercept=-log10(0.05), col="grey") +
      #geom_text_repel(max.overlaps = 30, color = "black", box.padding=0.4, size=4) +
      geom_label_repel(max.overlaps = 30, color = "black", box.padding=0.4, size=4, force=2) +
      geom_point(data = df[df$gene_label2 %in% goi_list,], color = "black") +
      labs(title=save_filename) +
      #guides(fill = guide_legend(override.aes = aes(label = ""))) +
      xlim(min(df$log2FoldChange), max(df$avg_log2FC)) + 
      ylim(0, max(-log10(df$padj))+3)

    # Save the volcano plot
    png(paste0(paste(save_dir, save_filename, save_filename, sep="/"),"_label.png"), width=1400, height=1400, res=300)
    print(vplot2)
    dev.off()

    # Return the 1st volcano plot
    return(vplot)
}

vplot1 <- plot_volcano_2(df, 'Tumor_MHCI_High', 'Tumor_MHCI_Low', '/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/15_DE_and_GSEA')
# Display the volcano plot
vplot1
```

    Warning message:
    “[1m[22mRemoved 17936 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).”


    [1] "plot_volcano complete"


    Warning message:
    “[1m[22mRemoved 17936 rows containing missing values or values outside the scale range
    (`geom_label_repel()`).”
    Warning message:
    “[1m[22mRemoved 17936 rows containing missing values or values outside the scale range
    (`geom_text_repel()`).”



    
![png](output_48_3.png)
    


# zoom in
# A function to make the volcano plot
plot_volcano_zoom <- function(df, group_name_1, group_name_2, save_dir) {
    
    # Make the volcano plot
    save_filename = paste(group_name_1, "vs", group_name_2, sep="_")

    # genes of interest
    old_goi_list <- c('HLA-A', 'HLA-B', 'HLA-C', 'TAP1', 'TAP2', 'B2M', 'CXCL9', 'CXCL10', 'AXL', 'PRRX1', 'VIM', 'SPARC', 'ZEB1', 'SNAI1', 'SNAI2', 'SNAI3')
    goi_list <- c('TAP1', 'TAP2', 'B2M', 'CXCL9', 'VIM', 'SPARC')
    
    #vplot <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val), col=diff_exp, label=gene_label2)) +
    vplot_zoom <- ggplot(data=df, aes(x=avg_log2FC, y=-log10(p_val), label=gene_label2)) +
      geom_point(data = df[!df$gene_label %in% goi_list,], color = "grey", alpha=0.65) + 
      theme_minimal() +
      geom_vline(xintercept=c(-0.8, 0.8), col="grey") +
      geom_hline(yintercept=-log10(0.05), col="grey") +
      geom_text_repel(max.overlaps = 30, color = "black", box.padding=0.4, size=4, force=2) +
      #geom_label_repel(max.overlaps = 30, color = "black") +
      geom_point(data = df[df$gene_label2 %in% goi_list,], color = "black") +
      labs(title=save_filename) +
      #guides(fill = guide_legend(override.aes = aes(label = ""))) +
      xlim(-0.8, 3) + 
      ylim(0, 15)
    
    # Save the volcano plot
    png(paste0(paste(save_dir, save_filename, save_filename, sep="/"),"_zoom.png"), width=1400, height=1400, res=300)
    print(vplot_zoom)
    dev.off()
    print('plot_volcano complete')

    # Return the 1st volcano plot
    return(vplot_zoom)
}

vplot_zoom <- plot_volcano_zoom(df, 'Tumor_MHCI_High', 'Tumor_MHCI_Low', 'DESEQ_GSEA_20240816_TumorHigh_MHCIHigh_vs_MHCILow_combinedMHCINegLow')
# Display the volcano plot
vplot_zoom

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
