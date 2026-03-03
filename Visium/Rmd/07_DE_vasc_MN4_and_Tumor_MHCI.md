# Find DE genes between comparisons of interest


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
data@meta.data[1:5,]
```


<table class="dataframe">
<caption>A data.frame: 5 × 57</caption>
<thead>
	<tr><th></th><th scope=col>Tissue_Slice</th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Macrophage</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>T_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>T_Reg</th><th scope=col>NK_Cell</th><th scope=col>B_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Cell</th><th scope=col>LEC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Angiogenesis</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>AyersIFNG</th><th scope=col>Exhaustion</th><th scope=col>Proliferation</th><th scope=col>ICB_Targets</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.4004691</td><td>-0.69531482</td><td>-0.3662634</td><td>-0.7699365</td><td>-0.4507923</td><td>-0.02904416</td><td>-0.5571237</td><td>-0.21851366</td><td>-0.5514423</td><td>-0.3546455</td><td>-0.4311886</td><td> 0.12079973</td><td>-0.8006772</td><td>-0.9052709</td><td>-0.39185941</td><td>-0.5925154</td><td>-0.3079612</td><td>-0.5925154</td><td>-0.5641991</td><td>-0.5501370</td><td>-0.3899313</td><td> 0.1970419</td><td>-0.433578691</td><td>-0.18834557</td><td>-0.2435803</td><td>-0.3726727</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.5186053</td><td>-0.74484775</td><td>-0.4261746</td><td>-0.7808247</td><td>-0.6767274</td><td>-0.24535449</td><td>-0.6095094</td><td>-0.28557023</td><td>-0.6946626</td><td>-0.5335451</td><td>-0.5121988</td><td>-0.01962366</td><td>-0.8379877</td><td>-0.9136749</td><td>-0.46665116</td><td>-0.6967574</td><td>-0.4306655</td><td>-0.6967574</td><td>-0.6308560</td><td>-0.6552919</td><td>-0.9176245</td><td>-0.1840407</td><td>-0.005009761</td><td>-0.28843952</td><td>-0.2999413</td><td>-0.3152315</td><td>Low </td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.4839912</td><td>-0.75372051</td><td>-0.4294818</td><td>-0.7867729</td><td>-0.6833256</td><td>-0.23086109</td><td>-0.1684702</td><td>-0.09813578</td><td>-0.6947400</td><td>-0.3313786</td><td>-0.5001427</td><td> 0.01283602</td><td>-0.8432302</td><td>-0.9216418</td><td>-0.40040051</td><td>-0.6930524</td><td>-0.3711015</td><td>-0.6930524</td><td>-0.6376792</td><td>-0.6627355</td><td>-0.2720599</td><td>-0.1513761</td><td> 0.024318103</td><td>-0.02579441</td><td>-0.3086330</td><td>-0.3965501</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1</td><td>VA1_AAACACCAATAACTGC-1</td><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td>4 </td><td>4 </td><td>VA1</td><td>0.907884</td><td>vasc_low </td><td> 9</td><td>Tumor_low </td><td>Tumor_low </td><td>1</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC+vasculature      </td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td>-0.4569040</td><td>-0.74726759</td><td>-0.1613635</td><td>-0.7629941</td><td>-0.6843917</td><td>-0.27864075</td><td>-0.6120296</td><td>-0.31814057</td><td>-0.6610014</td><td>-0.2153879</td><td>-0.5122670</td><td>-0.08010753</td><td>-0.4308938</td><td>-0.8679722</td><td>-0.42054109</td><td>-0.5740357</td><td>-0.4047324</td><td>-0.5740357</td><td>-0.4656548</td><td>-0.6234628</td><td>-0.7203440</td><td>-0.2440266</td><td>-0.164491205</td><td>-0.35470812</td><td>-0.3034304</td><td>-0.3416859</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td></tr>
	<tr><th scope=row>VA1_AAACAGCTTTCAGAAG-1</th><td>VA1</td><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>5132</td><td>3219</td><td>round1</td><td>4797</td><td>3218</td><td>4 </td><td>4 </td><td>VA1</td><td>3.751279</td><td>vasc_high</td><td>14</td><td>Tumor_low </td><td>Tumor_low </td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low          </td><td>vasc_high</td><td>SCLC+vasculature      </td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td>-0.1728081</td><td> 0.06870843</td><td>-0.4399879</td><td>-0.4423331</td><td>-0.1082299</td><td>-0.02566302</td><td>-0.1807694</td><td>-0.04870685</td><td>-0.4938486</td><td>-0.1624444</td><td>-0.1659391</td><td> 0.17713889</td><td>-0.8267151</td><td>-0.9143730</td><td>-0.06052303</td><td>-0.5525944</td><td>-0.3217238</td><td>-0.5525944</td><td>-0.3681237</td><td>-0.3736337</td><td>-0.2376679</td><td>-0.1671035</td><td> 0.088616684</td><td>-0.23390067</td><td>-0.2670847</td><td>-0.2853179</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td></tr>
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
<ol class=list-inline><li>'Tissue_Slice'</li><li>'Barcode'</li><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li><li>'Navin_annotations'</li><li>'Navin_annotations_simplified'</li><li>'Navin_Clusters1'</li><li>'Navin_Clusters'</li><li>'Navin_Clusters2'</li><li>'All_Immune'</li><li>'Macrophage'</li><li>'M1_Macrophage'</li><li>'M2_Macrophage'</li><li>'T_Cell'</li><li>'CD8_T_Cell'</li><li>'T_Reg'</li><li>'NK_Cell'</li><li>'B_Cell'</li><li>'DC'</li><li>'Endothelial_Cell'</li><li>'LEC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Chemokines'</li><li>'Angiogenesis'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'MN4_EC_Phenotype'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'AyersIFNG'</li><li>'Exhaustion'</li><li>'Proliferation'</li><li>'ICB_Targets'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'MN4_EC_Phenotype_Label'</li><li>'MN4_EC_Percentile_Label'</li><li>'vasc_MN4_label'</li><li>'vasc_MN4_percentile_label'</li></ol>




```R
Assays(data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Spatial'</li><li>'SCT'</li></ol>




```R
# check navin annotations
print(unique(data@meta.data$Navin_Clusters))
print(unique(data@meta.data$Navin_Clusters2))
```

     [1] 10_Respiratory_Epithelium_Lung 6_SCLC_Effector_Immune        
     [3] 4_SCLC                         1_SCLC                        
     [5] 2_SCLC_TSI                     3_SCLC_Stroma                 
     [7] 11_SCLC_Necrotic               0_Lymphocytes_B_Cells         
     [9] 8_Fibroblasts                  9_SCLC_Apoptotic              
    [11] 7_Normal_Lung_Stroma           5_SCLC                        
    12 Levels: 0_Lymphocytes_B_Cells 1_SCLC 2_SCLC_TSI 3_SCLC_Stroma ... 11_SCLC_Necrotic
    [1] 7_10_Normal_Lung       6_SCLC_Effector_Immune 1_4_5_11_SCLC         
    [4] 2_SCLC_TSI             3_SCLC_Stroma          0_Lymphocytes_B_Cells 
    [7] 8_Fibroblasts          9_SCLC_Apoptotic      
    8 Levels: 0_Lymphocytes_B_Cells 1_4_5_11_SCLC 2_SCLC_TSI ... 9_SCLC_Apoptotic



```R
unique(data@meta.data[c('Navin_Clusters','Navin_Clusters2')])
```


<table class="dataframe">
<caption>A data.frame: 12 × 2</caption>
<thead>
	<tr><th></th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td></tr>
	<tr><th scope=row>VA1_AAACCCGAACGAAATC-1</th><td>1_SCLC                        </td><td>1_4_5_11_SCLC         </td></tr>
	<tr><th scope=row>VA1_AAACCGGAAATGTTAA-1</th><td>2_SCLC_TSI                    </td><td>2_SCLC_TSI            </td></tr>
	<tr><th scope=row>VA1_AAACCGTTCGTCCAGG-1</th><td>3_SCLC_Stroma                 </td><td>3_SCLC_Stroma         </td></tr>
	<tr><th scope=row>VA1_AAAGGCTACGGACCAT-1</th><td>11_SCLC_Necrotic              </td><td>1_4_5_11_SCLC         </td></tr>
	<tr><th scope=row>VA1_AAAGGCTCTCGCGCCG-1</th><td>0_Lymphocytes_B_Cells         </td><td>0_Lymphocytes_B_Cells </td></tr>
	<tr><th scope=row>VA1_AAATAAGGTAGTGCCC-1</th><td>8_Fibroblasts                 </td><td>8_Fibroblasts         </td></tr>
	<tr><th scope=row>VA1_AAATGCTCGTTACGTT-1</th><td>9_SCLC_Apoptotic              </td><td>9_SCLC_Apoptotic      </td></tr>
	<tr><th scope=row>VA1_AACGCGAACGGCAACA-1</th><td>7_Normal_Lung_Stroma          </td><td>7_10_Normal_Lung      </td></tr>
	<tr><th scope=row>VA1_AATAGTCGCGAGTCGG-1</th><td>5_SCLC                        </td><td>1_4_5_11_SCLC         </td></tr>
</tbody>
</table>




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
data_VA <- subset(data, subset=Tissue_Slice=='VA1')
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



```R
# check umap
options(repr.plot.width=20, repr.plot.height=21)

navin_umap <- DimPlot(data_VA, reduction = "umap", group.by = c("Navin_Clusters", "Navin_annotations", 
                                                                "Navin_Clusters2", "Navin_annotations_simplified",
                                                                "Tumor_MHCI_label", "vasc_label"), 
                      ncol = 2,
                      cols = dittoColors()[1:15],
                      pt.size=2) +
    theme(text = element_text(size=20))
navin_umap
```


    
![png](output_13_0.png)
    


# Perform DE on vasc_high MN4_high vs MN4_low


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
Idents(data) <- "vasc_MN4_label"

de_vasc_MN4 <- FindMarkers(
  object = data,
  ident.1 = "vasc_high_MN4_high",  # First group to compare
  ident.2 = "vasc_high_MN4_low",    # Second group
  logfc.threshold = 0.25, # Minimum log2 fold change
  min.pct = 0.1,          # Gene detected in ≥10% of cells in either group
  test.use = "wilcox"     # Default Wilcoxon rank sum test
)

de_vasc_MN4$Gene <- rownames(de_vasc_MN4)
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
    



```R
BiocManager::install("EnhancedVolcano")
```

    'getOption("repos")' replaces Bioconductor standard repositories, see
    'help("repositories", package = "BiocManager")' for details.
    Replacement repositories:
        CRAN: https://cloud.r-project.org
    
    Bioconductor version 3.21 (BiocManager 1.30.25), R 4.5.0 (2025-04-11)
    
    Warning message:
    “package(s) not installed when version(s) same as or greater than current; use
      `force = TRUE` to re-install: 'EnhancedVolcano'”
    Old packages: 'BiocManager', 'BiocParallel', 'curl', 'duckdb', 'future',
      'future.apply', 'haven', 'magick', 'parallelly', 'pheatmap', 'rhdf5',
      'RSQLite', 'S4Arrays', 'tibble', 'utf8'
    



```R
library(EnhancedVolcano)
```

    Loading required package: ggrepel
    



```R
# function to label genes as significant up or down
# log2(1.5) is 0.58. A Fold change of 1.5 seems good, so let's use 0.58 as the cutoff.  
label_sig_genes <- function(df, p_val_cutoff = 0.05, avg_log2FC_cutoff = 0.58) {
    df$sig_DE <- "NO"
    df$sig_DE[(df$p_val < p_val_cutoff) & (df$avg_log2FC > avg_log2FC_cutoff)] = "UP"
    df$sig_DE[(df$p_val < p_val_cutoff) & (df$avg_log2FC < -avg_log2FC_cutoff)] = "DOWN"
    return(df)
}
```


```R
de_vasc_MN4 %>% head(15)
```


<table class="dataframe">
<caption>A data.frame: 15 × 6</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>Gene</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>S100A6</th><td>2.147527e-106</td><td> 0.9032005</td><td>0.989</td><td>0.905</td><td>3.853307e-102</td><td>S100A6</td></tr>
	<tr><th scope=row>LYZ</th><td>1.662859e-103</td><td> 1.4309588</td><td>0.907</td><td>0.686</td><td> 2.983668e-99</td><td>LYZ   </td></tr>
	<tr><th scope=row>SPARC</th><td> 5.066896e-99</td><td> 1.3106240</td><td>0.907</td><td>0.658</td><td> 9.091531e-95</td><td>SPARC </td></tr>
	<tr><th scope=row>IGFBP7</th><td> 1.825926e-98</td><td> 1.5237648</td><td>0.817</td><td>0.508</td><td> 3.276258e-94</td><td>IGFBP7</td></tr>
	<tr><th scope=row>DCN</th><td> 3.045898e-95</td><td> 1.3918593</td><td>0.864</td><td>0.609</td><td> 5.465256e-91</td><td>DCN   </td></tr>
	<tr><th scope=row>COL1A2</th><td> 6.042315e-95</td><td> 1.4412885</td><td>0.924</td><td>0.761</td><td> 1.084173e-90</td><td>COL1A2</td></tr>
	<tr><th scope=row>PSAP</th><td> 1.978932e-93</td><td> 1.1167877</td><td>0.971</td><td>0.877</td><td> 3.550798e-89</td><td>PSAP  </td></tr>
	<tr><th scope=row>SFTPB</th><td> 3.275107e-92</td><td> 2.0726144</td><td>0.792</td><td>0.501</td><td> 5.876524e-88</td><td>SFTPB </td></tr>
	<tr><th scope=row>COL3A1</th><td> 4.341552e-89</td><td> 1.4467838</td><td>0.898</td><td>0.719</td><td> 7.790046e-85</td><td>COL3A1</td></tr>
	<tr><th scope=row>CCN1</th><td> 3.323293e-88</td><td> 1.7609487</td><td>0.726</td><td>0.403</td><td> 5.962984e-84</td><td>CCN1  </td></tr>
	<tr><th scope=row>FN1</th><td> 1.072526e-87</td><td> 1.0037082</td><td>0.983</td><td>0.890</td><td> 1.924433e-83</td><td>FN1   </td></tr>
	<tr><th scope=row>CD74</th><td> 1.552077e-87</td><td> 0.7949820</td><td>0.998</td><td>0.964</td><td> 2.784891e-83</td><td>CD74  </td></tr>
	<tr><th scope=row>ICAM1</th><td> 3.989375e-84</td><td> 2.4969689</td><td>0.633</td><td>0.318</td><td> 7.158136e-80</td><td>ICAM1 </td></tr>
	<tr><th scope=row>PPDPF</th><td> 7.581220e-82</td><td>-0.6698764</td><td>0.980</td><td>0.986</td><td> 1.360298e-77</td><td>PPDPF </td></tr>
	<tr><th scope=row>HLA-E</th><td> 1.021765e-79</td><td> 0.9186353</td><td>0.940</td><td>0.729</td><td> 1.833352e-75</td><td>HLA-E </td></tr>
</tbody>
</table>




```R
Vasc_MN4_colors <- c('vasc_low' = '#9d9d9d',
                     'vasc_high_MN4_low' = '#FC8D59FF', 
                     'vasc_high_MN4_high' = '#B30000FF')

# spatial dimplot colored by navin annotations
Idents(data) <- "vasc_MN4_label"

options(repr.plot.width=28, repr.plot.height=9)

plot2 <- SpatialDimPlot(data, cols = Vasc_MN4_colors, pt.size.factor=2)
plot2
```


    
![png](output_21_0.png)
    



```R
Vasc_MN4_colors <- c('Tumor_low' = '#9d9d9d',
                     'Tumor_high_MHCI_low' = '#FC8D59FF', 
                     'Tumor_high_MHCI_high' = '#B30000FF')

# spatial dimplot colored by navin annotations
Idents(data) <- "Tumor_MHCI_label2"

options(repr.plot.width=28, repr.plot.height=9)

plot2 <- SpatialDimPlot(data, cols = Vasc_MN4_colors, pt.size.factor=2)
plot2
```


    
![png](output_22_0.png)
    



```R
options(repr.plot.width=20, repr.plot.height=28)
plot2 <- SpatialFeaturePlot(data, features=c('CD74', 'S100A6', 'VIM', 'DCN'), pt.size.factor=2)
plot2
```


    
![png](output_23_0.png)
    



```R
de_vasc_MN4_sig_label <- label_sig_genes(de_vasc_MN4)
de_vasc_MN4_sig_DE <- de_vasc_MN4_sig_label[de_vasc_MN4_sig_label$sig_DE != 'NO',]
print(dim(de_vasc_MN4))
print(dim(de_vasc_MN4_sig_DE[]))
de_vasc_MN4_sig_DE %>% head(5)
table(de_vasc_MN4_sig_DE$sig_DE)
```

    [1] 6234    6
    [1] 2127    7



<table class="dataframe">
<caption>A data.frame: 5 × 7</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>Gene</th><th scope=col>sig_DE</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>S100A6</th><td>2.147527e-106</td><td>0.9032005</td><td>0.989</td><td>0.905</td><td>3.853307e-102</td><td>S100A6</td><td>UP</td></tr>
	<tr><th scope=row>LYZ</th><td>1.662859e-103</td><td>1.4309588</td><td>0.907</td><td>0.686</td><td> 2.983668e-99</td><td>LYZ   </td><td>UP</td></tr>
	<tr><th scope=row>SPARC</th><td> 5.066896e-99</td><td>1.3106240</td><td>0.907</td><td>0.658</td><td> 9.091531e-95</td><td>SPARC </td><td>UP</td></tr>
	<tr><th scope=row>IGFBP7</th><td> 1.825926e-98</td><td>1.5237648</td><td>0.817</td><td>0.508</td><td> 3.276258e-94</td><td>IGFBP7</td><td>UP</td></tr>
	<tr><th scope=row>DCN</th><td> 3.045898e-95</td><td>1.3918593</td><td>0.864</td><td>0.609</td><td> 5.465256e-91</td><td>DCN   </td><td>UP</td></tr>
</tbody>
</table>




    
    DOWN   UP 
     972 1155 



```R
de_vasc_MN4_check_goi <- de_vasc_MN4 %>% rownames_to_column('gene') %>% subset(gene %in% c('ENG', 'VSIR', 'FASLG'))
de_vasc_MN4_check_goi
```


<table class="dataframe">
<caption>A data.frame: 2 × 7</caption>
<thead>
	<tr><th></th><th scope=col>gene</th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>Gene</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>438</th><td>ENG </td><td>1.139422e-17</td><td>0.7981987</td><td>0.438</td><td>0.288</td><td>2.044466e-13</td><td>ENG </td></tr>
	<tr><th scope=row>1259</th><td>VSIR</td><td>5.252948e-08</td><td>0.8004152</td><td>0.183</td><td>0.103</td><td>9.425364e-04</td><td>VSIR</td></tr>
</tbody>
</table>




```R
options(repr.plot.width=10, repr.plot.height=10)
de_plot <- EnhancedVolcano(de_vasc_MN4,
                lab = rownames(de_vasc_MN4),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'MN4_high vs MN4_low, in vasc_high spots')
de_plot

# save plot
png(file = "Output/07_DE_vasc_MN4_and_Tumor_MHCI/de_vasc_high_MN4_high_vs_low_volcano.png",
    width = 1000,
    height = 1000)
print(de_plot)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_26_1.png)
    



```R
options(repr.plot.width=10, repr.plot.height=12)
Idents(data) <- "vasc_MN4_label"
FeaturePlot(data, features = rownames(de_vasc_MN4)[1:6], reduction='umap')
```

    Warning message:
    “[1m[22mThe `slot` argument of `FetchData()` is deprecated as of SeuratObject 5.0.0.
    [36mℹ[39m Please use the `layer` argument instead.
    [36mℹ[39m The deprecated feature was likely used in the [34mSeurat[39m package.
      Please report the issue at [3m[34m<https://github.com/satijalab/seurat/issues>[39m[23m.”



    
![png](output_27_1.png)
    



```R
options(repr.plot.width=10, repr.plot.height=12)
Idents(data) <- "vasc_MN4_label"
VlnPlot(data, features = rownames(de_vasc_MN4)[1:6], split.by = 'vasc_MN4_label')
```

    The default behaviour of split.by has changed.
    Separate violin plots are now plotted side-by-side.
    To restore the old behaviour of a single split violin,
    set split.plot = TRUE.
          
    This message will be shown once per session.
    



    
![png](output_28_1.png)
    



```R
# Extract lists of genes UP and DOWN regulated
de_vasc_MN4_sig_DE_list <- list()
de_vasc_MN4_sig_DE_list$UP <- rownames(de_vasc_MN4_sig_DE[de_vasc_MN4_sig_DE$sig_DE=="UP",])
de_vasc_MN4_sig_DE_list$DOWN <- rownames(de_vasc_MN4_sig_DE[de_vasc_MN4_sig_DE$sig_DE=="DOWN",])
```


```R
# check list
de_vasc_MN4_sig_DE_list
```


<dl>
	<dt>$UP</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'S100A6'</li><li>'LYZ'</li><li>'SPARC'</li><li>'IGFBP7'</li><li>'DCN'</li><li>'COL1A2'</li><li>'PSAP'</li><li>'SFTPB'</li><li>'COL3A1'</li><li>'CCN1'</li><li>'FN1'</li><li>'CD74'</li><li>'ICAM1'</li><li>'HLA-E'</li><li>'PTGDS'</li><li>'APOE'</li><li>'THBS1'</li><li>'HLA-DMA'</li><li>'CTSD'</li><li>'SLC34A2'</li><li>'GADD45B'</li><li>'UBC'</li><li>'C1R'</li><li>'COL4A1'</li><li>'TPM2'</li><li>'S100A10'</li><li>'CCL18'</li><li>'C11orf96'</li><li>'TIMP1'</li><li>'TIMP3'</li><li>'FLNA'</li><li>'COL6A3'</li><li>'MCL1'</li><li>'VIM'</li><li>'CAVIN1'</li><li>'PLAUR'</li><li>'SERPINE1'</li><li>'JUNB'</li><li>'IGKC'</li><li>'C3'</li><li>'BGN'</li><li>'AQP1'</li><li>'LPCAT1'</li><li>'FOSB'</li><li>'SOCS3'</li><li>'ACTA2'</li><li>'SAT1'</li><li>'G0S2'</li><li>'IFI6'</li><li>'IGHG1'</li><li>'RHOB'</li><li>'ITM2B'</li><li>'ELL2'</li><li>'FBN1'</li><li>'IL32'</li><li>'LMNA'</li><li>'TYROBP'</li><li>'DUSP1'</li><li>'CEBPD'</li><li>'ZFP36'</li><li>'LAPTM5'</li><li>'EMP1'</li><li>'PFKFB3'</li><li>'HLA-DOA'</li><li>'FCGRT'</li><li>'NR4A1'</li><li>'ITGAX'</li><li>'CCN2'</li><li>'IGHJ6'</li><li>'SOD2'</li><li>'CD44'</li><li>'RCAN1'</li><li>'DDIT4'</li><li>'IER3'</li><li>'GSN'</li><li>'C4orf3'</li><li>'PPP1R14A'</li><li>'CFD'</li><li>'TCIM'</li><li>'PDGFRB'</li><li>'S100A4'</li><li>'NCOA7'</li><li>'FOS'</li><li>'JCHAIN'</li><li>'C1QC'</li><li>'RNF213'</li><li>'IGHG3'</li><li>'CRISPLD2'</li><li>'PPP1R15A'</li><li>'ADAMTS4'</li><li>'ATF3'</li><li>'STOM'</li><li>'TRIM22'</li><li>'ARRDC3'</li><li>'CTSZ'</li><li>'POSTN'</li><li>'CIITA'</li><li>'NAPSA'</li><li>'LSP1'</li><li>'COL1A1'</li><li>'TGM2'</li><li>'TCF21'</li><li>'HSPA5'</li><li>'PALM2-AKAP2'</li><li>'OLR1'</li><li>'TYMP'</li><li>'ANXA1'</li><li>'NABP1'</li><li>'PTPRE'</li><li>'ITGB1'</li><li>'HLA-DPB1'</li><li>'MYL9'</li><li>'MXRA7'</li><li>'BHLHE40'</li><li>'IFI30'</li><li>'NR4A3'</li><li>'SLC2A3'</li><li>'SGK1'</li><li>'LMO7'</li><li>'HLA-DPA1'</li><li>'SQSTM1'</li><li>'NAMPT'</li><li>'CXCR4'</li><li>'GJA1'</li><li>'HSPG2'</li><li>'SFTPD'</li><li>'PITPNA'</li><li>'SLC11A1'</li><li>'IL4R'</li><li>'GPRC5A'</li><li>'LBH'</li><li>'VWF'</li><li>'TPM1'</li><li>'MMP2'</li><li>'MYC'</li><li>'NUPR1'</li><li>'SMCHD1'</li><li>'LGALS1'</li><li>'ASAH1'</li><li>'ITGB2'</li><li>'AGER'</li><li>'LIF'</li><li>'PLXND1'</li><li>'ARPC1B'</li><li>'LAMA4'</li><li>'GABARAPL1'</li><li>'IL1R1'</li><li>'UPP1'</li><li>'IGFBP4'</li><li>'LUM'</li><li>'CD93'</li><li>'NFKB2'</li><li>'MMP14'</li><li>'PLAAT4'</li><li>'KLHL6'</li><li>'ZFP36L2'</li><li>'VAMP5'</li><li>'EHD2'</li><li>'EMILIN1'</li><li>'PTPRC'</li><li>'SQOR'</li><li>'PLEKHG2'</li><li>'NNMT'</li><li>'NBL1'</li><li>'SLC39A8'</li><li>'IL6'</li><li>'PLEC'</li><li>'COL4A2'</li><li>'LIMS2'</li><li>'ITGB6'</li><li>'ARF3'</li><li>'IRF1'</li><li>'C1S'</li><li>'ARL4C'</li><li>'GFPT2'</li><li>'EGR1'</li><li>'KDM6B'</li><li>'ZYX'</li><li>'LMCD1'</li><li>'CD37'</li><li>'ZFAND5'</li><li>'MUC1'</li><li>'SFTPC'</li><li>'HBEGF'</li><li>'TNC'</li><li>'S100A9'</li><li>'ERRFI1'</li><li>'PCOLCE'</li><li>'LRP1'</li><li>'ANKRD11'</li><li>'BST2'</li><li>'FAM43A'</li><li>'ARHGAP45'</li><li>'ARHGDIB'</li><li>'KRT7'</li><li>'SH3BP5'</li><li>'PRRX1'</li><li>'ADAMTS1'</li><li>'FNBP1'</li><li>'ADAMTS9'</li><li>⋯</li><li>'CTIF'</li><li>'PER1'</li><li>'SOGA1'</li><li>'PTPRT'</li><li>'CAPN3'</li><li>'EDNRA'</li><li>'CISH'</li><li>'RNF19A'</li><li>'CPXM2'</li><li>'IGHG4'</li><li>'TCF4'</li><li>'FBXL7'</li><li>'IFRD2'</li><li>'SFRP4'</li><li>'SNX21'</li><li>'CDYL2'</li><li>'CCL4'</li><li>'FAM53B'</li><li>'NDN'</li><li>'PPM1F'</li><li>'PLD2'</li><li>'SYNPO2'</li><li>'COL14A1'</li><li>'MAGI1'</li><li>'NFKBIE'</li><li>'COL12A1'</li><li>'NEXN'</li><li>'MADD'</li><li>'TMEM45A'</li><li>'APOC2'</li><li>'SEMA5A'</li><li>'MYRF'</li><li>'SNX30'</li><li>'TBC1D2'</li><li>'TK2'</li><li>'ADRB2'</li><li>'ANXA3'</li><li>'CYGB'</li><li>'ANKRD44'</li><li>'DLL4'</li><li>'PRR5L'</li><li>'MYO5A'</li><li>'MAFK'</li><li>'A4GALT'</li><li>'ARHGEF3'</li><li>'IGFLR1'</li><li>'CLIP4'</li><li>'PROCR'</li><li>'SDR16C5'</li><li>'ASPN'</li><li>'NCF4'</li><li>'TRIM21'</li><li>'PPRC1'</li><li>'VAV1'</li><li>'DHX58'</li><li>'ABI3BP'</li><li>'STIMATE'</li><li>'GPR107'</li><li>'PPARG'</li><li>'VASN'</li><li>'ZHX2'</li><li>'TRANK1'</li><li>'DCUN1D3'</li><li>'OCEL1'</li><li>'CST7'</li><li>'BNC2'</li><li>'SLC27A1'</li><li>'SH2B3'</li><li>'AP5Z1'</li><li>'ADCY4'</li><li>'DLG4'</li><li>'AVEN'</li><li>'SEMA3B'</li><li>'HSBP1L1'</li><li>'ITPR1'</li><li>'PODXL'</li><li>'RAPGEF5'</li><li>'CPEB2'</li><li>'TSPAN7'</li><li>'LZTS1'</li><li>'SHROOM2'</li><li>'TNFSF13B'</li><li>'UNC13D'</li><li>'SASH1'</li><li>'CCL19'</li><li>'LCK'</li><li>'GPR68'</li><li>'ABCA2'</li><li>'TRAF3IP3'</li><li>'MID1IP1'</li><li>'LPL'</li><li>'PALM'</li><li>'ACP2'</li><li>'APOBR'</li><li>'RELB'</li><li>'HECTD2'</li><li>'HIST1H2BM'</li><li>'GPSM1'</li><li>'FAT4'</li><li>'PARP10'</li><li>'PTGER4'</li><li>'ME3'</li><li>'ZFYVE28'</li><li>'FBXO44'</li><li>'BIN2'</li><li>'CAPS'</li><li>'SLC35F2'</li><li>'PRKCE'</li><li>'SLC16A3'</li><li>'MX2'</li><li>'IL15RA'</li><li>'RILP'</li><li>'GALNS'</li><li>'LOXL2'</li><li>'PLAT'</li><li>'IRAK3'</li><li>'ST3GAL1'</li><li>'ST3GAL2'</li><li>'LYN'</li><li>'ITGB4'</li><li>'KCNJ8'</li><li>'ST6GALNAC4'</li><li>'PROS1'</li><li>'SOCS2'</li><li>'RELT'</li><li>'IGLV6-57'</li><li>'MACROD2'</li><li>'LACTB'</li><li>'RABGEF1'</li><li>'GATA6'</li><li>'MAST3'</li><li>'ADAM12'</li><li>'ABCC4'</li><li>'GPR157'</li><li>'ANKDD1A'</li><li>'MAPKAPK3'</li><li>'PTPN13'</li><li>'SLAMF7'</li><li>'MCEE'</li><li>'TXNDC17'</li><li>'SMAD6'</li><li>'NMRK1'</li><li>'GUCY1A2'</li><li>'DUOX1'</li><li>'S1PR1'</li><li>'GALM'</li><li>'KITLG'</li><li>'ZSWIM4'</li><li>'JAM2'</li><li>'CD109'</li><li>'C8orf58'</li><li>'PTPRU'</li><li>'FCHSD1'</li><li>'HELZ2'</li><li>'REL'</li><li>'GAB3'</li><li>'DOK3'</li><li>'DGKQ'</li><li>'DACH1'</li><li>'SPRY2'</li><li>'MXD1'</li><li>'ABCA7'</li><li>'FAM241A'</li><li>'TLR4'</li><li>'VSTM4'</li><li>'ARAP3'</li><li>'STARD13'</li><li>'GLYCTK'</li><li>'CHST2'</li><li>'SLC17A9'</li><li>'ESM1'</li><li>'LAMB3'</li><li>'GIMAP8'</li><li>'FADS3'</li><li>'PITPNM1'</li><li>'ANTXR2'</li><li>'ARRDC2'</li><li>'VILL'</li><li>'SOX5'</li><li>'COPZ2'</li><li>'TRIP10'</li><li>'PLIN3'</li><li>'DUOXA1'</li><li>'S1PR3'</li><li>'LIX1L'</li><li>'DACT1'</li><li>'GOLGA6L9'</li><li>'HERC6'</li><li>'MET'</li><li>'DNAH1'</li><li>'SUSD1'</li><li>'OLFML3'</li><li>'ERMAP'</li><li>'IL17RA'</li><li>'LUZP1'</li><li>'DOCK6'</li><li>'STN1'</li><li>'ADD3'</li><li>'GSTM1'</li><li>'SH3PXD2B'</li></ol>
</dd>
	<dt>$DOWN</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'PPDPF'</li><li>'PCGF2'</li><li>'ERBB2'</li><li>'GRB7'</li><li>'PGAP3'</li><li>'CDK12'</li><li>'PSMB3'</li><li>'MLLT6'</li><li>'NUCKS1'</li><li>'SRCIN1'</li><li>'SYNE2'</li><li>'STARD3'</li><li>'HIST1H2BH'</li><li>'PIP4K2B'</li><li>'KRT15'</li><li>'HIST2H2AB'</li><li>'ALCAM'</li><li>'SOCS7'</li><li>'MIEN1'</li><li>'STMN1'</li><li>'ERV3-1'</li><li>'PCBP2'</li><li>'NKX2-1'</li><li>'TOP2A'</li><li>'LASP1'</li><li>'MDK'</li><li>'LMO3'</li><li>'H1F0'</li><li>'CKB'</li><li>'PAPOLA'</li><li>'RERGL'</li><li>'CXCL17'</li><li>'TUBA1B'</li><li>'AP1G2'</li><li>'ZNF117'</li><li>'PIK3C2G'</li><li>'HIST1H3H'</li><li>'HLA-A'</li><li>'SPP1'</li><li>'FOXA1'</li><li>'SMAD9'</li><li>'LSM14A'</li><li>'SFTA3'</li><li>'NFIB'</li><li>'FXR1'</li><li>'HDAC2'</li><li>'FAM111B'</li><li>'CACNB1'</li><li>'HIST1H2BL'</li><li>'GRIK3'</li><li>'HIST1H2BG'</li><li>'TYMS'</li><li>'UBE2C'</li><li>'TMEM123'</li><li>'MED1'</li><li>'SMC4'</li><li>'NSD2'</li><li>'EIF2A'</li><li>'BDH1'</li><li>'SP3'</li><li>'MLPH'</li><li>'CISD3'</li><li>'MDM4'</li><li>'ZNF146'</li><li>'CBX3'</li><li>'SLITRK6'</li><li>'HIST1H2AI'</li><li>'TOX3'</li><li>'PCP4'</li><li>'NFATC4'</li><li>'RNPS1'</li><li>'WDR72'</li><li>'ZNF507'</li><li>'SDK1'</li><li>'CEL'</li><li>'ARHGAP33'</li><li>'LINC00672'</li><li>'ACIN1'</li><li>'RIMKLA'</li><li>'ECT2'</li><li>'NDC80'</li><li>'SLC25A36'</li><li>'RGS17'</li><li>'HOXB3'</li><li>'CWC25'</li><li>'NSD3'</li><li>'RBP1'</li><li>'UMPS'</li><li>'NUDT21'</li><li>'MBOAT1'</li><li>'GRHL2'</li><li>'ASS1'</li><li>'CYP4X1'</li><li>'PPP1R1B'</li><li>'SYDE2'</li><li>'GSK3B'</li><li>'YEATS2'</li><li>'HELLS'</li><li>'IPO9'</li><li>'SOX2'</li><li>'GUCA1B'</li><li>'SRP9'</li><li>'SIX1'</li><li>'KHDC4'</li><li>'CHTOP'</li><li>'E2F1'</li><li>'GPBP1L1'</li><li>'TFDP2'</li><li>'PLXDC1'</li><li>'C19orf48'</li><li>'SAV1'</li><li>'PVALB'</li><li>'RHPN2'</li><li>'DNAH14'</li><li>'POF1B'</li><li>'B3GNT5'</li><li>'KLK11'</li><li>'CEBPG'</li><li>'TCAP'</li><li>'STIL'</li><li>'PERP'</li><li>'TRIM33'</li><li>'EXOC5'</li><li>'TTLL5'</li><li>'DLG1'</li><li>'TTK'</li><li>'ZBTB21'</li><li>'NET1'</li><li>'DHX36'</li><li>'STAG1'</li><li>'CIP2A'</li><li>'ARHGAP23'</li><li>'ATG3'</li><li>'RAD51B'</li><li>'POLE2'</li><li>'CACNA2D1'</li><li>'EIF1AX'</li><li>'RAB2B'</li><li>'TXNDC16'</li><li>'CDCA7'</li><li>'ATG12'</li><li>'CRABP2'</li><li>'CENPK'</li><li>'ZNF260'</li><li>'SNX4'</li><li>'RFC4'</li><li>'NAA50'</li><li>'GABRE'</li><li>'FBXO16'</li><li>'ZNF793'</li><li>'KIF18B'</li><li>'TTC6'</li><li>'SDR39U1'</li><li>'SEPTIN9'</li><li>'CCNC'</li><li>'POMT2'</li><li>'MED23'</li><li>'VANGL2'</li><li>'HOXB2'</li><li>'NSL1'</li><li>'SPTSSA'</li><li>'PELI2'</li><li>'U2SURP'</li><li>'EPHB1'</li><li>'ALG6'</li><li>'RECQL4'</li><li>'G2E3'</li><li>'DTL'</li><li>'R3HDM2'</li><li>'CCNE1'</li><li>'SIPA1L3'</li><li>'NLRP11'</li><li>'PLEKHA6'</li><li>'CP'</li><li>'SLC25A11'</li><li>'TRIP11'</li><li>'POLQ'</li><li>'PIAS3'</li><li>'ZBTB41'</li><li>'DYNLT1'</li><li>'SMG5'</li><li>'CNTNAP2'</li><li>'CASP8AP2'</li><li>'YLPM1'</li><li>'DNAJC19'</li><li>'WDR62'</li><li>'ZMAT4'</li><li>'TIMMDC1'</li><li>'KLHDC2'</li><li>'MFSD6'</li><li>'PHC3'</li><li>'SPAG5'</li><li>'BTF3L4'</li><li>'SPTB'</li><li>'MLF1'</li><li>'TBL1XR1'</li><li>'CRABP1'</li><li>'C1QL1'</li><li>'MIS18BP1'</li><li>'CDK19'</li><li>⋯</li><li>'USP48'</li><li>'ZNF354A'</li><li>'TRMT10C'</li><li>'REC8'</li><li>'HMCES'</li><li>'ZBTB5'</li><li>'SLC22A17'</li><li>'SLC4A8'</li><li>'SPEG'</li><li>'CRIPT'</li><li>'ADAM20'</li><li>'SAMD10'</li><li>'RFWD3'</li><li>'KLHL36'</li><li>'MARCH7'</li><li>'QRICH2'</li><li>'ATR'</li><li>'MMACHC'</li><li>'ZNF696'</li><li>'PLGRKT'</li><li>'IGSF9B'</li><li>'TMEM63C'</li><li>'AC015802.6'</li><li>'BGLAP'</li><li>'RANBP6'</li><li>'CDC14A'</li><li>'WDR24'</li><li>'HIF1AN'</li><li>'UNC5CL'</li><li>'TRMT13'</li><li>'GTF2H3'</li><li>'USP37'</li><li>'MEGF8'</li><li>'C5orf51'</li><li>'TC2N'</li><li>'MARC1'</li><li>'FBXO30'</li><li>'MPP2'</li><li>'KIF15'</li><li>'ARHGAP5'</li><li>'TMEM108'</li><li>'FRS2'</li><li>'ZNF214'</li><li>'C1orf112'</li><li>'NEK8'</li><li>'KLHL28'</li><li>'NDRG3'</li><li>'PGAP1'</li><li>'RAB33B'</li><li>'FAM160A1'</li><li>'ZNF713'</li><li>'NHS'</li><li>'BRCA2'</li><li>'DNASE1'</li><li>'EID2B'</li><li>'ZNF526'</li><li>'ST3GAL3'</li><li>'UBALD1'</li><li>'ZNF566'</li><li>'ZYG11B'</li><li>'CLDN11'</li><li>'GAS2L3'</li><li>'ZFHX2'</li><li>'NBPF20'</li><li>'HPN'</li><li>'GORAB'</li><li>'SRD5A3'</li><li>'USF3'</li><li>'PDK2'</li><li>'ZNF284'</li><li>'NAV2'</li><li>'NSUN4'</li><li>'TEAD2'</li><li>'WDR35'</li><li>'PPP1R1A'</li><li>'RB1CC1'</li><li>'OFD1'</li><li>'ERCC4'</li><li>'MBTPS2'</li><li>'PDK3'</li><li>'OXCT2'</li><li>'ZNF607'</li><li>'SPATA7'</li><li>'ALG1'</li><li>'TSEN15'</li><li>'ZBTB42'</li><li>'SLC30A5'</li><li>'KCTD16'</li><li>'DMKN'</li><li>'KIF2A'</li><li>'NUDT17'</li><li>'IGSF9'</li><li>'SLCO5A1'</li><li>'MTRF1'</li><li>'TOMM20'</li><li>'SLC16A10'</li><li>'ANGEL2'</li><li>'PPCS'</li><li>'TENT2'</li><li>'IRX3'</li><li>'NTF3'</li><li>'PABPC4L'</li><li>'ZFP3'</li><li>'ZNF570'</li><li>'ZNF34'</li><li>'SLF2'</li><li>'CASTOR3'</li><li>'HEXD'</li><li>'AREL1'</li><li>'DZIP1L'</li><li>'PGM2L1'</li><li>'PPP1R8'</li><li>'BMT2'</li><li>'EXOSC3'</li><li>'BAZ1A'</li><li>'ZNF514'</li><li>'SNRPB2'</li><li>'SCD5'</li><li>'DRD4'</li><li>'SGO1'</li><li>'CEACAM19'</li><li>'DUS4L'</li><li>'ZNF280C'</li><li>'ZNF181'</li><li>'KCNK1'</li><li>'TCEAL2'</li><li>'HIF3A'</li><li>'PTPRK'</li><li>'COA7'</li><li>'RASL10B'</li><li>'HENMT1'</li><li>'MPPE1'</li><li>'TAS2R5'</li><li>'FOXC1'</li><li>'MAP3K10'</li><li>'GPR161'</li><li>'RIDA'</li><li>'NSG1'</li><li>'MINPP1'</li><li>'FAIM'</li><li>'DAAM1'</li><li>'GPRASP1'</li><li>'PITPNM3'</li><li>'TAF3'</li><li>'USP10'</li><li>'TMEM54'</li><li>'SCNN1B'</li><li>'AGO4'</li><li>'STYX'</li><li>'ZNF135'</li><li>'NUS1'</li><li>'FBXO27'</li><li>'CLBA1'</li><li>'AMOTL2'</li><li>'BLOC1S4'</li><li>'EFNB1'</li><li>'DTX4'</li><li>'GABPB2'</li><li>'SLC38A9'</li><li>'NOL11'</li><li>'TFB2M'</li><li>'TMEM143'</li><li>'ERCC8'</li><li>'SRBD1'</li><li>'THBS3'</li><li>'ZNF461'</li><li>'MTO1'</li><li>'CASTOR2'</li><li>'RHBDF1'</li><li>'VBP1'</li><li>'MED7'</li><li>'ATAD5'</li><li>'CIDEB'</li><li>'TUB'</li><li>'NSMCE1'</li><li>'ZBTB49'</li><li>'SKA2'</li><li>'FAM102B'</li><li>'NSMCE2'</li><li>'TAF1B'</li><li>'ZNF790'</li><li>'RUBCN'</li><li>'ZNF90'</li><li>'ZNF302'</li><li>'FANCM'</li><li>'ANKRD12'</li><li>'UHRF1BP1L'</li><li>'ZNF595'</li><li>'TINAGL1'</li><li>'EFS'</li><li>'UBE2W'</li><li>'HOMER1'</li><li>'GRHL1'</li><li>'SULT1A1'</li><li>'BCL9'</li><li>'ZNF462'</li><li>'TFAP2C'</li><li>'NIFK'</li><li>'TRPV1'</li><li>'LLGL1'</li></ol>
</dd>
</dl>



# Subset data to spots not in Normal Lung or Fibroblast clusters (since there are lots of Normal Lung spots with vasc_high)


```R
print(unique(data@meta.data$Navin_Clusters2))
```

    [1] 7_10_Normal_Lung       6_SCLC_Effector_Immune 1_4_5_11_SCLC         
    [4] 2_SCLC_TSI             3_SCLC_Stroma          0_Lymphocytes_B_Cells 
    [7] 8_Fibroblasts          9_SCLC_Apoptotic      
    8 Levels: 0_Lymphocytes_B_Cells 1_4_5_11_SCLC 2_SCLC_TSI ... 9_SCLC_Apoptotic



```R
temp_data <- data
# remove pesky graphs
temp_data@graphs <- list()
data_not_normal <- subset(temp_data, subset=Navin_Clusters2!='7_10_Normal_Lung')
data_not_normal <- subset(data_not_normal, subset=Navin_Clusters2!='8_Fibroblasts')
print(unique(data_not_normal@meta.data$Navin_Clusters2))
data_not_normal
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


    [1] 6_SCLC_Effector_Immune 1_4_5_11_SCLC          2_SCLC_TSI            
    [4] 3_SCLC_Stroma          0_Lymphocytes_B_Cells  9_SCLC_Apoptotic      
    6 Levels: 0_Lymphocytes_B_Cells 1_4_5_11_SCLC 2_SCLC_TSI ... 9_SCLC_Apoptotic



    An object of class Seurat 
    35796 features across 8170 samples within 2 assays 
    Active assay: Spatial (17943 features, 2000 variable features)
     3 layers present: data, counts, scale.data
     1 other assay present: SCT
     3 dimensional reductions calculated: pca, integrated.cca, umap
     3 spatial fields of view present: slice1 slice1.2 slice1.3



```R
unique(data_not_normal@meta.data$vasc_MN4_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_low'</li><li>'vasc_high_MN4_high'</li><li>'vasc_high_MN4_low'</li></ol>




```R
Idents(data_not_normal) <- "vasc_MN4_label"

de_vasc_MN4_not_normal <- FindMarkers(
  object = data_not_normal,
  ident.1 = "vasc_high_MN4_high",  # First group to compare
  ident.2 = "vasc_high_MN4_low",    # Second group
  logfc.threshold = 0.25, # Minimum log2 fold change
  min.pct = 0.1,          # Gene detected in ≥10% of cells in either group
  test.use = "wilcox"     # Default Wilcoxon rank sum test
)

de_vasc_MN4_not_normal$Gene <- rownames(de_vasc_MN4_not_normal)
```


```R
de_vasc_MN4_not_normal %>% head(15)
```


<table class="dataframe">
<caption>A data.frame: 15 × 6</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>Gene</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>LYZ</th><td>9.023810e-69</td><td>1.3349424</td><td>0.895</td><td>0.679</td><td>1.619142e-64</td><td>LYZ   </td></tr>
	<tr><th scope=row>SPARC</th><td>8.976468e-65</td><td>1.1959381</td><td>0.883</td><td>0.650</td><td>1.610648e-60</td><td>SPARC </td></tr>
	<tr><th scope=row>COL1A2</th><td>1.509235e-60</td><td>1.3398002</td><td>0.902</td><td>0.753</td><td>2.708020e-56</td><td>COL1A2</td></tr>
	<tr><th scope=row>COL3A1</th><td>8.160517e-57</td><td>1.3482867</td><td>0.868</td><td>0.711</td><td>1.464242e-52</td><td>COL3A1</td></tr>
	<tr><th scope=row>FN1</th><td>6.245171e-56</td><td>0.9388052</td><td>0.982</td><td>0.886</td><td>1.120571e-51</td><td>FN1   </td></tr>
	<tr><th scope=row>S100A6</th><td>2.560301e-53</td><td>0.6524127</td><td>0.987</td><td>0.904</td><td>4.593948e-49</td><td>S100A6</td></tr>
	<tr><th scope=row>IGHG1</th><td>9.068149e-53</td><td>1.5647181</td><td>0.971</td><td>0.915</td><td>1.627098e-48</td><td>IGHG1 </td></tr>
	<tr><th scope=row>IGKC</th><td>3.354279e-52</td><td>1.4379913</td><td>1.000</td><td>0.996</td><td>6.018584e-48</td><td>IGKC  </td></tr>
	<tr><th scope=row>IGFBP7</th><td>3.644504e-52</td><td>1.2352471</td><td>0.757</td><td>0.505</td><td>6.539334e-48</td><td>IGFBP7</td></tr>
	<tr><th scope=row>B2M</th><td>5.373119e-52</td><td>0.4107423</td><td>0.997</td><td>0.993</td><td>9.640987e-48</td><td>B2M   </td></tr>
	<tr><th scope=row>PSAP</th><td>1.209349e-51</td><td>0.9757966</td><td>0.971</td><td>0.886</td><td>2.169934e-47</td><td>PSAP  </td></tr>
	<tr><th scope=row>CD74</th><td>1.525058e-49</td><td>0.6755053</td><td>0.996</td><td>0.965</td><td>2.736411e-45</td><td>CD74  </td></tr>
	<tr><th scope=row>CCN1</th><td>3.794567e-47</td><td>1.4428265</td><td>0.662</td><td>0.399</td><td>6.808591e-43</td><td>CCN1  </td></tr>
	<tr><th scope=row>TIMP1</th><td>6.504763e-47</td><td>0.8437818</td><td>0.926</td><td>0.806</td><td>1.167150e-42</td><td>TIMP1 </td></tr>
	<tr><th scope=row>ICAM1</th><td>9.666440e-47</td><td>1.9441492</td><td>0.576</td><td>0.321</td><td>1.734449e-42</td><td>ICAM1 </td></tr>
</tbody>
</table>




```R
# check significant genes
de_vasc_MN4_not_normal_sig_label <- label_sig_genes(de_vasc_MN4_not_normal)
de_vasc_MN4_not_normal_sig_DE <- de_vasc_MN4_not_normal_sig_label[de_vasc_MN4_not_normal_sig_label$sig_DE != 'NO',]
print(dim(de_vasc_MN4_not_normal))
print(dim(de_vasc_MN4_not_normal_sig_DE[]))
de_vasc_MN4_not_normal_sig_DE %>% head(5)
table(de_vasc_MN4_not_normal_sig_DE$sig_DE)
```

    [1] 5179    6
    [1] 1228    7



<table class="dataframe">
<caption>A data.frame: 5 × 7</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>Gene</th><th scope=col>sig_DE</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>LYZ</th><td>9.023810e-69</td><td>1.3349424</td><td>0.895</td><td>0.679</td><td>1.619142e-64</td><td>LYZ   </td><td>UP</td></tr>
	<tr><th scope=row>SPARC</th><td>8.976468e-65</td><td>1.1959381</td><td>0.883</td><td>0.650</td><td>1.610648e-60</td><td>SPARC </td><td>UP</td></tr>
	<tr><th scope=row>COL1A2</th><td>1.509235e-60</td><td>1.3398002</td><td>0.902</td><td>0.753</td><td>2.708020e-56</td><td>COL1A2</td><td>UP</td></tr>
	<tr><th scope=row>COL3A1</th><td>8.160517e-57</td><td>1.3482867</td><td>0.868</td><td>0.711</td><td>1.464242e-52</td><td>COL3A1</td><td>UP</td></tr>
	<tr><th scope=row>FN1</th><td>6.245171e-56</td><td>0.9388052</td><td>0.982</td><td>0.886</td><td>1.120571e-51</td><td>FN1   </td><td>UP</td></tr>
</tbody>
</table>




    
    DOWN   UP 
     328  900 



```R
de_vasc_MN4_not_normal_check_goi <- de_vasc_MN4_not_normal %>% rownames_to_column('gene') %>% subset(gene %in% c('ENG', 'VSIR', 'FASLG'))
rownames(de_vasc_MN4_not_normal_check_goi) <- de_vasc_MN4_not_normal_check_goi$gene
de_vasc_MN4_not_normal_check_goi
```


<table class="dataframe">
<caption>A data.frame: 2 × 7</caption>
<thead>
	<tr><th></th><th scope=col>gene</th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>Gene</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>ENG</th><td>ENG </td><td>7.195647e-15</td><td>0.6336688</td><td>0.448</td><td>0.290</td><td>1.291115e-10</td><td>ENG </td></tr>
	<tr><th scope=row>VSIR</th><td>VSIR</td><td>5.850418e-10</td><td>0.8136430</td><td>0.209</td><td>0.105</td><td>1.049741e-05</td><td>VSIR</td></tr>
</tbody>
</table>




```R
options(repr.plot.width=10, repr.plot.height=10)
de_plot <- EnhancedVolcano(de_vasc_MN4_not_normal,
                lab = rownames(de_vasc_MN4_not_normal),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'MN4_high vs MN4_low, in vasc_high spots, excluding normal spots')
de_plot

# save plot
png(file = "Output/07_DE_vasc_MN4_and_Tumor_MHCI/de_vasc_high_MN4_high_vs_low_not_normal_spots_volcano.png",
    width = 1000,
    height = 1000)
print(de_plot)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_39_1.png)
    



```R
options(repr.plot.width=10, repr.plot.height=12)
Idents(data) <- "vasc_MN4_label"
FeaturePlot(data, features = rownames(de_vasc_MN4_not_normal)[1:6], reduction='umap')
```


    
![png](output_40_0.png)
    



```R
options(repr.plot.width=10, repr.plot.height=12)
Idents(data) <- "vasc_MN4_label"
VlnPlot(data, features = rownames(de_vasc_MN4_not_normal)[1:6], split.by = 'vasc_MN4_label')
```


    
![png](output_41_0.png)
    



```R
# Extract lists of genes UP and DOWN regulated
de_vasc_MN4_not_normal_sig_DE_list <- list()
de_vasc_MN4_not_normal_sig_DE_list$UP <- rownames(de_vasc_MN4_not_normal_sig_DE[de_vasc_MN4_not_normal_sig_DE$sig_DE=="UP",])
de_vasc_MN4_not_normal_sig_DE_list$DOWN <- rownames(de_vasc_MN4_not_normal_sig_DE[de_vasc_MN4_not_normal_sig_DE$sig_DE=="DOWN",])
```


```R
# check list
de_vasc_MN4_not_normal_sig_DE_list
```


<dl>
	<dt>$UP</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'LYZ'</li><li>'SPARC'</li><li>'COL1A2'</li><li>'COL3A1'</li><li>'FN1'</li><li>'S100A6'</li><li>'IGHG1'</li><li>'IGKC'</li><li>'IGFBP7'</li><li>'PSAP'</li><li>'CD74'</li><li>'CCN1'</li><li>'TIMP1'</li><li>'ICAM1'</li><li>'HLA-E'</li><li>'DCN'</li><li>'COL6A3'</li><li>'COL4A1'</li><li>'HLA-DMA'</li><li>'THBS1'</li><li>'IL32'</li><li>'BGN'</li><li>'C1R'</li><li>'SLC34A2'</li><li>'ACTA2'</li><li>'S100A10'</li><li>'GADD45B'</li><li>'CTSD'</li><li>'APOE'</li><li>'IFI6'</li><li>'SFTPB'</li><li>'JCHAIN'</li><li>'TPM2'</li><li>'SERPINE1'</li><li>'CCL18'</li><li>'TIMP3'</li><li>'TYROBP'</li><li>'CAVIN1'</li><li>'HLA-DOA'</li><li>'PTGDS'</li><li>'C3'</li><li>'IGHG3'</li><li>'COL1A1'</li><li>'LAPTM5'</li><li>'IGHJ6'</li><li>'JUNB'</li><li>'SOCS3'</li><li>'AQP1'</li><li>'PLAUR'</li><li>'FLNA'</li><li>'POSTN'</li><li>'ITGAX'</li><li>'ELL2'</li><li>'C11orf96'</li><li>'FBN1'</li><li>'CIITA'</li><li>'CD44'</li><li>'LPCAT1'</li><li>'TRIM22'</li><li>'LSP1'</li><li>'EMP1'</li><li>'RNF213'</li><li>'MCL1'</li><li>'G0S2'</li><li>'CCN2'</li><li>'MYL9'</li><li>'LUM'</li><li>'FOSB'</li><li>'PLAAT4'</li><li>'CTSZ'</li><li>'SAT1'</li><li>'S100A4'</li><li>'MMP2'</li><li>'HLA-DPB1'</li><li>'LMNA'</li><li>'SOD2'</li><li>'EMILIN1'</li><li>'HSPG2'</li><li>'ITGB1'</li><li>'RHOB'</li><li>'C1QC'</li><li>'CD37'</li><li>'GSN'</li><li>'OLR1'</li><li>'LBH'</li><li>'ITGB2'</li><li>'NABP1'</li><li>'IFI30'</li><li>'ZFP36'</li><li>'MZB1'</li><li>'IER3'</li><li>'IGKV4-1'</li><li>'KLHL6'</li><li>'IKZF1'</li><li>'TYMP'</li><li>'CEBPD'</li><li>'MMP14'</li><li>'CYBB'</li><li>'STOM'</li><li>'C4orf3'</li><li>'PTPRE'</li><li>'PFKFB3'</li><li>'SULF1'</li><li>'PPP1R14A'</li><li>'VAMP5'</li><li>'PDGFRB'</li><li>'IGFBP4'</li><li>'NCOA7'</li><li>'DUSP1'</li><li>'CRISPLD2'</li><li>'TCIM'</li><li>'IGHM'</li><li>'S100A9'</li><li>'ARHGAP45'</li><li>'CCL2'</li><li>'IRF1'</li><li>'BST2'</li><li>'LGMN'</li><li>'ADAMTS4'</li><li>'GBP5'</li><li>'SFTPD'</li><li>'RCAN1'</li><li>'C1S'</li><li>'CERCAM'</li><li>'FCGRT'</li><li>'FCGR3A'</li><li>'HSPA5'</li><li>'TRAC'</li><li>'FNBP1'</li><li>'GMFG'</li><li>'SLC2A3'</li><li>'ATF3'</li><li>'IL4R'</li><li>'LRP1'</li><li>'DDIT4'</li><li>'COL5A1'</li><li>'NBL1'</li><li>'ZFP36L2'</li><li>'SAMD9L'</li><li>'CFD'</li><li>'NR4A1'</li><li>'PRRX1'</li><li>'SLC11A1'</li><li>'PPP1R15A'</li><li>'MXRA7'</li><li>'ARRDC3'</li><li>'IGLV3-1'</li><li>'PLXND1'</li><li>'CDH11'</li><li>'EPSTI1'</li><li>'CTSS'</li><li>'CTSC'</li><li>'SGK1'</li><li>'TRBC1'</li><li>'FOS'</li><li>'COTL1'</li><li>'PALM2-AKAP2'</li><li>'BHLHE40'</li><li>'PTPRC'</li><li>'RARRES2'</li><li>'TGFBR2'</li><li>'PITPNA'</li><li>'IL1R1'</li><li>'PDLIM3'</li><li>'SQOR'</li><li>'FAM49A'</li><li>'APOL1'</li><li>'IRF8'</li><li>'ITGB6'</li><li>'TMEM140'</li><li>'EHD2'</li><li>'GFPT2'</li><li>'CLEC14A'</li><li>'ZCCHC24'</li><li>'IRF4'</li><li>'KRT7'</li><li>'THBS2'</li><li>'SPI1'</li><li>'EMILIN2'</li><li>'CD93'</li><li>'TNC'</li><li>'IL7R'</li><li>'COL4A2'</li><li>'NR4A3'</li><li>'NFKB2'</li><li>'FGL2'</li><li>'WIPF1'</li><li>'PSMB10'</li><li>'PLEKHG2'</li><li>'CCDC80'</li><li>'XAF1'</li><li>'PIM2'</li><li>'NAPSA'</li><li>'ASAH1'</li><li>'CXCR4'</li><li>'COL16A1'</li><li>'IGHG2'</li><li>'EGR1'</li><li>'PGGHG'</li><li>'CNN2'</li><li>⋯</li><li>'VAV1'</li><li>'HIST1H3G'</li><li>'SAMD4A'</li><li>'C5AR1'</li><li>'ZC3H12A'</li><li>'ADGRG6'</li><li>'SH3TC1'</li><li>'CCL4'</li><li>'CMKLR1'</li><li>'PLCB2'</li><li>'ADAM12'</li><li>'ABI3'</li><li>'PPP3CC'</li><li>'CHI3L1'</li><li>'LTF'</li><li>'HMCN1'</li><li>'EDNRA'</li><li>'TNFAIP3'</li><li>'LMO2'</li><li>'ITGB3'</li><li>'MADD'</li><li>'HCST'</li><li>'CPXM2'</li><li>'EML3'</li><li>'GZMA'</li><li>'KSR1'</li><li>'TNFSF13'</li><li>'ITPKC'</li><li>'GNG11'</li><li>'CALCRL'</li><li>'LIPG'</li><li>'PXN'</li><li>'TRADD'</li><li>'MMP11'</li><li>'DNAJC1'</li><li>'APOC2'</li><li>'RABGEF1'</li><li>'ITPR1'</li><li>'PCDH17'</li><li>'CCL19'</li><li>'RELB'</li><li>'SHFL'</li><li>'FLI1'</li><li>'PDE4D'</li><li>'TRAF3IP3'</li><li>'WWC2'</li><li>'DDX3Y'</li><li>'NDNF'</li><li>'RFLNB'</li><li>'CD69'</li><li>'TBC1D2'</li><li>'ARHGEF3'</li><li>'TNFSF13B'</li><li>'CPM'</li><li>'KCNJ8'</li><li>'PRKCH'</li><li>'SPRY4'</li><li>'ARHGEF17'</li><li>'LY86'</li><li>'HECTD2'</li><li>'EVA1B'</li><li>'RNF19A'</li><li>'APBB3'</li><li>'ABLIM3'</li><li>'PLA2G2D'</li><li>'OCEL1'</li><li>'SLC35F2'</li><li>'PDPN'</li><li>'TMEM45A'</li><li>'ANKRD44'</li><li>'DENND5B'</li><li>'ABCA2'</li><li>'UNC13D'</li><li>'FBXO44'</li><li>'NFKB1'</li><li>'PPM1F'</li><li>'SLIT3'</li><li>'BNC2'</li><li>'STX12'</li><li>'ELL'</li><li>'IFFO1'</li><li>'ALDH1B1'</li><li>'SIGLEC1'</li><li>'NFAM1'</li><li>'PTPRT'</li><li>'SVEP1'</li><li>'PPRC1'</li><li>'TPGS1'</li><li>'ERG'</li><li>'SLC27A1'</li><li>'LZTS1'</li><li>'OSBPL10'</li><li>'PLD2'</li><li>'ING3'</li><li>'NFKBIE'</li><li>'CCDC68'</li><li>'RENBP'</li><li>'FKBP14'</li><li>'CAPN3'</li><li>'TRANK1'</li><li>'NEXN'</li><li>'AKAP12'</li><li>'RAPGEF5'</li><li>'GALM'</li><li>'AP5B1'</li><li>'MSC'</li><li>'CST7'</li><li>'LCK'</li><li>'SNX21'</li><li>'SHROOM2'</li><li>'SGPP1'</li><li>'NFIL3'</li><li>'DCUN1D3'</li><li>'MAPKAPK3'</li><li>'CISH'</li><li>'RIPOR2'</li><li>'MX2'</li><li>'PLCL2'</li><li>'CD34'</li><li>'DOK3'</li><li>'SYNE1'</li><li>'PODXL'</li><li>'NLRP1'</li><li>'AQP4'</li><li>'BIN2'</li><li>'ARSB'</li><li>'RARA'</li><li>'PDGFB'</li><li>'ADCY4'</li><li>'PPARG'</li><li>'GSTM1'</li><li>'ANTXR2'</li><li>'RILP'</li><li>'ANPEP'</li><li>'LACTB'</li><li>'MAST3'</li><li>'ST6GALNAC4'</li><li>'C1orf116'</li><li>'CPEB2'</li><li>'TLR4'</li><li>'DIPK1A'</li><li>'ABCA7'</li><li>'ABCC4'</li><li>'HSBP1L1'</li><li>'SEMA5A'</li><li>'MYRF'</li><li>'LIX1L'</li><li>'PLAT'</li><li>'DGKQ'</li><li>'ST3GAL2'</li><li>'PTGER4'</li><li>'SLC17A9'</li><li>'VSTM4'</li><li>'ST3GAL1'</li><li>'CCL14'</li><li>'DGKD'</li><li>'PLOD1'</li><li>'DACT1'</li><li>'S1PR1'</li><li>'PTPRU'</li><li>'PLIN3'</li><li>'GAB3'</li><li>'CYB5R2'</li><li>'SLC16A3'</li><li>'PPM1K'</li><li>'CPNE2'</li><li>'HERC5'</li><li>'LRRC53'</li><li>'JAG1'</li><li>'DLL4'</li><li>'MAGI1'</li><li>'TPPP3'</li><li>'C8orf58'</li><li>'INAFM2'</li><li>'POU6F1'</li><li>'FAT4'</li><li>'SUSD1'</li><li>'ABCA3'</li><li>'BMP2'</li><li>'FCHSD1'</li><li>'EXD3'</li><li>'ZNF264'</li><li>'S100A2'</li><li>'LUZP1'</li><li>'IL17RA'</li><li>'ANKDD1A'</li><li>'DACH1'</li><li>'JAM2'</li><li>'GIMAP8'</li><li>'GATA6'</li><li>'GPR157'</li><li>'KCNK6'</li><li>'FAM89A'</li><li>'GIPR'</li><li>'SNCA'</li><li>'RAP1GAP'</li><li>'GPX7'</li><li>'NPIPB15'</li><li>'SH3D21'</li><li>'PLA1A'</li></ol>
</dd>
	<dt>$DOWN</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'PCGF2'</li><li>'GRB7'</li><li>'ERBB2'</li><li>'PGAP3'</li><li>'CDK12'</li><li>'MLLT6'</li><li>'PSMB3'</li><li>'STARD3'</li><li>'SRCIN1'</li><li>'KRT15'</li><li>'PIP4K2B'</li><li>'ALCAM'</li><li>'NKX2-1'</li><li>'LASP1'</li><li>'SOCS7'</li><li>'LMO3'</li><li>'CKB'</li><li>'SFTA3'</li><li>'ERV3-1'</li><li>'MLPH'</li><li>'CXCL17'</li><li>'AP1G2'</li><li>'LSM14A'</li><li>'RERGL'</li><li>'FOXA1'</li><li>'TUBA1B'</li><li>'CACNB1'</li><li>'CWC25'</li><li>'SDK1'</li><li>'GRIK3'</li><li>'BDH1'</li><li>'WDR72'</li><li>'YEATS2'</li><li>'LINC00672'</li><li>'FAM111B'</li><li>'SYDE2'</li><li>'ECT2'</li><li>'MBOAT1'</li><li>'ZBTB21'</li><li>'NDC80'</li><li>'RIMKLA'</li><li>'NFATC4'</li><li>'PCP4'</li><li>'ALG6'</li><li>'MED1'</li><li>'SLC25A36'</li><li>'CBX3'</li><li>'MUC4'</li><li>'PPP1R1B'</li><li>'UMPS'</li><li>'ZC2HC1C'</li><li>'TCAP'</li><li>'CEL'</li><li>'CYP4X1'</li><li>'PLXDC1'</li><li>'SMC4'</li><li>'IPO9'</li><li>'E2F1'</li><li>'VANGL2'</li><li>'ZNF507'</li><li>'EPHB1'</li><li>'INPP5J'</li><li>'TFDP2'</li><li>'FBXO16'</li><li>'RPA2'</li><li>'TXNDC16'</li><li>'TRIM33'</li><li>'RHPN2'</li><li>'GUCA1B'</li><li>'RGS17'</li><li>'SLC6A20'</li><li>'GRHL2'</li><li>'MUC16'</li><li>'KDR'</li><li>'SIX1'</li><li>'MFSD6'</li><li>'PLEKHG6'</li><li>'ARHGAP23'</li><li>'CLDN1'</li><li>'TTC6'</li><li>'STAG1'</li><li>'GABRE'</li><li>'ATG3'</li><li>'KIF18B'</li><li>'CDK19'</li><li>'RAB2B'</li><li>'RET'</li><li>'PELI2'</li><li>'NET1'</li><li>'SAV1'</li><li>'POF1B'</li><li>'MME'</li><li>'TTK'</li><li>'SOX2'</li><li>'PLEKHA6'</li><li>'ATG12'</li><li>'MYO18A'</li><li>'NOS1AP'</li><li>'RAD51B'</li><li>'MYB'</li><li>'TRIP11'</li><li>'CRABP1'</li><li>'ZRANB3'</li><li>'ZNF217'</li><li>'ZMAT4'</li><li>'ZNF292'</li><li>'TRMT61A'</li><li>'G2E3'</li><li>'MIS18BP1'</li><li>'C1D'</li><li>'PERP'</li><li>'CENPK'</li><li>'KIAA0232'</li><li>'CRABP2'</li><li>'DEGS2'</li><li>'ALS2'</li><li>'SIPA1L3'</li><li>'LONRF2'</li><li>'EPM2AIP1'</li><li>'PANX2'</li><li>'DDX11'</li><li>'CREB3L4'</li><li>'POMT2'</li><li>'C1QL1'</li><li>'STIL'</li><li>'ZBTB41'</li><li>'CP'</li><li>'NLRP11'</li><li>'KCNK5'</li><li>'IRX5'</li><li>'BRD8'</li><li>'MTRF1L'</li><li>'GJB1'</li><li>'RECQL4'</li><li>'GTF3C4'</li><li>'EPOP'</li><li>'ESRRG'</li><li>'WNT11'</li><li>'ZNF260'</li><li>'SYT14'</li><li>'CEP78'</li><li>'IRAK1BP1'</li><li>'ZNF875'</li><li>'NEURL1B'</li><li>'MED23'</li><li>'STYK1'</li><li>'ADCY5'</li><li>'KLHL8'</li><li>'TUT4'</li><li>'ACSF3'</li><li>'COQ2'</li><li>'STEAP3'</li><li>'AHCTF1'</li><li>'POLE2'</li><li>'RNF6'</li><li>'ILDR1'</li><li>'C19orf54'</li><li>'COL21A1'</li><li>'TMC5'</li><li>'C19orf44'</li><li>'TKFC'</li><li>'ZNF84'</li><li>'NAA15'</li><li>'METTL27'</li><li>'SLC39A4'</li><li>'WDR62'</li><li>'PAH'</li><li>'MBIP'</li><li>'TMEM186'</li><li>'WDR49'</li><li>'KIAA0586'</li><li>'STRN3'</li><li>'TMEM74'</li><li>'CDC7'</li><li>'IP6K3'</li><li>'POLQ'</li><li>'IL17RE'</li><li>'CNIH2'</li><li>'PVALB'</li><li>'CCDC18'</li><li>'MFSD4A'</li><li>'MAPK7'</li><li>'SLC4A5'</li><li>'NSL1'</li><li>'LRP8'</li><li>'COX10'</li><li>'PPP1R13L'</li><li>'GGT2'</li><li>'FCF1'</li><li>'OSBP2'</li><li>'PLEKHA5'</li><li>'TRIM36'</li><li>'GDAP2'</li><li>'SIX4'</li><li>'YIPF6'</li><li>'CEP72'</li><li>'LRATD2'</li><li>'SCML2'</li><li>'ADCY9'</li><li>'CEP152'</li><li>'COCH'</li><li>'MND1'</li><li>'GGT1'</li><li>'COQ6'</li><li>'PAK6'</li><li>'PLAGL2'</li><li>'SUPV3L1'</li><li>'KALRN'</li><li>'IRF6'</li><li>'SPATA17'</li><li>'HAUS6'</li><li>'TFG'</li><li>'EEF1A2'</li><li>'LYSMD1'</li><li>'DNAJC22'</li><li>'APLP1'</li><li>'LATS1'</li><li>'TAB3'</li><li>'LRRC37A3'</li><li>'ZNF576'</li><li>'RNF150'</li><li>'ARHGAP8'</li><li>'ATP7B'</li><li>'PLPP6'</li><li>'PEX11A'</li><li>'PPP1R3E'</li><li>'ATP6V0E2'</li><li>'ZNF133'</li><li>'CEP55'</li><li>'PPP1R14D'</li><li>'SMPD3'</li><li>'UGT8'</li><li>'RDM1'</li><li>'ZFP82'</li><li>'TRMT12'</li><li>'FUNDC1'</li><li>'KIF19'</li><li>'AGPAT4'</li><li>'CCDC117'</li><li>'CNOT6'</li><li>'ZNF512B'</li><li>'NSUN7'</li><li>'FANCF'</li><li>'SMU1'</li><li>'RRNAD1'</li><li>'SAXO2'</li><li>'SAMD10'</li><li>'TRO'</li><li>'CUL3'</li><li>'ZNF585A'</li><li>'GEMIN5'</li><li>'E2F2'</li><li>'ABCE1'</li><li>'RAB22A'</li><li>'JAKMIP2'</li><li>'TRMT10B'</li><li>'DCAF6'</li><li>'MESP1'</li><li>'BORA'</li><li>'ACKR3'</li><li>'ADD2'</li><li>'RXYLT1'</li><li>'MLH3'</li><li>'CNFN'</li><li>'KCNIP3'</li><li>'MCOLN3'</li><li>'SLC12A2'</li><li>'APAF1'</li><li>'ZIC2'</li><li>'DEPDC1'</li><li>'CDC42BPG'</li><li>'PIAS2'</li><li>'TMEM251'</li><li>'NOXRED1'</li><li>'PLGRKT'</li><li>'WDR35'</li><li>'ARHGEF39'</li><li>'DDX18'</li><li>'ETV2'</li><li>'TMEM179'</li><li>'CEP83'</li><li>'METTL17'</li><li>'CDKN3'</li><li>'PAQR6'</li><li>'C1orf56'</li><li>'EHF'</li><li>'TET1'</li><li>'NUP37'</li><li>'RIOX1'</li><li>'BCKDHA'</li><li>'LRR1'</li><li>'PPP2R3B'</li><li>'WDR5B'</li><li>'TESMIN'</li><li>'GPC3'</li><li>'PRAF2'</li><li>'ABAT'</li><li>'HOXD13'</li><li>'DNAH11'</li><li>'IQCG'</li><li>'UAP1'</li><li>'KCNK17'</li><li>'OCLM'</li><li>'TSEN2'</li><li>'PLXNA2'</li><li>'TARBP1'</li><li>'CABLES2'</li><li>'C6orf223'</li><li>'DDX59'</li><li>'CALML5'</li><li>'C17orf100'</li><li>'BBS10'</li><li>'ZNF345'</li><li>'MID1'</li><li>'VTI1B'</li><li>'ZNF599'</li><li>'SP1'</li><li>'ERCC2'</li><li>'CTDSPL2'</li><li>'PCDH20'</li><li>'CDC14A'</li><li>'B4GALNT4'</li><li>'MCM9'</li><li>'GEMIN4'</li><li>'ZNF354A'</li><li>'PLCB1'</li><li>'MTIF2'</li><li>'RB1CC1'</li></ol>
</dd>
</dl>




```R

```

# DE for Tumor_high MHCI_high vs MHCI_Low


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
Idents(data) <- "Tumor_MHCI_label2"

de_Tumor_MHCI <- FindMarkers(
  object = data,
  ident.1 = "Tumor_high_MHCI_high",  # First group to compare
  ident.2 = "Tumor_high_MHCI_low",    # Second group
  logfc.threshold = 0.25, # Minimum log2 fold change
  min.pct = 0.1,          # Gene detected in ≥10% of cells in either group
  test.use = "wilcox"     # Default Wilcoxon rank sum test
)

de_Tumor_MHCI$Gene <- rownames(de_Tumor_MHCI)
```


```R
de_Tumor_MHCI %>% head(15)
```


<table class="dataframe">
<caption>A data.frame: 15 × 6</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>Gene</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B2M</th><td> 0.000000e+00</td><td>0.6443753</td><td>1.000</td><td>1.000</td><td> 0.000000e+00</td><td>B2M     </td></tr>
	<tr><th scope=row>TAP1</th><td>3.747407e-285</td><td>1.5546946</td><td>0.907</td><td>0.622</td><td>6.723973e-281</td><td>TAP1    </td></tr>
	<tr><th scope=row>TAP2</th><td>2.070688e-253</td><td>1.8187478</td><td>0.779</td><td>0.485</td><td>3.715436e-249</td><td>TAP2    </td></tr>
	<tr><th scope=row>CD74</th><td>1.321012e-155</td><td>0.7369142</td><td>0.995</td><td>0.981</td><td>2.370293e-151</td><td>CD74    </td></tr>
	<tr><th scope=row>HLA-DRA</th><td>4.071785e-136</td><td>0.7700200</td><td>0.970</td><td>0.948</td><td>7.306004e-132</td><td>HLA-DRA </td></tr>
	<tr><th scope=row>IFITM3</th><td>2.012365e-121</td><td>0.7131448</td><td>0.940</td><td>0.906</td><td>3.610787e-117</td><td>IFITM3  </td></tr>
	<tr><th scope=row>VIM</th><td>4.030422e-116</td><td>0.6374194</td><td>0.971</td><td>0.962</td><td>7.231787e-112</td><td>VIM     </td></tr>
	<tr><th scope=row>C1QB</th><td>8.725089e-105</td><td>0.7204158</td><td>0.919</td><td>0.872</td><td>1.565543e-100</td><td>C1QB    </td></tr>
	<tr><th scope=row>S100A6</th><td> 1.069319e-97</td><td>0.4796522</td><td>0.985</td><td>0.977</td><td> 1.918678e-93</td><td>S100A6  </td></tr>
	<tr><th scope=row>S100A10</th><td> 1.421101e-91</td><td>0.7294962</td><td>0.878</td><td>0.845</td><td> 2.549882e-87</td><td>S100A10 </td></tr>
	<tr><th scope=row>CTSB</th><td> 5.734065e-91</td><td>0.5737418</td><td>0.953</td><td>0.925</td><td> 1.028863e-86</td><td>CTSB    </td></tr>
	<tr><th scope=row>LYZ</th><td> 4.461857e-86</td><td>0.8472520</td><td>0.834</td><td>0.809</td><td> 8.005910e-82</td><td>LYZ     </td></tr>
	<tr><th scope=row>LAPTM5</th><td> 1.811341e-81</td><td>0.6645867</td><td>0.879</td><td>0.851</td><td> 3.250089e-77</td><td>LAPTM5  </td></tr>
	<tr><th scope=row>WARS</th><td> 1.436961e-78</td><td>0.6097824</td><td>0.985</td><td>0.968</td><td> 2.578339e-74</td><td>WARS    </td></tr>
	<tr><th scope=row>HLA-DPA1</th><td> 2.486971e-77</td><td>0.7884652</td><td>0.828</td><td>0.775</td><td> 4.462373e-73</td><td>HLA-DPA1</td></tr>
</tbody>
</table>




```R
# check significant genes
de_Tumor_MHCI_sig_label <- label_sig_genes(de_Tumor_MHCI)
de_Tumor_MHCI_sig_DE <- de_Tumor_MHCI_sig_label[de_Tumor_MHCI_sig_label$sig_DE != 'NO',]
print(dim(de_Tumor_MHCI))
print(dim(de_Tumor_MHCI_sig_DE[]))
de_Tumor_MHCI_sig_DE %>% head(5)
table(de_Tumor_MHCI_sig_DE$sig_DE)
```

    [1] 3215    6
    [1] 844   7



<table class="dataframe">
<caption>A data.frame: 5 × 7</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>Gene</th><th scope=col>sig_DE</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>B2M</th><td> 0.000000e+00</td><td>0.6443753</td><td>1.000</td><td>1.000</td><td> 0.000000e+00</td><td>B2M    </td><td>UP</td></tr>
	<tr><th scope=row>TAP1</th><td>3.747407e-285</td><td>1.5546946</td><td>0.907</td><td>0.622</td><td>6.723973e-281</td><td>TAP1   </td><td>UP</td></tr>
	<tr><th scope=row>TAP2</th><td>2.070688e-253</td><td>1.8187478</td><td>0.779</td><td>0.485</td><td>3.715436e-249</td><td>TAP2   </td><td>UP</td></tr>
	<tr><th scope=row>CD74</th><td>1.321012e-155</td><td>0.7369142</td><td>0.995</td><td>0.981</td><td>2.370293e-151</td><td>CD74   </td><td>UP</td></tr>
	<tr><th scope=row>HLA-DRA</th><td>4.071785e-136</td><td>0.7700200</td><td>0.970</td><td>0.948</td><td>7.306004e-132</td><td>HLA-DRA</td><td>UP</td></tr>
</tbody>
</table>




    
    DOWN   UP 
     309  535 



```R
# Extract lists of genes UP and DOWN regulated
de_Tumor_MHCI_sig_DE_list <- list()
de_Tumor_MHCI_sig_DE_list$UP <- rownames(de_Tumor_MHCI_sig_DE[de_Tumor_MHCI_sig_DE$sig_DE=="UP",])
de_Tumor_MHCI_sig_DE_list$DOWN <- rownames(de_Tumor_MHCI_sig_DE[de_Tumor_MHCI_sig_DE$sig_DE=="DOWN",])
```


```R
# check list
de_Tumor_MHCI_sig_DE_list
```


<dl>
	<dt>$UP</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'B2M'</li><li>'TAP1'</li><li>'TAP2'</li><li>'CD74'</li><li>'HLA-DRA'</li><li>'IFITM3'</li><li>'VIM'</li><li>'C1QB'</li><li>'S100A10'</li><li>'LYZ'</li><li>'LAPTM5'</li><li>'WARS'</li><li>'HLA-DPA1'</li><li>'BGN'</li><li>'APOE'</li><li>'TYROBP'</li><li>'HLA-DPB1'</li><li>'LUM'</li><li>'PLAAT4'</li><li>'A2M'</li><li>'IFI27'</li><li>'HLA-F'</li><li>'IFITM2'</li><li>'TIMP1'</li><li>'C1R'</li><li>'LBH'</li><li>'AEBP1'</li><li>'DCN'</li><li>'COL1A1'</li><li>'LGALS1'</li><li>'COL6A3'</li><li>'IL32'</li><li>'CD44'</li><li>'MYL9'</li><li>'GBP1'</li><li>'IGKC'</li><li>'FCER1G'</li><li>'C1QC'</li><li>'IGLC1'</li><li>'MMP2'</li><li>'BST2'</li><li>'TGM2'</li><li>'HLA-A'</li><li>'SPARC'</li><li>'RARRES2'</li><li>'COL1A2'</li><li>'LGALS3'</li><li>'HLA-DMB'</li><li>'TAGLN'</li><li>'NBL1'</li><li>'SFTPB'</li><li>'TGFBI'</li><li>'COL3A1'</li><li>'GBP2'</li><li>'CXCL9'</li><li>'GBP5'</li><li>'IGHA1'</li><li>'SULF1'</li><li>'FCGR3A'</li><li>'S100A9'</li><li>'C1S'</li><li>'RNASE1'</li><li>'IGHG1'</li><li>'SERPING1'</li><li>'LRP1'</li><li>'PSMB9'</li><li>'IFI44'</li><li>'TIMP3'</li><li>'CTSH'</li><li>'JCHAIN'</li><li>'C1QA'</li><li>'TRIM22'</li><li>'HLA-DMA'</li><li>'MS4A6A'</li><li>'S100A4'</li><li>'SRGN'</li><li>'CD163'</li><li>'IL1R1'</li><li>'CYBB'</li><li>'PARP14'</li><li>'LTBP2'</li><li>'RGCC'</li><li>'CIITA'</li><li>'COL5A1'</li><li>'SPI1'</li><li>'CD52'</li><li>'GBP4'</li><li>'MMP14'</li><li>'CD37'</li><li>'IFITM1'</li><li>'CTSS'</li><li>'IGFBP7'</li><li>'CCDC69'</li><li>'PTPRC'</li><li>'ANXA1'</li><li>'FGR'</li><li>'TRBC2'</li><li>'DPYD'</li><li>'CCN1'</li><li>'APOL1'</li><li>'GJB2'</li><li>'TSC22D3'</li><li>'AIF1'</li><li>'LSP1'</li><li>'GXYLT2'</li><li>'CCDC80'</li><li>'XAF1'</li><li>'CAVIN1'</li><li>'SAMD9L'</li><li>'IKZF1'</li><li>'FOS'</li><li>'ACTA2'</li><li>'IL4I1'</li><li>'S100A8'</li><li>'ENC1'</li><li>'CDH11'</li><li>'PRRX1'</li><li>'IL7R'</li><li>'HLA-DQA1'</li><li>'TYMP'</li><li>'TRAC'</li><li>'COL10A1'</li><li>'ZEB2'</li><li>'CALHM6'</li><li>'CYP1B1'</li><li>'GYPC'</li><li>'COL8A1'</li><li>'POSTN'</li><li>'NFIX'</li><li>'RAB31'</li><li>'RGS1'</li><li>'VSIG4'</li><li>'APOL3'</li><li>'CD53'</li><li>'CTSE'</li><li>'IFI44L'</li><li>'ITGB2'</li><li>'FPR3'</li><li>'PLEK'</li><li>'COL5A2'</li><li>'NEDD9'</li><li>'RAC2'</li><li>'CD14'</li><li>'PECAM1'</li><li>'SERPINE1'</li><li>'CXCL13'</li><li>'CNN2'</li><li>'TGFB1'</li><li>'IL10RA'</li><li>'CCL5'</li><li>'LAMA2'</li><li>'CD3E'</li><li>'VGLL3'</li><li>'PODN'</li><li>'CCL2'</li><li>'MYLK'</li><li>'FBLN2'</li><li>'EPSTI1'</li><li>'GAS7'</li><li>'ENG'</li><li>'SLCO2B1'</li><li>'CD3D'</li><li>'HLA-DOA'</li><li>'CERCAM'</li><li>'MSN'</li><li>'FAP'</li><li>'PLXDC2'</li><li>'CYTIP'</li><li>'MZB1'</li><li>'LCP2'</li><li>'LST1'</li><li>'C1orf162'</li><li>'THBS2'</li><li>'GPNMB'</li><li>'COL11A1'</li><li>'MYO1F'</li><li>'CSF1R'</li><li>'PTGDS'</li><li>'DUSP1'</li><li>'PRELP'</li><li>'CD4'</li><li>'OAS2'</li><li>'MARCO'</li><li>'CYP4B1'</li><li>'MAFB'</li><li>'SESN3'</li><li>'ADGRE5'</li><li>'GREM1'</li><li>'BICC1'</li><li>'VAMP5'</li><li>'ZCCHC24'</li><li>'CXCL12'</li><li>'MYO1G'</li><li>'IGHM'</li><li>'GLIPR1'</li><li>'SLC1A3'</li><li>'LAIR1'</li><li>'EMILIN1'</li><li>'PDLIM2'</li><li>'EVA1C'</li><li>⋯</li><li>'SLAMF8'</li><li>'CSF2RA'</li><li>'SFRP4'</li><li>'IDO1'</li><li>'PHLDA3'</li><li>'AFAP1L1'</li><li>'GPSM3'</li><li>'FAM20A'</li><li>'PREX1'</li><li>'ANGPTL2'</li><li>'LOXL1'</li><li>'PARVG'</li><li>'HTRA1'</li><li>'GIMAP1'</li><li>'RASAL3'</li><li>'ANXA3'</li><li>'MRVI1'</li><li>'SFRP1'</li><li>'TWIST2'</li><li>'ITGA11'</li><li>'SAMSN1'</li><li>'MYOF'</li><li>'CLEC11A'</li><li>'SCARA3'</li><li>'FCGR1B'</li><li>'AQP1'</li><li>'TFEC'</li><li>'FBLN5'</li><li>'CREB3L1'</li><li>'FAM49A'</li><li>'CCL4'</li><li>'TNFRSF1B'</li><li>'DOK2'</li><li>'LEF1'</li><li>'ANKRD22'</li><li>'TREM2'</li><li>'ARHGEF6'</li><li>'KYNU'</li><li>'PPP1R14A'</li><li>'CLEC14A'</li><li>'IRF4'</li><li>'RRAS'</li><li>'ACTG2'</li><li>'LRRK2'</li><li>'ENTPD1'</li><li>'SEPTIN6'</li><li>'NCF2'</li><li>'PRF1'</li><li>'CD96'</li><li>'CCL21'</li><li>'CD84'</li><li>'TRAJ23'</li><li>'SIGLEC10'</li><li>'RCSD1'</li><li>'EHD2'</li><li>'SFRP2'</li><li>'TNC'</li><li>'POU2AF1'</li><li>'TRAF3IP3'</li><li>'IGSF6'</li><li>'TAGAP'</li><li>'EFEMP2'</li><li>'GGT5'</li><li>'HIST1H4I'</li><li>'RNASE6'</li><li>'CMKLR1'</li><li>'ADGRA2'</li><li>'IGHG4'</li><li>'MSRB3'</li><li>'SSC5D'</li><li>'GPR65'</li><li>'TIGIT'</li><li>'ALOX5'</li><li>'ETV7'</li><li>'BNC2'</li><li>'AQP9'</li><li>'HCK'</li><li>'TSPAN2'</li><li>'EMILIN2'</li><li>'PLN'</li><li>'CYBRD1'</li><li>'HAVCR2'</li><li>'TRAF1'</li><li>'CD7'</li><li>'DDR2'</li><li>'RASGRP3'</li><li>'CTSO'</li><li>'CD93'</li><li>'OAF'</li><li>'OSMR'</li><li>'PGGHG'</li><li>'NCF4'</li><li>'ITM2A'</li><li>'CFB'</li><li>'FNDC1'</li><li>'CCN5'</li><li>'COL15A1'</li><li>'IL21R'</li><li>'MAP3K8'</li><li>'SFTPD'</li><li>'TGFB2'</li><li>'PTGIR'</li><li>'MOXD1'</li><li>'ARHGAP25'</li><li>'IGHJ2'</li><li>'TMEM140'</li><li>'RFLNB'</li><li>'FMO2'</li><li>'CTLA4'</li><li>'ATXN1'</li><li>'ZBP1'</li><li>'CYTH4'</li><li>'CCL18'</li><li>'LATS2'</li><li>'COL8A2'</li><li>'COL12A1'</li><li>'F2R'</li><li>'TNFRSF12A'</li><li>'PSTPIP2'</li><li>'CPED1'</li><li>'FOXF2'</li><li>'RHOJ'</li><li>'APOD'</li><li>'ADGRE2'</li><li>'MLKL'</li><li>'GEM'</li><li>'CLDN5'</li><li>'SLAMF7'</li><li>'GNG11'</li><li>'KCNE3'</li><li>'TMEM273'</li><li>'UNC5B'</li><li>'HCST'</li><li>'ACAP1'</li><li>'IGKV4-1'</li><li>'GZMA'</li><li>'NPL'</li><li>'KCNMB1'</li><li>'TMC8'</li><li>'ADAMTS14'</li><li>'GZMB'</li><li>'KIAA1755'</li><li>'INPP5D'</li><li>'PPP1R16B'</li><li>'CNN1'</li><li>'GJA1'</li><li>'ITGA4'</li><li>'PDPN'</li><li>'SHROOM4'</li><li>'ALDH1A3'</li><li>'SIGLEC1'</li><li>'CA13'</li><li>'ECM2'</li><li>'SH2D2A'</li><li>'SCARF1'</li><li>'TNFRSF4'</li><li>'ABI3'</li><li>'SAMD4A'</li><li>'TNFSF13B'</li><li>'PLEKHO2'</li><li>'PTGS2'</li><li>'ADARB1'</li><li>'PTPRU'</li><li>'NAIP'</li><li>'FZD4'</li><li>'SH2D4A'</li><li>'PIK3AP1'</li><li>'MARCH1'</li><li>'AKAP12'</li><li>'CCND2'</li><li>'SNAI2'</li><li>'CDK6'</li><li>'ASPN'</li><li>'IGHG3'</li><li>'ACVRL1'</li><li>'F3'</li><li>'FBXL7'</li><li>'GLRX'</li><li>'NFYB'</li><li>'CSGALNACT1'</li><li>'AASS'</li><li>'SERPINE2'</li><li>'ROBO4'</li><li>'NPR3'</li><li>'RASGRF2'</li><li>'EVA1B'</li><li>'NTN4'</li><li>'RSAD2'</li><li>'HEG1'</li><li>'GNLY'</li><li>'AOAH'</li><li>'EDN1'</li><li>'TRANK1'</li><li>'PLA2G15'</li><li>'SLC31A2'</li><li>'ABLIM3'</li><li>'ITGB7'</li><li>'SYDE1'</li><li>'SPHK1'</li><li>'CALCRL'</li></ol>
</dd>
	<dt>$DOWN</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CACNA1A'</li><li>'AZGP1'</li><li>'LAMP1'</li><li>'SCG3'</li><li>'AMACR'</li><li>'PCID2'</li><li>'SCG2'</li><li>'GRTP1'</li><li>'CHGB'</li><li>'AOC1'</li><li>'PKHD1L1'</li><li>'JSRP1'</li><li>'RETREG1'</li><li>'SV2B'</li><li>'ANXA13'</li><li>'COL22A1'</li><li>'SCGN'</li><li>'SLC7A14'</li><li>'LAMA1'</li><li>'PEX5L'</li><li>'NRXN1'</li><li>'SLC38A3'</li><li>'KCNH6'</li><li>'MUC13'</li><li>'RIMBP2'</li><li>'NKX2-2'</li><li>'AKR1C1'</li><li>'TNNT1'</li><li>'KCNH7'</li><li>'TMCO3'</li><li>'ADPRHL1'</li><li>'KCNJ3'</li><li>'GPX2'</li><li>'SCN3A'</li><li>'ELAVL4'</li><li>'RNF128'</li><li>'HOXD8'</li><li>'HEPACAM2'</li><li>'CLNK'</li><li>'SYT7'</li><li>'SNAP25'</li><li>'AKR1C2'</li><li>'AIFM3'</li><li>'MPP6'</li><li>'PLEKHB1'</li><li>'VIL1'</li><li>'RUNDC3A'</li><li>'CHGA'</li><li>'SYT1'</li><li>'UBE2QL1'</li><li>'RNF182'</li><li>'OGDHL'</li><li>'AKR1C3'</li><li>'SEZ6L2'</li><li>'SHISA2'</li><li>'PLS1'</li><li>'GJB1'</li><li>'OPCML'</li><li>'PTPRN2'</li><li>'MARCH4'</li><li>'TENT5A'</li><li>'CUL4A'</li><li>'SOGA3'</li><li>'SLC7A11'</li><li>'PCSK2'</li><li>'LTK'</li><li>'RFX6'</li><li>'NRN1'</li><li>'SPTBN4'</li><li>'CSAG3'</li><li>'UNC13A'</li><li>'NPC1L1'</li><li>'SATB2'</li><li>'OR51E1'</li><li>'STMN3'</li><li>'SLC35D3'</li><li>'HGD'</li><li>'GAS2'</li><li>'MARCH11'</li><li>'SGSM1'</li><li>'NOL4'</li><li>'PRR26'</li><li>'MS4A8'</li><li>'CXXC4'</li><li>'KCNH8'</li><li>'HNF4A'</li><li>'FAM167A'</li><li>'PKD1L3'</li><li>'PCDHB16'</li><li>'CHRNB2'</li><li>'ANKS4B'</li><li>'DENND4A'</li><li>'GRIA2'</li><li>'RPRML'</li><li>'BMP7'</li><li>'ONECUT2'</li><li>'NEURL1'</li><li>'HNF1A'</li><li>'TMEM178B'</li><li>'CBLC'</li><li>'CRYM'</li><li>'NPAS3'</li><li>'HCN1'</li><li>'XKR7'</li><li>'MAPK8IP2'</li><li>'CRYBA2'</li><li>'PHGR1'</li><li>'MLIP'</li><li>'KIRREL2'</li><li>'RELN'</li><li>'PPARGC1A'</li><li>'GABRA3'</li><li>'SPOCK3'</li><li>'A1CF'</li><li>'ANK1'</li><li>'CA9'</li><li>'PKD2L1'</li><li>'HSPA2'</li><li>'PCDHB14'</li><li>'GLYATL3'</li><li>'USH1C'</li><li>'HOXB9'</li><li>'PTPRN'</li><li>'CLDN18'</li><li>'CELF4'</li><li>'PPM1E'</li><li>'PGM5'</li><li>'PTGES3L'</li><li>'ELOVL2'</li><li>'IRS2'</li><li>'MOV10L1'</li><li>'TMEM163'</li><li>'DDC'</li><li>'UNC80'</li><li>'OR51F2'</li><li>'KIF5C'</li><li>'HOXD9'</li><li>'SEZ6L'</li><li>'ISL1'</li><li>'TEX101'</li><li>'MCF2L'</li><li>'FABP6'</li><li>'NEFL'</li><li>'TM4SF5'</li><li>'CCDC13'</li><li>'GFRA3'</li><li>'ZNF676'</li><li>'TMEM169'</li><li>'CHST9'</li><li>'KCNJ6'</li><li>'PCDHB5'</li><li>'PCDHB11'</li><li>'FRZB'</li><li>'LNX1'</li><li>'FOXA3'</li><li>'DTX1'</li><li>'NPIPB15'</li><li>'PDGFD'</li><li>'ST18'</li><li>'SKOR1'</li><li>'PROC'</li><li>'CYP4F2'</li><li>'F10'</li><li>'FGF14'</li><li>'RAB3C'</li><li>'DLL4'</li><li>'SRRM4'</li><li>'HS3ST4'</li><li>'DAPL1'</li><li>'SYT4'</li><li>'ENTPD8'</li><li>'CNTNAP5'</li><li>'PRAME'</li><li>'ZNF727'</li><li>'TENM1'</li><li>'JAKMIP2'</li><li>'CDH8'</li><li>'TMEM74'</li><li>'GPM6A'</li><li>'GREB1L'</li><li>'KLHL41'</li><li>'CACNA1H'</li><li>'SMIM6'</li><li>'HOXD4'</li><li>'KCTD16'</li><li>'SERPINA1'</li><li>'GOLT1A'</li><li>'CABP7'</li><li>'BSN'</li><li>'PCDHA12'</li><li>'FAM222A'</li><li>'GJC3'</li><li>'TAAR1'</li><li>'TLDC2'</li><li>'RTBDN'</li><li>'F7'</li><li>'HSPA4L'</li><li>'RNF186'</li><li>'RIPPLY2'</li><li>'ZNF304'</li><li>'MYO15A'</li><li>'PRRG2'</li><li>'KRT40'</li><li>'ACVR1C'</li><li>'PDIA2'</li><li>'TNNI3'</li><li>'DDO'</li><li>'ZNF208'</li><li>'FOXD1'</li><li>'TMED6'</li><li>'ZNF775'</li><li>'SOWAHA'</li><li>'AC008397.1'</li><li>'KCTD8'</li><li>'C2CD4A'</li><li>'SLC39A4'</li><li>'CA8'</li><li>'CYP3A5'</li><li>'SEMA3D'</li><li>'LYG2'</li><li>'LDLRAD1'</li><li>'PTPRH'</li><li>'LGSN'</li><li>'FBXL16'</li><li>'OR51C1P'</li><li>'CHRM4'</li><li>'PCDHGA10'</li><li>'GCNT3'</li><li>'NKAIN2'</li><li>'SLC39A5'</li><li>'AGMAT'</li><li>'MMP26'</li><li>'SLC7A4'</li><li>'GTSF1'</li><li>'HCN4'</li><li>'PROZ'</li><li>'LSAMP'</li><li>'CT83'</li><li>'PCDHGA2'</li><li>'GOLGA6L9'</li><li>'NDST4'</li><li>'PACRG'</li><li>'EPHA10'</li><li>'LIG4'</li><li>'NEUROD1'</li><li>'PLA2G12B'</li><li>'DCX'</li><li>'MYO7B'</li><li>'AMER3'</li><li>'MMP16'</li><li>'NR0B1'</li><li>'AGR3'</li><li>'CDH10'</li><li>'ZNF716'</li><li>'ZNF804A'</li><li>'SERPINA10'</li><li>'MISP3'</li><li>'XG'</li><li>'KIAA0319'</li><li>'KCNA4'</li><li>'NEFM'</li><li>'RDH12'</li><li>'ERBB4'</li><li>'SERPINA6'</li><li>'SPINK5'</li><li>'VNN3'</li><li>'GP2'</li><li>'TMEM61'</li><li>'AC115220.1'</li><li>'HOXD1'</li><li>'ITIH2'</li><li>'KCND2'</li><li>'B3GALT2'</li><li>'IZUMO2'</li><li>'FER1L6'</li><li>'ANKRD30B'</li><li>'TMEM121B'</li><li>'MAGEA1'</li><li>'ASB16'</li><li>'LHFPL4'</li><li>'OPRD1'</li><li>'TIGD3'</li><li>'ZNF679'</li><li>'MYT1'</li><li>'ARX'</li><li>'LINGO2'</li><li>'MGAT4C'</li><li>'WNK4'</li><li>'RETNLB'</li><li>'PCYT1B'</li><li>'GUCY2C'</li><li>'TEX15'</li><li>'IGSF1'</li><li>'P2RX6'</li><li>'ZP2'</li><li>'CERKL'</li><li>'GPR158'</li><li>'GSTO2'</li><li>'AADAC'</li><li>'SMIM24'</li><li>'ADAM2'</li><li>'KISS1R'</li><li>'TRIM72'</li><li>'PRSS1'</li><li>'MEP1B'</li><li>'KLHL17'</li><li>'CBLN2'</li><li>'OLFM3'</li><li>'STRA6'</li></ol>
</dd>
</dl>




```R
options(repr.plot.width=10, repr.plot.height=10)
de_plot <- EnhancedVolcano(de_Tumor_MHCI,
                lab = rownames(de_Tumor_MHCI),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                title = 'MHCI_high vs MHCI_low, in Tumor_high spots')
de_plot

# save plot
png(file = "Output/07_DE_vasc_MN4_and_Tumor_MHCI/de_tumor_high_MHCI_high_vs_low_volcano.png",
    width = 1000,
    height = 1000)
print(de_plot)
dev.off()
```

    Warning message:
    “One or more p-values is 0. Converting to 10^-1 * current lowest non-zero p-value...”



<strong>pdf:</strong> 2



    
![png](output_52_2.png)
    



```R
FeaturePlot(data, features = rownames(de_Tumor_MHCI)[1:6], reduction='umap')
```


    
![png](output_53_0.png)
    


# Compare DEGs to cell marker lists


```R
# Macrophages
Macrophage <- c("CD68", "CSF1R", "CD14", "ITGAM", "MARCO", "FCGR1A")
M1_Macrophage <- c("CD80", "CD86", "NOS2", "IRF5", "TNF", "IL1B", "IL6", "IL12A", "CXCL9", "CXCL10")
M2_Macrophage <- c("MRC1", "CD163", "MSR1", "ARG1", "IL10", "TGFBI", "RETNLA", "CHI3L3")
#Macrophage_Activation <- c("CD68", "MARCO", "FCGR1A", "IL1B", "AIF1", "C1QC", "C1QA", "C1QB")

# T Cells
T_Cell <- c('CD3E', 'CD3D', 'CD3G', 'TRBC1', 'TRBC2', 'TRAC', 'CD28', 'CD247', 'CD5', 'CD27')
CD8_T_Cell <-  c('CD8A','CD8B','CD3D','CD3E','GZMA','GZMB','GZMK','CD69','PRF1','NKG7','IL7R')
T_Reg <- c('CD3E', 'CD4','IL2RA','FOXP3')
# IL2RA is CD25

# NK Cells
NK_Cell <- c('NCR1','FCGR3A','KLRK1','KLRD1','NCAM1','NKG7','GNLY','KLRC1','KLRF1','KLRB1','GZMB','CD160','XCL1','NCR3','KLRC2','XCL2','PRF1','KIR2DL3')

# B Cells
B_Cell <- c('MS4A1', 'CD19', 'CD79A', 'IGKC', 'IGHG3', 'MZB1', 'CR2', 'CD79B', 'BLNK', 'JCHAIN', 'IGHG2', 'IGHG1', 'CD22', 'BANK1', 'IGHD', 'IGHA1', 'IGHM', 'PAX5', 'BLK')

# Dendritic Cells
DC <- c('LAMP3', 'ITGAX', 'ITGAM', 'HLA-DRA', 'CD80', 'CD86', 'BATF3', 'IRF4', 'IRF8','CD209')
#Plasmacytoid_DC <- c('CLEC4C', 'NRP1')

All_Immune <- c('PTPRC','CD4','FOXP3','IFNG',
                'CD68', 'CSF1R', 'CD14', 'ITGAM', 'MARCO', 'FCGR1A',
                'CD80', 'CD86', 'NOS2', 'IRF5', 'TNF', 'IL1B', 'IL6', 'IL12A', 'CXCL9', 'CXCL10',
                'MRC1', 'CD163', 'MSR1', 'ARG1', 'IL10', 'TGFBI', 'RETNLA', 'CHI3L3',
                'FCGR1A', 'IL1B', 'AIF1', 'C1QC', 'C1QA', 'C1QB',
               'CD3E', 'CD3D', 'CD3G', 'TRBC1', 'TRBC2', 'TRAC', 'CD28', 'CD247', 'CD5', 'CD27',
                'CD8A','CD8B','GZMA','GZMB','GZMK','CD69','PRF1','NKG7','IL7R',
                'CD3E', 'CD4','IL2RA','FOXP3',
                'NCR1','FCGR3A','KLRK1','KLRD1','NCAM1','NKG7','GNLY','KLRC1','KLRF1','KLRB1','GZMB','CD160','XCL1','NCR3','KLRC2','XCL2','PRF1','KIR2DL3',
                'MS4A1', 'CD19', 'CD79A', 'IGKC', 'IGHG3', 'MZB1', 'CR2', 'CD79B', 'BLNK', 'JCHAIN', 'IGHG2', 'IGHG1', 'CD22', 'BANK1', 'IGHD', 'IGHA1', 'IGHM', 'PAX5', 'BLK',
                'LAMP3', 'ITGAX', 'ITGAM', 'HLA-DRA', 'CD80', 'CD86', 'BATF3', 'IRF4', 'IRF8','CD209',
               'CEACAM8', 'ELANE', 'KIT', 'TPSAB1', 'CCR3', 'EPX', 'FCER1A', 'LYZ')
```


```R
cell_type_list <- list("Macrophage"=Macrophage, "M1_Macrophage"=M1_Macrophage, "M2_Macrophage"=M2_Macrophage, 
                       "T_Cell"=T_Cell, "CD8_T_Cell"=CD8_T_Cell, "T_Reg"=T_Reg, "NK_Cell"=NK_Cell, "B_Cell"=B_Cell, "DC"=DC)
```


```R
cell_type_list
```


<dl>
	<dt>$Macrophage</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CD68'</li><li>'CSF1R'</li><li>'CD14'</li><li>'ITGAM'</li><li>'MARCO'</li><li>'FCGR1A'</li></ol>
</dd>
	<dt>$M1_Macrophage</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CD80'</li><li>'CD86'</li><li>'NOS2'</li><li>'IRF5'</li><li>'TNF'</li><li>'IL1B'</li><li>'IL6'</li><li>'IL12A'</li><li>'CXCL9'</li><li>'CXCL10'</li></ol>
</dd>
	<dt>$M2_Macrophage</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'MRC1'</li><li>'CD163'</li><li>'MSR1'</li><li>'ARG1'</li><li>'IL10'</li><li>'TGFBI'</li><li>'RETNLA'</li><li>'CHI3L3'</li></ol>
</dd>
	<dt>$T_Cell</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CD3E'</li><li>'CD3D'</li><li>'CD3G'</li><li>'TRBC1'</li><li>'TRBC2'</li><li>'TRAC'</li><li>'CD28'</li><li>'CD247'</li><li>'CD5'</li><li>'CD27'</li></ol>
</dd>
	<dt>$CD8_T_Cell</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CD8A'</li><li>'CD8B'</li><li>'CD3D'</li><li>'CD3E'</li><li>'GZMA'</li><li>'GZMB'</li><li>'GZMK'</li><li>'CD69'</li><li>'PRF1'</li><li>'NKG7'</li><li>'IL7R'</li></ol>
</dd>
	<dt>$T_Reg</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CD3E'</li><li>'CD4'</li><li>'IL2RA'</li><li>'FOXP3'</li></ol>
</dd>
	<dt>$NK_Cell</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'NCR1'</li><li>'FCGR3A'</li><li>'KLRK1'</li><li>'KLRD1'</li><li>'NCAM1'</li><li>'NKG7'</li><li>'GNLY'</li><li>'KLRC1'</li><li>'KLRF1'</li><li>'KLRB1'</li><li>'GZMB'</li><li>'CD160'</li><li>'XCL1'</li><li>'NCR3'</li><li>'KLRC2'</li><li>'XCL2'</li><li>'PRF1'</li><li>'KIR2DL3'</li></ol>
</dd>
	<dt>$B_Cell</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'MS4A1'</li><li>'CD19'</li><li>'CD79A'</li><li>'IGKC'</li><li>'IGHG3'</li><li>'MZB1'</li><li>'CR2'</li><li>'CD79B'</li><li>'BLNK'</li><li>'JCHAIN'</li><li>'IGHG2'</li><li>'IGHG1'</li><li>'CD22'</li><li>'BANK1'</li><li>'IGHD'</li><li>'IGHA1'</li><li>'IGHM'</li><li>'PAX5'</li><li>'BLK'</li></ol>
</dd>
	<dt>$DC</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'LAMP3'</li><li>'ITGAX'</li><li>'ITGAM'</li><li>'HLA-DRA'</li><li>'CD80'</li><li>'CD86'</li><li>'BATF3'</li><li>'IRF4'</li><li>'IRF8'</li><li>'CD209'</li></ol>
</dd>
</dl>




```R
library(ComplexHeatmap)
library(circlize)
```

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
    
    


# Plot heatmap of genes in the cell marker lists


```R
### 1. Prepare gene annotations ###
# Convert list of markers to data frame with cell type annotations
gene_annotations <- stack(cell_type_list) %>% 
  rename(Gene = values, CellType = ind)
gene_annotations %>% head(8)
```


<table class="dataframe">
<caption>A data.frame: 8 × 2</caption>
<thead>
	<tr><th></th><th scope=col>Gene</th><th scope=col>CellType</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>CD68  </td><td>Macrophage   </td></tr>
	<tr><th scope=row>2</th><td>CSF1R </td><td>Macrophage   </td></tr>
	<tr><th scope=row>3</th><td>CD14  </td><td>Macrophage   </td></tr>
	<tr><th scope=row>4</th><td>ITGAM </td><td>Macrophage   </td></tr>
	<tr><th scope=row>5</th><td>MARCO </td><td>Macrophage   </td></tr>
	<tr><th scope=row>6</th><td>FCGR1A</td><td>Macrophage   </td></tr>
	<tr><th scope=row>7</th><td>CD80  </td><td>M1_Macrophage</td></tr>
	<tr><th scope=row>8</th><td>CD86  </td><td>M1_Macrophage</td></tr>
</tbody>
</table>




```R
# Merge cell type annotations into the sig DE dfs
de_vasc_MN4_sig_DE_celltypes <- merge(de_vasc_MN4_sig_DE, gene_annotations, by='Gene')
de_vasc_MN4_not_normal_sig_DE_celltypes <- merge(de_vasc_MN4_not_normal_sig_DE, gene_annotations, by='Gene')
de_Tumor_MHCI_sig_DE_celltypes <- merge(de_Tumor_MHCI_sig_DE, gene_annotations, by='Gene')
```


```R
# sort
de_vasc_MN4_sig_DE_celltypes <- de_vasc_MN4_sig_DE_celltypes %>% arrange(CellType)
de_vasc_MN4_not_normal_sig_DE_celltypes <- de_vasc_MN4_not_normal_sig_DE_celltypes %>% arrange(CellType)
de_Tumor_MHCI_sig_DE_celltypes <- de_Tumor_MHCI_sig_DE_celltypes %>% arrange(CellType)
```


```R
dim(de_vasc_MN4_sig_DE_celltypes)
dim(de_vasc_MN4_not_normal_sig_DE_celltypes)
dim(de_Tumor_MHCI_sig_DE_celltypes)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>37</li><li>8</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>37</li><li>8</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>42</li><li>8</li></ol>




```R
# Do the same for non significant, keeping all genes in the DE results
# Merge cell type annotations into all DE dfs
de_vasc_MN4_sig_label_celltypes <- left_join(de_vasc_MN4_sig_label, gene_annotations, by='Gene')
de_vasc_MN4_not_normal_sig_label_celltypes <- left_join(de_vasc_MN4_not_normal_sig_label, gene_annotations, by='Gene')
de_Tumor_MHCI_sig_label_celltypes <- left_join(de_Tumor_MHCI_sig_label, gene_annotations, by='Gene')
```


```R
dim(de_vasc_MN4_sig_label_celltypes)
dim(de_vasc_MN4_not_normal_sig_label_celltypes)
dim(de_Tumor_MHCI_sig_label_celltypes)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>6241</li><li>8</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>5186</li><li>8</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>3222</li><li>8</li></ol>




```R
de_vasc_MN4_sig_label_celltypes %>% head(20)
```


<table class="dataframe">
<caption>A data.frame: 20 × 8</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>Gene</th><th scope=col>sig_DE</th><th scope=col>CellType</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>2.147527e-106</td><td> 0.9032005</td><td>0.989</td><td>0.905</td><td>3.853307e-102</td><td>S100A6 </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>2</th><td>1.662859e-103</td><td> 1.4309588</td><td>0.907</td><td>0.686</td><td> 2.983668e-99</td><td>LYZ    </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>3</th><td> 5.066896e-99</td><td> 1.3106240</td><td>0.907</td><td>0.658</td><td> 9.091531e-95</td><td>SPARC  </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>4</th><td> 1.825926e-98</td><td> 1.5237648</td><td>0.817</td><td>0.508</td><td> 3.276258e-94</td><td>IGFBP7 </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>5</th><td> 3.045898e-95</td><td> 1.3918593</td><td>0.864</td><td>0.609</td><td> 5.465256e-91</td><td>DCN    </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>6</th><td> 6.042315e-95</td><td> 1.4412885</td><td>0.924</td><td>0.761</td><td> 1.084173e-90</td><td>COL1A2 </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>7</th><td> 1.978932e-93</td><td> 1.1167877</td><td>0.971</td><td>0.877</td><td> 3.550798e-89</td><td>PSAP   </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>8</th><td> 3.275107e-92</td><td> 2.0726144</td><td>0.792</td><td>0.501</td><td> 5.876524e-88</td><td>SFTPB  </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>9</th><td> 4.341552e-89</td><td> 1.4467838</td><td>0.898</td><td>0.719</td><td> 7.790046e-85</td><td>COL3A1 </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>10</th><td> 3.323293e-88</td><td> 1.7609487</td><td>0.726</td><td>0.403</td><td> 5.962984e-84</td><td>CCN1   </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>11</th><td> 1.072526e-87</td><td> 1.0037082</td><td>0.983</td><td>0.890</td><td> 1.924433e-83</td><td>FN1    </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>12</th><td> 1.552077e-87</td><td> 0.7949820</td><td>0.998</td><td>0.964</td><td> 2.784891e-83</td><td>CD74   </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>13</th><td> 3.989375e-84</td><td> 2.4969689</td><td>0.633</td><td>0.318</td><td> 7.158136e-80</td><td>ICAM1  </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>14</th><td> 7.581220e-82</td><td>-0.6698764</td><td>0.980</td><td>0.986</td><td> 1.360298e-77</td><td>PPDPF  </td><td>DOWN</td><td>NA</td></tr>
	<tr><th scope=row>15</th><td> 1.021765e-79</td><td> 0.9186353</td><td>0.940</td><td>0.729</td><td> 1.833352e-75</td><td>HLA-E  </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>16</th><td> 1.416476e-79</td><td> 2.6697985</td><td>0.567</td><td>0.212</td><td> 2.541583e-75</td><td>PTGDS  </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>17</th><td> 5.721378e-78</td><td> 1.6264280</td><td>0.910</td><td>0.742</td><td> 1.026587e-73</td><td>APOE   </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>18</th><td> 1.434357e-76</td><td> 1.4050323</td><td>0.795</td><td>0.529</td><td> 2.573666e-72</td><td>THBS1  </td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>19</th><td> 2.222201e-76</td><td> 1.7072714</td><td>0.747</td><td>0.460</td><td> 3.987296e-72</td><td>HLA-DMA</td><td>UP  </td><td>NA</td></tr>
	<tr><th scope=row>20</th><td> 6.502776e-74</td><td> 1.0678490</td><td>0.979</td><td>0.870</td><td> 1.166793e-69</td><td>CTSD   </td><td>UP  </td><td>NA</td></tr>
</tbody>
</table>




```R
# save overall DE results
# 07_DE_vasc_MN4_and_Tumor_MHCI
write.csv(de_vasc_MN4_sig_label_celltypes, 'Output/07_DE_vasc_MN4_and_Tumor_MHCI/de_vasc_high_MN4_high_vs_low.csv')
write.csv(de_vasc_MN4_not_normal_sig_label_celltypes, 'Output/07_DE_vasc_MN4_and_Tumor_MHCI/de_vasc_high_MN4_high_vs_low_not_normal_spots.csv')
write.csv(de_Tumor_MHCI_sig_label_celltypes, 'Output/07_DE_vasc_MN4_and_Tumor_MHCI/de_tumor_high_MHCI_high_vs_low.csv')
```

# Immune Cell Type Markers that are Significantly Differentially Expressed in spots with activated vasculature (vasc_high_MN4_high vs vasc_high_MN4_low)


```R
de_vasc_MN4_sig_DE_celltypes
```


<table class="dataframe">
<caption>A data.frame: 37 × 8</caption>
<thead>
	<tr><th scope=col>Gene</th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>sig_DE</th><th scope=col>CellType</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>CSF1R </td><td>9.940324e-15</td><td>0.8869399</td><td>0.363</td><td>0.228</td><td>1.783592e-10</td><td>UP</td><td>Macrophage   </td></tr>
	<tr><td>MARCO </td><td>2.696057e-10</td><td>1.3253457</td><td>0.146</td><td>0.061</td><td>4.837535e-06</td><td>UP</td><td>Macrophage   </td></tr>
	<tr><td>IL6   </td><td>1.246597e-27</td><td>3.2259883</td><td>0.260</td><td>0.088</td><td>2.236769e-23</td><td>UP</td><td>M1_Macrophage</td></tr>
	<tr><td>CD163 </td><td>1.192172e-18</td><td>0.6932890</td><td>0.444</td><td>0.271</td><td>2.139115e-14</td><td>UP</td><td>M2_Macrophage</td></tr>
	<tr><td>CD27  </td><td>1.138961e-08</td><td>1.0742770</td><td>0.134</td><td>0.061</td><td>2.043638e-04</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD28  </td><td>2.588757e-07</td><td>1.9330414</td><td>0.104</td><td>0.046</td><td>4.645007e-03</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3D  </td><td>3.636451e-13</td><td>0.9348484</td><td>0.262</td><td>0.144</td><td>6.524884e-09</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3E  </td><td>1.186132e-08</td><td>0.9373007</td><td>0.208</td><td>0.123</td><td>2.128276e-04</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3G  </td><td>4.528055e-06</td><td>0.5856099</td><td>0.105</td><td>0.051</td><td>8.124689e-02</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>TRAC  </td><td>2.328580e-18</td><td>0.8756531</td><td>0.312</td><td>0.161</td><td>4.178172e-14</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>TRBC1 </td><td>1.561717e-24</td><td>1.3220019</td><td>0.404</td><td>0.225</td><td>2.802189e-20</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3D  </td><td>3.636451e-13</td><td>0.9348484</td><td>0.262</td><td>0.144</td><td>6.524884e-09</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>CD3E  </td><td>1.186132e-08</td><td>0.9373007</td><td>0.208</td><td>0.123</td><td>2.128276e-04</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>CD69  </td><td>1.627922e-06</td><td>1.4329094</td><td>0.113</td><td>0.057</td><td>2.920980e-02</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>CD8A  </td><td>1.007274e-08</td><td>1.4269991</td><td>0.127</td><td>0.056</td><td>1.807352e-04</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>GZMB  </td><td>1.920854e-14</td><td>1.6663688</td><td>0.220</td><td>0.106</td><td>3.446589e-10</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>IL7R  </td><td>5.259619e-20</td><td>0.9973431</td><td>0.292</td><td>0.131</td><td>9.437335e-16</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>NKG7  </td><td>5.354142e-14</td><td>1.9118055</td><td>0.171</td><td>0.067</td><td>9.606938e-10</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>PRF1  </td><td>2.028522e-07</td><td>1.8378709</td><td>0.104</td><td>0.045</td><td>3.639778e-03</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>CD3E  </td><td>1.186132e-08</td><td>0.9373007</td><td>0.208</td><td>0.123</td><td>2.128276e-04</td><td>UP</td><td>T_Reg        </td></tr>
	<tr><td>CD4   </td><td>1.348500e-15</td><td>1.0861864</td><td>0.322</td><td>0.191</td><td>2.419614e-11</td><td>UP</td><td>T_Reg        </td></tr>
	<tr><td>FCGR3A</td><td>2.148464e-20</td><td>0.6568109</td><td>0.605</td><td>0.456</td><td>3.854989e-16</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>GNLY  </td><td>5.159571e-21</td><td>2.2525869</td><td>0.314</td><td>0.168</td><td>9.257819e-17</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>GZMB  </td><td>1.920854e-14</td><td>1.6663688</td><td>0.220</td><td>0.106</td><td>3.446589e-10</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>NKG7  </td><td>5.354142e-14</td><td>1.9118055</td><td>0.171</td><td>0.067</td><td>9.606938e-10</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>PRF1  </td><td>2.028522e-07</td><td>1.8378709</td><td>0.104</td><td>0.045</td><td>3.639778e-03</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>CD79A </td><td>3.066069e-09</td><td>0.8657483</td><td>0.184</td><td>0.098</td><td>5.501447e-05</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHG1 </td><td>1.380242e-58</td><td>1.1991099</td><td>0.977</td><td>0.914</td><td>2.476567e-54</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHG3 </td><td>1.216595e-42</td><td>1.6380036</td><td>0.692</td><td>0.487</td><td>2.182936e-38</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHM  </td><td>5.303325e-17</td><td>0.7870159</td><td>0.549</td><td>0.412</td><td>9.515757e-13</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGKC  </td><td>1.811489e-62</td><td>1.0863116</td><td>1.000</td><td>0.996</td><td>3.250354e-58</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>JCHAIN</td><td>3.724674e-43</td><td>1.1477494</td><td>0.712</td><td>0.502</td><td>6.683183e-39</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>MZB1  </td><td>1.303181e-13</td><td>1.1336022</td><td>0.310</td><td>0.192</td><td>2.338298e-09</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IRF4  </td><td>3.106965e-14</td><td>1.3699707</td><td>0.238</td><td>0.122</td><td>5.574827e-10</td><td>UP</td><td>DC           </td></tr>
	<tr><td>IRF8  </td><td>2.611800e-20</td><td>1.6474644</td><td>0.254</td><td>0.107</td><td>4.686352e-16</td><td>UP</td><td>DC           </td></tr>
	<tr><td>ITGAX </td><td>1.692581e-48</td><td>1.8771450</td><td>0.499</td><td>0.228</td><td>3.036998e-44</td><td>UP</td><td>DC           </td></tr>
	<tr><td>LAMP3 </td><td>1.522148e-17</td><td>2.0563685</td><td>0.265</td><td>0.133</td><td>2.731190e-13</td><td>UP</td><td>DC           </td></tr>
</tbody>
</table>



# Immune Cell Type Markers that are Significantly Differentially Expressed in spots with activated vasculature (vasc_high_MN4_high vs vasc_high_MN4_low)

## Excluding Normal Lung and Fibrotic tissue


```R
de_vasc_MN4_not_normal_sig_DE_celltypes
```


<table class="dataframe">
<caption>A data.frame: 37 × 8</caption>
<thead>
	<tr><th scope=col>Gene</th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>sig_DE</th><th scope=col>CellType</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>CSF1R </td><td>2.162110e-14</td><td>0.7843697</td><td>0.380</td><td>0.226</td><td>3.879475e-10</td><td>UP</td><td>Macrophage   </td></tr>
	<tr><td>MARCO </td><td>4.619122e-07</td><td>1.1303512</td><td>0.130</td><td>0.060</td><td>8.288090e-03</td><td>UP</td><td>Macrophage   </td></tr>
	<tr><td>IL6   </td><td>2.230266e-12</td><td>2.2965276</td><td>0.205</td><td>0.093</td><td>4.001767e-08</td><td>UP</td><td>M1_Macrophage</td></tr>
	<tr><td>CD27  </td><td>1.810858e-10</td><td>1.2809972</td><td>0.152</td><td>0.060</td><td>3.249223e-06</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD28  </td><td>6.283557e-06</td><td>1.5677298</td><td>0.104</td><td>0.049</td><td>1.127459e-01</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3D  </td><td>1.594157e-13</td><td>1.0461695</td><td>0.278</td><td>0.145</td><td>2.860396e-09</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3E  </td><td>1.305412e-08</td><td>0.9084409</td><td>0.220</td><td>0.125</td><td>2.342300e-04</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>TRAC  </td><td>7.568220e-20</td><td>0.9143036</td><td>0.336</td><td>0.159</td><td>1.357966e-15</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>TRBC1 </td><td>2.792941e-18</td><td>1.1335283</td><td>0.395</td><td>0.225</td><td>5.011374e-14</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3D  </td><td>1.594157e-13</td><td>1.0461695</td><td>0.278</td><td>0.145</td><td>2.860396e-09</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>CD3E  </td><td>1.305412e-08</td><td>0.9084409</td><td>0.220</td><td>0.125</td><td>2.342300e-04</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>CD69  </td><td>7.490997e-05</td><td>1.0522205</td><td>0.105</td><td>0.057</td><td>1.000000e+00</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>CD8A  </td><td>8.012461e-06</td><td>1.0364121</td><td>0.112</td><td>0.055</td><td>1.437676e-01</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>GZMA  </td><td>3.687254e-05</td><td>0.6312308</td><td>0.122</td><td>0.067</td><td>6.616040e-01</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>GZMB  </td><td>2.563580e-13</td><td>1.5691584</td><td>0.228</td><td>0.106</td><td>4.599831e-09</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>IL7R  </td><td>9.652168e-17</td><td>0.9905942</td><td>0.281</td><td>0.125</td><td>1.731889e-12</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>NKG7  </td><td>1.159740e-12</td><td>1.6939586</td><td>0.171</td><td>0.066</td><td>2.080922e-08</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>PRF1  </td><td>1.460099e-06</td><td>1.5575866</td><td>0.103</td><td>0.045</td><td>2.619855e-02</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>CD3E  </td><td>1.305412e-08</td><td>0.9084409</td><td>0.220</td><td>0.125</td><td>2.342300e-04</td><td>UP</td><td>T_Reg        </td></tr>
	<tr><td>CD4   </td><td>1.134605e-14</td><td>0.9954606</td><td>0.336</td><td>0.192</td><td>2.035821e-10</td><td>UP</td><td>T_Reg        </td></tr>
	<tr><td>FCGR3A</td><td>4.839338e-20</td><td>0.6978449</td><td>0.632</td><td>0.459</td><td>8.683223e-16</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>GNLY  </td><td>5.677795e-16</td><td>1.9597135</td><td>0.315</td><td>0.175</td><td>1.018767e-11</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>GZMB  </td><td>2.563580e-13</td><td>1.5691584</td><td>0.228</td><td>0.106</td><td>4.599831e-09</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>NKG7  </td><td>1.159740e-12</td><td>1.6939586</td><td>0.171</td><td>0.066</td><td>2.080922e-08</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>PRF1  </td><td>1.460099e-06</td><td>1.5575866</td><td>0.103</td><td>0.045</td><td>2.619855e-02</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>CD79A </td><td>9.401558e-14</td><td>0.9381482</td><td>0.228</td><td>0.101</td><td>1.686922e-09</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHG1 </td><td>9.068149e-53</td><td>1.5647181</td><td>0.971</td><td>0.915</td><td>1.627098e-48</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHG2 </td><td>3.789452e-16</td><td>0.6659853</td><td>0.661</td><td>0.519</td><td>6.799413e-12</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHG3 </td><td>1.319786e-33</td><td>1.9749428</td><td>0.667</td><td>0.493</td><td>2.368092e-29</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHM  </td><td>4.765999e-21</td><td>0.9816685</td><td>0.582</td><td>0.409</td><td>8.551632e-17</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGKC  </td><td>3.354279e-52</td><td>1.4379913</td><td>1.000</td><td>0.996</td><td>6.018584e-48</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>JCHAIN</td><td>6.372454e-37</td><td>1.2557126</td><td>0.705</td><td>0.505</td><td>1.143409e-32</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>MZB1  </td><td>2.776963e-23</td><td>1.4351588</td><td>0.381</td><td>0.193</td><td>4.982704e-19</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IRF4  </td><td>4.822215e-17</td><td>1.3665640</td><td>0.271</td><td>0.124</td><td>8.652501e-13</td><td>UP</td><td>DC           </td></tr>
	<tr><td>IRF8  </td><td>1.479582e-17</td><td>1.5821061</td><td>0.258</td><td>0.109</td><td>2.654814e-13</td><td>UP</td><td>DC           </td></tr>
	<tr><td>ITGAX </td><td>1.392352e-30</td><td>1.5410779</td><td>0.463</td><td>0.229</td><td>2.498298e-26</td><td>UP</td><td>DC           </td></tr>
	<tr><td>LAMP3 </td><td>2.971789e-08</td><td>1.4733039</td><td>0.226</td><td>0.135</td><td>5.332281e-04</td><td>UP</td><td>DC           </td></tr>
</tbody>
</table>



# Immune Cell Type Markers that are Significantly Differentially Expressed in spots with MHCI high Tumor (Tumor_high_MHCI_high vs Tumor_high_MHCI_low)


```R
de_Tumor_MHCI_sig_DE_celltypes
```


<table class="dataframe">
<caption>A data.frame: 42 × 8</caption>
<thead>
	<tr><th scope=col>Gene</th><th scope=col>p_val</th><th scope=col>avg_log2FC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>sig_DE</th><th scope=col>CellType</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><td>CD14   </td><td> 2.104167e-13</td><td>0.6248579</td><td>0.449</td><td>0.396</td><td> 3.775507e-09</td><td>UP</td><td>Macrophage   </td></tr>
	<tr><td>CSF1R  </td><td> 3.640817e-11</td><td>0.6495452</td><td>0.378</td><td>0.319</td><td> 6.532717e-07</td><td>UP</td><td>Macrophage   </td></tr>
	<tr><td>FCGR1A </td><td> 9.584505e-07</td><td>0.7531299</td><td>0.131</td><td>0.091</td><td> 1.719748e-02</td><td>UP</td><td>Macrophage   </td></tr>
	<tr><td>MARCO  </td><td> 6.460363e-11</td><td>1.0266053</td><td>0.129</td><td>0.077</td><td> 1.159183e-06</td><td>UP</td><td>Macrophage   </td></tr>
	<tr><td>CXCL10 </td><td> 7.037682e-07</td><td>0.9956671</td><td>0.139</td><td>0.099</td><td> 1.262771e-02</td><td>UP</td><td>M1_Macrophage</td></tr>
	<tr><td>CXCL9  </td><td> 4.168324e-28</td><td>0.8689559</td><td>0.574</td><td>0.522</td><td> 7.479225e-24</td><td>UP</td><td>M1_Macrophage</td></tr>
	<tr><td>CD163  </td><td> 2.690324e-23</td><td>0.7958491</td><td>0.446</td><td>0.342</td><td> 4.827249e-19</td><td>UP</td><td>M2_Macrophage</td></tr>
	<tr><td>TGFBI  </td><td> 1.634516e-28</td><td>0.7186537</td><td>0.656</td><td>0.647</td><td> 2.932812e-24</td><td>UP</td><td>M2_Macrophage</td></tr>
	<tr><td>CD27   </td><td> 8.712778e-10</td><td>1.1496332</td><td>0.133</td><td>0.084</td><td> 1.563334e-05</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3D   </td><td> 8.136714e-12</td><td>0.8081177</td><td>0.259</td><td>0.195</td><td> 1.459971e-07</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3E   </td><td> 1.794833e-12</td><td>1.0317713</td><td>0.214</td><td>0.150</td><td> 3.220469e-08</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3G   </td><td> 2.992562e-07</td><td>0.9123726</td><td>0.113</td><td>0.074</td><td> 5.369554e-03</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>TRAC   </td><td> 5.279117e-15</td><td>0.9909232</td><td>0.304</td><td>0.233</td><td> 9.472320e-11</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>TRBC1  </td><td> 1.157391e-09</td><td>0.7865319</td><td>0.334</td><td>0.288</td><td> 2.076707e-05</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>TRBC2  </td><td> 1.523052e-17</td><td>1.1361684</td><td>0.310</td><td>0.228</td><td> 2.732812e-13</td><td>UP</td><td>T_Cell       </td></tr>
	<tr><td>CD3D   </td><td> 8.136714e-12</td><td>0.8081177</td><td>0.259</td><td>0.195</td><td> 1.459971e-07</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>CD3E   </td><td> 1.794833e-12</td><td>1.0317713</td><td>0.214</td><td>0.150</td><td> 3.220469e-08</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>GZMA   </td><td> 8.722907e-04</td><td>0.6038671</td><td>0.137</td><td>0.111</td><td> 1.000000e+00</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>GZMB   </td><td> 1.146291e-03</td><td>0.7275884</td><td>0.164</td><td>0.139</td><td> 1.000000e+00</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>IL7R   </td><td> 4.454690e-15</td><td>1.0704423</td><td>0.241</td><td>0.163</td><td> 7.993050e-11</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>NKG7   </td><td> 9.793683e-08</td><td>1.0523836</td><td>0.140</td><td>0.097</td><td> 1.757281e-03</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>PRF1   </td><td> 1.828803e-05</td><td>0.7764318</td><td>0.100</td><td>0.069</td><td> 3.281421e-01</td><td>UP</td><td>CD8_T_Cell   </td></tr>
	<tr><td>CD3E   </td><td> 1.794833e-12</td><td>1.0317713</td><td>0.214</td><td>0.150</td><td> 3.220469e-08</td><td>UP</td><td>T_Reg        </td></tr>
	<tr><td>CD4    </td><td> 5.991989e-11</td><td>0.7338047</td><td>0.333</td><td>0.278</td><td> 1.075143e-06</td><td>UP</td><td>T_Reg        </td></tr>
	<tr><td>FCGR3A </td><td> 2.149556e-27</td><td>0.6744333</td><td>0.624</td><td>0.601</td><td> 3.856948e-23</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>GNLY   </td><td> 1.759905e-02</td><td>0.7414744</td><td>0.260</td><td>0.256</td><td> 1.000000e+00</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>GZMB   </td><td> 1.146291e-03</td><td>0.7275884</td><td>0.164</td><td>0.139</td><td> 1.000000e+00</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>NKG7   </td><td> 9.793683e-08</td><td>1.0523836</td><td>0.140</td><td>0.097</td><td> 1.757281e-03</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>PRF1   </td><td> 1.828803e-05</td><td>0.7764318</td><td>0.100</td><td>0.069</td><td> 3.281421e-01</td><td>UP</td><td>NK_Cell      </td></tr>
	<tr><td>CD79A  </td><td> 1.088140e-07</td><td>0.9721289</td><td>0.180</td><td>0.135</td><td> 1.952450e-03</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHA1  </td><td> 1.410477e-27</td><td>0.8121121</td><td>0.809</td><td>0.756</td><td> 2.530819e-23</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHG1  </td><td> 1.576773e-26</td><td>1.0698249</td><td>0.955</td><td>0.954</td><td> 2.829203e-22</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHG2  </td><td> 3.611903e-08</td><td>0.7400414</td><td>0.612</td><td>0.611</td><td> 6.480838e-04</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHG3  </td><td> 7.494074e-03</td><td>1.0763119</td><td>0.594</td><td>0.638</td><td> 1.000000e+00</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGHM   </td><td> 2.709671e-10</td><td>0.7264163</td><td>0.542</td><td>0.518</td><td> 4.861963e-06</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>IGKC   </td><td> 5.546598e-41</td><td>1.1595338</td><td>1.000</td><td>0.999</td><td> 9.952260e-37</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>JCHAIN </td><td> 6.602584e-25</td><td>0.9755628</td><td>0.631</td><td>0.607</td><td> 1.184702e-20</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>MZB1   </td><td> 1.574739e-11</td><td>1.1311611</td><td>0.298</td><td>0.244</td><td> 2.825554e-07</td><td>UP</td><td>B_Cell       </td></tr>
	<tr><td>HLA-DRA</td><td>4.071785e-136</td><td>0.7700200</td><td>0.970</td><td>0.948</td><td>7.306004e-132</td><td>UP</td><td>DC           </td></tr>
	<tr><td>IRF4   </td><td> 9.436127e-06</td><td>0.8678220</td><td>0.204</td><td>0.168</td><td> 1.693124e-01</td><td>UP</td><td>DC           </td></tr>
	<tr><td>IRF8   </td><td> 7.018075e-08</td><td>0.7688286</td><td>0.226</td><td>0.178</td><td> 1.259253e-03</td><td>UP</td><td>DC           </td></tr>
	<tr><td>ITGAX  </td><td> 1.750455e-08</td><td>0.6450652</td><td>0.389</td><td>0.346</td><td> 3.140842e-04</td><td>UP</td><td>DC           </td></tr>
</tbody>
</table>




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


    
![png](output_75_0.png)
    



```R
unique(data@meta.data$vasc_MN4_label)
unique(data@meta.data$Tumor_MHCI_label2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_high_MN4_high'</li><li>'vasc_low'</li><li>'vasc_high_MN4_low'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_high_MHCI_high'</li></ol>




```R
# define colors
vasc_MN4_colors = c("vasc_high_MN4_low" = dittoColors()[1],
                    "vasc_high_MN4_high" = dittoColors()[2])
Tumor_MHCI_colors = c("Tumor_high_MHCI_low" = dittoColors()[1],
                      "Tumor_high_MHCI_high" = dittoColors()[2])
Tissue_Slice_colors = c("VA1" = dittoColors()[13],
                        "VB1" = dittoColors()[11],
                        "VA2" = dittoColors()[15])
```

# Heatmap of Immune cell type genes: vasc high, MN4 high vs MN4 low


```R
options(repr.plot.width=12, repr.plot.height=13)

Idents(data) <- 'vasc_MN4_label'

# remove pesky graph
data_temp <- data
data_temp@graphs <- list()
# Sub to exclude vasc_low
data_sub <- subset(data_temp, subset = vasc_MN4_label != 'vasc_low')

### 2. Extract expression matrix ###
# Get normalized expression data (using RNA assay by default)
expr_matrix <- GetAssayData(data_sub, slot = "data")

# Subset to genes present in both object and marker lists
common_genes <- intersect(rownames(expr_matrix), gene_annotations$Gene)
expr_matrix <- expr_matrix[common_genes, ]
gene_annotations <- gene_annotations %>% filter(Gene %in% common_genes)

### 3. Order genes by cell type ###
# Create ordered factor for genes
gene_order <- gene_annotations %>% 
  arrange(CellType, Gene) %>% 
  pull(Gene)

# Subset and order expression matrix
expr_matrix <- expr_matrix[gene_order, ]

### 4. Prepare column annotations ###
# Get sample origins (assuming vasc_MN4_label is in metadata)
col_annot <- data_sub@meta.data %>% 
  arrange(vasc_MN4_label, Tissue_Slice) %>% 
  select(vasc_MN4_label, Tissue_Slice)

# Order expression matrix columns by orig.ident
expr_matrix <- expr_matrix[, rownames(col_annot)]

# Create annotation object
ha <- HeatmapAnnotation(
  Tissue = col_annot$Tissue_Slice,
  MITE = col_annot$vasc_MN4_label,
  col = list(MITE = vasc_MN4_colors,
             Tissue = Tissue_Slice_colors),
  annotation_name_side = "left"
)

# scale data for heatmap
mat <- t(scale(t(expr_matrix)))

### 5. Create heatmap ###
ht <- Heatmap(
  matrix = as.matrix(mat),
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  
  # Row configurations
  row_split = gene_annotations$CellType,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  
  # Column configurations
  column_split = col_annot[c("vasc_MN4_label", "Tissue_Slice")],
  cluster_columns = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  
  # Additional parameters
  use_raster = TRUE,  # For large datasets
  border_gp = gpar(col = "black", lty = 1),
  column_title = "DE vasc high, MN4 high vs low: Significant Immune Marker Genes"
)

# Draw the plot
draw(ht)

# save plot
png(file = "Output/07_DE_vasc_MN4_and_Tumor_MHCI/de_vasc_high_MN4_high_vs_low_sig_immune_genes_heatmap.png",
    width = 1200,
    height = 1300)
draw(ht)
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



<strong>pdf:</strong> 2



    
![png](output_79_2.png)
    


# Now only significantly DE immune cell type genes, vasc high, MN4 high vs MN4 low


```R
options(repr.plot.width=12, repr.plot.height=12)

### 2. Extract expression matrix ###
# Get normalized expression data (using RNA assay by default)
expr_matrix <- GetAssayData(data_sub, slot = "data")

# Subset to genes present in both object and marker lists
common_genes <- intersect(rownames(expr_matrix), gene_annotations$Gene)
# subset to only significant genes
common_genes <- intersect(common_genes,de_vasc_MN4_sig_DE_celltypes$Gene)
expr_matrix <- expr_matrix[common_genes, ]
gene_annotations <- gene_annotations %>% filter(Gene %in% common_genes)

### 3. Order genes by cell type ###
# Create ordered factor for genes
gene_order <- gene_annotations %>% 
  arrange(CellType, Gene) %>% 
  pull(Gene)

# Subset and order expression matrix
expr_matrix <- expr_matrix[gene_order, ]

### 4. Prepare column annotations ###
# Get sample origins (assuming vasc_MN4_label is in metadata)
col_annot <- data_sub@meta.data %>% 
  arrange(vasc_MN4_label, Tissue_Slice) %>% 
  select(vasc_MN4_label, Tissue_Slice)

# Order expression matrix columns by orig.ident
expr_matrix <- expr_matrix[, rownames(col_annot)]

# Create annotation object
ha <- HeatmapAnnotation(
  Tissue = col_annot$Tissue_Slice,
  MITE = col_annot$vasc_MN4_label,
  col = list(MITE = vasc_MN4_colors,
             Tissue = Tissue_Slice_colors),
  annotation_name_side = "left"
)

# scale data for heatmap
mat <- t(scale(t(expr_matrix)))

### 5. Create heatmap ###
ht <- Heatmap(
  matrix = as.matrix(mat),
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  
  # Row configurations
  row_split = gene_annotations$CellType,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  
  # Column configurations
  column_split = col_annot[c("vasc_MN4_label", "Tissue_Slice")],
  cluster_columns = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  
  # Additional parameters
  use_raster = TRUE,  # For large datasets
  border_gp = gpar(col = "black", lty = 1),
  column_title = "DE vasc high, MN4 high vs low: Significant Immune Marker Genes plot2"
)

# Draw the plot
draw(ht)

# save plot
png(file = "Output/07_DE_vasc_MN4_and_Tumor_MHCI/de_vasc_high_MN4_high_vs_low_sig_immune_genes_heatmap_plot2.png",
    width = 1200,
    height = 1200)
draw(ht)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_81_1.png)
    


# Heatmap of Immune cell type genes: Tumor high, MHCI high vs MHCI low


```R
options(repr.plot.width=12, repr.plot.height=12)

Idents(data) <- 'Tumor_MHCI_label2'

# remove pesky graph
data_temp <- data
data_temp@graphs <- list()
# Sub to exclude vasc_low
data_sub <- subset(data_temp, subset = Tumor_MHCI_label2 != 'Tumor_low')

### 2. Extract expression matrix ###
# Get normalized expression data (using RNA assay by default)
expr_matrix <- GetAssayData(data_sub, slot = "data")

# Subset to genes present in both object and marker lists
common_genes <- intersect(rownames(expr_matrix), gene_annotations$Gene)
expr_matrix <- expr_matrix[common_genes, ]
gene_annotations <- gene_annotations %>% filter(Gene %in% common_genes)

### 3. Order genes by cell type ###
# Create ordered factor for genes
gene_order <- gene_annotations %>% 
  arrange(CellType, Gene) %>% 
  pull(Gene)

# Subset and order expression matrix
expr_matrix <- expr_matrix[gene_order, ]

# Create annotation object
ha <- HeatmapAnnotation(
  Tissue = col_annot$Tissue_Slice,
  MHCI = col_annot$Tumor_MHCI_label2,
  col = list(MHCI = Tumor_MHCI_colors,
             Tissue = Tissue_Slice_colors),
  annotation_name_side = "left"
)
print(ha)
### 4. Prepare column annotations ###
# Get sample origins (assuming vasc_MN4_label is in metadata)
col_annot <- data_sub@meta.data %>% 
  arrange(Tumor_MHCI_label2, Tissue_Slice) %>% 
  select(Tumor_MHCI_label2, Tissue_Slice)

# Order expression matrix columns by orig.ident
expr_matrix <- expr_matrix[, rownames(col_annot)]

# scale data for heatmap
mat <- t(scale(t(expr_matrix)))

### 5. Create heatmap ###
ht <- Heatmap(
  matrix = as.matrix(mat),
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  
  # Row configurations
  row_split = gene_annotations$CellType,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  
  # Column configurations
  column_split = col_annot[c("Tumor_MHCI_label2", "Tissue_Slice")],
  cluster_columns = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  
  # Additional parameters
  use_raster = TRUE,  # For large datasets
  border_gp = gpar(col = "black", lty = 1),
  column_title = "DE tumor high, MHCI high vs low: Significant Immune Marker Genes"
)

# Draw the plot
draw(ht)

# save plot
png(file = "Output/07_DE_vasc_MN4_and_Tumor_MHCI/de_tumor_high_MHCI_high_vs_low_sig_immune_genes_heatmap.png",
    width = 1200,
    height = 1200)
draw(ht)
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


    A HeatmapAnnotation object with 2 annotations
      name: heatmap_annotation_3 
      position: column 
      items: 5483 
      width: 1npc 
      height: 10.3514598035146mm 
      this object is subsettable
      13.4671666666667mm extension on the left 
    
       name annotation_type color_mapping height
     Tissue discrete vector  user-defined    5mm
       MHCI discrete vector  user-defined    5mm



<strong>pdf:</strong> 2



    
![png](output_83_3.png)
    


# Now only significantly DE immune cell type genes, Tumor high, MHCI high vs MHCI low


```R
### 2. Extract expression matrix ###
# Get normalized expression data (using RNA assay by default)
expr_matrix <- GetAssayData(data_sub, slot = "data")

# Subset to genes present in both object and marker lists
common_genes <- intersect(rownames(expr_matrix), gene_annotations$Gene)
# subset to only significant genes
common_genes <- intersect(common_genes,de_Tumor_MHCI_sig_DE_celltypes$Gene)
expr_matrix <- expr_matrix[common_genes, ]
gene_annotations <- gene_annotations %>% filter(Gene %in% common_genes)

### 3. Order genes by cell type ###
# Create ordered factor for genes
gene_order <- gene_annotations %>% 
  arrange(CellType, Gene) %>% 
  pull(Gene)

# Subset and order expression matrix
expr_matrix <- expr_matrix[gene_order, ]

### 4. Prepare column annotations ###
# Get sample origins (assuming vasc_MN4_label is in metadata)
col_annot <- data_sub@meta.data %>% 
  arrange(Tumor_MHCI_label2, Tissue_Slice) %>% 
  select(Tumor_MHCI_label2, Tissue_Slice)

# Order expression matrix columns by orig.ident
expr_matrix <- expr_matrix[, rownames(col_annot)]

# Create annotation object
ha <- HeatmapAnnotation(
  Tissue = col_annot$Tissue_Slice,
  MHCI = col_annot$Tumor_MHCI_label2,
  col = list(MHCI = Tumor_MHCI_colors,
             Tissue = Tissue_Slice_colors),
  annotation_name_side = "left"
)

# scale data for heatmap
mat <- t(scale(t(expr_matrix)))

### 5. Create heatmap ###
ht <- Heatmap(
  matrix = as.matrix(mat),
  name = "Expression",
  col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  
  # Row configurations
  row_split = gene_annotations$CellType,
  cluster_rows = FALSE,
  show_row_names = TRUE,
  row_names_gp = gpar(fontsize = 8),
  
  # Column configurations
  column_split = col_annot[c("Tumor_MHCI_label2", "Tissue_Slice")],
  cluster_columns = FALSE,
  show_column_names = FALSE,
  top_annotation = ha,
  
  # Additional parameters
  use_raster = TRUE,  # For large datasets
  border_gp = gpar(col = "black", lty = 1),
  column_title = "DE tumor high, MHCI high vs low: Significant Immune Marker Genes plot2"
)

# Draw the plot
draw(ht)

# save plot
png(file = "Output/07_DE_vasc_MN4_and_Tumor_MHCI/de_tumor_high_MHCI_high_vs_low_sig_immune_genes_heatmap_plot2.png",
    width = 1200,
    height = 1200)
draw(ht)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_85_1.png)
    



```R

```


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
     [1] grid      stats4    parallel  stats     graphics  grDevices utils    
     [8] datasets  methods   base     
    
    other attached packages:
     [1] circlize_0.4.16        ComplexHeatmap_2.24.0  EnhancedVolcano_1.26.0
     [4] ggrepel_0.9.6          paletteer_1.6.0        ggbeeswarm_0.7.2      
     [7] GGally_2.2.1           scales_1.4.0           RColorBrewer_1.1-3    
    [10] pheatmap_1.0.12        dittoSeq_1.20.0        lubridate_1.9.4       
    [13] forcats_1.0.0          stringr_1.5.1          dplyr_1.1.4           
    [16] purrr_1.0.4            readr_2.1.5            tidyr_1.3.1           
    [19] tibble_3.2.1           tidyverse_2.0.0        msigdbr_24.1.0        
    [22] org.Hs.eg.db_3.21.0    AnnotationDbi_1.70.0   IRanges_2.42.0        
    [25] S4Vectors_0.46.0       Biobase_2.68.0         BiocGenerics_0.54.0   
    [28] generics_0.1.4         ggplot2_3.5.2          doParallel_1.0.17     
    [31] iterators_1.0.14       foreach_1.5.2          Matrix_1.7-3          
    [34] Seurat_5.3.0           SeuratObject_5.1.0     sp_2.2-0              
    
    loaded via a namespace (and not attached):
      [1] RcppAnnoy_0.0.22            splines_4.5.0              
      [3] later_1.4.2                 pbdZMQ_0.3-14              
      [5] polyclip_1.10-7             fastDummies_1.7.5          
      [7] lifecycle_1.0.4             globals_0.18.0             
      [9] lattice_0.22-7              MASS_7.3-65                
     [11] magrittr_2.0.3              limma_3.64.1               
     [13] plotly_4.10.4               httpuv_1.6.16              
     [15] sctransform_0.4.2           spam_2.11-1                
     [17] spatstat.sparse_3.1-0       reticulate_1.42.0          
     [19] cowplot_1.1.3               pbapply_1.7-2              
     [21] DBI_1.2.3                   abind_1.4-8                
     [23] Rtsne_0.17                  GenomicRanges_1.60.0       
     [25] GenomeInfoDbData_1.2.14     irlba_2.3.5.1              
     [27] listenv_0.9.1               spatstat.utils_3.1-4       
     [29] goftest_1.2-3               RSpectra_0.16-2            
     [31] spatstat.random_3.4-1       fitdistrplus_1.2-2         
     [33] parallelly_1.44.0           codetools_0.2-20           
     [35] DelayedArray_0.34.1         shape_1.4.6.1              
     [37] tidyselect_1.2.1            UCSC.utils_1.4.0           
     [39] farver_2.1.2                matrixStats_1.5.0          
     [41] base64enc_0.1-3             spatstat.explore_3.4-3     
     [43] jsonlite_2.0.0              GetoptLong_1.0.5           
     [45] progressr_0.15.1            ggridges_0.5.6             
     [47] survival_3.8-3              tools_4.5.0                
     [49] ica_1.0-3                   Rcpp_1.0.14                
     [51] glue_1.8.0                  gridExtra_2.3              
     [53] SparseArray_1.8.0           MatrixGenerics_1.20.0      
     [55] GenomeInfoDb_1.44.0         IRdisplay_1.1              
     [57] withr_3.0.2                 BiocManager_1.30.25        
     [59] fastmap_1.2.0               digest_0.6.37              
     [61] timechange_0.3.0            R6_2.6.1                   
     [63] mime_0.13                   colorspace_2.1-1           
     [65] Cairo_1.6-2                 scattermore_1.2            
     [67] tensor_1.5                  spatstat.data_3.1-6        
     [69] RSQLite_2.3.11              data.table_1.17.4          
     [71] httr_1.4.7                  htmlwidgets_1.6.4          
     [73] S4Arrays_1.8.0              ggstats_0.9.0              
     [75] uwot_0.2.3                  pkgconfig_2.0.3            
     [77] gtable_0.3.6                blob_1.2.4                 
     [79] lmtest_0.9-40               SingleCellExperiment_1.30.1
     [81] XVector_0.48.0              htmltools_0.5.8.1          
     [83] dotCall64_1.2               clue_0.3-66                
     [85] png_0.1-8                   spatstat.univar_3.1-3      
     [87] rjson_0.2.23                tzdb_0.5.0                 
     [89] reshape2_1.4.4              uuid_1.2-1                 
     [91] nlme_3.1-168                curl_6.2.3                 
     [93] GlobalOptions_0.1.2         repr_1.1.7                 
     [95] zoo_1.8-14                  cachem_1.1.0               
     [97] KernSmooth_2.23-26          vipor_0.4.7                
     [99] miniUI_0.1.2                ggrastr_1.0.2              
    [101] pillar_1.10.2               vctrs_0.6.5                
    [103] RANN_2.6.2                  promises_1.3.3             
    [105] xtable_1.8-4                cluster_2.1.8.1            
    [107] beeswarm_0.4.0              evaluate_1.0.3             
    [109] magick_2.8.6                cli_3.6.5                  
    [111] compiler_4.5.0              rlang_1.1.6                
    [113] crayon_1.5.3                future.apply_1.11.3        
    [115] labeling_0.4.3              rematch2_2.1.2             
    [117] plyr_1.8.9                  stringi_1.8.7              
    [119] viridisLite_0.4.2           deldir_2.0-4               
    [121] assertthat_0.2.1            babelgene_22.9             
    [123] Biostrings_2.76.0           lazyeval_0.2.2             
    [125] spatstat.geom_3.4-1         IRkernel_1.3.2             
    [127] RcppHNSW_0.6.0              hms_1.1.3                  
    [129] patchwork_1.3.0             bit64_4.6.0-1              
    [131] future_1.49.0               statmod_1.5.0              
    [133] KEGGREST_1.48.0             shiny_1.10.0               
    [135] SummarizedExperiment_1.38.1 ROCR_1.0-11                
    [137] igraph_2.1.4                memoise_2.0.1              
    [139] bit_4.6.0                  



```R

```


```R

```
