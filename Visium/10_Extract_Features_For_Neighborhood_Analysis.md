# Extracting Features for Neighborhood Analysis of the Visium Spots

1. Around a Tumor MHCI High spot, do we see vascular spots with high activation score (MN4_EC)? Compared to Tumor MHCI Low
2. Around a Tumor MHCI High spot, do we see more NK/T cells compared with Tumor MHCI Low?
3. If you’re a vasc_high spot with high vs low activation (MN4_EC), do you tend to have NK/T cells near you?


## 1. Around a Tumor MHCI High spot, do we see vascular spots with high activation score (MN4_EC)? Compared to Tumor MHCI Low
- Tumor_MHCI_Label: TumorHigh_MHCIHigh and TumorHigh_MHCILow
- Average MN4_EC enrichment in the neighborhood

## 2. Around a Tumor MHCI High spot, do we see more NK/T cells compared with Tumor MHCI Low?
- Tumor_MHCI_Label: TumorHigh_MHCIHigh and TumorHigh_MHCILow
- Average NK Cell enrichment in the neighborhood
- Average T Cell enrichment in the neighborhood

## 3. If you’re a vasc_high spot with high vs low activation (MN4_EC), do you tend to have NK/T cells near you?
- vasc_label: vasc_high (can also take vasc_low and vasc_neg and subset it later)
- MN4_EC_label: MN4_EC_high and MN4_EC_low
- Average NK Cell enrichment in the neighborhood
- Average T Cell enrichment in the neighborhood

# Overall things to measure in each neighborhood, keeping all labels:
- Average MN4_EC enrichment
- Average NK Cell enrichment
- Average T Cell enrichment
- Average other cell type enrichments


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
library(ggbeeswarm, quietly = T)
#library(hdf5r)

options(repr.matrix.max.cols=1000, repr.matrix.max.rows=100)
```

# Load data


```R
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
AB1_A2@meta.data[1:3,]
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
unique(AB1_A2@meta.data$orig.ident)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'VA1'</li><li>'VA2'</li><li>'VB1'</li></ol>



# Unfortunately, GetTissueCoordinates only grabs coordinates for one of the tissue slices... I'll bring in the raw data for VA and VB, get each of their coords, combine them, and merge those in. 


```R
install.packages("hdf5r")
```

    Installing package into ‘/usr/local/lib/R/site-library’
    (as ‘lib’ is unspecified)
    



```R
VA1 <- Load10X_Spatial("/immuno/ian/Projects/2023_02_SCLC_10X_Visium/VA_count/outs/")
VB1 <- Load10X_Spatial("/immuno/ian/Projects/2023_02_SCLC_10X_Visium/VB_count/outs/")
VA2 <- Load10X_Spatial("/immuno/ian/Projects/2024_09_Barbie_Visium_Round2/spaceranger_count/VisA_count/outs/")
VA1$orig.ident <- "VA1"
VB1$orig.ident <- "VB1"
VA2$orig.ident <- "VA2"

va1_coords <- GetTissueCoordinates(VA1)
va1_coords$Barcode <- rownames(va1_coords)
vb1_coords <- GetTissueCoordinates(VB1)
vb1_coords$Barcode <- rownames(vb1_coords)
va2_coords <- GetTissueCoordinates(VA2)
va2_coords$Barcode <- rownames(va2_coords)

va1_coords[1:3,]
vb1_coords[1:3,]
va2_coords[1:3,]
```


<table class="dataframe">
<caption>A data.frame: 3 × 4</caption>
<thead>
	<tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>Barcode</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACAACGAATAGTTC-1</th><td>30730</td><td>26443</td><td>AAACAACGAATAGTTC-1</td><td>AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>AAACAAGTATCTCCCA-1</th><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td><td>AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>AAACAATCTACTAGCA-1</th><td>29633</td><td>20870</td><td>AAACAATCTACTAGCA-1</td><td>AAACAATCTACTAGCA-1</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 3 × 4</caption>
<thead>
	<tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>Barcode</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACAACGAATAGTTC-1</th><td>31283</td><td>26488</td><td>AAACAACGAATAGTTC-1</td><td>AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>AAACAAGTATCTCCCA-1</th><td>13243</td><td> 8788</td><td>AAACAAGTATCTCCCA-1</td><td>AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>AAACACCAATAACTGC-1</th><td>10072</td><td>25947</td><td>AAACACCAATAACTGC-1</td><td>AAACACCAATAACTGC-1</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 3 × 4</caption>
<thead>
	<tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>Barcode</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACAAGTATCTCCCA-1</th><td> 7764.086</td><td>20534.804</td><td>AAACAAGTATCTCCCA-1</td><td>AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>AAACACCAATAACTGC-1</th><td>24672.634</td><td>23846.305</td><td>AAACACCAATAACTGC-1</td><td>AAACACCAATAACTGC-1</td></tr>
	<tr><th scope=row>AAACAGAGCGACTCCT-1</th><td> 9484.343</td><td> 7768.555</td><td>AAACAGAGCGACTCCT-1</td><td>AAACAGAGCGACTCCT-1</td></tr>
</tbody>
</table>




```R
#test...
va1_coords_test1 <- GetTissueCoordinates(VA1, scale = 'hires')
va1_coords_test2 <- GetTissueCoordinates(VA1, scale = 'lowres')
va1_coords_test3 <- GetTissueCoordinates(VA1, scale = NULL)
va1_coords_test1[1:3,]
va1_coords_test2[1:3,]
va1_coords_test3[1:3,]
```


<table class="dataframe">
<caption>A data.frame: 3 × 3</caption>
<thead>
	<tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACAACGAATAGTTC-1</th><td>1778.356</td><td>1530.2662</td><td>AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>AAACAAGTATCTCCCA-1</th><td> 734.838</td><td> 505.9606</td><td>AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>AAACAATCTACTAGCA-1</th><td>1714.873</td><td>1207.7546</td><td>AAACAATCTACTAGCA-1</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 3 × 3</caption>
<thead>
	<tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACAACGAATAGTTC-1</th><td>533.5070</td><td>459.0799</td><td>AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>AAACAAGTATCTCCCA-1</th><td>220.4514</td><td>151.7882</td><td>AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>AAACAATCTACTAGCA-1</th><td>514.4618</td><td>362.3264</td><td>AAACAATCTACTAGCA-1</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 3 × 3</caption>
<thead>
	<tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACAACGAATAGTTC-1</th><td>30730</td><td>26443</td><td>AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>AAACAAGTATCTCCCA-1</th><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>AAACAATCTACTAGCA-1</th><td>29633</td><td>20870</td><td>AAACAATCTACTAGCA-1</td></tr>
</tbody>
</table>




```R
# reset
#va_coords <- va_temp
#vb_coords <- vb_temp

# copy
#va_temp <- va_coords
#vb_temp <- vb_coords

va1_coords$orig.ident <- 'VA1'
vb1_coords$orig.ident <- 'VB1'
va2_coords$orig.ident <- 'VA2'

va1_coords$Barcode2 <- va1_coords$Barcode
vb1_coords$Barcode2 <- vb1_coords$Barcode
va2_coords$Barcode2 <- va2_coords$Barcode

va1_coords$Barcode <- paste(va1_coords$orig.ident, va1_coords$Barcode2, sep='_')
vb1_coords$Barcode <- paste(vb1_coords$orig.ident, vb1_coords$Barcode2, sep='_')
va2_coords$Barcode <- paste(va2_coords$orig.ident, va2_coords$Barcode2, sep='_')

va1_coords$Barcode3 <- va1_coords$Barcode2
vb1_coords$Barcode3 <- vb1_coords$Barcode2
va2_coords$Barcode3 <- va2_coords$Barcode2

va1_coords$Barcode2 <- va1_coords$Barcode
vb1_coords$Barcode2 <- vb1_coords$Barcode
va2_coords$Barcode2 <- va2_coords$Barcode

va1_coords <- va1_coords %>% remove_rownames() %>% column_to_rownames('Barcode2') %>% select(-c('Barcode3','orig.ident'))
vb1_coords <- vb1_coords %>% remove_rownames() %>% column_to_rownames('Barcode2') %>% select(-c('Barcode3','orig.ident'))
va2_coords <- va2_coords %>% remove_rownames() %>% column_to_rownames('Barcode2') %>% select(-c('Barcode3','orig.ident'))

va1_coords[1:3,]
vb1_coords[1:3,]
va2_coords[1:3,]
```


<table class="dataframe">
<caption>A data.frame: 3 × 4</caption>
<thead>
	<tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>Barcode</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>30730</td><td>26443</td><td>AAACAACGAATAGTTC-1</td><td>VA1_AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>29633</td><td>20870</td><td>AAACAATCTACTAGCA-1</td><td>VA1_AAACAATCTACTAGCA-1</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 3 × 4</caption>
<thead>
	<tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>Barcode</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VB1_AAACAACGAATAGTTC-1</th><td>31283</td><td>26488</td><td>AAACAACGAATAGTTC-1</td><td>VB1_AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>VB1_AAACAAGTATCTCCCA-1</th><td>13243</td><td> 8788</td><td>AAACAAGTATCTCCCA-1</td><td>VB1_AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>VB1_AAACACCAATAACTGC-1</th><td>10072</td><td>25947</td><td>AAACACCAATAACTGC-1</td><td>VB1_AAACACCAATAACTGC-1</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 3 × 4</caption>
<thead>
	<tr><th></th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>Barcode</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA2_AAACAAGTATCTCCCA-1</th><td> 7764.086</td><td>20534.804</td><td>AAACAAGTATCTCCCA-1</td><td>VA2_AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>VA2_AAACACCAATAACTGC-1</th><td>24672.634</td><td>23846.305</td><td>AAACACCAATAACTGC-1</td><td>VA2_AAACACCAATAACTGC-1</td></tr>
	<tr><th scope=row>VA2_AAACAGAGCGACTCCT-1</th><td> 9484.343</td><td> 7768.555</td><td>AAACAGAGCGACTCCT-1</td><td>VA2_AAACAGAGCGACTCCT-1</td></tr>
</tbody>
</table>



# merge with metadata


```R
temp_meta_AB1_A2 <- AB1_A2@meta.data
temp_meta_va1 <- temp_meta_AB1_A2[temp_meta_AB1_A2$orig.ident == 'VA1',]
temp_meta_vb1 <- temp_meta_AB1_A2[temp_meta_AB1_A2$orig.ident == 'VB1',]
temp_meta_va2 <- temp_meta_AB1_A2[temp_meta_AB1_A2$orig.ident == 'VA2',]
```


```R
temp_meta_va1 <- merge(temp_meta_va1, va1_coords, by.x='Barcode', by.y='Barcode')
rownames(temp_meta_va1) <- temp_meta_va1$Barcode
temp_meta_va1[1:3,]
```


<table class="dataframe">
<caption>A data.frame: 3 × 60</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td> 0.6906627</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>30730</td><td>26443</td><td>AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.2739867</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>29633</td><td>20870</td><td>AAACAATCTACTAGCA-1</td></tr>
</tbody>
</table>




```R
temp_meta_vb1 <- merge(temp_meta_vb1, vb1_coords, by.x='Barcode', by.y='Barcode')
rownames(temp_meta_vb1) <- temp_meta_vb1$Barcode
temp_meta_vb1[1:3,]
```


<table class="dataframe">
<caption>A data.frame: 3 × 60</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VB1_AAACAACGAATAGTTC-1</th><td>VB1_AAACAACGAATAGTTC-1</td><td>VB1</td><td>VB1</td><td>7174</td><td>3783</td><td>round1</td><td>5068</td><td>3601</td><td>6</td><td>6</td><td>VB1</td><td>0.0000000</td><td>vasc_low </td><td>20</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low </td><td>Tumor_high_MHCI_low </td><td>vasc_neg </td><td>NA</td><td>NA</td><td>SCLC + Effector Immune Cells (T+NK)   </td><td>6_SCLC_Effector_Immune</td><td>6_SCLC_Effector_Immune</td><td>-0.4091901</td><td>-0.2934343</td><td> 0.01579317</td><td>-0.6232956</td><td>-0.2943380</td><td>-0.5389701</td><td> 0.1609967</td><td>-0.44028788</td><td>-0.04994617</td><td>-0.2194454</td><td>-0.18496858</td><td>0.01012297</td><td>0.00174647</td><td>0.52323585</td><td>-0.3835384</td><td> 0.2639777</td><td>-0.3175423</td><td>-0.18950681</td><td> 0.1715214</td><td>-0.5652826</td><td>-0.4481054</td><td>-0.3094859</td><td>-0.3636408</td><td>-0.3739748</td><td>-0.18496858</td><td>0.09003211</td><td>Low </td><td>Low              </td><td>vasc_low         </td><td>vasc_low</td><td>31283</td><td>26488</td><td>AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>VB1_AAACAAGTATCTCCCA-1</th><td>VB1_AAACAAGTATCTCCCA-1</td><td>VB1</td><td>VB1</td><td>3817</td><td>2574</td><td>round1</td><td>3824</td><td>2572</td><td>1</td><td>1</td><td>VB1</td><td>3.1173698</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high </td><td>Tumor_low           </td><td>vasc_high</td><td>NA</td><td>NA</td><td>SCLC                                  </td><td>1_SCLC                </td><td>1_4_5_11_SCLC         </td><td> 0.1248427</td><td> 0.5082822</td><td>-0.37090396</td><td>-0.2081911</td><td>-0.1350080</td><td> 0.2548888</td><td> 0.3546991</td><td> 0.61804117</td><td> 0.17432327</td><td> 0.6448965</td><td> 0.22323805</td><td>0.06597860</td><td>0.01068036</td><td>0.64566927</td><td> 0.9555956</td><td>-0.3258627</td><td> 0.2594864</td><td>-0.09579307</td><td>-0.2807481</td><td>-0.6112452</td><td> 0.4179383</td><td>-0.3445947</td><td>-0.3024936</td><td>-0.2643529</td><td> 0.22323805</td><td>0.69416330</td><td>Low </td><td>Low              </td><td>vasc_high_MN4_low</td><td>vasc_low</td><td>13243</td><td> 8788</td><td>AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>VB1_AAACACCAATAACTGC-1</th><td>VB1_AAACACCAATAACTGC-1</td><td>VB1</td><td>VB1</td><td>6137</td><td>3622</td><td>round1</td><td>4990</td><td>3585</td><td>0</td><td>0</td><td>VB1</td><td>0.9667787</td><td>vasc_low </td><td>17</td><td>Tumor_high</td><td>Tumor_high</td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_high_MHCI_high</td><td>Tumor_high_MHCI_high</td><td>vasc_low </td><td>NA</td><td>NA</td><td>Lymphocytes B Cells (and Plasma Cells)</td><td>0_Lymphocytes_B_Cells </td><td>0_Lymphocytes_B_Cells </td><td>-0.2232842</td><td> 0.1240651</td><td> 0.14827892</td><td>-0.4825369</td><td>-0.2947074</td><td>-0.3440408</td><td>-0.2352929</td><td> 0.08386593</td><td>-0.56773770</td><td>-0.2392742</td><td>-0.01450242</td><td>0.04427482</td><td>0.10580519</td><td>0.06204871</td><td>-0.3869387</td><td>-0.4880361</td><td>-0.4476507</td><td> 0.07267653</td><td> 0.3163856</td><td> 0.4089414</td><td> 0.2907777</td><td>-0.2117291</td><td>-0.3645462</td><td>-0.3794759</td><td>-0.01450242</td><td>0.14351748</td><td>High</td><td>Above80Percentile</td><td>vasc_low         </td><td>vasc_low</td><td>10072</td><td>25947</td><td>AAACACCAATAACTGC-1</td></tr>
</tbody>
</table>




```R
temp_meta_va2 <- merge(temp_meta_va2, va2_coords, by.x='Barcode', by.y='Barcode')
rownames(temp_meta_va2) <- temp_meta_va2$Barcode
temp_meta_va2[1:3,]
```


<table class="dataframe">
<caption>A data.frame: 3 × 60</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA2_AAACAAGTATCTCCCA-1</th><td>VA2_AAACAAGTATCTCCCA-1</td><td>VA2</td><td>VA2</td><td>45867</td><td>8515</td><td>round2</td><td>24897</td><td>7041</td><td>5</td><td>5</td><td>VA2</td><td>1.2963355</td><td>vasc_low</td><td>23</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low</td><td>NA</td><td>NA</td><td>SCLC                                  </td><td>5_SCLC               </td><td>1_4_5_11_SCLC        </td><td>-0.4509680</td><td>-0.4483304</td><td>-0.6914766</td><td>-0.3743552</td><td>-0.6934548</td><td>-0.8061806</td><td>-0.4244671</td><td>-0.2619432</td><td>-0.6846815</td><td>-0.7654531</td><td>-0.1999233</td><td>-0.1339729</td><td>-0.2137573</td><td>-0.2902687</td><td>0.02930414</td><td>-0.2612448</td><td>-0.7413616</td><td>-0.3953259</td><td>-0.2975777</td><td> 0.1823718</td><td>-0.5831718</td><td> 0.484040526</td><td>-0.3325986</td><td>-0.7049820</td><td>-0.1999233</td><td>-0.1565646</td><td>Low </td><td>Below20Percentile</td><td>vasc_low</td><td>vasc_low</td><td> 7764.086</td><td>20534.804</td><td>AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>VA2_AAACACCAATAACTGC-1</th><td>VA2_AAACACCAATAACTGC-1</td><td>VA2</td><td>VA2</td><td>11849</td><td>5082</td><td>round2</td><td>22823</td><td>5881</td><td>0</td><td>0</td><td>VA2</td><td>1.2238234</td><td>vasc_low</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>6</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low          </td><td>vasc_low</td><td>NA</td><td>NA</td><td>Lymphocytes B Cells (and Plasma Cells)</td><td>0_Lymphocytes_B_Cells</td><td>0_Lymphocytes_B_Cells</td><td> 0.4657553</td><td> 0.4768693</td><td> 0.8017860</td><td> 0.6856195</td><td> 0.7281010</td><td>-0.4407441</td><td>-0.3554181</td><td> 0.4510906</td><td> 0.6792392</td><td> 0.5655572</td><td> 0.3014952</td><td> 0.0995948</td><td> 0.1687630</td><td> 0.3944400</td><td>0.27930871</td><td> 0.3967318</td><td> 0.5340971</td><td> 0.3684659</td><td> 0.6188102</td><td> 0.6471133</td><td> 0.1671480</td><td>-0.005985246</td><td> 0.4556002</td><td> 0.6415331</td><td> 0.3014952</td><td> 0.1536071</td><td>High</td><td>Above80Percentile</td><td>vasc_low</td><td>vasc_low</td><td>24672.634</td><td>23846.305</td><td>AAACACCAATAACTGC-1</td></tr>
	<tr><th scope=row>VA2_AAACAGAGCGACTCCT-1</th><td>VA2_AAACAGAGCGACTCCT-1</td><td>VA2</td><td>VA2</td><td>43834</td><td>8662</td><td>round2</td><td>24885</td><td>7403</td><td>5</td><td>5</td><td>VA2</td><td>0.2054955</td><td>vasc_low</td><td>24</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low</td><td>NA</td><td>NA</td><td>SCLC                                  </td><td>5_SCLC               </td><td>1_4_5_11_SCLC        </td><td>-0.5878032</td><td>-0.6471249</td><td>-0.2034801</td><td>-0.6089440</td><td>-0.6035510</td><td>-0.8401840</td><td>-0.3283846</td><td>-0.5893354</td><td>-0.6599931</td><td>-0.7992599</td><td>-0.4414986</td><td>-0.2771161</td><td>-0.3490007</td><td>-0.6546928</td><td>0.15123322</td><td>-0.5942471</td><td>-0.6668816</td><td>-0.4782615</td><td>-0.6769066</td><td>-0.8178636</td><td>-0.6608763</td><td> 0.140619235</td><td>-0.7686918</td><td>-0.7469988</td><td>-0.4414986</td><td>-0.1267453</td><td>Low </td><td>Below20Percentile</td><td>vasc_low</td><td>vasc_low</td><td> 9484.343</td><td> 7768.555</td><td>AAACAGAGCGACTCCT-1</td></tr>
</tbody>
</table>



# save for neighborhood analysis in python


```R
write.csv(temp_meta_va1, 'Processing/10_Neighborhood_Analysis/va1_nn_pre.csv', row.names=FALSE)
write.csv(temp_meta_vb1, 'Processing/10_Neighborhood_Analysis/vb1_nn_pre.csv', row.names=FALSE)
write.csv(temp_meta_va2, 'Processing/10_Neighborhood_Analysis/va2_nn_pre.csv', row.names=FALSE)
```


```R

```
