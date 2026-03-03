# Visualize Neighborhoods - for the Tumor_MHCI spots


```R
# tidyverse includes: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats. 
# https://www.tidyverse.org/
library(tidyverse)
library(Seurat, quietly = T)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
library(ggthemes)
#library(DESeq2)
#library(tximport)
#library(tximeta)
#library(EnhancedVolcano)
library(pheatmap)
library(RColorBrewer)
#library(devtools)
#library(ashr)
#library(enrichplot)
library(fgsea)
library(msigdbr)
#library(org.Hs.eg.db)
#library(AnnotationDbi)
library(openxlsx)
#library(optparse)

library(dittoSeq, quietly = T)
library(scales, quietly = T)
library(GGally, quietly = T)

options(repr.matrix.max.cols=1000, repr.matrix.max.rows=100)
```


```R
# Load seurat objects
#AB_norm <- readRDS('Processing_2024/data_AB_TumorMHCI_TumorMHCILabel_20240419.rds')
#VA_norm <- readRDS('Processing_2024/data_VA_TumorMHCI_TumorMHCILabel_20240419.rds')

data <- readRDS('Processing/AB1_A2_Annotated_05_do_GSVA_by_Tissue_Slice_20250425.rds')
```


```R
# Load neighborhood data
nns_all <- read_csv('Processing/10_Neighborhood_Analysis/AB1_A2_Neighborhoods.csv')
nns_6adj <- read_csv('Processing/10_Neighborhood_Analysis/AB1_A2_Neighborhoods_All_Six_Adjacent.csv')
```

    [1m[22mNew names:
    [36m•[39m `` -> `...1`
    [1mRows: [22m[34m9711[39m [1mColumns: [22m[34m142[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m  (40): Barcode, Tissue_Slice, orig.ident, visium_round, slice_ident, vas...
    [32mdbl[39m (102): ...1, nCount_Spatial, nFeature_Spatial, nCount_SCT, nFeature_SCT,...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.
    [1m[22mNew names:
    [36m•[39m `` -> `...1`
    [1mRows: [22m[34m8522[39m [1mColumns: [22m[34m142[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m  (40): Barcode, Tissue_Slice, orig.ident, visium_round, slice_ident, vas...
    [32mdbl[39m (102): ...1, nCount_Spatial, nFeature_Spatial, nCount_SCT, nFeature_SCT,...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.



```R
# drop unnecessary first column and set rownames
nns_all <- nns_all %>% select(-c('...1')) %>% as.data.frame()
nns_6adj <- nns_6adj %>% select(-c('...1')) %>% as.data.frame()

rownames(nns_all) <- nns_all$Barcode
rownames(nns_6adj) <- nns_6adj$Barcode
```


```R
data@meta.data %>% head(3)
nns_all %>% head(3)
nns_6adj %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 53</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td> 0.6906627</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.2739867</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 3 × 141</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_mean</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_mean</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_mean</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_mean</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_mean</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_mean</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_mean</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_mean</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_mean</th><th scope=col>MN4_EC_Phenotype_nn_vals</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td> 0.6906627</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>30730</td><td>26443</td><td>AAACAACGAATAGTTC-1</td><td>VA1_0</td><td>VA1_854 </td><td>VA1_1863</td><td>VA1_2089</td><td>VA1_517 </td><td>VA1_2692</td><td>VA1_2472</td><td>3</td><td>['VA1_854', 'VA1_1863', 'VA1_2089', 'VA1_517', 'VA1_2692', 'VA1_2472'] </td><td>-0.1403735</td><td>-0.68346143</td><td>-0.6824274</td><td>-0.02534835</td><td>-0.1830202</td><td>-0.6898109</td><td>-0.4007403</td><td>[-0.140373503285981, -0.683461428931692, -0.682427363760641, -0.0253483482945374, -0.183020236470598, -0.689810896482086] </td><td> 0.4787432</td><td>-0.11469175</td><td>-0.05304243</td><td>-0.07726181</td><td>0.05871469</td><td> 0.1634258</td><td> 0.07598128</td><td>[0.478743199137973, -0.114691753402848, -0.053042433947209, -0.0772618094476313, 0.058714690391722, 0.163425797244824]  </td><td>0.4985932</td><td>-0.1559936</td><td>-0.09335601</td><td>-0.1178707</td><td>-0.2441883</td><td>-0.2408739</td><td>-0.05894823</td><td>[0.498593160799533, -0.155993596157572, -0.0933560136080898, -0.11787072243338, -0.244188290075738, -0.240873901729021]</td><td>0.04912155</td><td> 0.46869062</td><td>0.1456998</td><td> 0.4895150</td><td>-0.05807054</td><td>0.42619319</td><td>0.25352493</td><td>[0.0491215471122242, 0.468690624713802, 0.145699754382207, 0.489514984067094, -0.058070538457045, 0.42619318624764]  </td><td>0.18087518</td><td>0.4527944</td><td> 0.21247564</td><td> 0.42031179</td><td> 0.2671539</td><td> 0.5003215</td><td> 0.3389888</td><td>[0.180875180364392, 0.452794442072525, 0.212475639001522, 0.420311792896499, 0.267153935558, 0.500321548802953]         </td><td>-0.0428909507</td><td> 0.1220257</td><td>-0.5131311</td><td>-0.53971745</td><td>-0.38422937</td><td> 0.003417859</td><td>-0.2257542</td><td>[-0.0428909507251363, 0.122025706475451, -0.513131142574183, -0.539717446893455, -0.384229368789918, 0.0034178586983256] </td><td>-0.1024639</td><td> 0.6480175</td><td>-0.2458688</td><td> 0.39703261</td><td>-0.2214054</td><td> 0.62847798</td><td> 0.1839650</td><td>[-0.102463880569944, 0.648017492909985, -0.245868803204754, 0.397032612718896, -0.221405399324788, 0.62847797680395]  </td><td>-0.1681126</td><td>-0.1728810</td><td>-0.07891578</td><td> 0.9548053</td><td>-0.2141428</td><td>-0.2022404</td><td> 0.0197521</td><td>[-0.168112618939413, -0.172881026928877, -0.0789157831565245, 0.954805319001111, -0.214142828565577, -0.202240448089484]</td><td> 0.08796249</td><td> 0.04550749</td><td>-0.1515810</td><td>-0.06665673</td><td> 0.08863924</td><td> 0.1215510</td><td> 0.02090375</td><td>[0.0879624949143438, 0.0455074885754431, -0.151580958192362, -0.0666567297841213, 0.0886392350786018, 0.121550954065541] </td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284 </td><td>VA1_771 </td><td>VA1_352 </td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']  </td><td>-0.3313744</td><td>-0.25345755</td><td>-0.3396933</td><td>-0.70457534</td><td>-0.1946746</td><td>-0.6988578</td><td>-0.4204388</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]   </td><td>-0.3655433</td><td> 0.08738646</td><td>-0.01253212</td><td>-0.19925824</td><td>0.07594644</td><td>-0.1103883</td><td>-0.08739819</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.33109866</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>-0.17753110</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.54673283</td><td>-0.11670278</td><td>0.2619499</td><td> 0.4510973</td><td> 0.38571992</td><td>0.07931599</td><td>0.26801886</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073] </td><td>0.55563396</td><td>0.2090081</td><td>-0.14647071</td><td> 0.07915905</td><td> 0.6203140</td><td>-0.3870134</td><td> 0.1551052</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]   </td><td>-0.1416772948</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.52396701</td><td>-0.548042038</td><td>-0.3423706</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td> 0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td> 0.14113379</td><td>-0.4395605</td><td>-0.06779542</td><td>-0.1273297</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.32306461</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>-0.2725464</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.14347061</td><td>-0.3029000</td><td>-0.24465014</td><td>-0.35970213</td><td>-0.2758603</td><td>-0.22987461</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td> 6</td><td> 6</td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.2739867</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>29633</td><td>20870</td><td>AAACAATCTACTAGCA-1</td><td>VA1_2</td><td>VA1_2207</td><td>VA1_1663</td><td>VA1_1116</td><td>VA1_1210</td><td>VA1_920 </td><td>VA1_1124</td><td>4</td><td>['VA1_2207', 'VA1_1663', 'VA1_1116', 'VA1_1210', 'VA1_920', 'VA1_1124']</td><td>-0.1021611</td><td>-0.09404125</td><td>-0.7143306</td><td>-0.09243213</td><td>-0.1110161</td><td>-0.2277834</td><td>-0.2236274</td><td>[-0.102161065707172, -0.0940412510144283, -0.714330573226588, -0.0924321318329605, -0.111016118158986, -0.227783356007723]</td><td>-0.1750505</td><td> 0.16211636</td><td> 0.37237171</td><td> 0.51536446</td><td>0.14442029</td><td>-0.3628057</td><td> 0.10940276</td><td>[-0.175050506937197, 0.162116360694768, 0.372371708563783, 0.515364456905993, 0.144420288347746, -0.362805745678248]    </td><td>0.5555935</td><td> 0.5465806</td><td> 0.32440788</td><td>-0.2791675</td><td> 0.5603602</td><td>-0.2980788</td><td> 0.23494930</td><td>[0.555593475471297, 0.546580604053197, 0.324407882945447, -0.279167500500132, 0.560360187021303, -0.298078847308211]   </td><td>0.02470210</td><td> 0.02862245</td><td>0.3301684</td><td>-0.1450918</td><td>-0.05224927</td><td>0.39084903</td><td>0.09616682</td><td>[0.02470210181234, 0.0286224547429464, 0.330168428714374, -0.145091843512385, -0.0522492652549711, 0.390849032936518]</td><td>0.07100433</td><td>0.2218239</td><td>-0.06996729</td><td>-0.54888330</td><td>-0.5257730</td><td>-0.4789928</td><td>-0.2217980</td><td>[0.0710043326473108, 0.221823916630173, -0.0699672876994677, -0.548883300106269, -0.525772957160213, -0.478992784712071]</td><td> 0.0002582967</td><td>-0.5967590</td><td> 0.3099307</td><td>-0.37202321</td><td>-0.03360161</td><td> 0.083418241</td><td>-0.1014628</td><td>[0.0002582966608133, -0.596759018177968, 0.309930663526349, -0.372023213928184, -0.0336016075053905, 0.0834182412782388] </td><td>-0.1388668</td><td> 0.1185314</td><td>-0.4495049</td><td>-0.07106897</td><td>-0.4148401</td><td> 0.31776205</td><td>-0.1063312</td><td>[-0.138866804624008, 0.118531428276318, -0.44950493543207, -0.0710689737394324, -0.414840143712403, 0.317762048735336]</td><td> 0.9090881</td><td>-0.1755351</td><td> 0.78090959</td><td>-0.2698540</td><td>-0.2240448</td><td>-0.2925585</td><td> 0.1213342</td><td>[0.90908809327848, -0.175535107021276, 0.78090959121468, -0.26985397079401, -0.224044808961654, -0.292558511702187]     </td><td>-0.20592704</td><td>-0.18814347</td><td>-0.1059103</td><td>-0.01350095</td><td>-0.09839253</td><td> 0.1376595</td><td>-0.07903581</td><td>[-0.205927037792345, -0.188143473994638, -0.105910293614525, -0.0135009514713148, -0.0983925267325698, 0.137659452545575]</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 3 × 141</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_mean</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_mean</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_mean</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_mean</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_mean</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_mean</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_mean</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_mean</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_mean</th><th scope=col>MN4_EC_Phenotype_nn_vals</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6</td><td>6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low</td><td>vasc_low </td><td>SCLC            </td><td>SCLC</td><td>SCLC + Effector Immune Cells (T+NK)</td><td>6_SCLC_Effector_Immune</td><td>6_SCLC_Effector_Immune</td><td>-0.4320402</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.3802492</td><td>-0.2729638</td><td>-0.6165697</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.53207656</td><td>-0.2841903</td><td>-0.4421974</td><td>-0.20627210</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.4421974</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284 </td><td>VA1_771 </td><td>VA1_352 </td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']  </td><td>-0.3313744</td><td>-0.2534575</td><td>-0.3396933</td><td>-0.7045753</td><td>-0.1946746</td><td>-0.6988578</td><td>-0.4204388</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165] </td><td>-0.3655433</td><td> 0.08738646</td><td>-0.01253212</td><td>-0.1992582</td><td> 0.07594644</td><td>-0.1103883</td><td>-0.08739819</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td> 0.1532445</td><td>-0.2668601</td><td>-0.3310987</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>-0.17753110</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td> 0.54673283</td><td>-0.1167028</td><td> 0.2619499</td><td>0.4510973</td><td> 0.38571992</td><td> 0.07931599</td><td>0.2680189</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073]  </td><td> 0.5556340</td><td> 0.2090081</td><td>-0.1464707</td><td> 0.07915905</td><td> 0.6203139711</td><td>-0.3870134</td><td> 0.1551052</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]   </td><td>-0.1416773</td><td>-0.607476142</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.5239670</td><td>-0.548042038</td><td>-0.3423706</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td> 0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td> 0.14113379</td><td>-0.4395605</td><td>-0.06779542</td><td>-0.1273297</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]   </td><td>-0.3917784</td><td>-0.2591518</td><td>-0.3230646</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>-0.27254643</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.14347061</td><td>-0.302899967</td><td>-0.24465014</td><td>-0.3597021</td><td>-0.27586034</td><td>-0.2298746</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1_AAACACCAATAACTGC-1</td><td>VA1</td><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td>4</td><td>4</td><td>VA1</td><td>0.907884</td><td>vasc_low </td><td> 9</td><td>Tumor_low</td><td>Tumor_low</td><td>1</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low</td><td>vasc_low </td><td>SCLC+vasculature</td><td>SCLC</td><td>SCLC                               </td><td>4_SCLC                </td><td>1_4_5_11_SCLC         </td><td>-0.2949967</td><td>-0.1592542</td><td>-0.6246499</td><td>-0.7091880</td><td>-0.3126876</td><td>-0.1285324</td><td> 0.5479671</td><td>-0.2873825</td><td> 0.05564459</td><td>-0.3094928</td><td>-0.2340697</td><td>-0.03504189</td><td>-0.1584708</td><td>-0.2611306</td><td>-0.3051610</td><td> 0.3381445</td><td>-0.3300223</td><td>-0.1418120</td><td>-0.1902565</td><td>-0.7206696</td><td>-0.2190269</td><td>-0.45491178</td><td>-0.3728931</td><td>-0.4603460</td><td>-0.2340697</td><td> 0.07776597</td><td>Low </td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td> 9524</td><td>25898</td><td>AAACACCAATAACTGC-1</td><td>VA1_3</td><td>VA1_1747</td><td>VA1_425 </td><td>VA1_1286</td><td>VA1_1385</td><td>VA1_2332</td><td>VA1_2567</td><td>6</td><td>['VA1_1747', 'VA1_425', 'VA1_1286', 'VA1_1385', 'VA1_2332', 'VA1_2567']</td><td>-0.7080698</td><td>-0.1727255</td><td>-0.3059416</td><td>-0.7100941</td><td>-0.7212701</td><td>-0.7231039</td><td>-0.5568675</td><td>[-0.708069838814006, -0.172725483080612, -0.305941553061541, -0.710094101019417, -0.721270084194434, -0.723103920735993]</td><td>-0.2331998</td><td> 0.07857171</td><td> 0.47593594</td><td>-0.3569468</td><td>-0.37384910</td><td>-0.3794862</td><td>-0.13149570</td><td>[-0.23319977872425, 0.0785717073459754, 0.475935940981545, -0.356946754827358, -0.373849101912343, -0.379486233065752]  </td><td>-0.2393405</td><td> 0.3552654</td><td> 0.2957953</td><td>-0.2810686</td><td> 0.3142241</td><td>-0.3293976</td><td> 0.01924632</td><td>[-0.239340532575186, 0.355265387687536, 0.295795295662033, -0.281068641184542, 0.314224061582033, -0.329397638582968]  </td><td>-0.01724305</td><td> 0.5798389</td><td> 0.5448870</td><td>0.3572084</td><td> 0.06235329</td><td>-0.25844957</td><td>0.2114325</td><td>[-0.0172430474396086, 0.579838912374754, 0.544886975743813, 0.357208398080687, 0.0623532857737059, -0.258449565629918]</td><td>-0.6919329</td><td>-0.8139959</td><td>-0.8053769</td><td>-0.59253595</td><td>-0.6808669993</td><td>-0.8014868</td><td>-0.7310326</td><td>[-0.69193292191416, -0.813995861950977, -0.805376914068887, -0.592535945325797, -0.680866999335305, -0.801486804482884] </td><td>-0.5362091</td><td> 0.002455842</td><td> 0.3126581</td><td>-0.15292795</td><td>-0.6217053</td><td>-0.625378821</td><td>-0.2701845</td><td>[-0.536209058824282, 0.0024558419296655, 0.31265811387002, -0.152927946116749, -0.621705278122557, -0.625378820848625]   </td><td>-0.4001988</td><td>-0.2774540</td><td> 0.1709900</td><td>-0.28708766</td><td>-0.4553370</td><td>-0.33424428</td><td>-0.2638886</td><td>[-0.400198776367015, -0.277454023705374, 0.170989973541542, -0.287087656671685, -0.455337042372491, -0.334244282995581]  </td><td>-0.2081416</td><td>-0.2570514</td><td>-0.2973595</td><td>-0.2734547</td><td>-0.2881576</td><td>-0.3242649</td><td>-0.27473828</td><td>[-0.20814162832553, -0.25705141028191, -0.297359471894224, -0.273454690938038, -0.288157631526152, -0.324264852970433]  </td><td>-0.24604744</td><td>-0.08491335</td><td>-0.007453716</td><td> 0.02875985</td><td>-0.3204698</td><td>-0.31904223</td><td>-0.1581944</td><td>[-0.246047435280206, -0.084913354050493, -0.0074537156101513, 0.028759847156343, -0.320469802379713, -0.319042226520557] </td></tr>
	<tr><th scope=row>VA1_AAACAGCTTTCAGAAG-1</th><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>VA1</td><td>5132</td><td>3219</td><td>round1</td><td>4797</td><td>3218</td><td>4</td><td>4</td><td>VA1</td><td>3.751279</td><td>vasc_high</td><td>14</td><td>Tumor_low</td><td>Tumor_low</td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low</td><td>vasc_high</td><td>SCLC+vasculature</td><td>SCLC</td><td>SCLC                               </td><td>4_SCLC                </td><td>1_4_5_11_SCLC         </td><td> 0.1486568</td><td> 0.6701079</td><td> 0.1810290</td><td>-0.6826419</td><td> 0.5112075</td><td> 0.3512003</td><td> 0.2988667</td><td> 0.5369628</td><td>-0.22259553</td><td>-0.2723653</td><td>-0.2440354</td><td>-0.02474443</td><td>-0.0993218</td><td> 0.5413331</td><td> 0.7253772</td><td>-0.4164351</td><td> 0.1290473</td><td> 0.1598434</td><td> 0.4337860</td><td> 0.6304333</td><td> 0.5986082</td><td>-0.09760237</td><td> 0.8719730</td><td> 0.9183918</td><td>-0.2440354</td><td>-0.15730040</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>15282</td><td>27943</td><td>AAACAGCTTTCAGAAG-1</td><td>VA1_4</td><td>VA1_964 </td><td>VA1_1350</td><td>VA1_2187</td><td>VA1_320 </td><td>VA1_885 </td><td>VA1_1441</td><td>6</td><td>['VA1_964', 'VA1_1350', 'VA1_2187', 'VA1_320', 'VA1_885', 'VA1_1441']  </td><td>-0.2747007</td><td>-0.3598579</td><td>-0.7267706</td><td>-0.2614986</td><td>-0.3453775</td><td>-0.1007939</td><td>-0.3448332</td><td>[-0.274700715316639, -0.35985792380019, -0.726770556323605, -0.261498631836216, -0.345377508433023, -0.1007939422355]   </td><td> 0.6469571</td><td>-0.37709477</td><td>-0.38767968</td><td>-0.3679927</td><td>-0.38363109</td><td>-0.1733069</td><td>-0.17379134</td><td>[0.646957089052535, -0.377094771708025, -0.387679680538107, -0.367992732636768, -0.383631090919492, -0.173306877574485] </td><td> 0.7398643</td><td> 0.2864618</td><td>-0.3462077</td><td> 0.3057052</td><td>-0.3378027</td><td> 0.5554457</td><td> 0.20057777</td><td>[0.739864288138109, 0.286461826481311, -0.346207724634595, 0.305705163870573, -0.337802681608782, 0.555445748014179]   </td><td> 0.56832168</td><td>-0.3177625</td><td>-0.2912681</td><td>0.5579031</td><td>-0.27516198</td><td> 0.51495656</td><td>0.1261648</td><td>[0.568321679326659, -0.317762474933888, -0.291268134566025, 0.557903093968171, -0.275161976240716, 0.514956561063913] </td><td>-0.7979204</td><td>-0.7024213</td><td>-0.7951722</td><td>-0.49489641</td><td>-0.0002819505</td><td>-0.1563592</td><td>-0.4911752</td><td>[-0.79792035031997, -0.702421272421103, -0.795172162092652, -0.494896413487619, -0.0002819505313625, -0.156359175593681]</td><td>-0.5093777</td><td>-0.158298382</td><td>-0.6316153</td><td> 0.53755664</td><td>-0.1075707</td><td> 0.003694195</td><td>-0.1442685</td><td>[-0.50937767833692, -0.158298382085461, -0.631615307337922, 0.537556639866296, -0.107570690365855, 0.0036941951266885]   </td><td>-0.3089724</td><td>-0.1347193</td><td>-0.1689729</td><td>-0.05015917</td><td>-0.3393764</td><td>-0.14524716</td><td>-0.1912412</td><td>[-0.308972369596385, -0.134719262393709, -0.168972855426038, -0.0501591707515464, -0.339376448811238, -0.145247155185077]</td><td> 0.7439191</td><td>-0.3464693</td><td>-0.3373675</td><td>-0.2891578</td><td> 0.7462068</td><td>-0.1784357</td><td> 0.05644928</td><td>[0.743919132154805, -0.346469293858606, -0.337367473494535, -0.28915783156616, 0.746206813940613, -0.178435687137299]   </td><td>-0.01687242</td><td>-0.12972981</td><td>-0.239731861</td><td>-0.09475526</td><td>-0.1666009</td><td> 0.02857231</td><td>-0.1031863</td><td>[-0.0168724150989142, -0.129729812963318, -0.239731860944855, -0.0947552583169846, -0.166600945093969, 0.02857231209327] </td></tr>
</tbody>
</table>




```R
#vasc_MN4_label
unique(nns_all$Tumor_label)
unique(nns_all$Tumor_MHCI_label)
unique(nns_all$Tumor_MHCI_label2)
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
<ol class=list-inline><li>'Tumor_low_MHCI_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_low_MHCI_high'</li><li>'Tumor_high_MHCI_high'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_high_MHCI_high'</li></ol>




```R
nns_all %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 141</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_mean</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_mean</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_mean</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_mean</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_mean</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_mean</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_mean</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_mean</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_mean</th><th scope=col>MN4_EC_Phenotype_nn_vals</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td> 0.6906627</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>30730</td><td>26443</td><td>AAACAACGAATAGTTC-1</td><td>VA1_0</td><td>VA1_854 </td><td>VA1_1863</td><td>VA1_2089</td><td>VA1_517 </td><td>VA1_2692</td><td>VA1_2472</td><td>3</td><td>['VA1_854', 'VA1_1863', 'VA1_2089', 'VA1_517', 'VA1_2692', 'VA1_2472'] </td><td>-0.1403735</td><td>-0.68346143</td><td>-0.6824274</td><td>-0.02534835</td><td>-0.1830202</td><td>-0.6898109</td><td>-0.4007403</td><td>[-0.140373503285981, -0.683461428931692, -0.682427363760641, -0.0253483482945374, -0.183020236470598, -0.689810896482086] </td><td> 0.4787432</td><td>-0.11469175</td><td>-0.05304243</td><td>-0.07726181</td><td>0.05871469</td><td> 0.1634258</td><td> 0.07598128</td><td>[0.478743199137973, -0.114691753402848, -0.053042433947209, -0.0772618094476313, 0.058714690391722, 0.163425797244824]  </td><td>0.4985932</td><td>-0.1559936</td><td>-0.09335601</td><td>-0.1178707</td><td>-0.2441883</td><td>-0.2408739</td><td>-0.05894823</td><td>[0.498593160799533, -0.155993596157572, -0.0933560136080898, -0.11787072243338, -0.244188290075738, -0.240873901729021]</td><td>0.04912155</td><td> 0.46869062</td><td>0.1456998</td><td> 0.4895150</td><td>-0.05807054</td><td>0.42619319</td><td>0.25352493</td><td>[0.0491215471122242, 0.468690624713802, 0.145699754382207, 0.489514984067094, -0.058070538457045, 0.42619318624764]  </td><td>0.18087518</td><td>0.4527944</td><td> 0.21247564</td><td> 0.42031179</td><td> 0.2671539</td><td> 0.5003215</td><td> 0.3389888</td><td>[0.180875180364392, 0.452794442072525, 0.212475639001522, 0.420311792896499, 0.267153935558, 0.500321548802953]         </td><td>-0.0428909507</td><td> 0.1220257</td><td>-0.5131311</td><td>-0.53971745</td><td>-0.38422937</td><td> 0.003417859</td><td>-0.2257542</td><td>[-0.0428909507251363, 0.122025706475451, -0.513131142574183, -0.539717446893455, -0.384229368789918, 0.0034178586983256] </td><td>-0.1024639</td><td> 0.6480175</td><td>-0.2458688</td><td> 0.39703261</td><td>-0.2214054</td><td> 0.62847798</td><td> 0.1839650</td><td>[-0.102463880569944, 0.648017492909985, -0.245868803204754, 0.397032612718896, -0.221405399324788, 0.62847797680395]  </td><td>-0.1681126</td><td>-0.1728810</td><td>-0.07891578</td><td> 0.9548053</td><td>-0.2141428</td><td>-0.2022404</td><td> 0.0197521</td><td>[-0.168112618939413, -0.172881026928877, -0.0789157831565245, 0.954805319001111, -0.214142828565577, -0.202240448089484]</td><td> 0.08796249</td><td> 0.04550749</td><td>-0.1515810</td><td>-0.06665673</td><td> 0.08863924</td><td> 0.1215510</td><td> 0.02090375</td><td>[0.0879624949143438, 0.0455074885754431, -0.151580958192362, -0.0666567297841213, 0.0886392350786018, 0.121550954065541] </td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284 </td><td>VA1_771 </td><td>VA1_352 </td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']  </td><td>-0.3313744</td><td>-0.25345755</td><td>-0.3396933</td><td>-0.70457534</td><td>-0.1946746</td><td>-0.6988578</td><td>-0.4204388</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]   </td><td>-0.3655433</td><td> 0.08738646</td><td>-0.01253212</td><td>-0.19925824</td><td>0.07594644</td><td>-0.1103883</td><td>-0.08739819</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.33109866</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>-0.17753110</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.54673283</td><td>-0.11670278</td><td>0.2619499</td><td> 0.4510973</td><td> 0.38571992</td><td>0.07931599</td><td>0.26801886</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073] </td><td>0.55563396</td><td>0.2090081</td><td>-0.14647071</td><td> 0.07915905</td><td> 0.6203140</td><td>-0.3870134</td><td> 0.1551052</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]   </td><td>-0.1416772948</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.52396701</td><td>-0.548042038</td><td>-0.3423706</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td> 0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td> 0.14113379</td><td>-0.4395605</td><td>-0.06779542</td><td>-0.1273297</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.32306461</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>-0.2725464</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.14347061</td><td>-0.3029000</td><td>-0.24465014</td><td>-0.35970213</td><td>-0.2758603</td><td>-0.22987461</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td> 6</td><td> 6</td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.2739867</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>29633</td><td>20870</td><td>AAACAATCTACTAGCA-1</td><td>VA1_2</td><td>VA1_2207</td><td>VA1_1663</td><td>VA1_1116</td><td>VA1_1210</td><td>VA1_920 </td><td>VA1_1124</td><td>4</td><td>['VA1_2207', 'VA1_1663', 'VA1_1116', 'VA1_1210', 'VA1_920', 'VA1_1124']</td><td>-0.1021611</td><td>-0.09404125</td><td>-0.7143306</td><td>-0.09243213</td><td>-0.1110161</td><td>-0.2277834</td><td>-0.2236274</td><td>[-0.102161065707172, -0.0940412510144283, -0.714330573226588, -0.0924321318329605, -0.111016118158986, -0.227783356007723]</td><td>-0.1750505</td><td> 0.16211636</td><td> 0.37237171</td><td> 0.51536446</td><td>0.14442029</td><td>-0.3628057</td><td> 0.10940276</td><td>[-0.175050506937197, 0.162116360694768, 0.372371708563783, 0.515364456905993, 0.144420288347746, -0.362805745678248]    </td><td>0.5555935</td><td> 0.5465806</td><td> 0.32440788</td><td>-0.2791675</td><td> 0.5603602</td><td>-0.2980788</td><td> 0.23494930</td><td>[0.555593475471297, 0.546580604053197, 0.324407882945447, -0.279167500500132, 0.560360187021303, -0.298078847308211]   </td><td>0.02470210</td><td> 0.02862245</td><td>0.3301684</td><td>-0.1450918</td><td>-0.05224927</td><td>0.39084903</td><td>0.09616682</td><td>[0.02470210181234, 0.0286224547429464, 0.330168428714374, -0.145091843512385, -0.0522492652549711, 0.390849032936518]</td><td>0.07100433</td><td>0.2218239</td><td>-0.06996729</td><td>-0.54888330</td><td>-0.5257730</td><td>-0.4789928</td><td>-0.2217980</td><td>[0.0710043326473108, 0.221823916630173, -0.0699672876994677, -0.548883300106269, -0.525772957160213, -0.478992784712071]</td><td> 0.0002582967</td><td>-0.5967590</td><td> 0.3099307</td><td>-0.37202321</td><td>-0.03360161</td><td> 0.083418241</td><td>-0.1014628</td><td>[0.0002582966608133, -0.596759018177968, 0.309930663526349, -0.372023213928184, -0.0336016075053905, 0.0834182412782388] </td><td>-0.1388668</td><td> 0.1185314</td><td>-0.4495049</td><td>-0.07106897</td><td>-0.4148401</td><td> 0.31776205</td><td>-0.1063312</td><td>[-0.138866804624008, 0.118531428276318, -0.44950493543207, -0.0710689737394324, -0.414840143712403, 0.317762048735336]</td><td> 0.9090881</td><td>-0.1755351</td><td> 0.78090959</td><td>-0.2698540</td><td>-0.2240448</td><td>-0.2925585</td><td> 0.1213342</td><td>[0.90908809327848, -0.175535107021276, 0.78090959121468, -0.26985397079401, -0.224044808961654, -0.292558511702187]     </td><td>-0.20592704</td><td>-0.18814347</td><td>-0.1059103</td><td>-0.01350095</td><td>-0.09839253</td><td> 0.1376595</td><td>-0.07903581</td><td>[-0.205927037792345, -0.188143473994638, -0.105910293614525, -0.0135009514713148, -0.0983925267325698, 0.137659452545575]</td></tr>
</tbody>
</table>



# Add labels for plotting


```R
# First layer: all negative, grey
nns_all$tumor_plot_label <- 'Tumor_Neg_or_Low'

# set VascHigh_MN4Low_Neighbor and VascHigh_MN4High_Neighbor columns
nns_all$TumorHigh_MHCILow_Neighbor <- 'Not_Neighbor'
nns_all$TumorHigh_MHCIHigh_Neighbor <- 'Not_Neighbor'

# Second layer: neighbors, colored by GSVA Enrichment
# get list of vasc_high spots. Use "id" column
tumor_high_id_list <- unique(nns_all$id[nns_all$Tumor_label == 'Tumor_high'])
# for each id in vasc_high_id_list, look at its neighbor ids. then, go and set vasc_plot_label to VascHigh_Neighbor
# and set Neighbor or Not_Neighbor in two other columns, VascHigh_MN4Low_Neighbor and VascHigh_MN4High_Neighbor

for (curr_id in tumor_high_id_list) {
    # skip if the center spot doesn't have 6 adjacent nearest neighbors
    if(nns_all$num_adjacent_nns[nns_all$id == curr_id] < 6) next
    
    # get neighbor ids fur the current center spot id
    curr_nn1_id <- nns_all$nn1_id[nns_all$id == curr_id]
    curr_nn2_id <- nns_all$nn2_id[nns_all$id == curr_id]
    curr_nn3_id <- nns_all$nn3_id[nns_all$id == curr_id]
    curr_nn4_id <- nns_all$nn4_id[nns_all$id == curr_id]
    curr_nn5_id <- nns_all$nn5_id[nns_all$id == curr_id]
    curr_nn6_id <- nns_all$nn6_id[nns_all$id == curr_id]

    curr_neighbor_id_list <- c(curr_nn1_id, curr_nn2_id, curr_nn3_id, curr_nn4_id, curr_nn5_id, curr_nn6_id)
    # set the tumor_plot_label to Neighbor for the spot ids which are neighbors of the current center spot
    # Also set TumorHigh_MN4Low_Neighbor and TumorHigh_MN4High_Neighbor columns
    for (curr_neighbor_id in curr_neighbor_id_list) {
        nns_all$tumor_plot_label[nns_all$id == curr_neighbor_id] <- 'TumorHigh_Neighbor'

        # if curr_id is a VascHigh_MN4Low or VascHigh_MN4High spot, set accordingly...
        if (nns_all$Tumor_MHCI_label[nns_all$id == curr_id] == 'Tumor_high_MHCI_low') {
            nns_all$TumorHigh_MHCILow_Neighbor[nns_all$id == curr_neighbor_id] <- 'Neighbor'
        }
        if (nns_all$Tumor_MHCI_label[nns_all$id == curr_id] == 'Tumor_high_MHCI_high') {
            nns_all$TumorHigh_MHCIHigh_Neighbor[nns_all$id == curr_neighbor_id] <- 'Neighbor'
        }
    }
}

# re-loop for the center spots so they don't get covered up if they occurred earlier
for (curr_id in tumor_high_id_list) {
    # this is another loop cycling through all vasc_high ids...
    # Third layer: center spot... Vasc High, MN4 Low or High. white
    # Even if one of these spots is already labeled as a neighbor, the fact that it's
    # vasc_high takes priority!

    # set tumor_plot_label
    # if curr_id is a VascHigh_MN4Low or VascHigh_MN4High spot, set accordingly...
    if (nns_all$Tumor_MHCI_label[nns_all$id == curr_id] == 'Tumor_high_MHCI_low') {
        nns_all$tumor_plot_label[nns_all$id == curr_id] <- 'TumorHigh_MHCILow'
    }
    if (nns_all$Tumor_MHCI_label[nns_all$id == curr_id] == 'Tumor_high_MHCI_high') {
        nns_all$tumor_plot_label[nns_all$id == curr_id] <- 'TumorHigh_MHCIHigh'
    }
}


```


```R
nns_all %>% head(5)
```


<table class="dataframe">
<caption>A data.frame: 5 × 144</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_mean</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_mean</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_mean</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_mean</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_mean</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_mean</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_mean</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_mean</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_mean</th><th scope=col>MN4_EC_Phenotype_nn_vals</th><th scope=col>tumor_plot_label</th><th scope=col>TumorHigh_MHCILow_Neighbor</th><th scope=col>TumorHigh_MHCIHigh_Neighbor</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.15846527</td><td> 0.6906627</td><td>-0.04649791</td><td> 0.11080007</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>30730</td><td>26443</td><td>AAACAACGAATAGTTC-1</td><td>VA1_0</td><td>VA1_854 </td><td>VA1_1863</td><td>VA1_2089</td><td>VA1_517 </td><td>VA1_2692</td><td>VA1_2472</td><td>3</td><td>['VA1_854', 'VA1_1863', 'VA1_2089', 'VA1_517', 'VA1_2692', 'VA1_2472'] </td><td>-0.1403735</td><td>-0.68346143</td><td>-0.6824274</td><td>-0.02534835</td><td>-0.1830202</td><td>-0.6898109</td><td>-0.4007403</td><td>[-0.140373503285981, -0.683461428931692, -0.682427363760641, -0.0253483482945374, -0.183020236470598, -0.689810896482086] </td><td> 0.4787432</td><td>-0.11469175</td><td>-0.05304243</td><td>-0.07726181</td><td> 0.05871469</td><td> 0.1634258</td><td> 0.07598128</td><td>[0.478743199137973, -0.114691753402848, -0.053042433947209, -0.0772618094476313, 0.058714690391722, 0.163425797244824]  </td><td> 0.4985932</td><td>-0.1559936</td><td>-0.09335601</td><td>-0.1178707</td><td>-0.2441883</td><td>-0.2408739</td><td>-0.05894823</td><td>[0.498593160799533, -0.155993596157572, -0.0933560136080898, -0.11787072243338, -0.244188290075738, -0.240873901729021]</td><td> 0.04912155</td><td> 0.46869062</td><td> 0.1456998</td><td> 0.4895150</td><td>-0.05807054</td><td> 0.42619319</td><td>0.25352493</td><td>[0.0491215471122242, 0.468690624713802, 0.145699754382207, 0.489514984067094, -0.058070538457045, 0.42619318624764]   </td><td> 0.18087518</td><td> 0.4527944</td><td> 0.21247564</td><td> 0.42031179</td><td> 0.2671539356</td><td> 0.5003215</td><td> 0.3389888</td><td>[0.180875180364392, 0.452794442072525, 0.212475639001522, 0.420311792896499, 0.267153935558, 0.500321548802953]         </td><td>-0.0428909507</td><td> 0.122025706</td><td>-0.5131311</td><td>-0.53971745</td><td>-0.38422937</td><td> 0.003417859</td><td>-0.2257542</td><td>[-0.0428909507251363, 0.122025706475451, -0.513131142574183, -0.539717446893455, -0.384229368789918, 0.0034178586983256] </td><td>-0.1024639</td><td> 0.6480175</td><td>-0.2458688</td><td> 0.39703261</td><td>-0.2214054</td><td> 0.62847798</td><td> 0.1839650</td><td>[-0.102463880569944, 0.648017492909985, -0.245868803204754, 0.397032612718896, -0.221405399324788, 0.62847797680395]     </td><td>-0.1681126</td><td>-0.1728810</td><td>-0.07891578</td><td> 0.9548053</td><td>-0.2141428</td><td>-0.2022404</td><td> 0.01975210</td><td>[-0.168112618939413, -0.172881026928877, -0.0789157831565245, 0.954805319001111, -0.214142828565577, -0.202240448089484]</td><td> 0.08796249</td><td> 0.04550749</td><td>-0.151580958</td><td>-0.06665673</td><td> 0.08863924</td><td> 0.12155095</td><td> 0.02090375</td><td>[0.0879624949143438, 0.0455074885754431, -0.151580958192362, -0.0666567297841213, 0.0886392350786018, 0.121550954065541] </td><td>Tumor_Neg_or_Low  </td><td>Not_Neighbor</td><td>Not_Neighbor</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.53207656</td><td>-0.2841903</td><td>-0.44219743</td><td>-0.20627210</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284 </td><td>VA1_771 </td><td>VA1_352 </td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']  </td><td>-0.3313744</td><td>-0.25345755</td><td>-0.3396933</td><td>-0.70457534</td><td>-0.1946746</td><td>-0.6988578</td><td>-0.4204388</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]   </td><td>-0.3655433</td><td> 0.08738646</td><td>-0.01253212</td><td>-0.19925824</td><td> 0.07594644</td><td>-0.1103883</td><td>-0.08739819</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td> 0.1532445</td><td>-0.2668601</td><td>-0.33109866</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>-0.17753110</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td> 0.54673283</td><td>-0.11670278</td><td> 0.2619499</td><td> 0.4510973</td><td> 0.38571992</td><td> 0.07931599</td><td>0.26801886</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073]  </td><td> 0.55563396</td><td> 0.2090081</td><td>-0.14647071</td><td> 0.07915905</td><td> 0.6203139711</td><td>-0.3870134</td><td> 0.1551052</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]   </td><td>-0.1416772948</td><td>-0.607476142</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.52396701</td><td>-0.548042038</td><td>-0.3423706</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td> 0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td> 0.14113379</td><td>-0.4395605</td><td>-0.06779542</td><td>-0.1273297</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]   </td><td>-0.3917784</td><td>-0.2591518</td><td>-0.32306461</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>-0.27254643</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.14347061</td><td>-0.302899967</td><td>-0.24465014</td><td>-0.35970213</td><td>-0.27586034</td><td>-0.22987461</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td><td>TumorHigh_Neighbor</td><td>Neighbor    </td><td>Not_Neighbor</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td> 6</td><td> 6</td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.51946615</td><td>-0.2739867</td><td>-0.29691635</td><td>-0.11225417</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>29633</td><td>20870</td><td>AAACAATCTACTAGCA-1</td><td>VA1_2</td><td>VA1_2207</td><td>VA1_1663</td><td>VA1_1116</td><td>VA1_1210</td><td>VA1_920 </td><td>VA1_1124</td><td>4</td><td>['VA1_2207', 'VA1_1663', 'VA1_1116', 'VA1_1210', 'VA1_920', 'VA1_1124']</td><td>-0.1021611</td><td>-0.09404125</td><td>-0.7143306</td><td>-0.09243213</td><td>-0.1110161</td><td>-0.2277834</td><td>-0.2236274</td><td>[-0.102161065707172, -0.0940412510144283, -0.714330573226588, -0.0924321318329605, -0.111016118158986, -0.227783356007723]</td><td>-0.1750505</td><td> 0.16211636</td><td> 0.37237171</td><td> 0.51536446</td><td> 0.14442029</td><td>-0.3628057</td><td> 0.10940276</td><td>[-0.175050506937197, 0.162116360694768, 0.372371708563783, 0.515364456905993, 0.144420288347746, -0.362805745678248]    </td><td> 0.5555935</td><td> 0.5465806</td><td> 0.32440788</td><td>-0.2791675</td><td> 0.5603602</td><td>-0.2980788</td><td> 0.23494930</td><td>[0.555593475471297, 0.546580604053197, 0.324407882945447, -0.279167500500132, 0.560360187021303, -0.298078847308211]   </td><td> 0.02470210</td><td> 0.02862245</td><td> 0.3301684</td><td>-0.1450918</td><td>-0.05224927</td><td> 0.39084903</td><td>0.09616682</td><td>[0.02470210181234, 0.0286224547429464, 0.330168428714374, -0.145091843512385, -0.0522492652549711, 0.390849032936518] </td><td> 0.07100433</td><td> 0.2218239</td><td>-0.06996729</td><td>-0.54888330</td><td>-0.5257729572</td><td>-0.4789928</td><td>-0.2217980</td><td>[0.0710043326473108, 0.221823916630173, -0.0699672876994677, -0.548883300106269, -0.525772957160213, -0.478992784712071]</td><td> 0.0002582967</td><td>-0.596759018</td><td> 0.3099307</td><td>-0.37202321</td><td>-0.03360161</td><td> 0.083418241</td><td>-0.1014628</td><td>[0.0002582966608133, -0.596759018177968, 0.309930663526349, -0.372023213928184, -0.0336016075053905, 0.0834182412782388] </td><td>-0.1388668</td><td> 0.1185314</td><td>-0.4495049</td><td>-0.07106897</td><td>-0.4148401</td><td> 0.31776205</td><td>-0.1063312</td><td>[-0.138866804624008, 0.118531428276318, -0.44950493543207, -0.0710689737394324, -0.414840143712403, 0.317762048735336]   </td><td> 0.9090881</td><td>-0.1755351</td><td> 0.78090959</td><td>-0.2698540</td><td>-0.2240448</td><td>-0.2925585</td><td> 0.12133421</td><td>[0.90908809327848, -0.175535107021276, 0.78090959121468, -0.26985397079401, -0.224044808961654, -0.292558511702187]     </td><td>-0.20592704</td><td>-0.18814347</td><td>-0.105910294</td><td>-0.01350095</td><td>-0.09839253</td><td> 0.13765945</td><td>-0.07903581</td><td>[-0.205927037792345, -0.188143473994638, -0.105910293614525, -0.0135009514713148, -0.0983925267325698, 0.137659452545575]</td><td>TumorHigh_MHCILow </td><td>Not_Neighbor</td><td>Neighbor    </td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1_AAACACCAATAACTGC-1</td><td>VA1</td><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td> 4</td><td> 4</td><td>VA1</td><td>0.907884</td><td>vasc_low </td><td> 9</td><td>Tumor_low </td><td>Tumor_low </td><td>1</td><td>MHCI_low </td><td>MHCI_low </td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC+vasculature      </td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td>-0.29499671</td><td>-0.1592542</td><td>-0.6246499</td><td>-0.70918797</td><td>-0.3126876</td><td>-0.128532387</td><td> 0.54796712</td><td>-0.28738253</td><td> 0.05564459</td><td>-0.3094928</td><td>-0.23406975</td><td>-0.03504189</td><td>-0.1584708</td><td>-0.2611306</td><td>-0.3051610</td><td> 0.3381445</td><td>-0.3300223</td><td>-0.14181196</td><td>-0.19025652</td><td>-0.7206696</td><td>-0.21902693</td><td>-0.45491178</td><td>-0.3728931</td><td>-0.4603460</td><td>-0.23406975</td><td> 0.07776597</td><td>Low </td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td> 9524</td><td>25898</td><td>AAACACCAATAACTGC-1</td><td>VA1_3</td><td>VA1_1747</td><td>VA1_425 </td><td>VA1_1286</td><td>VA1_1385</td><td>VA1_2332</td><td>VA1_2567</td><td>6</td><td>['VA1_1747', 'VA1_425', 'VA1_1286', 'VA1_1385', 'VA1_2332', 'VA1_2567']</td><td>-0.7080698</td><td>-0.17272548</td><td>-0.3059416</td><td>-0.71009410</td><td>-0.7212701</td><td>-0.7231039</td><td>-0.5568675</td><td>[-0.708069838814006, -0.172725483080612, -0.305941553061541, -0.710094101019417, -0.721270084194434, -0.723103920735993]  </td><td>-0.2331998</td><td> 0.07857171</td><td> 0.47593594</td><td>-0.35694675</td><td>-0.37384910</td><td>-0.3794862</td><td>-0.13149570</td><td>[-0.23319977872425, 0.0785717073459754, 0.475935940981545, -0.356946754827358, -0.373849101912343, -0.379486233065752]  </td><td>-0.2393405</td><td> 0.3552654</td><td> 0.29579530</td><td>-0.2810686</td><td> 0.3142241</td><td>-0.3293976</td><td> 0.01924632</td><td>[-0.239340532575186, 0.355265387687536, 0.295795295662033, -0.281068641184542, 0.314224061582033, -0.329397638582968]  </td><td>-0.01724305</td><td> 0.57983891</td><td> 0.5448870</td><td> 0.3572084</td><td> 0.06235329</td><td>-0.25844957</td><td>0.21143249</td><td>[-0.0172430474396086, 0.579838912374754, 0.544886975743813, 0.357208398080687, 0.0623532857737059, -0.258449565629918]</td><td>-0.69193292</td><td>-0.8139959</td><td>-0.80537691</td><td>-0.59253595</td><td>-0.6808669993</td><td>-0.8014868</td><td>-0.7310326</td><td>[-0.69193292191416, -0.813995861950977, -0.805376914068887, -0.592535945325797, -0.680866999335305, -0.801486804482884] </td><td>-0.5362090588</td><td> 0.002455842</td><td> 0.3126581</td><td>-0.15292795</td><td>-0.62170528</td><td>-0.625378821</td><td>-0.2701845</td><td>[-0.536209058824282, 0.0024558419296655, 0.31265811387002, -0.152927946116749, -0.621705278122557, -0.625378820848625]   </td><td>-0.4001988</td><td>-0.2774540</td><td> 0.1709900</td><td>-0.28708766</td><td>-0.4553370</td><td>-0.33424428</td><td>-0.2638886</td><td>[-0.400198776367015, -0.277454023705374, 0.170989973541542, -0.287087656671685, -0.455337042372491, -0.334244282995581]  </td><td>-0.2081416</td><td>-0.2570514</td><td>-0.29735947</td><td>-0.2734547</td><td>-0.2881576</td><td>-0.3242649</td><td>-0.27473828</td><td>[-0.20814162832553, -0.25705141028191, -0.297359471894224, -0.273454690938038, -0.288157631526152, -0.324264852970433]  </td><td>-0.24604744</td><td>-0.08491335</td><td>-0.007453716</td><td> 0.02875985</td><td>-0.32046980</td><td>-0.31904223</td><td>-0.15819445</td><td>[-0.246047435280206, -0.084913354050493, -0.0074537156101513, 0.028759847156343, -0.320469802379713, -0.319042226520557] </td><td>TumorHigh_Neighbor</td><td>Neighbor    </td><td>Not_Neighbor</td></tr>
	<tr><th scope=row>VA1_AAACAGCTTTCAGAAG-1</th><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>VA1</td><td>5132</td><td>3219</td><td>round1</td><td>4797</td><td>3218</td><td> 4</td><td> 4</td><td>VA1</td><td>3.751279</td><td>vasc_high</td><td>14</td><td>Tumor_low </td><td>Tumor_low </td><td>4</td><td>MHCI_high</td><td>MHCI_high</td><td>Tumor_low_MHCI_high</td><td>Tumor_low          </td><td>vasc_high</td><td>SCLC+vasculature      </td><td>SCLC                  </td><td>SCLC                                </td><td>4_SCLC                        </td><td>1_4_5_11_SCLC         </td><td> 0.14865678</td><td> 0.6701079</td><td> 0.1810290</td><td>-0.68264193</td><td> 0.5112075</td><td> 0.351200349</td><td> 0.29886669</td><td> 0.53696284</td><td>-0.22259553</td><td>-0.2723653</td><td>-0.24403536</td><td>-0.02474443</td><td>-0.0993218</td><td> 0.5413331</td><td> 0.7253772</td><td>-0.4164351</td><td> 0.1290473</td><td> 0.15984342</td><td> 0.43378600</td><td> 0.6304333</td><td> 0.59860816</td><td>-0.09760237</td><td> 0.8719730</td><td> 0.9183918</td><td>-0.24403536</td><td>-0.15730040</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>15282</td><td>27943</td><td>AAACAGCTTTCAGAAG-1</td><td>VA1_4</td><td>VA1_964 </td><td>VA1_1350</td><td>VA1_2187</td><td>VA1_320 </td><td>VA1_885 </td><td>VA1_1441</td><td>6</td><td>['VA1_964', 'VA1_1350', 'VA1_2187', 'VA1_320', 'VA1_885', 'VA1_1441']  </td><td>-0.2747007</td><td>-0.35985792</td><td>-0.7267706</td><td>-0.26149863</td><td>-0.3453775</td><td>-0.1007939</td><td>-0.3448332</td><td>[-0.274700715316639, -0.35985792380019, -0.726770556323605, -0.261498631836216, -0.345377508433023, -0.1007939422355]     </td><td> 0.6469571</td><td>-0.37709477</td><td>-0.38767968</td><td>-0.36799273</td><td>-0.38363109</td><td>-0.1733069</td><td>-0.17379134</td><td>[0.646957089052535, -0.377094771708025, -0.387679680538107, -0.367992732636768, -0.383631090919492, -0.173306877574485] </td><td> 0.7398643</td><td> 0.2864618</td><td>-0.34620772</td><td> 0.3057052</td><td>-0.3378027</td><td> 0.5554457</td><td> 0.20057777</td><td>[0.739864288138109, 0.286461826481311, -0.346207724634595, 0.305705163870573, -0.337802681608782, 0.555445748014179]   </td><td> 0.56832168</td><td>-0.31776247</td><td>-0.2912681</td><td> 0.5579031</td><td>-0.27516198</td><td> 0.51495656</td><td>0.12616479</td><td>[0.568321679326659, -0.317762474933888, -0.291268134566025, 0.557903093968171, -0.275161976240716, 0.514956561063913] </td><td>-0.79792035</td><td>-0.7024213</td><td>-0.79517216</td><td>-0.49489641</td><td>-0.0002819505</td><td>-0.1563592</td><td>-0.4911752</td><td>[-0.79792035031997, -0.702421272421103, -0.795172162092652, -0.494896413487619, -0.0002819505313625, -0.156359175593681]</td><td>-0.5093776783</td><td>-0.158298382</td><td>-0.6316153</td><td> 0.53755664</td><td>-0.10757069</td><td> 0.003694195</td><td>-0.1442685</td><td>[-0.50937767833692, -0.158298382085461, -0.631615307337922, 0.537556639866296, -0.107570690365855, 0.0036941951266885]   </td><td>-0.3089724</td><td>-0.1347193</td><td>-0.1689729</td><td>-0.05015917</td><td>-0.3393764</td><td>-0.14524716</td><td>-0.1912412</td><td>[-0.308972369596385, -0.134719262393709, -0.168972855426038, -0.0501591707515464, -0.339376448811238, -0.145247155185077]</td><td> 0.7439191</td><td>-0.3464693</td><td>-0.33736747</td><td>-0.2891578</td><td> 0.7462068</td><td>-0.1784357</td><td> 0.05644928</td><td>[0.743919132154805, -0.346469293858606, -0.337367473494535, -0.28915783156616, 0.746206813940613, -0.178435687137299]   </td><td>-0.01687242</td><td>-0.12972981</td><td>-0.239731861</td><td>-0.09475526</td><td>-0.16660095</td><td> 0.02857231</td><td>-0.10318633</td><td>[-0.0168724150989142, -0.129729812963318, -0.239731860944855, -0.0947552583169846, -0.166600945093969, 0.02857231209327] </td><td>TumorHigh_Neighbor</td><td>Neighbor    </td><td>Neighbor    </td></tr>
</tbody>
</table>




```R
# counts...
dplyr::count(nns_all, tumor_plot_label)
dplyr::count(nns_all, TumorHigh_MHCILow_Neighbor)
dplyr::count(nns_all, TumorHigh_MHCIHigh_Neighbor)
```


<table class="dataframe">
<caption>A data.frame: 4 × 2</caption>
<thead>
	<tr><th scope=col>tumor_plot_label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>TumorHigh_MHCIHigh</td><td>1903</td></tr>
	<tr><td>TumorHigh_MHCILow </td><td>3580</td></tr>
	<tr><td>TumorHigh_Neighbor</td><td>2365</td></tr>
	<tr><td>Tumor_Neg_or_Low  </td><td>1863</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 2 × 2</caption>
<thead>
	<tr><th scope=col>TumorHigh_MHCILow_Neighbor</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Neighbor    </td><td>6765</td></tr>
	<tr><td>Not_Neighbor</td><td>2946</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 2 × 2</caption>
<thead>
	<tr><th scope=col>TumorHigh_MHCIHigh_Neighbor</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>Neighbor    </td><td>5369</td></tr>
	<tr><td>Not_Neighbor</td><td>4342</td></tr>
</tbody>
</table>



# Now we have labels for each spot. Now we need to collect the GSVA Enrichment of interest and plot that for the TumorHigh_Neighbor spots. 


```R
dim(data@meta.data)
dim(nns_all)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>53</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>144</li></ol>




```R
#set as metadata
data@meta.data <- nns_all
```


```R
options(repr.plot.width=15, repr.plot.height=12)
Idents(data) <- 'tumor_plot_label'
SpatialPlot(data, pt.size.factor=2.2, image.alpha = 0.9, stroke = 0.3)[[1]] + ggplot2::theme_minimal()
```


    
![png](output_15_0.png)
    



```R
options(repr.plot.width=15, repr.plot.height=12)
Idents(data) <- 'tumor_plot_label'
SpatialPlot(data, pt.size.factor=2.2, image.alpha = 0.9, stroke = 0.3)[[2]] + ggplot2::theme_minimal()
```


    
![png](output_16_0.png)
    



```R
options(repr.plot.width=15, repr.plot.height=12)
Idents(data) <- 'tumor_plot_label'
SpatialPlot(data, pt.size.factor=2.2, image.alpha = 0.9, stroke = 0.3)[[3]] + ggplot2::theme_minimal()
```


    
![png](output_17_0.png)
    



```R
data@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 144</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_mean</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_mean</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_mean</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_mean</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_mean</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_mean</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_mean</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_mean</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_mean</th><th scope=col>MN4_EC_Phenotype_nn_vals</th><th scope=col>tumor_plot_label</th><th scope=col>TumorHigh_MHCILow_Neighbor</th><th scope=col>TumorHigh_MHCIHigh_Neighbor</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td> 0.6906627</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>30730</td><td>26443</td><td>AAACAACGAATAGTTC-1</td><td>VA1_0</td><td>VA1_854 </td><td>VA1_1863</td><td>VA1_2089</td><td>VA1_517 </td><td>VA1_2692</td><td>VA1_2472</td><td>3</td><td>['VA1_854', 'VA1_1863', 'VA1_2089', 'VA1_517', 'VA1_2692', 'VA1_2472'] </td><td>-0.1403735</td><td>-0.68346143</td><td>-0.6824274</td><td>-0.02534835</td><td>-0.1830202</td><td>-0.6898109</td><td>-0.4007403</td><td>[-0.140373503285981, -0.683461428931692, -0.682427363760641, -0.0253483482945374, -0.183020236470598, -0.689810896482086] </td><td> 0.4787432</td><td>-0.11469175</td><td>-0.05304243</td><td>-0.07726181</td><td>0.05871469</td><td> 0.1634258</td><td> 0.07598128</td><td>[0.478743199137973, -0.114691753402848, -0.053042433947209, -0.0772618094476313, 0.058714690391722, 0.163425797244824]  </td><td>0.4985932</td><td>-0.1559936</td><td>-0.09335601</td><td>-0.1178707</td><td>-0.2441883</td><td>-0.2408739</td><td>-0.05894823</td><td>[0.498593160799533, -0.155993596157572, -0.0933560136080898, -0.11787072243338, -0.244188290075738, -0.240873901729021]</td><td>0.04912155</td><td> 0.46869062</td><td>0.1456998</td><td> 0.4895150</td><td>-0.05807054</td><td>0.42619319</td><td>0.25352493</td><td>[0.0491215471122242, 0.468690624713802, 0.145699754382207, 0.489514984067094, -0.058070538457045, 0.42619318624764]  </td><td>0.18087518</td><td>0.4527944</td><td> 0.21247564</td><td> 0.42031179</td><td> 0.2671539</td><td> 0.5003215</td><td> 0.3389888</td><td>[0.180875180364392, 0.452794442072525, 0.212475639001522, 0.420311792896499, 0.267153935558, 0.500321548802953]         </td><td>-0.0428909507</td><td> 0.1220257</td><td>-0.5131311</td><td>-0.53971745</td><td>-0.38422937</td><td> 0.003417859</td><td>-0.2257542</td><td>[-0.0428909507251363, 0.122025706475451, -0.513131142574183, -0.539717446893455, -0.384229368789918, 0.0034178586983256] </td><td>-0.1024639</td><td> 0.6480175</td><td>-0.2458688</td><td> 0.39703261</td><td>-0.2214054</td><td> 0.62847798</td><td> 0.1839650</td><td>[-0.102463880569944, 0.648017492909985, -0.245868803204754, 0.397032612718896, -0.221405399324788, 0.62847797680395]  </td><td>-0.1681126</td><td>-0.1728810</td><td>-0.07891578</td><td> 0.9548053</td><td>-0.2141428</td><td>-0.2022404</td><td> 0.0197521</td><td>[-0.168112618939413, -0.172881026928877, -0.0789157831565245, 0.954805319001111, -0.214142828565577, -0.202240448089484]</td><td> 0.08796249</td><td> 0.04550749</td><td>-0.1515810</td><td>-0.06665673</td><td> 0.08863924</td><td> 0.1215510</td><td> 0.02090375</td><td>[0.0879624949143438, 0.0455074885754431, -0.151580958192362, -0.0666567297841213, 0.0886392350786018, 0.121550954065541] </td><td>Tumor_Neg_or_Low  </td><td>Not_Neighbor</td><td>Not_Neighbor</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284 </td><td>VA1_771 </td><td>VA1_352 </td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']  </td><td>-0.3313744</td><td>-0.25345755</td><td>-0.3396933</td><td>-0.70457534</td><td>-0.1946746</td><td>-0.6988578</td><td>-0.4204388</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]   </td><td>-0.3655433</td><td> 0.08738646</td><td>-0.01253212</td><td>-0.19925824</td><td>0.07594644</td><td>-0.1103883</td><td>-0.08739819</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.33109866</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>-0.17753110</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.54673283</td><td>-0.11670278</td><td>0.2619499</td><td> 0.4510973</td><td> 0.38571992</td><td>0.07931599</td><td>0.26801886</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073] </td><td>0.55563396</td><td>0.2090081</td><td>-0.14647071</td><td> 0.07915905</td><td> 0.6203140</td><td>-0.3870134</td><td> 0.1551052</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]   </td><td>-0.1416772948</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.52396701</td><td>-0.548042038</td><td>-0.3423706</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td> 0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td> 0.14113379</td><td>-0.4395605</td><td>-0.06779542</td><td>-0.1273297</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.32306461</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>-0.2725464</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.14347061</td><td>-0.3029000</td><td>-0.24465014</td><td>-0.35970213</td><td>-0.2758603</td><td>-0.22987461</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td><td>TumorHigh_Neighbor</td><td>Neighbor    </td><td>Not_Neighbor</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td> 6</td><td> 6</td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.2739867</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>29633</td><td>20870</td><td>AAACAATCTACTAGCA-1</td><td>VA1_2</td><td>VA1_2207</td><td>VA1_1663</td><td>VA1_1116</td><td>VA1_1210</td><td>VA1_920 </td><td>VA1_1124</td><td>4</td><td>['VA1_2207', 'VA1_1663', 'VA1_1116', 'VA1_1210', 'VA1_920', 'VA1_1124']</td><td>-0.1021611</td><td>-0.09404125</td><td>-0.7143306</td><td>-0.09243213</td><td>-0.1110161</td><td>-0.2277834</td><td>-0.2236274</td><td>[-0.102161065707172, -0.0940412510144283, -0.714330573226588, -0.0924321318329605, -0.111016118158986, -0.227783356007723]</td><td>-0.1750505</td><td> 0.16211636</td><td> 0.37237171</td><td> 0.51536446</td><td>0.14442029</td><td>-0.3628057</td><td> 0.10940276</td><td>[-0.175050506937197, 0.162116360694768, 0.372371708563783, 0.515364456905993, 0.144420288347746, -0.362805745678248]    </td><td>0.5555935</td><td> 0.5465806</td><td> 0.32440788</td><td>-0.2791675</td><td> 0.5603602</td><td>-0.2980788</td><td> 0.23494930</td><td>[0.555593475471297, 0.546580604053197, 0.324407882945447, -0.279167500500132, 0.560360187021303, -0.298078847308211]   </td><td>0.02470210</td><td> 0.02862245</td><td>0.3301684</td><td>-0.1450918</td><td>-0.05224927</td><td>0.39084903</td><td>0.09616682</td><td>[0.02470210181234, 0.0286224547429464, 0.330168428714374, -0.145091843512385, -0.0522492652549711, 0.390849032936518]</td><td>0.07100433</td><td>0.2218239</td><td>-0.06996729</td><td>-0.54888330</td><td>-0.5257730</td><td>-0.4789928</td><td>-0.2217980</td><td>[0.0710043326473108, 0.221823916630173, -0.0699672876994677, -0.548883300106269, -0.525772957160213, -0.478992784712071]</td><td> 0.0002582967</td><td>-0.5967590</td><td> 0.3099307</td><td>-0.37202321</td><td>-0.03360161</td><td> 0.083418241</td><td>-0.1014628</td><td>[0.0002582966608133, -0.596759018177968, 0.309930663526349, -0.372023213928184, -0.0336016075053905, 0.0834182412782388] </td><td>-0.1388668</td><td> 0.1185314</td><td>-0.4495049</td><td>-0.07106897</td><td>-0.4148401</td><td> 0.31776205</td><td>-0.1063312</td><td>[-0.138866804624008, 0.118531428276318, -0.44950493543207, -0.0710689737394324, -0.414840143712403, 0.317762048735336]</td><td> 0.9090881</td><td>-0.1755351</td><td> 0.78090959</td><td>-0.2698540</td><td>-0.2240448</td><td>-0.2925585</td><td> 0.1213342</td><td>[0.90908809327848, -0.175535107021276, 0.78090959121468, -0.26985397079401, -0.224044808961654, -0.292558511702187]     </td><td>-0.20592704</td><td>-0.18814347</td><td>-0.1059103</td><td>-0.01350095</td><td>-0.09839253</td><td> 0.1376595</td><td>-0.07903581</td><td>[-0.205927037792345, -0.188143473994638, -0.105910293614525, -0.0135009514713148, -0.0983925267325698, 0.137659452545575]</td><td>TumorHigh_MHCILow </td><td>Not_Neighbor</td><td>Neighbor    </td></tr>
</tbody>
</table>



# Check on some colors....


```R
clrs <- rev(paletteer::paletteer_c("ggthemes::Red-Blue Diverging",n=100))
show_col(clrs)
```


    
![png](output_20_0.png)
    



```R
clrs[1]
clrs[50]
clrs[100]
```


    <colors>
    [37m[48;5;67m#2E5A87FF[49m[39m 



    <colors>
    [30m[48;5;224m#DCD3D1FF[49m[39m 



    <colors>
    [37m[48;5;125m#A90C38FF[49m[39m 



```R
c(clrs)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'#2E5A87FF'</li><li>'#305C89FF'</li><li>'#335F8CFF'</li><li>'#35618EFF'</li><li>'#376491FF'</li><li>'#396693FF'</li><li>'#3C6995FF'</li><li>'#3E6B98FF'</li><li>'#406E9AFF'</li><li>'#42709DFF'</li><li>'#45739FFF'</li><li>'#4775A1FF'</li><li>'#4978A4FF'</li><li>'#4B7AA6FF'</li><li>'#4D7DA9FF'</li><li>'#507FABFF'</li><li>'#5282AEFF'</li><li>'#5484B0FF'</li><li>'#5787B2FF'</li><li>'#5A89B4FF'</li><li>'#5C8CB7FF'</li><li>'#5F8EB9FF'</li><li>'#6191BBFF'</li><li>'#6494BDFF'</li><li>'#6796BFFF'</li><li>'#6999C1FF'</li><li>'#6C9BC4FF'</li><li>'#6E9EC6FF'</li><li>'#71A1C8FF'</li><li>'#74A3CAFF'</li><li>'#76A6CCFF'</li><li>'#79A9CFFF'</li><li>'#7BABD1FF'</li><li>'#7EAED3FF'</li><li>'#85B0D3FF'</li><li>'#8CB3D3FF'</li><li>'#92B5D3FF'</li><li>'#98B7D3FF'</li><li>'#9EB9D3FF'</li><li>'#A4BCD2FF'</li><li>'#AABED2FF'</li><li>'#B0C0D2FF'</li><li>'#B6C3D2FF'</li><li>'#BCC5D2FF'</li><li>'#C1C7D2FF'</li><li>'#C7CAD2FF'</li><li>'#CCCCD2FF'</li><li>'#D2CED1FF'</li><li>'#D7D0D1FF'</li><li>'#DCD3D1FF'</li><li>'#E0D2CEFF'</li><li>'#E3CDC7FF'</li><li>'#E6C8C1FF'</li><li>'#E8C3BAFF'</li><li>'#EABEB4FF'</li><li>'#ECB9ADFF'</li><li>'#EEB4A7FF'</li><li>'#F0AFA1FF'</li><li>'#F1AA9AFF'</li><li>'#F3A494FF'</li><li>'#F49F8EFF'</li><li>'#F59A88FF'</li><li>'#F69581FF'</li><li>'#F68F7BFF'</li><li>'#F78A75FF'</li><li>'#F8856FFF'</li><li>'#F87F69FF'</li><li>'#F77B67FF'</li><li>'#F57864FF'</li><li>'#F47462FF'</li><li>'#F3705FFF'</li><li>'#F16C5DFF'</li><li>'#F0695AFF'</li><li>'#EE6558FF'</li><li>'#ED6156FF'</li><li>'#EB5D53FF'</li><li>'#EA5951FF'</li><li>'#E8544FFF'</li><li>'#E7504CFF'</li><li>'#E54C4AFF'</li><li>'#E44748FF'</li><li>'#E24245FF'</li><li>'#E13E43FF'</li><li>'#DE3A42FF'</li><li>'#DB3741FF'</li><li>'#D83541FF'</li><li>'#D43340FF'</li><li>'#D1303FFF'</li><li>'#CD2E3FFF'</li><li>'#CA2B3EFF'</li><li>'#C7283EFF'</li><li>'#C3263DFF'</li><li>'#C0233CFF'</li><li>'#BD203CFF'</li><li>'#B91D3BFF'</li><li>'#B61A3BFF'</li><li>'#B3173AFF'</li><li>'#B01439FF'</li><li>'#AC1039FF'</li><li>'#A90C38FF'</li></ol>




```R
# set colors...
clrs <- rev(paletteer::paletteer_c("ggthemes::Red-Blue Diverging",n=100))

# gsva enrichment colors...
my_keys <- 1:100
mylist <- c(clrs)
names(mylist) <- my_keys

#other spots (Neg, center spots)
# 111 is grey, Neg spots
# 222 is light green, VascHigh_MN4Low spots
# 333 is darker green, VascHigh_MN4High spots
#other_spots_list <- c('#666568', '#79ff74', '#598d1d')
other_spots_list <- c('#666568', '#93d519', '#cdc60a')
other_spots_keys <- c(111, 222, 333)
names(other_spots_list) <- other_spots_keys

# combine the GSVA Enrichment colors with the three categoricals
test_colors <- c(mylist, other_spots_list)
test_colors
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>1</dt><dd>'#2E5A87FF'</dd><dt>2</dt><dd>'#305C89FF'</dd><dt>3</dt><dd>'#335F8CFF'</dd><dt>4</dt><dd>'#35618EFF'</dd><dt>5</dt><dd>'#376491FF'</dd><dt>6</dt><dd>'#396693FF'</dd><dt>7</dt><dd>'#3C6995FF'</dd><dt>8</dt><dd>'#3E6B98FF'</dd><dt>9</dt><dd>'#406E9AFF'</dd><dt>10</dt><dd>'#42709DFF'</dd><dt>11</dt><dd>'#45739FFF'</dd><dt>12</dt><dd>'#4775A1FF'</dd><dt>13</dt><dd>'#4978A4FF'</dd><dt>14</dt><dd>'#4B7AA6FF'</dd><dt>15</dt><dd>'#4D7DA9FF'</dd><dt>16</dt><dd>'#507FABFF'</dd><dt>17</dt><dd>'#5282AEFF'</dd><dt>18</dt><dd>'#5484B0FF'</dd><dt>19</dt><dd>'#5787B2FF'</dd><dt>20</dt><dd>'#5A89B4FF'</dd><dt>21</dt><dd>'#5C8CB7FF'</dd><dt>22</dt><dd>'#5F8EB9FF'</dd><dt>23</dt><dd>'#6191BBFF'</dd><dt>24</dt><dd>'#6494BDFF'</dd><dt>25</dt><dd>'#6796BFFF'</dd><dt>26</dt><dd>'#6999C1FF'</dd><dt>27</dt><dd>'#6C9BC4FF'</dd><dt>28</dt><dd>'#6E9EC6FF'</dd><dt>29</dt><dd>'#71A1C8FF'</dd><dt>30</dt><dd>'#74A3CAFF'</dd><dt>31</dt><dd>'#76A6CCFF'</dd><dt>32</dt><dd>'#79A9CFFF'</dd><dt>33</dt><dd>'#7BABD1FF'</dd><dt>34</dt><dd>'#7EAED3FF'</dd><dt>35</dt><dd>'#85B0D3FF'</dd><dt>36</dt><dd>'#8CB3D3FF'</dd><dt>37</dt><dd>'#92B5D3FF'</dd><dt>38</dt><dd>'#98B7D3FF'</dd><dt>39</dt><dd>'#9EB9D3FF'</dd><dt>40</dt><dd>'#A4BCD2FF'</dd><dt>41</dt><dd>'#AABED2FF'</dd><dt>42</dt><dd>'#B0C0D2FF'</dd><dt>43</dt><dd>'#B6C3D2FF'</dd><dt>44</dt><dd>'#BCC5D2FF'</dd><dt>45</dt><dd>'#C1C7D2FF'</dd><dt>46</dt><dd>'#C7CAD2FF'</dd><dt>47</dt><dd>'#CCCCD2FF'</dd><dt>48</dt><dd>'#D2CED1FF'</dd><dt>49</dt><dd>'#D7D0D1FF'</dd><dt>50</dt><dd>'#DCD3D1FF'</dd><dt>51</dt><dd>'#E0D2CEFF'</dd><dt>52</dt><dd>'#E3CDC7FF'</dd><dt>53</dt><dd>'#E6C8C1FF'</dd><dt>54</dt><dd>'#E8C3BAFF'</dd><dt>55</dt><dd>'#EABEB4FF'</dd><dt>56</dt><dd>'#ECB9ADFF'</dd><dt>57</dt><dd>'#EEB4A7FF'</dd><dt>58</dt><dd>'#F0AFA1FF'</dd><dt>59</dt><dd>'#F1AA9AFF'</dd><dt>60</dt><dd>'#F3A494FF'</dd><dt>61</dt><dd>'#F49F8EFF'</dd><dt>62</dt><dd>'#F59A88FF'</dd><dt>63</dt><dd>'#F69581FF'</dd><dt>64</dt><dd>'#F68F7BFF'</dd><dt>65</dt><dd>'#F78A75FF'</dd><dt>66</dt><dd>'#F8856FFF'</dd><dt>67</dt><dd>'#F87F69FF'</dd><dt>68</dt><dd>'#F77B67FF'</dd><dt>69</dt><dd>'#F57864FF'</dd><dt>70</dt><dd>'#F47462FF'</dd><dt>71</dt><dd>'#F3705FFF'</dd><dt>72</dt><dd>'#F16C5DFF'</dd><dt>73</dt><dd>'#F0695AFF'</dd><dt>74</dt><dd>'#EE6558FF'</dd><dt>75</dt><dd>'#ED6156FF'</dd><dt>76</dt><dd>'#EB5D53FF'</dd><dt>77</dt><dd>'#EA5951FF'</dd><dt>78</dt><dd>'#E8544FFF'</dd><dt>79</dt><dd>'#E7504CFF'</dd><dt>80</dt><dd>'#E54C4AFF'</dd><dt>81</dt><dd>'#E44748FF'</dd><dt>82</dt><dd>'#E24245FF'</dd><dt>83</dt><dd>'#E13E43FF'</dd><dt>84</dt><dd>'#DE3A42FF'</dd><dt>85</dt><dd>'#DB3741FF'</dd><dt>86</dt><dd>'#D83541FF'</dd><dt>87</dt><dd>'#D43340FF'</dd><dt>88</dt><dd>'#D1303FFF'</dd><dt>89</dt><dd>'#CD2E3FFF'</dd><dt>90</dt><dd>'#CA2B3EFF'</dd><dt>91</dt><dd>'#C7283EFF'</dd><dt>92</dt><dd>'#C3263DFF'</dd><dt>93</dt><dd>'#C0233CFF'</dd><dt>94</dt><dd>'#BD203CFF'</dd><dt>95</dt><dd>'#B91D3BFF'</dd><dt>96</dt><dd>'#B61A3BFF'</dd><dt>97</dt><dd>'#B3173AFF'</dd><dt>98</dt><dd>'#B01439FF'</dd><dt>99</dt><dd>'#AC1039FF'</dd><dt>100</dt><dd>'#A90C38FF'</dd><dt>111</dt><dd>'#666568'</dd><dt>222</dt><dd>'#93d519'</dd><dt>333</dt><dd>'#cdc60a'</dd></dl>




```R
other_spots_list
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>111</dt><dd>'#666568'</dd><dt>222</dt><dd>'#93d519'</dd><dt>333</dt><dd>'#cdc60a'</dd></dl>




```R
show_col(other_spots_list)
```


    
![png](output_25_0.png)
    

#install.packages('ch')# use package ch to use show_color function
# https://rdrr.io/cran/ch/man/ch-package.html
library(ch)col_legend1 <- show_color(other_spots_list, ncol=1, label=TRUE) +
theme_bw() +
labs(title='Center Spot Legend: grey = Vasc_Neg or Vasc_Low, Green = VascHigh_MN4Low, Yellow = VascHighMN4High')
col_legend1
#ggsave('/immuno/ian/Projects/2023_02_SCLC_10X_Visium/Seurat_Analysis/Neighborhood_Analysis_20240808/color_legend.png')
# Function to melt each GSVA Enrichment of interest, split into quantiles for easier color paletting, and then plotting


```R
# function to get colors for the categorical variables and the quantiled GSVA Enrichment scores
get_plotting_colors <- function() {
    # set colors...
    clrs <- rev(paletteer::paletteer_c("ggthemes::Red-Blue Diverging",n=100))
    # gsva enrichment colors...
    my_keys <- 1:100
    mylist <- c(clrs)
    names(mylist) <- my_keys
    #other spots (Neg, center spots)
    # 111 is grey, Neg spots
    # 222 is light greenish yellow, VascHigh_MN4Low spots #93d519
    # 333 is brighter yellow, VascHigh_MN4High spots #eee54f
   #other_spots_list <- c('#666568', '#7b9b32', '#cdc60a')
    other_spots_list <- c('#666568', '#93d519', '#cdc60a')
    other_spots_keys <- c(111, 222, 333)
    names(other_spots_list) <- other_spots_keys
    # combine the GSVA Enrichment colors with the three categoricals
    test_colors <- c(mylist, other_spots_list)
    return(test_colors)
}


# function to make the plot
plot_reshaped_data <- function(temp_seu_obj, curr_enrichment) {

    # get plotting colors using our function from above
    plot_colors <- get_plotting_colors()
    plot_colors
    
    options(repr.plot.width=15, repr.plot.height=12)
    Idents(temp_seu_obj) <- 'GSVA_Enrichment2'
    #pre_plot <- SpatialDimPlot(temp_seu_obj, cols = plot_colors)
    pre_plot <- SpatialPlot(temp_seu_obj, group.by = 'GSVA_Enrichment2', cols = plot_colors, pt.size.factor=2.6, image.alpha = 0.9, stroke = 0.3)
    plot_VA2 <- pre_plot[[1]] + 
                labs(title=curr_enrichment) + NoLegend()
    plot_VA1 <- pre_plot[[2]] + 
                labs(title=curr_enrichment) + NoLegend()
    plot_VB1 <- pre_plot[[3]] + 
                labs(title=curr_enrichment) + NoLegend()

    # define filename
    png_filename_VA1 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/Overall_Plots_Tumor_MHCI/png/Overall_Plot_VA1_",curr_enrichment)
    pdf_filename_VA1 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/Overall_Plots_Tumor_MHCI/pdf/Overall_Plot_VA1_",curr_enrichment)

    # save png
    png(file = paste0(png_filename_VA1,".png"),
        width = 1000,
        height = 800)
    print(plot_VA1)
    dev.off()

    # save pdf
    pdf(file = paste0(pdf_filename_VA1,".pdf"),
        width = 10,
        height = 8)
    print(plot_VA1)
    dev.off()

    # VB
    # define filename
    png_filename_VB1 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/Overall_Plots_Tumor_MHCI/png/Overall_Plot_VB1_",curr_enrichment)
    pdf_filename_VB1 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/Overall_Plots_Tumor_MHCI/pdf/Overall_Plot_VB1_",curr_enrichment)

    # save png
    png(file = paste0(png_filename_VB1,".png"),
        width = 1000,
        height = 800)
    print(plot_VB1)
    dev.off()

    # save pdf
    pdf(file = paste0(pdf_filename_VB1,".pdf"),
        width = 10,
        height = 8)
    print(plot_VB1)
    dev.off()

    # define filename
    png_filename_VA2 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/Overall_Plots_Tumor_MHCI/png/Overall_Plot_VA2_",curr_enrichment)
    pdf_filename_VA2 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/Overall_Plots_Tumor_MHCI/pdf/Overall_Plot_VA2_",curr_enrichment)

    # save png
    png(file = paste0(png_filename_VA2,".png"),
        width = 1000,
        height = 800)
    print(plot_VA2)
    dev.off()

    # save pdf
    pdf(file = paste0(pdf_filename_VA2,".pdf"),
        width = 10,
        height = 8)
    print(plot_VA2)
    dev.off()
}


reshape_quantile_and_plot <- function(seu_obj, enrichment_list) {
    
    # extract metadata from seurat object for reshaping...
    df <- seu_obj@meta.data
    
    # First, reshape the data by pivoting longer on the GSVA Enrichments. 
    # We will subset to focus on each individual GSVA Enrichment one at a time later. 
    df_melt <- df %>% pivot_longer(enrichment_list, names_to = "Enrichment_Label", values_to = "GSVA_Enrichment")
    
    # cycle through each GSVA Enrichment of interest
    for (curr_enrichment in enrichment_list) {
        print('------------------------------------------')
        print(curr_enrichment)
        print('------------------------------------------')
        # subset the melted df to focus on the current GSVA Enrichment
        df_melt_sub <- df_melt[df_melt$Enrichment_Label==curr_enrichment,]
        # remove duplicated rows if there are any
        df_melt_sub <- df_melt_sub[!duplicated(df_melt_sub),]
        # Measure Quantiles for the GSVA Enrichment so we can plot them with easier color scales...
        # Essentially converting the quantitative variable into a categorical one so we can also plot with 
        # the other categorical variables (Neg, center spots)
        df_melt_sub <- df_melt_sub %>% mutate(GSVA_Enrichment2 = ntile(GSVA_Enrichment, 100))
        # Label the Neg (111), VascHigh_MN4Low (222) and VascHigh_MN4High (333) spots as well...
        df_melt_sub$GSVA_Enrichment2[df_melt_sub$tumor_plot_label == 'Tumor_Neg_or_Low'] <- 111
        df_melt_sub$GSVA_Enrichment2[df_melt_sub$tumor_plot_label == 'TumorHigh_MHCILow'] <- 222
        df_melt_sub$GSVA_Enrichment2[df_melt_sub$tumor_plot_label == 'TumorHigh_MHCIHigh'] <- 333

        # Make the GSVA_Enrichment2 column a factor so we can use it to plot...
        df_melt_sub$GSVA_Enrichment2 <- as.factor(df_melt_sub$GSVA_Enrichment2)
        
        # Make copy of seurat object for plotting, using the reshaped metadata
        temp_seu_obj <- seu_obj
        df_melt_sub <- as.data.frame(df_melt_sub)
        rownames(df_melt_sub) <- df_melt_sub$Barcode
        temp_seu_obj@meta.data <- df_melt_sub

        #print(df_melt_sub %>% select(c(Barcode, Vasc_MN4_Label, vasc_plot_label, GSVA_Enrichment, GSVA_Enrichment2)) %>%head(3))
        
        #print(temp_seu_obj)

        #return(temp_seu_obj)
        
        # Now we can proceed with the plotting using our function from above
        plot_reshaped_data(temp_seu_obj, curr_enrichment)
        
    }
}
```


```R
#enrichment_list_to_use <- c('Macrophage','T_Cell','CD8_T_Cell','NK_Cell','B_Cell','DC', 'MN4_EC_Phenotype')
enrichment_list_to_use <- c('NK_Cell')
reshape_quantile_and_plot(data, enrichment_list_to_use)
#temp_seu_obj <- reshape_quantile_and_plot(AB_norm, enrichment_list_to_use)
```

    [1] "------------------------------------------"
    [1] "NK_Cell"
    [1] "------------------------------------------"


# Now let's do the same, but plot TumorHigh_MHCILow and TumorHigh_MHCIHigh separately


```R
data@meta.data %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 144</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_mean</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_mean</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_mean</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_mean</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_mean</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_mean</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_mean</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_mean</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_mean</th><th scope=col>MN4_EC_Phenotype_nn_vals</th><th scope=col>tumor_plot_label</th><th scope=col>TumorHigh_MHCILow_Neighbor</th><th scope=col>TumorHigh_MHCIHigh_Neighbor</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td> 0.6906627</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>30730</td><td>26443</td><td>AAACAACGAATAGTTC-1</td><td>VA1_0</td><td>VA1_854 </td><td>VA1_1863</td><td>VA1_2089</td><td>VA1_517 </td><td>VA1_2692</td><td>VA1_2472</td><td>3</td><td>['VA1_854', 'VA1_1863', 'VA1_2089', 'VA1_517', 'VA1_2692', 'VA1_2472'] </td><td>-0.1403735</td><td>-0.68346143</td><td>-0.6824274</td><td>-0.02534835</td><td>-0.1830202</td><td>-0.6898109</td><td>-0.4007403</td><td>[-0.140373503285981, -0.683461428931692, -0.682427363760641, -0.0253483482945374, -0.183020236470598, -0.689810896482086] </td><td> 0.4787432</td><td>-0.11469175</td><td>-0.05304243</td><td>-0.07726181</td><td>0.05871469</td><td> 0.1634258</td><td> 0.07598128</td><td>[0.478743199137973, -0.114691753402848, -0.053042433947209, -0.0772618094476313, 0.058714690391722, 0.163425797244824]  </td><td>0.4985932</td><td>-0.1559936</td><td>-0.09335601</td><td>-0.1178707</td><td>-0.2441883</td><td>-0.2408739</td><td>-0.05894823</td><td>[0.498593160799533, -0.155993596157572, -0.0933560136080898, -0.11787072243338, -0.244188290075738, -0.240873901729021]</td><td>0.04912155</td><td> 0.46869062</td><td>0.1456998</td><td> 0.4895150</td><td>-0.05807054</td><td>0.42619319</td><td>0.25352493</td><td>[0.0491215471122242, 0.468690624713802, 0.145699754382207, 0.489514984067094, -0.058070538457045, 0.42619318624764]  </td><td>0.18087518</td><td>0.4527944</td><td> 0.21247564</td><td> 0.42031179</td><td> 0.2671539</td><td> 0.5003215</td><td> 0.3389888</td><td>[0.180875180364392, 0.452794442072525, 0.212475639001522, 0.420311792896499, 0.267153935558, 0.500321548802953]         </td><td>-0.0428909507</td><td> 0.1220257</td><td>-0.5131311</td><td>-0.53971745</td><td>-0.38422937</td><td> 0.003417859</td><td>-0.2257542</td><td>[-0.0428909507251363, 0.122025706475451, -0.513131142574183, -0.539717446893455, -0.384229368789918, 0.0034178586983256] </td><td>-0.1024639</td><td> 0.6480175</td><td>-0.2458688</td><td> 0.39703261</td><td>-0.2214054</td><td> 0.62847798</td><td> 0.1839650</td><td>[-0.102463880569944, 0.648017492909985, -0.245868803204754, 0.397032612718896, -0.221405399324788, 0.62847797680395]  </td><td>-0.1681126</td><td>-0.1728810</td><td>-0.07891578</td><td> 0.9548053</td><td>-0.2141428</td><td>-0.2022404</td><td> 0.0197521</td><td>[-0.168112618939413, -0.172881026928877, -0.0789157831565245, 0.954805319001111, -0.214142828565577, -0.202240448089484]</td><td> 0.08796249</td><td> 0.04550749</td><td>-0.1515810</td><td>-0.06665673</td><td> 0.08863924</td><td> 0.1215510</td><td> 0.02090375</td><td>[0.0879624949143438, 0.0455074885754431, -0.151580958192362, -0.0666567297841213, 0.0886392350786018, 0.121550954065541] </td><td>Tumor_Neg_or_Low  </td><td>Not_Neighbor</td><td>Not_Neighbor</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284 </td><td>VA1_771 </td><td>VA1_352 </td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']  </td><td>-0.3313744</td><td>-0.25345755</td><td>-0.3396933</td><td>-0.70457534</td><td>-0.1946746</td><td>-0.6988578</td><td>-0.4204388</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]   </td><td>-0.3655433</td><td> 0.08738646</td><td>-0.01253212</td><td>-0.19925824</td><td>0.07594644</td><td>-0.1103883</td><td>-0.08739819</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.33109866</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>-0.17753110</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.54673283</td><td>-0.11670278</td><td>0.2619499</td><td> 0.4510973</td><td> 0.38571992</td><td>0.07931599</td><td>0.26801886</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073] </td><td>0.55563396</td><td>0.2090081</td><td>-0.14647071</td><td> 0.07915905</td><td> 0.6203140</td><td>-0.3870134</td><td> 0.1551052</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]   </td><td>-0.1416772948</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.52396701</td><td>-0.548042038</td><td>-0.3423706</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td> 0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td> 0.14113379</td><td>-0.4395605</td><td>-0.06779542</td><td>-0.1273297</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.32306461</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>-0.2725464</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.14347061</td><td>-0.3029000</td><td>-0.24465014</td><td>-0.35970213</td><td>-0.2758603</td><td>-0.22987461</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td><td>TumorHigh_Neighbor</td><td>Neighbor    </td><td>Not_Neighbor</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td> 6</td><td> 6</td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.2739867</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.3895103</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>High</td><td>Low              </td><td>vasc_low          </td><td>vasc_low             </td><td>29633</td><td>20870</td><td>AAACAATCTACTAGCA-1</td><td>VA1_2</td><td>VA1_2207</td><td>VA1_1663</td><td>VA1_1116</td><td>VA1_1210</td><td>VA1_920 </td><td>VA1_1124</td><td>4</td><td>['VA1_2207', 'VA1_1663', 'VA1_1116', 'VA1_1210', 'VA1_920', 'VA1_1124']</td><td>-0.1021611</td><td>-0.09404125</td><td>-0.7143306</td><td>-0.09243213</td><td>-0.1110161</td><td>-0.2277834</td><td>-0.2236274</td><td>[-0.102161065707172, -0.0940412510144283, -0.714330573226588, -0.0924321318329605, -0.111016118158986, -0.227783356007723]</td><td>-0.1750505</td><td> 0.16211636</td><td> 0.37237171</td><td> 0.51536446</td><td>0.14442029</td><td>-0.3628057</td><td> 0.10940276</td><td>[-0.175050506937197, 0.162116360694768, 0.372371708563783, 0.515364456905993, 0.144420288347746, -0.362805745678248]    </td><td>0.5555935</td><td> 0.5465806</td><td> 0.32440788</td><td>-0.2791675</td><td> 0.5603602</td><td>-0.2980788</td><td> 0.23494930</td><td>[0.555593475471297, 0.546580604053197, 0.324407882945447, -0.279167500500132, 0.560360187021303, -0.298078847308211]   </td><td>0.02470210</td><td> 0.02862245</td><td>0.3301684</td><td>-0.1450918</td><td>-0.05224927</td><td>0.39084903</td><td>0.09616682</td><td>[0.02470210181234, 0.0286224547429464, 0.330168428714374, -0.145091843512385, -0.0522492652549711, 0.390849032936518]</td><td>0.07100433</td><td>0.2218239</td><td>-0.06996729</td><td>-0.54888330</td><td>-0.5257730</td><td>-0.4789928</td><td>-0.2217980</td><td>[0.0710043326473108, 0.221823916630173, -0.0699672876994677, -0.548883300106269, -0.525772957160213, -0.478992784712071]</td><td> 0.0002582967</td><td>-0.5967590</td><td> 0.3099307</td><td>-0.37202321</td><td>-0.03360161</td><td> 0.083418241</td><td>-0.1014628</td><td>[0.0002582966608133, -0.596759018177968, 0.309930663526349, -0.372023213928184, -0.0336016075053905, 0.0834182412782388] </td><td>-0.1388668</td><td> 0.1185314</td><td>-0.4495049</td><td>-0.07106897</td><td>-0.4148401</td><td> 0.31776205</td><td>-0.1063312</td><td>[-0.138866804624008, 0.118531428276318, -0.44950493543207, -0.0710689737394324, -0.414840143712403, 0.317762048735336]</td><td> 0.9090881</td><td>-0.1755351</td><td> 0.78090959</td><td>-0.2698540</td><td>-0.2240448</td><td>-0.2925585</td><td> 0.1213342</td><td>[0.90908809327848, -0.175535107021276, 0.78090959121468, -0.26985397079401, -0.224044808961654, -0.292558511702187]     </td><td>-0.20592704</td><td>-0.18814347</td><td>-0.1059103</td><td>-0.01350095</td><td>-0.09839253</td><td> 0.1376595</td><td>-0.07903581</td><td>[-0.205927037792345, -0.188143473994638, -0.105910293614525, -0.0135009514713148, -0.0983925267325698, 0.137659452545575]</td><td>TumorHigh_MHCILow </td><td>Not_Neighbor</td><td>Neighbor    </td></tr>
</tbody>
</table>




```R
unique(data@meta.data$tumor_plot_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_Neg_or_Low'</li><li>'TumorHigh_Neighbor'</li><li>'TumorHigh_MHCILow'</li><li>'TumorHigh_MHCIHigh'</li></ol>




```R
# function to make the plot
plot_reshaped_data_2 <- function(temp_seu_obj, curr_enrichment, plot_group_name) {

    # get plotting colors using our function from above
    plot_colors <- get_plotting_colors()
    plot_colors
    
    options(repr.plot.width=15, repr.plot.height=12)
    Idents(temp_seu_obj) <- 'GSVA_Enrichment2'
    #pre_plot <- SpatialDimPlot(temp_seu_obj, cols = plot_colors)
    pre_plot <- SpatialPlot(temp_seu_obj, group.by = 'GSVA_Enrichment2', cols = plot_colors, pt.size.factor=2.6, image.alpha = 0.9, stroke = 0.3)
    plot_VA2 <- pre_plot[[1]] + 
                labs(title=curr_enrichment) + NoLegend()
    plot_VA1 <- pre_plot[[2]] + 
                labs(title=curr_enrichment) + NoLegend()
    plot_VB1 <- pre_plot[[3]] + 
                labs(title=curr_enrichment) + NoLegend()


    # define filename
    png_filename_VA1 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/",plot_group_name,"/png/Overall_Plot_VA1_",curr_enrichment)
    pdf_filename_VA1 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/",plot_group_name,"/pdf/Overall_Plot_VA1_",curr_enrichment)

    # save png
    png(file = paste0(png_filename_VA1,".png"),
        width = 1000,
        height = 800)
    print(plot_VA1)
    dev.off()

    # save pdf
    pdf(file = paste0(pdf_filename_VA1,".pdf"),
        width = 10,
        height = 8)
    print(plot_VA1)
    dev.off()

    # VB
    # define filename
    png_filename_VB1 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/",plot_group_name,"/png/Overall_Plot_VB1_",curr_enrichment)
    pdf_filename_VB1 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/",plot_group_name,"/pdf/Overall_Plot_VB1_",curr_enrichment)

    # save png
    png(file = paste0(png_filename_VB1,".png"),
        width = 1000,
        height = 800)
    print(plot_VB1)
    dev.off()

    # save pdf
    pdf(file = paste0(pdf_filename_VB1,".pdf"),
        width = 10,
        height = 8)
    print(plot_VB1)
    dev.off()

    # define filename
    png_filename_VA2 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/",plot_group_name,"/png/Overall_Plot_VA2_",curr_enrichment)
    pdf_filename_VA2 <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/",plot_group_name,"/pdf/Overall_Plot_VA2_",curr_enrichment)

    # save png
    png(file = paste0(png_filename_VA2,".png"),
        width = 1000,
        height = 800)
    print(plot_VA2)
    dev.off()

    # save pdf
    pdf(file = paste0(pdf_filename_VA2,".pdf"),
        width = 10,
        height = 8)
    print(plot_VA2)
    dev.off()
}


reshape_quantile_and_plot_2 <- function(seu_obj, enrichment_list) {
    
    # extract metadata from seurat object for reshaping...
    df <- seu_obj@meta.data
    
    # First, reshape the data by pivoting longer on the GSVA Enrichments. 
    # We will subset to focus on each individual GSVA Enrichment one at a time later. 
    df_melt <- df %>% pivot_longer(enrichment_list, names_to = "Enrichment_Label", values_to = "GSVA_Enrichment")
    
    # cycle through each GSVA Enrichment of interest
    for (curr_enrichment in enrichment_list) {
        print('------------------------------------------')
        print(curr_enrichment)
        print('------------------------------------------')
        # subset the melted df to focus on the current GSVA Enrichment
        df_melt_sub <- df_melt[df_melt$Enrichment_Label==curr_enrichment,]
        # remove duplicated rows if there are any
        df_melt_sub <- df_melt_sub[!duplicated(df_melt_sub),]
        # Measure Quantiles for the GSVA Enrichment so we can plot them with easier color scales...
        # Essentially converting the quantitative variable into a categorical one so we can also plot with 
        # the other categorical variables (Neg, center spots)
        df_melt_sub <- df_melt_sub %>% mutate(GSVA_Enrichment2 = ntile(GSVA_Enrichment, 100))
        # Label the Neg (111), TumorHigh_MHCILow (222) and TumorHigh_MHCIHigh (333) spots as well...
        df_melt_sub$GSVA_Enrichment2[df_melt_sub$tumor_plot_label == 'Tumor_Neg_or_Low'] <- 111

        # Make separate dfs for MN4Low and MN4High focus
        df_melt_sub_MHCILow <- df_melt_sub
        df_melt_sub_MHCIHigh <- df_melt_sub

        
        # MN4Low focus
        # Grey out MHCIHigh spots if they aren't MHCILow neighbors
        df_melt_sub_MHCILow$GSVA_Enrichment2[(df_melt_sub_MHCILow$tumor_plot_label == 'TumorHigh_MCIHigh') &
                                            (df_melt_sub_MHCILow$TumorHigh_MHCILow_Neighbor == 'Not_Neighbor')] <- 111
        # Grey out TumorHigh_Neighbor spots if they aren't MHCILow Neighbors
        df_melt_sub_MHCILow$GSVA_Enrichment2[(df_melt_sub_MHCILow$tumor_plot_label == 'TumorHigh_Neighbor') &
                                             (df_melt_sub_MHCILow$TumorHigh_MHCILow_Neighbor == 'Not_Neighbor')] <- 111
        # Keep MHCILow spots
        df_melt_sub_MHCILow$GSVA_Enrichment2[df_melt_sub_MHCILow$tumor_plot_label == 'TumorHigh_MHCILow'] <- 222

        
        # MHCIHigh focus
        # Grey out MHCILow spots if they aren't MHCIHigh neighbors
        df_melt_sub_MHCIHigh$GSVA_Enrichment2[(df_melt_sub_MHCIHigh$tumor_plot_label == 'TumorHigh_MHCILow') &
                                             (df_melt_sub_MHCIHigh$TumorHigh_MHCIHigh_Neighbor == 'Not_Neighbor')] <- 111
        # Grey out TumorHigh_Neighbor spots if they aren't MHCIHigh Neighbors
        df_melt_sub_MHCIHigh$GSVA_Enrichment2[(df_melt_sub_MHCIHigh$tumor_plot_label == 'TumorHigh_Neighbor') &
                                             (df_melt_sub_MHCIHigh$TumorHigh_MHCIHigh_Neighbor == 'Not_Neighbor')] <- 111
        # Keep MN4High spots
        df_melt_sub_MHCIHigh$GSVA_Enrichment2[df_melt_sub_MHCIHigh$tumor_plot_label == 'TumorHigh_MHCIHigh'] <- 333

        # Make the GSVA_Enrichment2 column a factor so we can use it to plot...
        df_melt_sub_MHCILow$GSVA_Enrichment2 <- as.factor(df_melt_sub_MHCILow$GSVA_Enrichment2)
        df_melt_sub_MHCIHigh$GSVA_Enrichment2 <- as.factor(df_melt_sub_MHCIHigh$GSVA_Enrichment2)
        
        # Make copy of seurat object for plotting, using the reshaped metadata
        temp_seu_obj_MHCILow <- seu_obj
        temp_seu_obj_MHCIHigh <- seu_obj
        df_melt_sub_MHCILow <- as.data.frame(df_melt_sub_MHCILow)
        df_melt_sub_MHCIHigh <- as.data.frame(df_melt_sub_MHCIHigh)
        rownames(df_melt_sub_MHCILow) <- df_melt_sub_MHCILow$Barcode
        rownames(df_melt_sub_MHCIHigh) <- df_melt_sub_MHCIHigh$Barcode
        temp_seu_obj_MHCILow@meta.data <- df_melt_sub_MHCILow
        temp_seu_obj_MHCIHigh@meta.data <- df_melt_sub_MHCIHigh

        
        # Now we can proceed with the plotting using our function from above
        plot_reshaped_data_2(temp_seu_obj_MHCILow, curr_enrichment, plot_group_name='TumorHigh_MHCILow_only')
        plot_reshaped_data_2(temp_seu_obj_MHCIHigh, curr_enrichment, plot_group_name='TumorHigh_MHCIHigh_only')
        
    }
}
```


```R
#enrichment_list_to_use <- c('Macrophage','T_Cell','CD8_T_Cell','NK_Cell','B_Cell','DC', 'MN4_EC_Phenotype')
enrichment_list_to_use <- c('NK_Cell')
reshape_quantile_and_plot_2(data, enrichment_list_to_use)
```

    [1] "------------------------------------------"
    [1] "NK_Cell"
    [1] "------------------------------------------"



```R
unique(data@meta.data$Tumor_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high'</li></ol>



# Now let's plot each individual spot/neighborhood on its own


```R
# function to make the plot
plot_reshaped_data_3 <- function(temp_seu_obj, curr_enrichment, curr_id) {

    # Only plot VA or VB depending on where the curr spot is...
    if (strsplit(curr_id, "_")[[1]][1] == 'VA2') {
        img_to_use = 1
        }
    if (strsplit(curr_id, "_")[[1]][1] == 'VA1') {
        img_to_use = 2
        }
    if (strsplit(curr_id, "_")[[1]][1] == 'VB1') {
        img_to_use = 3
        }
    
    # get plotting colors using our function from above
    plot_colors <- get_plotting_colors()
    plot_colors
    
    options(repr.plot.width=15, repr.plot.height=12)
    Idents(temp_seu_obj) <- 'GSVA_Enrichment2'
    #pre_plot <- SpatialDimPlot(temp_seu_obj, cols = plot_colors)
    pre_plot <- SpatialPlot(temp_seu_obj, group.by = 'GSVA_Enrichment2', cols = plot_colors, pt.size.factor=2.6, image.alpha = 0.9, stroke = 0.3)
    plot <- pre_plot[[img_to_use]] + 
                labs(title=curr_enrichment) + NoLegend()
    plot

    # define filename
    png_filename <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/Individual_Spots_Tumor_MHCI/png/",curr_id,"_",curr_enrichment)
    pdf_filename <- paste0("/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/14_Visualize_Neighborhoods_tumor_MHCI/Individual_Spots_Tumor_MHCI/pdf/",curr_id,"_",curr_enrichment)

    # save png
    png(file = paste0(png_filename,".png"),
        width = 1000,
        height = 800)
    print(plot)
    dev.off()

    # save pdf
    pdf(file = paste0(pdf_filename,".pdf"),
        width = 10,
        height = 8)
    print(plot)
    dev.off()
}


reshape_quantile_and_plot_3 <- function(seu_obj, enrichment_list) {
    
    # extract metadata from seurat object for reshaping...
    df <- seu_obj@meta.data
    
    # First, reshape the data by pivoting longer on the GSVA Enrichments. 
    # We will subset to focus on each individual GSVA Enrichment one at a time later. 
    df_melt <- df %>% pivot_longer(enrichment_list, names_to = "Enrichment_Label", values_to = "GSVA_Enrichment")
    
    # cycle through each GSVA Enrichment of interest
    for (curr_enrichment in enrichment_list) {
        print('------------------------------------------')
        print(curr_enrichment)
        print('------------------------------------------')
        # subset the melted df to focus on the current GSVA Enrichment
        df_melt_sub <- df_melt[df_melt$Enrichment_Label==curr_enrichment,]
        # remove duplicated rows if there are any
        df_melt_sub <- df_melt_sub[!duplicated(df_melt_sub),]
        # Measure Quantiles for the GSVA Enrichment so we can plot them with easier color scales...
        # Essentially converting the quantitative variable into a categorical one so we can also plot with 
        # the other categorical variables (Neg, center spots)
        df_melt_sub <- df_melt_sub %>% mutate(GSVA_Enrichment2 = ntile(GSVA_Enrichment, 100))

        # Store a copy of these quantiles.... in column GSVA_Enrichment3
        df_melt_sub$GSVA_Enrichment3 <- df_melt_sub$GSVA_Enrichment2

        # Now cycle through all vasc_high spots, and plot them individually
        # -- reused code from before --
        # get list of vasc_high spots. Use "id" column
        tumor_high_id_list <- unique(df_melt_sub$id[df_melt_sub$Tumor_label == 'Tumor_high'])
        # for each id in tumor_high_id_list, look at its neighbor ids. then, go and set tumor_plot_label to TumorHigh_Neighbor
        # and set Neighbor or Not_Neighbor in two other columns, TumorHigh_MHCILow_Neighbor and TumorHigh_MHCIHigh_Neighbor
        
        for (curr_id in tumor_high_id_list) {
            
            # get neighbor ids fur the current center spot id
            curr_nn1_id <- df_melt_sub$nn1_id[df_melt_sub$id == curr_id]
            curr_nn2_id <- df_melt_sub$nn2_id[df_melt_sub$id == curr_id]
            curr_nn3_id <- df_melt_sub$nn3_id[df_melt_sub$id == curr_id]
            curr_nn4_id <- df_melt_sub$nn4_id[df_melt_sub$id == curr_id]
            curr_nn5_id <- df_melt_sub$nn5_id[df_melt_sub$id == curr_id]
            curr_nn6_id <- df_melt_sub$nn6_id[df_melt_sub$id == curr_id]
        
            curr_neighbor_id_list <- c(curr_nn1_id, curr_nn2_id, curr_nn3_id, curr_nn4_id, curr_nn5_id, curr_nn6_id)

            # initialize all spots to grey at first
            df_melt_sub$GSVA_Enrichment2 <- 111

            # fill in GSVA Enrichments (quantized) for each neighbor of the curr_id
            for (curr_neighbor_id in curr_neighbor_id_list) {
                # access the copied quantized enrichment in GSVA_Enrichment3
                df_melt_sub$GSVA_Enrichment2[df_melt_sub$id == curr_neighbor_id] <- df_melt_sub$GSVA_Enrichment3[df_melt_sub$id == curr_neighbor_id]
            }


            # Now fill in the center spot itself with the categorical label color
            df_melt_sub$GSVA_Enrichment2[(df_melt_sub$tumor_plot_label == 'TumorHigh_MHCILow') &
                                         (df_melt_sub$id == curr_id)] <- 222
            df_melt_sub$GSVA_Enrichment2[(df_melt_sub$tumor_plot_label == 'TumorHigh_MHCIHigh') &
                                         (df_melt_sub$id == curr_id)] <- 333
            
            
            # Make the GSVA_Enrichment2 column a factor so we can use it to plot...
            df_melt_sub$GSVA_Enrichment2 <- as.factor(df_melt_sub$GSVA_Enrichment2)
            
            # Make copy of seurat object for plotting, using the reshaped metadata
            temp_seu_obj <- seu_obj
            df_melt_sub <- as.data.frame(df_melt_sub)
            rownames(df_melt_sub) <- df_melt_sub$Barcode
            temp_seu_obj@meta.data <- df_melt_sub
            
            # Now we can proceed with the plotting using our function from above
            plot_reshaped_data_3(temp_seu_obj, curr_enrichment, curr_id)
        }
    }
}
```


```R
#enrichment_list_to_use <- c('Macrophage','T_Cell','CD8_T_Cell','NK_Cell','B_Cell','DC', 'MN4_EC_Phenotype')
enrichment_list_to_use <- c('NK_Cell')
reshape_quantile_and_plot_3(data, enrichment_list_to_use)
```

    [1] "------------------------------------------"
    [1] "NK_Cell"
    [1] "------------------------------------------"



```R
print('done')
```

    [1] "done"



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


```R

```


```R

```


```R

```


```R

```
