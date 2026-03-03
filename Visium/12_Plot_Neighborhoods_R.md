# Measured neighborhoods in python, now plot some results in R

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
# tidyverse includes: ggplot2, dplyr, tidyr, readr, purrr, tibble, stringr, forcats. 
# https://www.tidyverse.org/
library(tidyverse)
#library(Seurat, quietly = T)
library(ggplot2)
library(ggrepel)
library(ggbeeswarm)
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

options(repr.matrix.max.cols=1000, repr.matrix.max.rows=100)
```

# load the neighborhoods data


```R
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
nns_all <- as.data.frame(nns_all) %>% select(-c('...1'))
nns_6adj <- as.data.frame(nns_6adj) %>% select(-c('...1'))

rownames(nns_all) <- nns_all$Barcode
rownames(nns_6adj) <- nns_6adj$Barcode
```


```R
nns_all %>% head(2)
nns_6adj %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 141</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_mean</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_mean</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_mean</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_mean</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_mean</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_mean</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_mean</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_mean</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_mean</th><th scope=col>MN4_EC_Phenotype_nn_vals</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td> 0.6906627</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.3217252</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>High</td><td>Above80Percentile</td><td>vasc_high_MN4_high</td><td>vasc_high_MN4_upper80</td><td>30730</td><td>26443</td><td>AAACAACGAATAGTTC-1</td><td>VA1_0</td><td>VA1_854</td><td>VA1_1863</td><td>VA1_2089</td><td>VA1_517 </td><td>VA1_2692</td><td>VA1_2472</td><td>3</td><td>['VA1_854', 'VA1_1863', 'VA1_2089', 'VA1_517', 'VA1_2692', 'VA1_2472']</td><td>-0.1403735</td><td>-0.6834614</td><td>-0.6824274</td><td>-0.02534835</td><td>-0.1830202</td><td>-0.6898109</td><td>-0.4007403</td><td>[-0.140373503285981, -0.683461428931692, -0.682427363760641, -0.0253483482945374, -0.183020236470598, -0.689810896482086]</td><td> 0.4787432</td><td>-0.11469175</td><td>-0.05304243</td><td>-0.07726181</td><td>0.05871469</td><td> 0.1634258</td><td> 0.07598128</td><td>[0.478743199137973, -0.114691753402848, -0.053042433947209, -0.0772618094476313, 0.058714690391722, 0.163425797244824]  </td><td>0.4985932</td><td>-0.1559936</td><td>-0.09335601</td><td>-0.1178707</td><td>-0.2441883</td><td>-0.2408739</td><td>-0.05894823</td><td>[0.498593160799533, -0.155993596157572, -0.0933560136080898, -0.11787072243338, -0.244188290075738, -0.240873901729021]</td><td>0.04912155</td><td> 0.4686906</td><td>0.1456998</td><td>0.4895150</td><td>-0.05807054</td><td>0.42619319</td><td>0.2535249</td><td>[0.0491215471122242, 0.468690624713802, 0.145699754382207, 0.489514984067094, -0.058070538457045, 0.42619318624764] </td><td>0.1808752</td><td>0.4527944</td><td> 0.2124756</td><td>0.42031179</td><td>0.2671539</td><td> 0.5003215</td><td>0.3389888</td><td>[0.180875180364392, 0.452794442072525, 0.212475639001522, 0.420311792896499, 0.267153935558, 0.500321548802953]      </td><td>-0.04289095</td><td> 0.1220257</td><td>-0.5131311</td><td>-0.53971745</td><td>-0.3842294</td><td> 0.003417859</td><td>-0.2257542</td><td>[-0.0428909507251363, 0.122025706475451, -0.513131142574183, -0.539717446893455, -0.384229368789918, 0.0034178586983256] </td><td>-0.1024639</td><td> 0.6480175</td><td>-0.2458688</td><td>0.3970326</td><td>-0.2214054</td><td> 0.62847798</td><td> 0.1839650</td><td>[-0.102463880569944, 0.648017492909985, -0.245868803204754, 0.397032612718896, -0.221405399324788, 0.62847797680395]  </td><td>-0.1681126</td><td>-0.1728810</td><td>-0.07891578</td><td> 0.9548053</td><td>-0.2141428</td><td>-0.2022404</td><td> 0.0197521</td><td>[-0.168112618939413, -0.172881026928877, -0.0789157831565245, 0.954805319001111, -0.214142828565577, -0.202240448089484]</td><td> 0.08796249</td><td> 0.04550749</td><td>-0.151581</td><td>-0.06665673</td><td> 0.08863924</td><td> 0.1215510</td><td> 0.02090375</td><td>[0.0879624949143438, 0.0455074885754431, -0.151580958192362, -0.0666567297841213, 0.0886392350786018, 0.121550954065541] </td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td> 6</td><td> 6</td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td><td>SCLC                  </td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>Low </td><td>Below20Percentile</td><td>vasc_low          </td><td>vasc_low             </td><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284</td><td>VA1_771 </td><td>VA1_352 </td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604'] </td><td>-0.3313744</td><td>-0.2534575</td><td>-0.3396933</td><td>-0.70457534</td><td>-0.1946746</td><td>-0.6988578</td><td>-0.4204388</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]  </td><td>-0.3655433</td><td> 0.08738646</td><td>-0.01253212</td><td>-0.19925824</td><td>0.07594644</td><td>-0.1103883</td><td>-0.08739819</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.33109866</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>-0.17753110</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.54673283</td><td>-0.1167028</td><td>0.2619499</td><td>0.4510973</td><td> 0.38571992</td><td>0.07931599</td><td>0.2680189</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073]</td><td>0.5556340</td><td>0.2090081</td><td>-0.1464707</td><td>0.07915905</td><td>0.6203140</td><td>-0.3870134</td><td>0.1551052</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]</td><td>-0.14167729</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.5239670</td><td>-0.548042038</td><td>-0.3423706</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td> 0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td>0.1411338</td><td>-0.4395605</td><td>-0.06779542</td><td>-0.1273297</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.32306461</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>-0.2725464</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.14347061</td><td>-0.302900</td><td>-0.24465014</td><td>-0.35970213</td><td>-0.2758603</td><td>-0.22987461</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 2 × 141</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_mean</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_mean</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_mean</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_mean</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_mean</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_mean</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_mean</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_mean</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_mean</th><th scope=col>MN4_EC_Phenotype_nn_vals</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6</td><td>6</td><td>VA1</td><td>1.001716</td><td>vasc_low</td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low</td><td>SCLC            </td><td>SCLC</td><td>SCLC + Effector Immune Cells (T+NK)</td><td>6_SCLC_Effector_Immune</td><td>6_SCLC_Effector_Immune</td><td>-0.4320402</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.3802492</td><td>-0.2729638</td><td>-0.6165697</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.53207656</td><td>-0.2841903</td><td>-0.4421974</td><td>-0.20627210</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.323644</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.4421974</td><td>-0.48410628</td><td>Low</td><td>Below20Percentile</td><td>vasc_low</td><td>vasc_low</td><td>12698</td><td> 8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284 </td><td>VA1_771</td><td>VA1_352 </td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']  </td><td>-0.3313744</td><td>-0.2534575</td><td>-0.3396933</td><td>-0.7045753</td><td>-0.1946746</td><td>-0.6988578</td><td>-0.4204388</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165] </td><td>-0.3655433</td><td>0.08738646</td><td>-0.01253212</td><td>-0.1992582</td><td> 0.07594644</td><td>-0.1103883</td><td>-0.08739819</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td> 0.1532445</td><td>-0.2668601</td><td>-0.3310987</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>-0.17753110</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td> 0.54673283</td><td>-0.1167028</td><td>0.2619499</td><td>0.4510973</td><td>0.38571992</td><td> 0.07931599</td><td>0.2680189</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073]  </td><td> 0.5556340</td><td> 0.2090081</td><td>-0.1464707</td><td> 0.07915905</td><td> 0.620314</td><td>-0.3870134</td><td> 0.1551052</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]  </td><td>-0.1416773</td><td>-0.607476142</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.5239670</td><td>-0.5480420</td><td>-0.3423706</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td> 0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td> 0.1411338</td><td>-0.4395605</td><td>-0.06779542</td><td>-0.1273297</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941] </td><td>-0.3917784</td><td>-0.2591518</td><td>-0.3230646</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>-0.2725464</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.14347061</td><td>-0.302899967</td><td>-0.24465014</td><td>-0.3597021</td><td>-0.2758603</td><td>-0.2298746</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td></tr>
	<tr><th scope=row>VA1_AAACACCAATAACTGC-1</th><td>VA1_AAACACCAATAACTGC-1</td><td>VA1</td><td>VA1</td><td>6761</td><td>3762</td><td>round1</td><td>4936</td><td>3641</td><td>4</td><td>4</td><td>VA1</td><td>0.907884</td><td>vasc_low</td><td> 9</td><td>Tumor_low</td><td>Tumor_low</td><td>1</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low</td><td>SCLC+vasculature</td><td>SCLC</td><td>SCLC                               </td><td>4_SCLC                </td><td>1_4_5_11_SCLC         </td><td>-0.2949967</td><td>-0.1592542</td><td>-0.6246499</td><td>-0.7091880</td><td>-0.3126876</td><td>-0.1285324</td><td> 0.5479671</td><td>-0.2873825</td><td> 0.05564459</td><td>-0.3094928</td><td>-0.2340697</td><td>-0.03504189</td><td>-0.1584708</td><td>-0.2611306</td><td>-0.3051610</td><td> 0.3381445</td><td>-0.3300223</td><td>-0.141812</td><td>-0.1902565</td><td>-0.7206696</td><td>-0.2190269</td><td>-0.4549118</td><td>-0.3728931</td><td>-0.4603460</td><td>-0.2340697</td><td> 0.07776597</td><td>Low</td><td>Low              </td><td>vasc_low</td><td>vasc_low</td><td> 9524</td><td>25898</td><td>AAACACCAATAACTGC-1</td><td>VA1_3</td><td>VA1_1747</td><td>VA1_425</td><td>VA1_1286</td><td>VA1_1385</td><td>VA1_2332</td><td>VA1_2567</td><td>6</td><td>['VA1_1747', 'VA1_425', 'VA1_1286', 'VA1_1385', 'VA1_2332', 'VA1_2567']</td><td>-0.7080698</td><td>-0.1727255</td><td>-0.3059416</td><td>-0.7100941</td><td>-0.7212701</td><td>-0.7231039</td><td>-0.5568675</td><td>[-0.708069838814006, -0.172725483080612, -0.305941553061541, -0.710094101019417, -0.721270084194434, -0.723103920735993]</td><td>-0.2331998</td><td>0.07857171</td><td> 0.47593594</td><td>-0.3569468</td><td>-0.37384910</td><td>-0.3794862</td><td>-0.13149570</td><td>[-0.23319977872425, 0.0785717073459754, 0.475935940981545, -0.356946754827358, -0.373849101912343, -0.379486233065752]  </td><td>-0.2393405</td><td> 0.3552654</td><td> 0.2957953</td><td>-0.2810686</td><td> 0.3142241</td><td>-0.3293976</td><td> 0.01924632</td><td>[-0.239340532575186, 0.355265387687536, 0.295795295662033, -0.281068641184542, 0.314224061582033, -0.329397638582968]  </td><td>-0.01724305</td><td> 0.5798389</td><td>0.5448870</td><td>0.3572084</td><td>0.06235329</td><td>-0.25844957</td><td>0.2114325</td><td>[-0.0172430474396086, 0.579838912374754, 0.544886975743813, 0.357208398080687, 0.0623532857737059, -0.258449565629918]</td><td>-0.6919329</td><td>-0.8139959</td><td>-0.8053769</td><td>-0.59253595</td><td>-0.680867</td><td>-0.8014868</td><td>-0.7310326</td><td>[-0.69193292191416, -0.813995861950977, -0.805376914068887, -0.592535945325797, -0.680866999335305, -0.801486804482884]</td><td>-0.5362091</td><td> 0.002455842</td><td> 0.3126581</td><td>-0.15292795</td><td>-0.6217053</td><td>-0.6253788</td><td>-0.2701845</td><td>[-0.536209058824282, 0.0024558419296655, 0.31265811387002, -0.152927946116749, -0.621705278122557, -0.625378820848625]   </td><td>-0.4001988</td><td>-0.2774540</td><td> 0.1709900</td><td>-0.2870877</td><td>-0.4553370</td><td>-0.33424428</td><td>-0.2638886</td><td>[-0.400198776367015, -0.277454023705374, 0.170989973541542, -0.287087656671685, -0.455337042372491, -0.334244282995581]</td><td>-0.2081416</td><td>-0.2570514</td><td>-0.2973595</td><td>-0.2734547</td><td>-0.2881576</td><td>-0.3242649</td><td>-0.2747383</td><td>[-0.20814162832553, -0.25705141028191, -0.297359471894224, -0.273454690938038, -0.288157631526152, -0.324264852970433]  </td><td>-0.24604744</td><td>-0.08491335</td><td>-0.007453716</td><td> 0.02875985</td><td>-0.3204698</td><td>-0.3190422</td><td>-0.1581944</td><td>[-0.246047435280206, -0.084913354050493, -0.0074537156101513, 0.028759847156343, -0.320469802379713, -0.319042226520557] </td></tr>
</tbody>
</table>




```R
unique(nns_all$orig.ident)
unique(nns_6adj$orig.ident)
unique(nns_6adj$Tumor_MHCI_label2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'VA1'</li><li>'VB1'</li><li>'VA2'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'VA1'</li><li>'VB1'</li><li>'VA2'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high_MHCI_low'</li><li>'Tumor_high_MHCI_high'</li></ol>



# Melt to consolidate average neighbor scores into one column


```R
nns_6adj_melt <- nns_6adj %>% pivot_longer(c('Macrophage_nn_mean','T_Cell_nn_mean','CD8_T_Cell_nn_mean','NK_Cell_nn_mean','B_Cell_nn_mean',
                                             'DC_nn_mean','Endothelial_Cell_nn_mean', 'LEC_nn_mean','MN4_EC_Phenotype_nn_mean'), 
                                           names_to = "nn_mean_Cell_Type_Label", values_to = "nn_mean_Cell_Type_Enrichment")

```

## 1. Around a Tumor MHCI High spot, do we see vascular spots with high activation score (MN4_EC)? Compared to Tumor MHCI Low
### and...
## 2. Around a Tumor MHCI High spot, do we see more NK/T cells compared with Tumor MHCI Low?


```R
options(repr.plot.width=20, repr.plot.height=16)

NN_celltype_vs_TumorMHCILabel_plot2 <- ggplot(nns_6adj_melt[nns_6adj_melt$Tumor_MHCI_label2 %in% c('Tumor_high_MHCI_high','Tumor_high_MHCI_low'),], 
                                          aes(x=Tumor_MHCI_label2, y=nn_mean_Cell_Type_Enrichment, color=orig.ident)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey') + 
  geom_beeswarm(cex = 0.25, size=1) +
  #geom_point() + 
  facet_wrap(vars(nn_mean_Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Nearest Neighbors Cell Type GSVA Enrichment vs. Tumor_MHCI_Label")
             
NN_celltype_vs_TumorMHCILabel_plot2

# save plot
pdf(file = "Output/12_Plot_Neighborhoods/VAB1_VA2_NN_cell_type_GSVA_per_Tumor_MHCI_Label_HiVsLow.pdf",
    width = 1800,
    height = 1400)
print(NN_celltype_vs_TumorMHCILabel_plot2)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_12_1.png)
    


# subset to VA1 and VB1, of VA2


```R
nns_6adj_melt_VAB1 <- nns_6adj_melt[nns_6adj_melt$orig.ident %in% c('VA1','VB1'),]
nns_6adj_melt_VA2 <- nns_6adj_melt[nns_6adj_melt$orig.ident %in% c('VA2'),]
dim(nns_6adj_melt_VAB1)
dim(nns_6adj_melt_VA2)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>40509</li><li>134</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>36189</li><li>134</li></ol>



# VA1 and VB1 only


```R
options(repr.plot.width=20, repr.plot.height=16)

NN_celltype_vs_TumorMHCILabel_VAB1_plot2 <- ggplot(nns_6adj_melt_VAB1[nns_6adj_melt_VAB1$Tumor_MHCI_label2 %in% c('Tumor_high_MHCI_high','Tumor_high_MHCI_low'),], 
                                          aes(x=Tumor_MHCI_label2, y=nn_mean_Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey') + 
  geom_beeswarm(cex = 0.35, size=1) +
  #geom_point() + 
  facet_wrap(vars(nn_mean_Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "VAB1 Nearest Neighbors Cell Type GSVA Enrichment vs. Tumor_MHCI_Label")
             
NN_celltype_vs_TumorMHCILabel_VAB1_plot2

# save plot
pdf(file = "Output/12_Plot_Neighborhoods/VAB1_NN_cell_type_GSVA_per_Tumor_MHCI_Label_HiVsLow.pdf",
    width = 1800,
    height = 1400)
print(NN_celltype_vs_TumorMHCILabel_VAB1_plot2)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_16_1.png)
    


# VA2 only


```R
options(repr.plot.width=20, repr.plot.height=16)

NN_celltype_vs_TumorMHCILabel_VA2_plot2 <- ggplot(nns_6adj_melt_VA2[nns_6adj_melt_VA2$Tumor_MHCI_label2 %in% c('Tumor_high_MHCI_high','Tumor_high_MHCI_low'),], 
                                          aes(x=Tumor_MHCI_label2, y=nn_mean_Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey') + 
  geom_beeswarm(cex = 0.35, size=1) +
  #geom_point() + 
  facet_wrap(vars(nn_mean_Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "VA2 Nearest Neighbors Cell Type GSVA Enrichment vs. Tumor_MHCI_Label")
             
NN_celltype_vs_TumorMHCILabel_VA2_plot2

# save plot
pdf(file = "Output/12_Plot_Neighborhoods/VA2_NN_cell_type_GSVA_per_Tumor_MHCI_Label_HiVsLow.pdf",
    width = 1800,
    height = 1400)
print(NN_celltype_vs_TumorMHCILabel_VA2_plot2)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_18_1.png)
    



```R
# save data...
write.csv(nns_6adj_melt, 'Processing/10_Neighborhood_Analysis/nns_6adj_melt.csv', row.names = FALSE)
write.csv(nns_6adj_melt[nns_6adj_melt$Tumor_MHCI_label2 %in% c('Tumor_high_MHCI_high','Tumor_high_MHCI_low'),], 'Processing/10_Neighborhood_Analysis/nns_6adj_melt_TumorMHCI_Hi_vs_Low.csv', row.names = FALSE)
```


```R
nns_6adj_melt %>% head(2)
```


<table class="dataframe">
<caption>A tibble: 2 × 134</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_vals</th><th scope=col>nn_mean_Cell_Type_Label</th><th scope=col>nn_mean_Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6</td><td>6</td><td>VA1</td><td>1.001716</td><td>vasc_low</td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low</td><td>SCLC</td><td>SCLC</td><td>SCLC + Effector Immune Cells (T+NK)</td><td>6_SCLC_Effector_Immune</td><td>6_SCLC_Effector_Immune</td><td>-0.4320402</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.3802492</td><td>-0.2729638</td><td>-0.6165697</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.4421974</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td>0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.323644</td><td>-0.4311868</td><td>-0.713442</td><td>-0.1225925</td><td>0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.4421974</td><td>-0.4841063</td><td>Low</td><td>Below20Percentile</td><td>vasc_low</td><td>vasc_low</td><td>12698</td><td>8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284</td><td>VA1_771</td><td>VA1_352</td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']</td><td>-0.3313744</td><td>-0.2534575</td><td>-0.3396933</td><td>-0.7045753</td><td>-0.1946746</td><td>-0.6988578</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]</td><td>-0.3655433</td><td>0.08738646</td><td>-0.01253212</td><td>-0.1992582</td><td>0.07594644</td><td>-0.1103883</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.3310987</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.5467328</td><td>-0.1167028</td><td>0.2619499</td><td>0.4510973</td><td>0.3857199</td><td>0.07931599</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073]</td><td>0.555634</td><td>0.2090081</td><td>-0.1464707</td><td>0.07915905</td><td>0.620314</td><td>-0.3870134</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]</td><td>-0.1416773</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.523967</td><td>-0.548042</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td>0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td>0.1411338</td><td>-0.4395605</td><td>-0.06779542</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.3230646</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.1434706</td><td>-0.3029</td><td>-0.2446501</td><td>-0.3597021</td><td>-0.2758603</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td><td>Macrophage_nn_mean</td><td>-0.42043885</td></tr>
	<tr><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6</td><td>6</td><td>VA1</td><td>1.001716</td><td>vasc_low</td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low</td><td>SCLC</td><td>SCLC</td><td>SCLC + Effector Immune Cells (T+NK)</td><td>6_SCLC_Effector_Immune</td><td>6_SCLC_Effector_Immune</td><td>-0.4320402</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.3802492</td><td>-0.2729638</td><td>-0.6165697</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.4421974</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td>0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.323644</td><td>-0.4311868</td><td>-0.713442</td><td>-0.1225925</td><td>0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.4421974</td><td>-0.4841063</td><td>Low</td><td>Below20Percentile</td><td>vasc_low</td><td>vasc_low</td><td>12698</td><td>8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284</td><td>VA1_771</td><td>VA1_352</td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']</td><td>-0.3313744</td><td>-0.2534575</td><td>-0.3396933</td><td>-0.7045753</td><td>-0.1946746</td><td>-0.6988578</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]</td><td>-0.3655433</td><td>0.08738646</td><td>-0.01253212</td><td>-0.1992582</td><td>0.07594644</td><td>-0.1103883</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.3310987</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.5467328</td><td>-0.1167028</td><td>0.2619499</td><td>0.4510973</td><td>0.3857199</td><td>0.07931599</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073]</td><td>0.555634</td><td>0.2090081</td><td>-0.1464707</td><td>0.07915905</td><td>0.620314</td><td>-0.3870134</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]</td><td>-0.1416773</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.523967</td><td>-0.548042</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td>0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td>0.1411338</td><td>-0.4395605</td><td>-0.06779542</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.3230646</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.1434706</td><td>-0.3029</td><td>-0.2446501</td><td>-0.3597021</td><td>-0.2758603</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td><td>T_Cell_nn_mean    </td><td>-0.08739819</td></tr>
</tbody>
</table>



# save subset to only relevant columns. 


```R
tumor_high_data <- nns_6adj_melt[(nns_6adj_melt$Tumor_label == 'Tumor_high'),]

tumor_high_data_tosave = tumor_high_data %>% select(Barcode, orig.ident, Tumor_label, MHCI_label, Tumor_MHCI_label, nn_mean_Cell_Type_Label, nn_mean_Cell_Type_Enrichment)

# drop duplicates
tumor_high_data_tosave <- tumor_high_data_tosave[!duplicated(tumor_high_data_tosave), ]

tumor_high_data_tosave[1:3,]
```


<table class="dataframe">
<caption>A tibble: 3 × 7</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>Tumor_label</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>nn_mean_Cell_Type_Label</th><th scope=col>nn_mean_Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAGGGTCTATATT-1</td><td>VA1</td><td>Tumor_high</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Macrophage_nn_mean</td><td>-0.3936181</td></tr>
	<tr><td>VA1_AAACAGGGTCTATATT-1</td><td>VA1</td><td>Tumor_high</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>T_Cell_nn_mean    </td><td>-0.2079122</td></tr>
	<tr><td>VA1_AAACAGGGTCTATATT-1</td><td>VA1</td><td>Tumor_high</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>CD8_T_Cell_nn_mean</td><td>-0.2222651</td></tr>
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
	<tr><td>VA1</td><td>Tumor_high_MHCI_high</td><td> 4725</td></tr>
	<tr><td>VA1</td><td>Tumor_high_MHCI_low </td><td> 6597</td></tr>
	<tr><td>VA2</td><td>Tumor_high_MHCI_high</td><td> 6480</td></tr>
	<tr><td>VA2</td><td>Tumor_high_MHCI_low </td><td>16182</td></tr>
	<tr><td>VB1</td><td>Tumor_high_MHCI_high</td><td> 4383</td></tr>
	<tr><td>VB1</td><td>Tumor_high_MHCI_low </td><td> 6732</td></tr>
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
<ol class=list-inline><li>45099</li><li>7</li></ol>




```R
# save
write.csv(tumor_high_data_tosave, 'Output/12_Plot_Neighborhoods/nns_6adj_tumor_high_cell_type_enrichment_by_MHCI.csv', row.names = FALSE)

```

# save a version that's unmelted


```R
tumor_high_data_tosave_unmelted <- tumor_high_data_tosave %>% pivot_wider(names_from = nn_mean_Cell_Type_Label, values_from = nn_mean_Cell_Type_Enrichment)
```


```R
# save
write.csv(tumor_high_data_tosave_unmelted, 'Output/12_Plot_Neighborhoods/nns_6adj_tumor_high_cell_type_enrichment_by_MHCI_wideform.csv', row.names = FALSE)

```

# stats


```R
nns_6adj_melt %>% head(2)
```


<table class="dataframe">
<caption>A tibble: 2 × 134</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_vals</th><th scope=col>nn_mean_Cell_Type_Label</th><th scope=col>nn_mean_Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6</td><td>6</td><td>VA1</td><td>1.001716</td><td>vasc_low</td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low</td><td>SCLC</td><td>SCLC</td><td>SCLC + Effector Immune Cells (T+NK)</td><td>6_SCLC_Effector_Immune</td><td>6_SCLC_Effector_Immune</td><td>-0.4320402</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.3802492</td><td>-0.2729638</td><td>-0.6165697</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.4421974</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td>0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.323644</td><td>-0.4311868</td><td>-0.713442</td><td>-0.1225925</td><td>0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.4421974</td><td>-0.4841063</td><td>Low</td><td>Below20Percentile</td><td>vasc_low</td><td>vasc_low</td><td>12698</td><td>8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284</td><td>VA1_771</td><td>VA1_352</td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']</td><td>-0.3313744</td><td>-0.2534575</td><td>-0.3396933</td><td>-0.7045753</td><td>-0.1946746</td><td>-0.6988578</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]</td><td>-0.3655433</td><td>0.08738646</td><td>-0.01253212</td><td>-0.1992582</td><td>0.07594644</td><td>-0.1103883</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.3310987</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.5467328</td><td>-0.1167028</td><td>0.2619499</td><td>0.4510973</td><td>0.3857199</td><td>0.07931599</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073]</td><td>0.555634</td><td>0.2090081</td><td>-0.1464707</td><td>0.07915905</td><td>0.620314</td><td>-0.3870134</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]</td><td>-0.1416773</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.523967</td><td>-0.548042</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td>0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td>0.1411338</td><td>-0.4395605</td><td>-0.06779542</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.3230646</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.1434706</td><td>-0.3029</td><td>-0.2446501</td><td>-0.3597021</td><td>-0.2758603</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td><td>Macrophage_nn_mean</td><td>-0.42043885</td></tr>
	<tr><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6</td><td>6</td><td>VA1</td><td>1.001716</td><td>vasc_low</td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low</td><td>SCLC</td><td>SCLC</td><td>SCLC + Effector Immune Cells (T+NK)</td><td>6_SCLC_Effector_Immune</td><td>6_SCLC_Effector_Immune</td><td>-0.4320402</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.3802492</td><td>-0.2729638</td><td>-0.6165697</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.4421974</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td>0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.323644</td><td>-0.4311868</td><td>-0.713442</td><td>-0.1225925</td><td>0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.4421974</td><td>-0.4841063</td><td>Low</td><td>Below20Percentile</td><td>vasc_low</td><td>vasc_low</td><td>12698</td><td>8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284</td><td>VA1_771</td><td>VA1_352</td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']</td><td>-0.3313744</td><td>-0.2534575</td><td>-0.3396933</td><td>-0.7045753</td><td>-0.1946746</td><td>-0.6988578</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]</td><td>-0.3655433</td><td>0.08738646</td><td>-0.01253212</td><td>-0.1992582</td><td>0.07594644</td><td>-0.1103883</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.3310987</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.5467328</td><td>-0.1167028</td><td>0.2619499</td><td>0.4510973</td><td>0.3857199</td><td>0.07931599</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073]</td><td>0.555634</td><td>0.2090081</td><td>-0.1464707</td><td>0.07915905</td><td>0.620314</td><td>-0.3870134</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]</td><td>-0.1416773</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.523967</td><td>-0.548042</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td>0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td>0.1411338</td><td>-0.4395605</td><td>-0.06779542</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.3230646</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.1434706</td><td>-0.3029</td><td>-0.2446501</td><td>-0.3597021</td><td>-0.2758603</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td><td>T_Cell_nn_mean    </td><td>-0.08739819</td></tr>
</tbody>
</table>




```R
unique(nns_6adj_melt$nn_mean_Cell_Type_Label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Macrophage_nn_mean'</li><li>'T_Cell_nn_mean'</li><li>'CD8_T_Cell_nn_mean'</li><li>'NK_Cell_nn_mean'</li><li>'B_Cell_nn_mean'</li><li>'DC_nn_mean'</li><li>'Endothelial_Cell_nn_mean'</li><li>'LEC_nn_mean'</li><li>'MN4_EC_Phenotype_nn_mean'</li></ol>




```R
# get only groups of interest
temp <- nns_6adj_melt[nns_6adj_melt$Tumor_MHCI_label2 %in% c('Tumor_high_MHCI_high','Tumor_high_MHCI_low'),]
temp_VAB1 <- temp[temp$orig.ident %in% c('VA1','VB1'),]
temp_VA2 <- temp[temp$orig.ident %in% c('VA2'),]
# do t-test on groups of interest
celltype_vs_temp_stats1 <- temp %>% group_by(nn_mean_Cell_Type_Label) %>% summarise(p_TumorHigh_MHCI_HiVsLo=t.test(nn_mean_Cell_Type_Enrichment[Tumor_MHCI_label2=='Tumor_high_MHCI_high'], nn_mean_Cell_Type_Enrichment[Tumor_MHCI_label2=='Tumor_high_MHCI_low'], paired=FALSE)$p.value)
celltype_vs_temp_stats1_VA2 <- temp_VA2 %>% group_by(nn_mean_Cell_Type_Label) %>% summarise(p_TumorHigh_MHCI_HiVsLo_VA2=t.test(nn_mean_Cell_Type_Enrichment[Tumor_MHCI_label2=='Tumor_high_MHCI_high'], nn_mean_Cell_Type_Enrichment[Tumor_MHCI_label2=='Tumor_high_MHCI_low'], paired=FALSE)$p.value)
celltype_vs_temp_stats1_VAB1 <- temp_VAB1 %>% group_by(nn_mean_Cell_Type_Label) %>% summarise(p_TumorHigh_MHCI_HiVsLo_VAB1=t.test(nn_mean_Cell_Type_Enrichment[Tumor_MHCI_label2=='Tumor_high_MHCI_high'], nn_mean_Cell_Type_Enrichment[Tumor_MHCI_label2=='Tumor_high_MHCI_low'], paired=FALSE)$p.value)
# add multiple testing corrections
celltype_vs_temp_stats1 <- celltype_vs_temp_stats1 %>% mutate(padj_TumorHigh_MHCI_HiVsLo = p.adjust(p_TumorHigh_MHCI_HiVsLo, method = "BH"))
celltype_vs_temp_stats1_VA2 <- celltype_vs_temp_stats1_VA2 %>% mutate(padj_TumorHigh_MHCI_HiVsLo_VA2 = p.adjust(p_TumorHigh_MHCI_HiVsLo_VA2, method = "BH"))
celltype_vs_temp_stats1_VAB1 <- celltype_vs_temp_stats1_VAB1 %>% mutate(padj_TumorHigh_MHCI_HiVsLo_VAB1 = p.adjust(p_TumorHigh_MHCI_HiVsLo_VAB1, method = "BH"))
# combine results
stats_celltypesigs3_all <- merge(celltype_vs_temp_stats1, 
                                 celltype_vs_temp_stats1_VAB1, by.x='nn_mean_Cell_Type_Label', by.y='nn_mean_Cell_Type_Label')
stats_celltypesigs3_all <- merge(stats_celltypesigs3_all, 
                                 celltype_vs_temp_stats1_VA2, by.x='nn_mean_Cell_Type_Label', by.y='nn_mean_Cell_Type_Label')

print('t-test: MHCI high vs MHCI low, for Tumor High spots')
stats_celltypesigs3_all
```

    [1] "t-test: MHCI high vs MHCI low, for Tumor High spots"



<table class="dataframe">
<caption>A data.frame: 9 × 7</caption>
<thead>
	<tr><th scope=col>nn_mean_Cell_Type_Label</th><th scope=col>p_TumorHigh_MHCI_HiVsLo</th><th scope=col>padj_TumorHigh_MHCI_HiVsLo</th><th scope=col>p_TumorHigh_MHCI_HiVsLo_VAB1</th><th scope=col>padj_TumorHigh_MHCI_HiVsLo_VAB1</th><th scope=col>p_TumorHigh_MHCI_HiVsLo_VA2</th><th scope=col>padj_TumorHigh_MHCI_HiVsLo_VA2</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>B_Cell_nn_mean          </td><td> 1.093495e-32</td><td> 1.230181e-32</td><td>1.320055e-03</td><td>1.485062e-03</td><td> 1.701877e-42</td><td> 2.188128e-42</td></tr>
	<tr><td>CD8_T_Cell_nn_mean      </td><td>4.679451e-109</td><td>1.403835e-108</td><td>4.533938e-06</td><td>6.800907e-06</td><td> 6.619590e-83</td><td> 1.985877e-82</td></tr>
	<tr><td>DC_nn_mean              </td><td> 2.025534e-63</td><td> 2.604258e-63</td><td>3.068931e-11</td><td>1.381019e-10</td><td> 4.631159e-38</td><td> 4.631159e-38</td></tr>
	<tr><td>Endothelial_Cell_nn_mean</td><td> 1.804236e-76</td><td> 3.247625e-76</td><td>1.588998e-06</td><td>2.860196e-06</td><td> 9.381710e-52</td><td> 1.407257e-51</td></tr>
	<tr><td>LEC_nn_mean             </td><td> 3.534735e-24</td><td> 3.534735e-24</td><td>3.308574e-01</td><td>3.308574e-01</td><td> 3.576106e-39</td><td> 4.023119e-39</td></tr>
	<tr><td>Macrophage_nn_mean      </td><td> 2.720064e-65</td><td> 4.080095e-65</td><td>6.680895e-07</td><td>1.503201e-06</td><td> 4.962532e-71</td><td> 1.116570e-70</td></tr>
	<tr><td>MN4_EC_Phenotype_nn_mean</td><td>2.720994e-171</td><td>2.448894e-170</td><td>6.919226e-38</td><td>6.227304e-37</td><td>3.341982e-107</td><td>3.007784e-106</td></tr>
	<tr><td>NK_Cell_nn_mean         </td><td> 8.057661e-83</td><td> 1.812974e-82</td><td>1.701274e-04</td><td>2.187352e-04</td><td> 3.949403e-70</td><td> 7.108925e-70</td></tr>
	<tr><td>T_Cell_nn_mean          </td><td>7.283240e-128</td><td>3.277458e-127</td><td>1.880704e-09</td><td>5.642113e-09</td><td> 2.742600e-96</td><td> 1.234170e-95</td></tr>
</tbody>
</table>




```R
write.csv(stats_celltypesigs3_all, 'Output/12_Plot_Neighborhoods/STATS_nns_6adj_melt_TumorMHCI_Hi_vs_Low.csv', row.names = FALSE)
```

## 3. If you’re a vasc_high spot with high vs low activation (MN4_EC), do you tend to have NK/T cells near you?

# plot


```R
options(repr.plot.width=20, repr.plot.height=16)

NN_celltype_vs_TumorMHCILabel_plot3 <- ggplot(nns_6adj_melt[nns_6adj_melt$vasc_label == 'vasc_high',], 
                                          aes(x=MN4_EC_Phenotype_Label, y=nn_mean_Cell_Type_Enrichment, color=orig.ident)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey') + 
  geom_beeswarm(cex = 0.25, size=1) +
  #geom_point() + 
  facet_wrap(vars(nn_mean_Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Nearest Neighbors Cell Type GSVA Enrichment vs. center MN4_EC")
             
NN_celltype_vs_TumorMHCILabel_plot3

# save plot
pdf(file = "Output/12_Plot_Neighborhoods/AB_NN_cell_type_GSVA_per_MN4_EC_Label.pdf",
    width = 1800,
    height = 1400)
print(NN_celltype_vs_TumorMHCILabel_plot3)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_36_1.png)
    


# VAB1 and VA2 separately


```R
options(repr.plot.width=20, repr.plot.height=16)

NN_celltype_vs_TumorMHCILabel_plot3 <- ggplot(nns_6adj_melt_VAB1[nns_6adj_melt_VAB1$vasc_label == 'vasc_high',], 
                                          aes(x=MN4_EC_Phenotype_Label, y=nn_mean_Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey') + 
  geom_beeswarm(cex = 0.4, size=1) +
  #geom_point() + 
  facet_wrap(vars(nn_mean_Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "VAB1 Nearest Neighbors Cell Type GSVA Enrichment vs. center MN4_EC")
             
NN_celltype_vs_TumorMHCILabel_plot3

# save plot
pdf(file = "Output/12_Plot_Neighborhoods/AB_NN_cell_type_GSVA_per_MN4_EC_Label_VAB1.pdf",
    width = 1800,
    height = 1400)
print(NN_celltype_vs_TumorMHCILabel_plot3)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_38_1.png)
    



```R
options(repr.plot.width=20, repr.plot.height=16)

NN_celltype_vs_TumorMHCILabel_plot3 <- ggplot(nns_6adj_melt_VA2[nns_6adj_melt_VA2$vasc_label == 'vasc_high',], 
                                          aes(x=MN4_EC_Phenotype_Label, y=nn_mean_Cell_Type_Enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5), fill='light grey') + 
  geom_beeswarm(cex = 0.4, size=1) +
  #geom_point() + 
  facet_wrap(vars(nn_mean_Cell_Type_Label)) +
  theme_bw() +
  theme(text=element_text(size=20),
       plot.title = element_text(size=24),
       strip.text = element_text(size=18)) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "VA2 Nearest Neighbors Cell Type GSVA Enrichment vs. center MN4_EC")
             
NN_celltype_vs_TumorMHCILabel_plot3

# save plot
pdf(file = "Output/12_Plot_Neighborhoods/AB_NN_cell_type_GSVA_per_MN4_EC_Label_VA2.pdf",
    width = 1800,
    height = 1400)
print(NN_celltype_vs_TumorMHCILabel_plot3)
dev.off()
```


<strong>pdf:</strong> 2



    
![png](output_39_1.png)
    



```R
# save data
write.csv(nns_6adj_melt[nns_6adj_melt$vasc_label == 'vasc_high',], 'Processing/10_Neighborhood_Analysis/nns_6adj_melt_VascHigh_MN4_Hi_vs_Low.csv', row.names = FALSE)
```


```R
nns_6adj_melt %>% head(2)
```


<table class="dataframe">
<caption>A tibble: 2 × 134</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>Tissue_Slice</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Exhaustion</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M2_Macrophage</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>MN4_EC_Percentile_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>vasc_MN4_percentile_label</th><th scope=col>x</th><th scope=col>y</th><th scope=col>cell</th><th scope=col>id</th><th scope=col>nn1_id</th><th scope=col>nn2_id</th><th scope=col>nn3_id</th><th scope=col>nn4_id</th><th scope=col>nn5_id</th><th scope=col>nn6_id</th><th scope=col>num_adjacent_nns</th><th scope=col>nn_list</th><th scope=col>Macrophage_nn_1</th><th scope=col>Macrophage_nn_2</th><th scope=col>Macrophage_nn_3</th><th scope=col>Macrophage_nn_4</th><th scope=col>Macrophage_nn_5</th><th scope=col>Macrophage_nn_6</th><th scope=col>Macrophage_nn_vals</th><th scope=col>T_Cell_nn_1</th><th scope=col>T_Cell_nn_2</th><th scope=col>T_Cell_nn_3</th><th scope=col>T_Cell_nn_4</th><th scope=col>T_Cell_nn_5</th><th scope=col>T_Cell_nn_6</th><th scope=col>T_Cell_nn_vals</th><th scope=col>CD8_T_Cell_nn_1</th><th scope=col>CD8_T_Cell_nn_2</th><th scope=col>CD8_T_Cell_nn_3</th><th scope=col>CD8_T_Cell_nn_4</th><th scope=col>CD8_T_Cell_nn_5</th><th scope=col>CD8_T_Cell_nn_6</th><th scope=col>CD8_T_Cell_nn_vals</th><th scope=col>NK_Cell_nn_1</th><th scope=col>NK_Cell_nn_2</th><th scope=col>NK_Cell_nn_3</th><th scope=col>NK_Cell_nn_4</th><th scope=col>NK_Cell_nn_5</th><th scope=col>NK_Cell_nn_6</th><th scope=col>NK_Cell_nn_vals</th><th scope=col>B_Cell_nn_1</th><th scope=col>B_Cell_nn_2</th><th scope=col>B_Cell_nn_3</th><th scope=col>B_Cell_nn_4</th><th scope=col>B_Cell_nn_5</th><th scope=col>B_Cell_nn_6</th><th scope=col>B_Cell_nn_vals</th><th scope=col>DC_nn_1</th><th scope=col>DC_nn_2</th><th scope=col>DC_nn_3</th><th scope=col>DC_nn_4</th><th scope=col>DC_nn_5</th><th scope=col>DC_nn_6</th><th scope=col>DC_nn_vals</th><th scope=col>Endothelial_Cell_nn_1</th><th scope=col>Endothelial_Cell_nn_2</th><th scope=col>Endothelial_Cell_nn_3</th><th scope=col>Endothelial_Cell_nn_4</th><th scope=col>Endothelial_Cell_nn_5</th><th scope=col>Endothelial_Cell_nn_6</th><th scope=col>Endothelial_Cell_nn_vals</th><th scope=col>LEC_nn_1</th><th scope=col>LEC_nn_2</th><th scope=col>LEC_nn_3</th><th scope=col>LEC_nn_4</th><th scope=col>LEC_nn_5</th><th scope=col>LEC_nn_6</th><th scope=col>LEC_nn_vals</th><th scope=col>MN4_EC_Phenotype_nn_1</th><th scope=col>MN4_EC_Phenotype_nn_2</th><th scope=col>MN4_EC_Phenotype_nn_3</th><th scope=col>MN4_EC_Phenotype_nn_4</th><th scope=col>MN4_EC_Phenotype_nn_5</th><th scope=col>MN4_EC_Phenotype_nn_6</th><th scope=col>MN4_EC_Phenotype_nn_vals</th><th scope=col>nn_mean_Cell_Type_Label</th><th scope=col>nn_mean_Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6</td><td>6</td><td>VA1</td><td>1.001716</td><td>vasc_low</td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low</td><td>SCLC</td><td>SCLC</td><td>SCLC + Effector Immune Cells (T+NK)</td><td>6_SCLC_Effector_Immune</td><td>6_SCLC_Effector_Immune</td><td>-0.4320402</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.3802492</td><td>-0.2729638</td><td>-0.6165697</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.4421974</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td>0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.323644</td><td>-0.4311868</td><td>-0.713442</td><td>-0.1225925</td><td>0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.4421974</td><td>-0.4841063</td><td>Low</td><td>Below20Percentile</td><td>vasc_low</td><td>vasc_low</td><td>12698</td><td>8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284</td><td>VA1_771</td><td>VA1_352</td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']</td><td>-0.3313744</td><td>-0.2534575</td><td>-0.3396933</td><td>-0.7045753</td><td>-0.1946746</td><td>-0.6988578</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]</td><td>-0.3655433</td><td>0.08738646</td><td>-0.01253212</td><td>-0.1992582</td><td>0.07594644</td><td>-0.1103883</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.3310987</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.5467328</td><td>-0.1167028</td><td>0.2619499</td><td>0.4510973</td><td>0.3857199</td><td>0.07931599</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073]</td><td>0.555634</td><td>0.2090081</td><td>-0.1464707</td><td>0.07915905</td><td>0.620314</td><td>-0.3870134</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]</td><td>-0.1416773</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.523967</td><td>-0.548042</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td>0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td>0.1411338</td><td>-0.4395605</td><td>-0.06779542</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.3230646</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.1434706</td><td>-0.3029</td><td>-0.2446501</td><td>-0.3597021</td><td>-0.2758603</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td><td>Macrophage_nn_mean</td><td>-0.42043885</td></tr>
	<tr><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6</td><td>6</td><td>VA1</td><td>1.001716</td><td>vasc_low</td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low</td><td>SCLC</td><td>SCLC</td><td>SCLC + Effector Immune Cells (T+NK)</td><td>6_SCLC_Effector_Immune</td><td>6_SCLC_Effector_Immune</td><td>-0.4320402</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.3802492</td><td>-0.2729638</td><td>-0.6165697</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.2841903</td><td>-0.4421974</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td>0.8164109</td><td>-0.4222185</td><td>-0.5429759</td><td>-0.323644</td><td>-0.4311868</td><td>-0.713442</td><td>-0.1225925</td><td>0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.4421974</td><td>-0.4841063</td><td>Low</td><td>Below20Percentile</td><td>vasc_low</td><td>vasc_low</td><td>12698</td><td>8743</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_1</td><td>VA1_284</td><td>VA1_771</td><td>VA1_352</td><td>VA1_1783</td><td>VA1_2738</td><td>VA1_1604</td><td>6</td><td>['VA1_284', 'VA1_771', 'VA1_352', 'VA1_1783', 'VA1_2738', 'VA1_1604']</td><td>-0.3313744</td><td>-0.2534575</td><td>-0.3396933</td><td>-0.7045753</td><td>-0.1946746</td><td>-0.6988578</td><td>[-0.33137442702708, -0.253457547402507, -0.339693309699658, -0.704575342834968, -0.194674612981697, -0.698857839643165]</td><td>-0.3655433</td><td>0.08738646</td><td>-0.01253212</td><td>-0.1992582</td><td>0.07594644</td><td>-0.1103883</td><td>[-0.365543343516556, 0.0873864569432376, -0.0125321158506752, -0.199258240763675, 0.075946441965502, -0.110388310648632]</td><td>0.1532445</td><td>-0.2668601</td><td>-0.3310987</td><td>-0.2024215</td><td>-0.2665599</td><td>-0.1514909</td><td>[0.153244482371191, -0.266860116069477, -0.331098659195336, -0.202421452871575, -0.266559935961412, -0.151490894536615]</td><td>0.5467328</td><td>-0.1167028</td><td>0.2619499</td><td>0.4510973</td><td>0.3857199</td><td>0.07931599</td><td>[0.546732833950551, -0.116702778523135, 0.261949922898662, 0.451097284530692, 0.385719915781331, 0.0793159900391073]</td><td>0.555634</td><td>0.2090081</td><td>-0.1464707</td><td>0.07915905</td><td>0.620314</td><td>-0.3870134</td><td>[0.555633961432688, 0.209008116797947, -0.146470705958925, 0.0791590535930388, 0.620313971068745, -0.387013420681998]</td><td>-0.1416773</td><td>-0.6074761</td><td>-0.1887733</td><td>-0.04428771</td><td>-0.523967</td><td>-0.548042</td><td>[-0.141677294770587, -0.607476141948519, -0.188773278840846, -0.0442877125238972, -0.523967008528942, -0.548042038055844]</td><td>0.1655037</td><td>-0.2721193</td><td>-0.2911405</td><td>0.1411338</td><td>-0.4395605</td><td>-0.06779542</td><td>[0.165503689409276, -0.272119318571468, -0.291140519099978, 0.141133788722132, -0.439560478618792, -0.067795415267941]</td><td>-0.3917784</td><td>-0.2591518</td><td>-0.3230646</td><td>-0.1908382</td><td>-0.2596519</td><td>-0.2107937</td><td>[-0.391778355670959, -0.259151830365927, -0.323064612922424, -0.190838167633395, -0.259651930385931, -0.210793700590192]</td><td>-0.05266449</td><td>-0.1434706</td><td>-0.3029</td><td>-0.2446501</td><td>-0.3597021</td><td>-0.2758603</td><td>[-0.0526644931955577, -0.143470613692345, -0.302899967001359, -0.244650143833782, -0.359702126487288, -0.275860343619002]</td><td>T_Cell_nn_mean    </td><td>-0.08739819</td></tr>
</tbody>
</table>




```R
unique(nns_6adj_melt$MN4_EC_Phenotype_Label)
unique(nns_6adj_melt$vasc_MN4_label)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Low'</li><li>'High'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'vasc_low'</li><li>'vasc_high_MN4_high'</li><li>'vasc_high_MN4_low'</li></ol>



# save subset to only relevant columns. 


```R
vasc_high_data <- nns_6adj_melt[(nns_6adj_melt$vasc_label == 'vasc_high'),]

vasc_high_data_tosave = vasc_high_data %>% select(Barcode, orig.ident, vasc_label, MN4_EC_Phenotype_Label, vasc_MN4_label, nn_mean_Cell_Type_Label, nn_mean_Cell_Type_Enrichment)

# drop duplicates
vasc_high_data_tosave <- vasc_high_data_tosave[!duplicated(vasc_high_data_tosave), ]

vasc_high_data_tosave[1:3,]
```


<table class="dataframe">
<caption>A tibble: 3 × 7</caption>
<thead>
	<tr><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>vasc_label</th><th scope=col>MN4_EC_Phenotype_Label</th><th scope=col>vasc_MN4_label</th><th scope=col>nn_mean_Cell_Type_Label</th><th scope=col>nn_mean_Cell_Type_Enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>vasc_high_MN4_high</td><td>Macrophage_nn_mean</td><td>-0.3448332</td></tr>
	<tr><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>vasc_high_MN4_high</td><td>T_Cell_nn_mean    </td><td>-0.1737913</td></tr>
	<tr><td>VA1_AAACAGCTTTCAGAAG-1</td><td>VA1</td><td>vasc_high</td><td>High</td><td>vasc_high_MN4_high</td><td>CD8_T_Cell_nn_mean</td><td> 0.2005778</td></tr>
</tbody>
</table>




```R
dplyr::count(vasc_high_data_tosave, orig.ident, vasc_MN4_label)
```


<table class="dataframe">
<caption>A tibble: 6 × 3</caption>
<thead>
	<tr><th scope=col>orig.ident</th><th scope=col>vasc_MN4_label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1</td><td>vasc_high_MN4_high</td><td> 3960</td></tr>
	<tr><td>VA1</td><td>vasc_high_MN4_low </td><td> 1872</td></tr>
	<tr><td>VA2</td><td>vasc_high_MN4_high</td><td>12456</td></tr>
	<tr><td>VA2</td><td>vasc_high_MN4_low </td><td> 2061</td></tr>
	<tr><td>VB1</td><td>vasc_high_MN4_high</td><td> 3141</td></tr>
	<tr><td>VB1</td><td>vasc_high_MN4_low </td><td> 1620</td></tr>
</tbody>
</table>




```R
dim(vasc_high_data_tosave)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>25110</li><li>7</li></ol>




```R
# save
write.csv(vasc_high_data_tosave, 'Output/12_Plot_Neighborhoods/nns_6adj_vasc_high_cell_type_enrichment_by_MN4.csv', row.names = FALSE)

```

# save a version that's unmelted


```R
vasc_high_data_tosave_unmelted <- vasc_high_data_tosave %>% pivot_wider(names_from = nn_mean_Cell_Type_Label, values_from = nn_mean_Cell_Type_Enrichment)
```


```R
# save
write.csv(vasc_high_data_tosave_unmelted, 'Output/12_Plot_Neighborhoods/nns_6adj_vasc_high_cell_type_enrichment_by_MN4_wideform.csv', row.names = FALSE)

```

# stats

stats_celltypesigs3_all <- nns_6adj_melt_lab[nns_6adj_melt_lab$vasc_label %in% c('vasc_high'),] %>% group_by(nn_mean_Cell_Type_Label) %>% summarise(p_NN_GSVA_for_VascHigh_MN4_EC_HiVsLo_AB=wilcox.test(nn_mean_Cell_Type_Enrichment[MN4_EC_Label=='MN4_EC_High'], nn_mean_Cell_Type_Enrichment[MN4_EC_Label=='MN4_EC_Low'], paired=FALSE)$p.value)

print('Mann-Whitney U Test: Nearest Neighbor Cell Type Enrichments for MN4_EC high vs MN4_EC low center spots, vasc_high only')
stats_celltypesigs3_all

```R
# get only groups of interest
temp <- nns_6adj_melt[nns_6adj_melt$vasc_label %in% c('vasc_high'),]
temp_VAB1 <- temp[temp$orig.ident %in% c('VA1','VB1'),]
temp_VA2 <- temp[temp$orig.ident %in% c('VA2'),]
# do t-test on groups of interest
celltype_vs_temp_stats1 <- temp %>% group_by(nn_mean_Cell_Type_Label) %>% summarise(p_VascHigh_MN4_HiVsLo=t.test(nn_mean_Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='High'], nn_mean_Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='Low'], paired=FALSE)$p.value)
celltype_vs_temp_stats1_VA2 <- temp_VA2 %>% group_by(nn_mean_Cell_Type_Label) %>% summarise(p_VascHigh_MN4_HiVsLo_VA2=t.test(nn_mean_Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='High'], nn_mean_Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='Low'], paired=FALSE)$p.value)
celltype_vs_temp_stats1_VAB1 <- temp_VAB1 %>% group_by(nn_mean_Cell_Type_Label) %>% summarise(p_VascHigh_MN4_HiVsLo_VAB1=t.test(nn_mean_Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='High'], nn_mean_Cell_Type_Enrichment[MN4_EC_Phenotype_Label=='Low'], paired=FALSE)$p.value)
# add multiple testing corrections
celltype_vs_temp_stats1 <- celltype_vs_temp_stats1 %>% mutate(padj_VascHigh_MN4_HiVsLo = p.adjust(p_VascHigh_MN4_HiVsLo, method = "BH"))
celltype_vs_temp_stats1_VA2 <- celltype_vs_temp_stats1_VA2 %>% mutate(padj_VascHigh_MN4_HiVsLo_VA2 = p.adjust(p_VascHigh_MN4_HiVsLo_VA2, method = "BH"))
celltype_vs_temp_stats1_VAB1 <- celltype_vs_temp_stats1_VAB1 %>% mutate(padj_VascHigh_MN4_HiVsLo_VAB1 = p.adjust(p_VascHigh_MN4_HiVsLo_VAB1, method = "BH"))
# combine results
stats_celltypesigs3_all <- merge(celltype_vs_temp_stats1, 
                                 celltype_vs_temp_stats1_VAB1, by.x='nn_mean_Cell_Type_Label', by.y='nn_mean_Cell_Type_Label')
stats_celltypesigs3_all <- merge(stats_celltypesigs3_all, 
                                 celltype_vs_temp_stats1_VA2, by.x='nn_mean_Cell_Type_Label', by.y='nn_mean_Cell_Type_Label')

print('t-test: MN4_EC high vs MN4_EC low, for vasc_high spots')
stats_celltypesigs3_all
```

    [1] "t-test: MN4_EC high vs MN4_EC low, for vasc_high spots"



<table class="dataframe">
<caption>A data.frame: 9 × 7</caption>
<thead>
	<tr><th scope=col>nn_mean_Cell_Type_Label</th><th scope=col>p_VascHigh_MN4_HiVsLo</th><th scope=col>padj_VascHigh_MN4_HiVsLo</th><th scope=col>p_VascHigh_MN4_HiVsLo_VAB1</th><th scope=col>padj_VascHigh_MN4_HiVsLo_VAB1</th><th scope=col>p_VascHigh_MN4_HiVsLo_VA2</th><th scope=col>padj_VascHigh_MN4_HiVsLo_VA2</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>B_Cell_nn_mean          </td><td> 4.161439e-75</td><td> 7.490590e-75</td><td>1.577460e-14</td><td>3.549285e-14</td><td>3.563609e-31</td><td>5.345413e-31</td></tr>
	<tr><td>CD8_T_Cell_nn_mean      </td><td> 9.768843e-76</td><td> 2.197990e-75</td><td>2.068516e-12</td><td>2.659521e-12</td><td>2.404106e-39</td><td>5.409239e-39</td></tr>
	<tr><td>DC_nn_mean              </td><td> 1.125289e-38</td><td> 1.446800e-38</td><td>4.368569e-17</td><td>1.965856e-16</td><td>2.619493e-25</td><td>2.946929e-25</td></tr>
	<tr><td>Endothelial_Cell_nn_mean</td><td>2.240616e-108</td><td>1.008277e-107</td><td>2.529793e-16</td><td>7.589378e-16</td><td>4.606428e-46</td><td>1.381928e-45</td></tr>
	<tr><td>LEC_nn_mean             </td><td> 5.786673e-26</td><td> 6.510007e-26</td><td>8.995761e-01</td><td>8.995761e-01</td><td>5.264176e-17</td><td>5.264176e-17</td></tr>
	<tr><td>Macrophage_nn_mean      </td><td> 2.357136e-56</td><td> 3.535704e-56</td><td>7.789765e-14</td><td>1.402158e-13</td><td>7.714422e-29</td><td>9.918543e-29</td></tr>
	<tr><td>MN4_EC_Phenotype_nn_mean</td><td>1.637291e-185</td><td>1.473562e-184</td><td>2.744431e-69</td><td>2.469988e-68</td><td>4.815232e-62</td><td>4.333709e-61</td></tr>
	<tr><td>NK_Cell_nn_mean         </td><td> 5.703303e-17</td><td> 5.703303e-17</td><td>5.824012e-02</td><td>6.552013e-02</td><td>2.292753e-37</td><td>4.126956e-37</td></tr>
	<tr><td>T_Cell_nn_mean          </td><td> 2.325440e-93</td><td> 6.976320e-93</td><td>2.618700e-13</td><td>3.928050e-13</td><td>5.101913e-48</td><td>2.295861e-47</td></tr>
</tbody>
</table>




```R
write.csv(stats_celltypesigs3_all, 'Output/12_Plot_Neighborhoods/STATS_nns_6adj_melt_VascHigh_MN4_Hi_vs_Low.csv', row.names = FALSE)
```


```R
# save data

```


```R

```


```R

```


```R

```
