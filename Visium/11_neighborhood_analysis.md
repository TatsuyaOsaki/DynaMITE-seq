# Neighborhood analysis on visium spots

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


```python
import pandas as pd
import numpy as np
import h5py, os, re, math, sys
import numpy as np
from scipy import spatial

#np.set_printoptions(threshold=sys.maxsize)
#np.set_printoptions(threshold=10000)

pd.set_option('display.max_columns', None)
```


```python
# import data
#ab_nn_pre = pd.read_csv('Neighborhood_Analysis_20240711/ab_nn_pre.csv')
va1_nn_pre = pd.read_csv('Processing/10_Neighborhood_Analysis/va1_nn_pre.csv')
vb1_nn_pre = pd.read_csv('Processing/10_Neighborhood_Analysis/vb1_nn_pre.csv')
va2_nn_pre = pd.read_csv('Processing/10_Neighborhood_Analysis/va2_nn_pre.csv')
```


```python
va1_nn_pre.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VA1_AAACAACGAATAGTTC-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>3625</td>
      <td>2571</td>
      <td>round1</td>
      <td>3649</td>
      <td>2571</td>
      <td>10</td>
      <td>10</td>
      <td>VA1</td>
      <td>2.648104</td>
      <td>vasc_high</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>2</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_high</td>
      <td>Respiratory epithelium</td>
      <td>Respiratory epithelium</td>
      <td>Respiratory Epithelium (Normal Lung)</td>
      <td>10_Respiratory_Epithelium_Lung</td>
      <td>7_10_Normal_Lung</td>
      <td>-0.056807</td>
      <td>-0.259847</td>
      <td>0.296748</td>
      <td>-0.038431</td>
      <td>0.513527</td>
      <td>-0.007136</td>
      <td>0.412196</td>
      <td>0.404657</td>
      <td>-0.158465</td>
      <td>0.690663</td>
      <td>-0.046498</td>
      <td>0.110800</td>
      <td>0.102126</td>
      <td>0.654483</td>
      <td>-0.176035</td>
      <td>-0.321725</td>
      <td>-0.469550</td>
      <td>0.189499</td>
      <td>0.139214</td>
      <td>-0.165127</td>
      <td>0.015559</td>
      <td>-0.755217</td>
      <td>0.633032</td>
      <td>-0.326533</td>
      <td>-0.046498</td>
      <td>0.031719</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_high_MN4_high</td>
      <td>vasc_high_MN4_upper80</td>
      <td>30730</td>
      <td>26443</td>
      <td>AAACAACGAATAGTTC-1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VA1_AAACAAGTATCTCCCA-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>5804</td>
      <td>3329</td>
      <td>round1</td>
      <td>4956</td>
      <td>3308</td>
      <td>6</td>
      <td>6</td>
      <td>VA1</td>
      <td>1.001716</td>
      <td>vasc_low</td>
      <td>13</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>2</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_low</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.432040</td>
      <td>-0.363696</td>
      <td>-0.665610</td>
      <td>-0.380249</td>
      <td>-0.272964</td>
      <td>-0.616570</td>
      <td>-0.133588</td>
      <td>-0.277191</td>
      <td>-0.532077</td>
      <td>-0.284190</td>
      <td>-0.442197</td>
      <td>-0.206272</td>
      <td>-0.292946</td>
      <td>-0.220110</td>
      <td>0.816411</td>
      <td>-0.422218</td>
      <td>-0.542976</td>
      <td>-0.323644</td>
      <td>-0.431187</td>
      <td>-0.713442</td>
      <td>-0.122593</td>
      <td>0.079799</td>
      <td>-0.348075</td>
      <td>-0.420942</td>
      <td>-0.442197</td>
      <td>-0.484106</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>12698</td>
      <td>8743</td>
      <td>AAACAAGTATCTCCCA-1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VA1_AAACAATCTACTAGCA-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>4980</td>
      <td>3146</td>
      <td>round1</td>
      <td>4711</td>
      <td>3146</td>
      <td>6</td>
      <td>6</td>
      <td>VA1</td>
      <td>1.101286</td>
      <td>vasc_low</td>
      <td>19</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_low</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.312254</td>
      <td>-0.142246</td>
      <td>0.215974</td>
      <td>-0.522873</td>
      <td>-0.254499</td>
      <td>0.344170</td>
      <td>-0.096425</td>
      <td>0.015784</td>
      <td>-0.519466</td>
      <td>-0.273987</td>
      <td>-0.296916</td>
      <td>-0.112254</td>
      <td>-0.139708</td>
      <td>0.574189</td>
      <td>-0.242849</td>
      <td>-0.389510</td>
      <td>-0.529344</td>
      <td>-0.039327</td>
      <td>0.055921</td>
      <td>-0.707380</td>
      <td>-0.091447</td>
      <td>0.049511</td>
      <td>-0.304126</td>
      <td>0.924392</td>
      <td>-0.296916</td>
      <td>-0.466480</td>
      <td>High</td>
      <td>Low</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>29633</td>
      <td>20870</td>
      <td>AAACAATCTACTAGCA-1</td>
    </tr>
  </tbody>
</table>
</div>




```python
vb1_nn_pre.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VB1_AAACAACGAATAGTTC-1</td>
      <td>VB1</td>
      <td>VB1</td>
      <td>7174</td>
      <td>3783</td>
      <td>round1</td>
      <td>5068</td>
      <td>3601</td>
      <td>6</td>
      <td>6</td>
      <td>VB1</td>
      <td>0.000000</td>
      <td>vasc_low</td>
      <td>20</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_neg</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.409190</td>
      <td>-0.293434</td>
      <td>0.015793</td>
      <td>-0.623296</td>
      <td>-0.294338</td>
      <td>-0.538970</td>
      <td>0.160997</td>
      <td>-0.440288</td>
      <td>-0.049946</td>
      <td>-0.219445</td>
      <td>-0.184969</td>
      <td>0.010123</td>
      <td>0.001746</td>
      <td>0.523236</td>
      <td>-0.383538</td>
      <td>0.263978</td>
      <td>-0.317542</td>
      <td>-0.189507</td>
      <td>0.171521</td>
      <td>-0.565283</td>
      <td>-0.448105</td>
      <td>-0.309486</td>
      <td>-0.363641</td>
      <td>-0.373975</td>
      <td>-0.184969</td>
      <td>0.090032</td>
      <td>Low</td>
      <td>Low</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>31283</td>
      <td>26488</td>
      <td>AAACAACGAATAGTTC-1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VB1_AAACAAGTATCTCCCA-1</td>
      <td>VB1</td>
      <td>VB1</td>
      <td>3817</td>
      <td>2574</td>
      <td>round1</td>
      <td>3824</td>
      <td>2572</td>
      <td>1</td>
      <td>1</td>
      <td>VB1</td>
      <td>3.117370</td>
      <td>vasc_high</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>4</td>
      <td>MHCI_high</td>
      <td>MHCI_high</td>
      <td>Tumor_low_MHCI_high</td>
      <td>Tumor_low</td>
      <td>vasc_high</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC</td>
      <td>1_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>0.124843</td>
      <td>0.508282</td>
      <td>-0.370904</td>
      <td>-0.208191</td>
      <td>-0.135008</td>
      <td>0.254889</td>
      <td>0.354699</td>
      <td>0.618041</td>
      <td>0.174323</td>
      <td>0.644896</td>
      <td>0.223238</td>
      <td>0.065979</td>
      <td>0.010680</td>
      <td>0.645669</td>
      <td>0.955596</td>
      <td>-0.325863</td>
      <td>0.259486</td>
      <td>-0.095793</td>
      <td>-0.280748</td>
      <td>-0.611245</td>
      <td>0.417938</td>
      <td>-0.344595</td>
      <td>-0.302494</td>
      <td>-0.264353</td>
      <td>0.223238</td>
      <td>0.694163</td>
      <td>Low</td>
      <td>Low</td>
      <td>vasc_high_MN4_low</td>
      <td>vasc_low</td>
      <td>13243</td>
      <td>8788</td>
      <td>AAACAAGTATCTCCCA-1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VB1_AAACACCAATAACTGC-1</td>
      <td>VB1</td>
      <td>VB1</td>
      <td>6137</td>
      <td>3622</td>
      <td>round1</td>
      <td>4990</td>
      <td>3585</td>
      <td>0</td>
      <td>0</td>
      <td>VB1</td>
      <td>0.966779</td>
      <td>vasc_low</td>
      <td>17</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>4</td>
      <td>MHCI_high</td>
      <td>MHCI_high</td>
      <td>Tumor_high_MHCI_high</td>
      <td>Tumor_high_MHCI_high</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Lymphocytes B Cells (and Plasma Cells)</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>-0.223284</td>
      <td>0.124065</td>
      <td>0.148279</td>
      <td>-0.482537</td>
      <td>-0.294707</td>
      <td>-0.344041</td>
      <td>-0.235293</td>
      <td>0.083866</td>
      <td>-0.567738</td>
      <td>-0.239274</td>
      <td>-0.014502</td>
      <td>0.044275</td>
      <td>0.105805</td>
      <td>0.062049</td>
      <td>-0.386939</td>
      <td>-0.488036</td>
      <td>-0.447651</td>
      <td>0.072677</td>
      <td>0.316386</td>
      <td>0.408941</td>
      <td>0.290778</td>
      <td>-0.211729</td>
      <td>-0.364546</td>
      <td>-0.379476</td>
      <td>-0.014502</td>
      <td>0.143517</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>10072</td>
      <td>25947</td>
      <td>AAACACCAATAACTGC-1</td>
    </tr>
  </tbody>
</table>
</div>




```python
va2_nn_pre.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VA2_AAACAAGTATCTCCCA-1</td>
      <td>VA2</td>
      <td>VA2</td>
      <td>45867</td>
      <td>8515</td>
      <td>round2</td>
      <td>24897</td>
      <td>7041</td>
      <td>5</td>
      <td>5</td>
      <td>VA2</td>
      <td>1.296336</td>
      <td>vasc_low</td>
      <td>23</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC</td>
      <td>5_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>-0.450968</td>
      <td>-0.448330</td>
      <td>-0.691477</td>
      <td>-0.374355</td>
      <td>-0.693455</td>
      <td>-0.806181</td>
      <td>-0.424467</td>
      <td>-0.261943</td>
      <td>-0.684682</td>
      <td>-0.765453</td>
      <td>-0.199923</td>
      <td>-0.133973</td>
      <td>-0.213757</td>
      <td>-0.290269</td>
      <td>0.029304</td>
      <td>-0.261245</td>
      <td>-0.741362</td>
      <td>-0.395326</td>
      <td>-0.297578</td>
      <td>0.182372</td>
      <td>-0.583172</td>
      <td>0.484041</td>
      <td>-0.332599</td>
      <td>-0.704982</td>
      <td>-0.199923</td>
      <td>-0.156565</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>7764.085688</td>
      <td>20534.804310</td>
      <td>AAACAAGTATCTCCCA-1</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VA2_AAACACCAATAACTGC-1</td>
      <td>VA2</td>
      <td>VA2</td>
      <td>11849</td>
      <td>5082</td>
      <td>round2</td>
      <td>22823</td>
      <td>5881</td>
      <td>0</td>
      <td>0</td>
      <td>VA2</td>
      <td>1.223823</td>
      <td>vasc_low</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>6</td>
      <td>MHCI_high</td>
      <td>MHCI_high</td>
      <td>Tumor_low_MHCI_high</td>
      <td>Tumor_low</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Lymphocytes B Cells (and Plasma Cells)</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>0.465755</td>
      <td>0.476869</td>
      <td>0.801786</td>
      <td>0.685619</td>
      <td>0.728101</td>
      <td>-0.440744</td>
      <td>-0.355418</td>
      <td>0.451091</td>
      <td>0.679239</td>
      <td>0.565557</td>
      <td>0.301495</td>
      <td>0.099595</td>
      <td>0.168763</td>
      <td>0.394440</td>
      <td>0.279309</td>
      <td>0.396732</td>
      <td>0.534097</td>
      <td>0.368466</td>
      <td>0.618810</td>
      <td>0.647113</td>
      <td>0.167148</td>
      <td>-0.005985</td>
      <td>0.455600</td>
      <td>0.641533</td>
      <td>0.301495</td>
      <td>0.153607</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>24672.634008</td>
      <td>23846.305076</td>
      <td>AAACACCAATAACTGC-1</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VA2_AAACAGAGCGACTCCT-1</td>
      <td>VA2</td>
      <td>VA2</td>
      <td>43834</td>
      <td>8662</td>
      <td>round2</td>
      <td>24885</td>
      <td>7403</td>
      <td>5</td>
      <td>5</td>
      <td>VA2</td>
      <td>0.205495</td>
      <td>vasc_low</td>
      <td>24</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC</td>
      <td>5_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>-0.587803</td>
      <td>-0.647125</td>
      <td>-0.203480</td>
      <td>-0.608944</td>
      <td>-0.603551</td>
      <td>-0.840184</td>
      <td>-0.328385</td>
      <td>-0.589335</td>
      <td>-0.659993</td>
      <td>-0.799260</td>
      <td>-0.441499</td>
      <td>-0.277116</td>
      <td>-0.349001</td>
      <td>-0.654693</td>
      <td>0.151233</td>
      <td>-0.594247</td>
      <td>-0.666882</td>
      <td>-0.478262</td>
      <td>-0.676907</td>
      <td>-0.817864</td>
      <td>-0.660876</td>
      <td>0.140619</td>
      <td>-0.768692</td>
      <td>-0.746999</td>
      <td>-0.441499</td>
      <td>-0.126745</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>9484.342696</td>
      <td>7768.555088</td>
      <td>AAACAGAGCGACTCCT-1</td>
    </tr>
  </tbody>
</table>
</div>




```python
print(va1_nn_pre['orig.ident'].value_counts())
print(vb1_nn_pre['orig.ident'].value_counts())
print(va2_nn_pre['orig.ident'].value_counts())
```

    orig.ident
    VA1    2806
    Name: count, dtype: int64
    orig.ident
    VB1    2589
    Name: count, dtype: int64
    orig.ident
    VA2    4316
    Name: count, dtype: int64



```python

```


```python
# exctract x and y coords as a numpy array
va1_coordinates = va1_nn_pre[['x','y']].to_numpy()
print(va1_coordinates)
```

    [[30730 26443]
     [12698  8743]
     [29633 20870]
     ...
     [ 9867 21146]
     [ 9157 23627]
     [14550 24228]]



```python
# exctract x and y coords as a numpy array
vb1_coordinates = vb1_nn_pre[['x','y']].to_numpy()
print(vb1_coordinates)
```

    [[31283 26488]
     [13243  8788]
     [10072 25947]
     ...
     [10413 21194]
     [ 9704 23676]
     [15098 24276]]



```python
# exctract x and y coords as a numpy array
va2_coordinates = va2_nn_pre[['x','y']].to_numpy()
print(va2_coordinates)
```

    [[ 7764.08568781 20534.80431007]
     [24672.63400836 23846.30507581]
     [ 9484.34269641  7768.55508804]
     ...
     [23075.13754142 18865.97437049]
     [20150.62015553 28784.30455988]
     [18272.82729726  5344.70864869]]



```python
# find nearest neighbors

va1_tree = spatial.KDTree(va1_coordinates)
va1_nearest_dist, va1_nearest_ind = va1_tree.query(va1_coordinates, k=7)  # k=7, six nearest neighbors where k1 = identity

#print(nearest_dist[:, 1])    # drop id; assumes sorted -> see args!
#print(nearest_ind[:, 1])     # drop id 
print(va1_nearest_dist[0:10])    # drop id; assumes sorted -> see args!
print(va1_nearest_ind[0:10])     # drop id 
```

    [[  0.         413.00121065 413.00484259 413.53718092 716.30300851
      716.57658349 718.00278551]
     [  0.         413.00121065 413.00484259 414.03864554 414.27647773
      414.77222665 414.90360326]
     [  0.         413.00121065 414.03864554 414.27647773 414.90360326
      715.71293128 716.80192522]
     [  0.         413.00121065 413.00484259 414.03864554 414.27647773
      414.27647773 414.90360326]
     [  0.         413.00121065 413.00484259 414.03864554 414.27647773
      414.27647773 414.90360326]
     [  0.         413.00121065 413.00484259 414.03864554 414.27647773
      414.27647773 414.90360326]
     [  0.         413.00121065 413.00484259 414.03864554 414.27647773
      414.27647773 414.90360326]
     [  0.         413.00121065 413.00484259 413.53718092 414.77222665
      414.90360326 415.14575754]
     [  0.         413.00121065 414.00483089 414.03864554 414.27647773
      414.40318532 414.77222665]
     [  0.         413.00484259 413.53718092 414.00120773 414.27647773
      414.77222665 414.90360326]]
    [[   0  854 1863 2089  517 2692 2472]
     [   1  284  771  352 1783 2738 1604]
     [   2 2207 1663 1116 1210  920 1124]
     [   3 1747  425 1286 1385 2332 2567]
     [   4  964 1350 2187  320  885 1441]
     [   5 1045  976 1053 2706 1504  190]
     [   6 2428 1102 2758 2766 2012 1580]
     [   7 2124 2026 2623  665 1780  335]
     [   8 1128  375  567 1573 2777 1101]
     [   9  865   77 1703 2772  773 1671]]



```python
print(va1_nearest_dist[0:10, 1:])    # drop id; assumes sorted -> see args!
print(va1_nearest_ind[0:10, 1:])     # drop id 
```

    [[413.00121065 413.00484259 413.53718092 716.30300851 716.57658349
      718.00278551]
     [413.00121065 413.00484259 414.03864554 414.27647773 414.77222665
      414.90360326]
     [413.00121065 414.03864554 414.27647773 414.90360326 715.71293128
      716.80192522]
     [413.00121065 413.00484259 414.03864554 414.27647773 414.27647773
      414.90360326]
     [413.00121065 413.00484259 414.03864554 414.27647773 414.27647773
      414.90360326]
     [413.00121065 413.00484259 414.03864554 414.27647773 414.27647773
      414.90360326]
     [413.00121065 413.00484259 414.03864554 414.27647773 414.27647773
      414.90360326]
     [413.00121065 413.00484259 413.53718092 414.77222665 414.90360326
      415.14575754]
     [413.00121065 414.00483089 414.03864554 414.27647773 414.40318532
      414.77222665]
     [413.00484259 413.53718092 414.00120773 414.27647773 414.77222665
      414.90360326]]
    [[ 854 1863 2089  517 2692 2472]
     [ 284  771  352 1783 2738 1604]
     [2207 1663 1116 1210  920 1124]
     [1747  425 1286 1385 2332 2567]
     [ 964 1350 2187  320  885 1441]
     [1045  976 1053 2706 1504  190]
     [2428 1102 2758 2766 2012 1580]
     [2124 2026 2623  665 1780  335]
     [1128  375  567 1573 2777 1101]
     [ 865   77 1703 2772  773 1671]]



```python
print(va1_nn_pre.shape)
print(va1_nearest_dist.shape)
print(va1_nearest_ind.shape)
```

    (2806, 60)
    (2806, 7)
    (2806, 7)



```python
# find nearest neighbors

vb1_tree = spatial.KDTree(vb1_coordinates)
vb1_nearest_dist, vb1_nearest_ind = vb1_tree.query(vb1_coordinates, k=7)  # k=7, six nearest neighbors where k1 = identity

#print(nearest_dist[:, 1])    # drop id; assumes sorted -> see args!
#print(nearest_ind[:, 1])     # drop id 
print(vb1_nearest_dist[0:10])    # drop id; assumes sorted -> see args!
print(vb1_nearest_ind[0:10])     # drop id 
```

    [[  0.         413.00121065 413.00484259 414.27647773 414.90360326
      716.30300851 717.08088805]
     [  0.         413.00121065 414.00483089 414.03864554 414.27647773
      414.40318532 414.77222665]
     [  0.         413.00121065 413.00484259 414.03864554 414.77222665
      414.90360326 415.14575754]
     [  0.         413.00484259 414.03864554 414.27647773 415.14575754
      716.21784396 826.00544792]
     [  0.         413.00121065 413.00484259 414.27647773 414.77222665
      414.90360326 414.90360326]
     [  0.         413.00121065 413.00484259 414.27647773 414.77222665
      414.90360326 414.90360326]
     [  0.         414.00483089 414.03864554 415.64046964 716.30300851
      717.08088805 719.00278164]
     [  0.         413.00121065 413.00484259 414.27647773 414.27647773
      414.90360326 414.90360326]
     [  0.         413.00121065 413.00484259 414.03864554 414.27647773
      414.90360326 415.14575754]
     [  0.         413.00121065 413.00484259 414.27647773 414.27647773
      414.90360326 414.90360326]]
    [[   0  804 1719  654 1936  492 2478]
     [   1  724  270 1481 2524  336 1646]
     [   2  410 1610 2358 2145 1193 1279]
     [   3 1455 2110  727 1576  157 1115]
     [   4  899 1253  833  306 2021 1329]
     [   5  912  976 2494 1394  189  986]
     [   6 2115 1130  535 1689  901 1614]
     [   7 1035 2229 2554 1870 1463 2546]
     [   8  360 1050 2564 1034  537 1456]
     [   9 1566  815  728 2561 1539   77]]



```python
print(vb1_nearest_dist[0:10, 1:])    # drop id; assumes sorted -> see args!
print(vb1_nearest_ind[0:10, 1:])     # drop id 
```

    [[413.00121065 413.00484259 414.27647773 414.90360326 716.30300851
      717.08088805]
     [413.00121065 414.00483089 414.03864554 414.27647773 414.40318532
      414.77222665]
     [413.00121065 413.00484259 414.03864554 414.77222665 414.90360326
      415.14575754]
     [413.00484259 414.03864554 414.27647773 415.14575754 716.21784396
      826.00544792]
     [413.00121065 413.00484259 414.27647773 414.77222665 414.90360326
      414.90360326]
     [413.00121065 413.00484259 414.27647773 414.77222665 414.90360326
      414.90360326]
     [414.00483089 414.03864554 415.64046964 716.30300851 717.08088805
      719.00278164]
     [413.00121065 413.00484259 414.27647773 414.27647773 414.90360326
      414.90360326]
     [413.00121065 413.00484259 414.03864554 414.27647773 414.90360326
      415.14575754]
     [413.00121065 413.00484259 414.27647773 414.27647773 414.90360326
      414.90360326]]
    [[ 804 1719  654 1936  492 2478]
     [ 724  270 1481 2524  336 1646]
     [ 410 1610 2358 2145 1193 1279]
     [1455 2110  727 1576  157 1115]
     [ 899 1253  833  306 2021 1329]
     [ 912  976 2494 1394  189  986]
     [2115 1130  535 1689  901 1614]
     [1035 2229 2554 1870 1463 2546]
     [ 360 1050 2564 1034  537 1456]
     [1566  815  728 2561 1539   77]]



```python
print(vb1_nn_pre.shape)
print(vb1_nearest_dist.shape)
print(vb1_nearest_ind.shape)
```

    (2589, 60)
    (2589, 7)
    (2589, 7)



```python
# find nearest neighbors

va2_tree = spatial.KDTree(va2_coordinates)
va2_nearest_dist, va2_nearest_ind = va2_tree.query(va2_coordinates, k=7)  # k=7, six nearest neighbors where k1 = identity

#print(nearest_dist[:, 1])    # drop id; assumes sorted -> see args!
#print(nearest_ind[:, 1])     # drop id 
print(va2_nearest_dist[0:10])    # drop id; assumes sorted -> see args!
print(va2_nearest_ind[0:10])     # drop id 
```

    [[  0.         407.9771153  407.9771153  409.38221671 409.38221671
      409.38221671 409.38221671]
     [  0.         407.9771153  407.9771153  409.38221671 409.38221671
      409.38221671 409.38221671]
     [  0.         407.9771153  407.9771153  409.38221671 409.38221671
      409.38221671 409.38221671]
     [  0.         407.9771153  407.9771153  409.38221671 409.38221671
      409.38221671 409.38221671]
     [  0.         407.9771153  407.9771153  409.38221671 409.38221671
      409.38221671 409.38221671]
     [  0.         407.9771153  407.9771153  409.38221671 409.38221671
      409.38221671 409.38221671]
     [  0.         407.9771153  409.38221671 409.38221671 707.44925795
      707.44925795 709.88018061]
     [  0.         407.9771153  409.38221671 409.38221671 409.38221671
      409.38221671 707.44925795]
     [  0.         407.9771153  407.9771153  409.38221671 409.38221671
      409.38221671 409.38221671]
     [  0.         407.9771153  407.9771153  409.38221671 409.38221671
      409.38221671 409.38221671]]
    [[   0 1217  467 2483 4204  565 2751]
     [   1  681 2687 3589 2004 3930 2163]
     [   2 2444 4259 3524 1222 2630 3521]
     [   3 1502 2113 3380 1392 2248  523]
     [   4 1630 1520 2349 1647  312 4157]
     [   5 4314 1720 2459 3933  890 4167]
     [   6 3534  893 1897 2826 1504 2696]
     [   7 3507 1575  845 1531 2538 2155]
     [   8 4232 2837 3213   15 3862 3487]
     [   9  608 1754  895 2445 1723 4269]]



```python
print(va2_nearest_dist[0:10, 1:])    # drop id; assumes sorted -> see args!
print(va2_nearest_ind[0:10, 1:])     # drop id 
```

    [[407.9771153  407.9771153  409.38221671 409.38221671 409.38221671
      409.38221671]
     [407.9771153  407.9771153  409.38221671 409.38221671 409.38221671
      409.38221671]
     [407.9771153  407.9771153  409.38221671 409.38221671 409.38221671
      409.38221671]
     [407.9771153  407.9771153  409.38221671 409.38221671 409.38221671
      409.38221671]
     [407.9771153  407.9771153  409.38221671 409.38221671 409.38221671
      409.38221671]
     [407.9771153  407.9771153  409.38221671 409.38221671 409.38221671
      409.38221671]
     [407.9771153  409.38221671 409.38221671 707.44925795 707.44925795
      709.88018061]
     [407.9771153  409.38221671 409.38221671 409.38221671 409.38221671
      707.44925795]
     [407.9771153  407.9771153  409.38221671 409.38221671 409.38221671
      409.38221671]
     [407.9771153  407.9771153  409.38221671 409.38221671 409.38221671
      409.38221671]]
    [[1217  467 2483 4204  565 2751]
     [ 681 2687 3589 2004 3930 2163]
     [2444 4259 3524 1222 2630 3521]
     [1502 2113 3380 1392 2248  523]
     [1630 1520 2349 1647  312 4157]
     [4314 1720 2459 3933  890 4167]
     [3534  893 1897 2826 1504 2696]
     [3507 1575  845 1531 2538 2155]
     [4232 2837 3213   15 3862 3487]
     [ 608 1754  895 2445 1723 4269]]



```python
print(va2_nn_pre.shape)
print(va2_nearest_dist.shape)
print(va2_nearest_ind.shape)
```

    (4316, 60)
    (4316, 7)
    (4316, 7)


# - - - - - - - - - -
# VA1 VA1 VA1 VA1 VA1 VA1
# - - - - - - - - - -

# Some neighborhoods include spots that are farther than the immediate ring... ( distances > 8 ) these must be spots near the edge. I should remove these from the analysis. 


```python
# make into dfs
va1_id_df = pd.DataFrame(va1_nearest_ind, columns=['id','nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id'])
va1_dist_df = pd.DataFrame(va1_nearest_dist[:, 1:], columns=['nn1_dist','nn2_dist','nn3_dist','nn4_dist','nn5_dist','nn6_dist'])
print(va1_id_df.shape)
print(va1_dist_df.shape)
va1_id_df = va1_id_df.merge(va1_dist_df, left_index=True, right_index=True)
print(va1_id_df.shape)
va1_id_df.head(3)
```

    (2806, 7)
    (2806, 6)
    (2806, 13)





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>nn1_dist</th>
      <th>nn2_dist</th>
      <th>nn3_dist</th>
      <th>nn4_dist</th>
      <th>nn5_dist</th>
      <th>nn6_dist</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>854</td>
      <td>1863</td>
      <td>2089</td>
      <td>517</td>
      <td>2692</td>
      <td>2472</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>413.537181</td>
      <td>716.303009</td>
      <td>716.576583</td>
      <td>718.002786</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>284</td>
      <td>771</td>
      <td>352</td>
      <td>1783</td>
      <td>2738</td>
      <td>1604</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>414.038646</td>
      <td>414.276478</td>
      <td>414.772227</td>
      <td>414.903603</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>2207</td>
      <td>1663</td>
      <td>1116</td>
      <td>1210</td>
      <td>920</td>
      <td>1124</td>
      <td>413.001211</td>
      <td>414.038646</td>
      <td>414.276478</td>
      <td>414.903603</td>
      <td>715.712931</td>
      <td>716.801925</td>
    </tr>
  </tbody>
</table>
</div>




```python
# count the number of nearest neighbors directly adjacent
va1_id_df['num_adjacent_nns'] = va1_id_df[['nn1_dist','nn2_dist','nn3_dist','nn4_dist','nn5_dist','nn6_dist']].apply(lambda x: x[x < 420].count(), axis=1)
va1_id_df.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>nn1_dist</th>
      <th>nn2_dist</th>
      <th>nn3_dist</th>
      <th>nn4_dist</th>
      <th>nn5_dist</th>
      <th>nn6_dist</th>
      <th>num_adjacent_nns</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>854</td>
      <td>1863</td>
      <td>2089</td>
      <td>517</td>
      <td>2692</td>
      <td>2472</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>413.537181</td>
      <td>716.303009</td>
      <td>716.576583</td>
      <td>718.002786</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>284</td>
      <td>771</td>
      <td>352</td>
      <td>1783</td>
      <td>2738</td>
      <td>1604</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>414.038646</td>
      <td>414.276478</td>
      <td>414.772227</td>
      <td>414.903603</td>
      <td>6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>2207</td>
      <td>1663</td>
      <td>1116</td>
      <td>1210</td>
      <td>920</td>
      <td>1124</td>
      <td>413.001211</td>
      <td>414.038646</td>
      <td>414.276478</td>
      <td>414.903603</td>
      <td>715.712931</td>
      <td>716.801925</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>




```python
# function to add tissue slice identifier to an id column
def add_tissue_slice_id(df, tissue_slice_id, col_list):
    for col in col_list:
        df[col] = tissue_slice_id + '_' + df[col].astype(str)
    return(df)
```


```python
# Add tissue slice identifier VA to all ids
va1_temp = va1_id_df.copy()

col_list = ['id','nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id']
va1_temp = add_tissue_slice_id(va1_temp, 'VA1', col_list)
display(va1_temp.head(3))

va1_id_df = va1_temp.copy()
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>nn1_dist</th>
      <th>nn2_dist</th>
      <th>nn3_dist</th>
      <th>nn4_dist</th>
      <th>nn5_dist</th>
      <th>nn6_dist</th>
      <th>num_adjacent_nns</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VA1_0</td>
      <td>VA1_854</td>
      <td>VA1_1863</td>
      <td>VA1_2089</td>
      <td>VA1_517</td>
      <td>VA1_2692</td>
      <td>VA1_2472</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>413.537181</td>
      <td>716.303009</td>
      <td>716.576583</td>
      <td>718.002786</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VA1_1</td>
      <td>VA1_284</td>
      <td>VA1_771</td>
      <td>VA1_352</td>
      <td>VA1_1783</td>
      <td>VA1_2738</td>
      <td>VA1_1604</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>414.038646</td>
      <td>414.276478</td>
      <td>414.772227</td>
      <td>414.903603</td>
      <td>6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VA1_2</td>
      <td>VA1_2207</td>
      <td>VA1_1663</td>
      <td>VA1_1116</td>
      <td>VA1_1210</td>
      <td>VA1_920</td>
      <td>VA1_1124</td>
      <td>413.001211</td>
      <td>414.038646</td>
      <td>414.276478</td>
      <td>414.903603</td>
      <td>715.712931</td>
      <td>716.801925</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>



```python
# drop neighborhoods with spots farther than 420 in nn6_dist
va1_id_df_sub1 = va1_id_df[va1_id_df['nn6_dist']<420]
print(va1_id_df_sub1.shape)

# allow one farther (so 5 neighbors immediately next to center spot)
va1_id_df_sub2 = va1_id_df[va1_id_df['nn5_dist']<420]
print(va1_id_df_sub2.shape)

# allow two farther (so 4 neighbors immediately next to center spot)
va1_id_df_sub3 = va1_id_df[va1_id_df['nn4_dist']<420]
print(va1_id_df_sub3.shape)
```

    (2392, 14)
    (2519, 14)
    (2678, 14)



```python
# let's require all 6 neighbors
va1_nn = va1_nn_pre.merge(va1_id_df[['id','nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id','num_adjacent_nns']], left_index=True, right_index=True)
print(va1_nn.shape)
va1_nn.head(3)
```

    (2806, 68)





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>num_adjacent_nns</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VA1_AAACAACGAATAGTTC-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>3625</td>
      <td>2571</td>
      <td>round1</td>
      <td>3649</td>
      <td>2571</td>
      <td>10</td>
      <td>10</td>
      <td>VA1</td>
      <td>2.648104</td>
      <td>vasc_high</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>2</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_high</td>
      <td>Respiratory epithelium</td>
      <td>Respiratory epithelium</td>
      <td>Respiratory Epithelium (Normal Lung)</td>
      <td>10_Respiratory_Epithelium_Lung</td>
      <td>7_10_Normal_Lung</td>
      <td>-0.056807</td>
      <td>-0.259847</td>
      <td>0.296748</td>
      <td>-0.038431</td>
      <td>0.513527</td>
      <td>-0.007136</td>
      <td>0.412196</td>
      <td>0.404657</td>
      <td>-0.158465</td>
      <td>0.690663</td>
      <td>-0.046498</td>
      <td>0.110800</td>
      <td>0.102126</td>
      <td>0.654483</td>
      <td>-0.176035</td>
      <td>-0.321725</td>
      <td>-0.469550</td>
      <td>0.189499</td>
      <td>0.139214</td>
      <td>-0.165127</td>
      <td>0.015559</td>
      <td>-0.755217</td>
      <td>0.633032</td>
      <td>-0.326533</td>
      <td>-0.046498</td>
      <td>0.031719</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_high_MN4_high</td>
      <td>vasc_high_MN4_upper80</td>
      <td>30730</td>
      <td>26443</td>
      <td>AAACAACGAATAGTTC-1</td>
      <td>VA1_0</td>
      <td>VA1_854</td>
      <td>VA1_1863</td>
      <td>VA1_2089</td>
      <td>VA1_517</td>
      <td>VA1_2692</td>
      <td>VA1_2472</td>
      <td>3</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VA1_AAACAAGTATCTCCCA-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>5804</td>
      <td>3329</td>
      <td>round1</td>
      <td>4956</td>
      <td>3308</td>
      <td>6</td>
      <td>6</td>
      <td>VA1</td>
      <td>1.001716</td>
      <td>vasc_low</td>
      <td>13</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>2</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_low</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.432040</td>
      <td>-0.363696</td>
      <td>-0.665610</td>
      <td>-0.380249</td>
      <td>-0.272964</td>
      <td>-0.616570</td>
      <td>-0.133588</td>
      <td>-0.277191</td>
      <td>-0.532077</td>
      <td>-0.284190</td>
      <td>-0.442197</td>
      <td>-0.206272</td>
      <td>-0.292946</td>
      <td>-0.220110</td>
      <td>0.816411</td>
      <td>-0.422218</td>
      <td>-0.542976</td>
      <td>-0.323644</td>
      <td>-0.431187</td>
      <td>-0.713442</td>
      <td>-0.122593</td>
      <td>0.079799</td>
      <td>-0.348075</td>
      <td>-0.420942</td>
      <td>-0.442197</td>
      <td>-0.484106</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>12698</td>
      <td>8743</td>
      <td>AAACAAGTATCTCCCA-1</td>
      <td>VA1_1</td>
      <td>VA1_284</td>
      <td>VA1_771</td>
      <td>VA1_352</td>
      <td>VA1_1783</td>
      <td>VA1_2738</td>
      <td>VA1_1604</td>
      <td>6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VA1_AAACAATCTACTAGCA-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>4980</td>
      <td>3146</td>
      <td>round1</td>
      <td>4711</td>
      <td>3146</td>
      <td>6</td>
      <td>6</td>
      <td>VA1</td>
      <td>1.101286</td>
      <td>vasc_low</td>
      <td>19</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_low</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.312254</td>
      <td>-0.142246</td>
      <td>0.215974</td>
      <td>-0.522873</td>
      <td>-0.254499</td>
      <td>0.344170</td>
      <td>-0.096425</td>
      <td>0.015784</td>
      <td>-0.519466</td>
      <td>-0.273987</td>
      <td>-0.296916</td>
      <td>-0.112254</td>
      <td>-0.139708</td>
      <td>0.574189</td>
      <td>-0.242849</td>
      <td>-0.389510</td>
      <td>-0.529344</td>
      <td>-0.039327</td>
      <td>0.055921</td>
      <td>-0.707380</td>
      <td>-0.091447</td>
      <td>0.049511</td>
      <td>-0.304126</td>
      <td>0.924392</td>
      <td>-0.296916</td>
      <td>-0.466480</td>
      <td>High</td>
      <td>Low</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>29633</td>
      <td>20870</td>
      <td>AAACAATCTACTAGCA-1</td>
      <td>VA1_2</td>
      <td>VA1_2207</td>
      <td>VA1_1663</td>
      <td>VA1_1116</td>
      <td>VA1_1210</td>
      <td>VA1_920</td>
      <td>VA1_1124</td>
      <td>4</td>
    </tr>
  </tbody>
</table>
</div>




```python
# turn nn_ids into a list in a new column
va1_nn_list_temp = ['nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id']
va1_nn['nn_list'] = va1_nn[va1_nn_list_temp].values.tolist()
va1_nn.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>num_adjacent_nns</th>
      <th>nn_list</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VA1_AAACAACGAATAGTTC-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>3625</td>
      <td>2571</td>
      <td>round1</td>
      <td>3649</td>
      <td>2571</td>
      <td>10</td>
      <td>10</td>
      <td>VA1</td>
      <td>2.648104</td>
      <td>vasc_high</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>2</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_high</td>
      <td>Respiratory epithelium</td>
      <td>Respiratory epithelium</td>
      <td>Respiratory Epithelium (Normal Lung)</td>
      <td>10_Respiratory_Epithelium_Lung</td>
      <td>7_10_Normal_Lung</td>
      <td>-0.056807</td>
      <td>-0.259847</td>
      <td>0.296748</td>
      <td>-0.038431</td>
      <td>0.513527</td>
      <td>-0.007136</td>
      <td>0.412196</td>
      <td>0.404657</td>
      <td>-0.158465</td>
      <td>0.690663</td>
      <td>-0.046498</td>
      <td>0.110800</td>
      <td>0.102126</td>
      <td>0.654483</td>
      <td>-0.176035</td>
      <td>-0.321725</td>
      <td>-0.469550</td>
      <td>0.189499</td>
      <td>0.139214</td>
      <td>-0.165127</td>
      <td>0.015559</td>
      <td>-0.755217</td>
      <td>0.633032</td>
      <td>-0.326533</td>
      <td>-0.046498</td>
      <td>0.031719</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_high_MN4_high</td>
      <td>vasc_high_MN4_upper80</td>
      <td>30730</td>
      <td>26443</td>
      <td>AAACAACGAATAGTTC-1</td>
      <td>VA1_0</td>
      <td>VA1_854</td>
      <td>VA1_1863</td>
      <td>VA1_2089</td>
      <td>VA1_517</td>
      <td>VA1_2692</td>
      <td>VA1_2472</td>
      <td>3</td>
      <td>[VA1_854, VA1_1863, VA1_2089, VA1_517, VA1_269...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VA1_AAACAAGTATCTCCCA-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>5804</td>
      <td>3329</td>
      <td>round1</td>
      <td>4956</td>
      <td>3308</td>
      <td>6</td>
      <td>6</td>
      <td>VA1</td>
      <td>1.001716</td>
      <td>vasc_low</td>
      <td>13</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>2</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_low</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.432040</td>
      <td>-0.363696</td>
      <td>-0.665610</td>
      <td>-0.380249</td>
      <td>-0.272964</td>
      <td>-0.616570</td>
      <td>-0.133588</td>
      <td>-0.277191</td>
      <td>-0.532077</td>
      <td>-0.284190</td>
      <td>-0.442197</td>
      <td>-0.206272</td>
      <td>-0.292946</td>
      <td>-0.220110</td>
      <td>0.816411</td>
      <td>-0.422218</td>
      <td>-0.542976</td>
      <td>-0.323644</td>
      <td>-0.431187</td>
      <td>-0.713442</td>
      <td>-0.122593</td>
      <td>0.079799</td>
      <td>-0.348075</td>
      <td>-0.420942</td>
      <td>-0.442197</td>
      <td>-0.484106</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>12698</td>
      <td>8743</td>
      <td>AAACAAGTATCTCCCA-1</td>
      <td>VA1_1</td>
      <td>VA1_284</td>
      <td>VA1_771</td>
      <td>VA1_352</td>
      <td>VA1_1783</td>
      <td>VA1_2738</td>
      <td>VA1_1604</td>
      <td>6</td>
      <td>[VA1_284, VA1_771, VA1_352, VA1_1783, VA1_2738...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VA1_AAACAATCTACTAGCA-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>4980</td>
      <td>3146</td>
      <td>round1</td>
      <td>4711</td>
      <td>3146</td>
      <td>6</td>
      <td>6</td>
      <td>VA1</td>
      <td>1.101286</td>
      <td>vasc_low</td>
      <td>19</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_low</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.312254</td>
      <td>-0.142246</td>
      <td>0.215974</td>
      <td>-0.522873</td>
      <td>-0.254499</td>
      <td>0.344170</td>
      <td>-0.096425</td>
      <td>0.015784</td>
      <td>-0.519466</td>
      <td>-0.273987</td>
      <td>-0.296916</td>
      <td>-0.112254</td>
      <td>-0.139708</td>
      <td>0.574189</td>
      <td>-0.242849</td>
      <td>-0.389510</td>
      <td>-0.529344</td>
      <td>-0.039327</td>
      <td>0.055921</td>
      <td>-0.707380</td>
      <td>-0.091447</td>
      <td>0.049511</td>
      <td>-0.304126</td>
      <td>0.924392</td>
      <td>-0.296916</td>
      <td>-0.466480</td>
      <td>High</td>
      <td>Low</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>29633</td>
      <td>20870</td>
      <td>AAACAATCTACTAGCA-1</td>
      <td>VA1_2</td>
      <td>VA1_2207</td>
      <td>VA1_1663</td>
      <td>VA1_1116</td>
      <td>VA1_1210</td>
      <td>VA1_920</td>
      <td>VA1_1124</td>
      <td>4</td>
      <td>[VA1_2207, VA1_1663, VA1_1116, VA1_1210, VA1_9...</td>
    </tr>
  </tbody>
</table>
</div>



# - - - - - - - - - -
# VB1 VB1 VB1 VB1 VB1 VB1
# - - - - - - - - - -

# Some neighborhoods include spots that are farther than the immediate ring... ( distances > 8 ) these must be spots near the edge. I should remove these from the analysis. 


```python
# make into dfs
vb1_id_df = pd.DataFrame(vb1_nearest_ind, columns=['id','nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id'])
vb1_dist_df = pd.DataFrame(vb1_nearest_dist[:, 1:], columns=['nn1_dist','nn2_dist','nn3_dist','nn4_dist','nn5_dist','nn6_dist'])
print(vb1_id_df.shape)
print(vb1_dist_df.shape)
vb1_id_df = vb1_id_df.merge(vb1_dist_df, left_index=True, right_index=True)
print(vb1_id_df.shape)
vb1_id_df.head(3)
```

    (2589, 7)
    (2589, 6)
    (2589, 13)





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>nn1_dist</th>
      <th>nn2_dist</th>
      <th>nn3_dist</th>
      <th>nn4_dist</th>
      <th>nn5_dist</th>
      <th>nn6_dist</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>804</td>
      <td>1719</td>
      <td>654</td>
      <td>1936</td>
      <td>492</td>
      <td>2478</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>414.276478</td>
      <td>414.903603</td>
      <td>716.303009</td>
      <td>717.080888</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>724</td>
      <td>270</td>
      <td>1481</td>
      <td>2524</td>
      <td>336</td>
      <td>1646</td>
      <td>413.001211</td>
      <td>414.004831</td>
      <td>414.038646</td>
      <td>414.276478</td>
      <td>414.403185</td>
      <td>414.772227</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>410</td>
      <td>1610</td>
      <td>2358</td>
      <td>2145</td>
      <td>1193</td>
      <td>1279</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>414.038646</td>
      <td>414.772227</td>
      <td>414.903603</td>
      <td>415.145758</td>
    </tr>
  </tbody>
</table>
</div>




```python
# count the number of nearest neighbors directly adjacent
vb1_id_df['num_adjacent_nns'] = vb1_id_df[['nn1_dist','nn2_dist','nn3_dist','nn4_dist','nn5_dist','nn6_dist']].apply(lambda x: x[x < 420].count(), axis=1)
vb1_id_df.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>nn1_dist</th>
      <th>nn2_dist</th>
      <th>nn3_dist</th>
      <th>nn4_dist</th>
      <th>nn5_dist</th>
      <th>nn6_dist</th>
      <th>num_adjacent_nns</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>804</td>
      <td>1719</td>
      <td>654</td>
      <td>1936</td>
      <td>492</td>
      <td>2478</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>414.276478</td>
      <td>414.903603</td>
      <td>716.303009</td>
      <td>717.080888</td>
      <td>4</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>724</td>
      <td>270</td>
      <td>1481</td>
      <td>2524</td>
      <td>336</td>
      <td>1646</td>
      <td>413.001211</td>
      <td>414.004831</td>
      <td>414.038646</td>
      <td>414.276478</td>
      <td>414.403185</td>
      <td>414.772227</td>
      <td>6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>410</td>
      <td>1610</td>
      <td>2358</td>
      <td>2145</td>
      <td>1193</td>
      <td>1279</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>414.038646</td>
      <td>414.772227</td>
      <td>414.903603</td>
      <td>415.145758</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
</div>




```python
# function to add tissue slice identifier to an id column
def add_tissue_slice_id(df, tissue_slice_id, col_list):
    for col in col_list:
        df[col] = tissue_slice_id + '_' + df[col].astype(str)
    return(df)
```


```python
# Add tissue slice identifier VB to all ids
vb1_temp = vb1_id_df.copy()

col_list = ['id','nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id']
vb1_temp = add_tissue_slice_id(vb1_temp, 'VB1', col_list)
display(vb1_temp.head(3))

vb1_id_df = vb1_temp.copy()
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>nn1_dist</th>
      <th>nn2_dist</th>
      <th>nn3_dist</th>
      <th>nn4_dist</th>
      <th>nn5_dist</th>
      <th>nn6_dist</th>
      <th>num_adjacent_nns</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VB1_0</td>
      <td>VB1_804</td>
      <td>VB1_1719</td>
      <td>VB1_654</td>
      <td>VB1_1936</td>
      <td>VB1_492</td>
      <td>VB1_2478</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>414.276478</td>
      <td>414.903603</td>
      <td>716.303009</td>
      <td>717.080888</td>
      <td>4</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VB1_1</td>
      <td>VB1_724</td>
      <td>VB1_270</td>
      <td>VB1_1481</td>
      <td>VB1_2524</td>
      <td>VB1_336</td>
      <td>VB1_1646</td>
      <td>413.001211</td>
      <td>414.004831</td>
      <td>414.038646</td>
      <td>414.276478</td>
      <td>414.403185</td>
      <td>414.772227</td>
      <td>6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VB1_2</td>
      <td>VB1_410</td>
      <td>VB1_1610</td>
      <td>VB1_2358</td>
      <td>VB1_2145</td>
      <td>VB1_1193</td>
      <td>VB1_1279</td>
      <td>413.001211</td>
      <td>413.004843</td>
      <td>414.038646</td>
      <td>414.772227</td>
      <td>414.903603</td>
      <td>415.145758</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
</div>



```python
# drop neighborhoods with spots farther than 420 in nn6_dist
vb1_id_df_sub1 = vb1_id_df[vb1_id_df['nn6_dist']<420]
print(vb1_id_df_sub1.shape)

# allow one farther (so 5 neighbors immediately next to center spot)
vb1_id_df_sub2 = vb1_id_df[vb1_id_df['nn5_dist']<420]
print(vb1_id_df_sub2.shape)

# allow two farther (so 4 neighbors immediately next to center spot)
vb1_id_df_sub3 = vb1_id_df[vb1_id_df['nn4_dist']<420]
print(vb1_id_df_sub3.shape)
```

    (2109, 14)
    (2278, 14)
    (2427, 14)



```python
# let's require all 6 neighbors
vb1_nn = vb1_nn_pre.merge(vb1_id_df[['id','nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id','num_adjacent_nns']], left_index=True, right_index=True)
print(vb1_nn.shape)
vb1_nn.head(3)
```

    (2589, 68)





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>num_adjacent_nns</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VB1_AAACAACGAATAGTTC-1</td>
      <td>VB1</td>
      <td>VB1</td>
      <td>7174</td>
      <td>3783</td>
      <td>round1</td>
      <td>5068</td>
      <td>3601</td>
      <td>6</td>
      <td>6</td>
      <td>VB1</td>
      <td>0.000000</td>
      <td>vasc_low</td>
      <td>20</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_neg</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.409190</td>
      <td>-0.293434</td>
      <td>0.015793</td>
      <td>-0.623296</td>
      <td>-0.294338</td>
      <td>-0.538970</td>
      <td>0.160997</td>
      <td>-0.440288</td>
      <td>-0.049946</td>
      <td>-0.219445</td>
      <td>-0.184969</td>
      <td>0.010123</td>
      <td>0.001746</td>
      <td>0.523236</td>
      <td>-0.383538</td>
      <td>0.263978</td>
      <td>-0.317542</td>
      <td>-0.189507</td>
      <td>0.171521</td>
      <td>-0.565283</td>
      <td>-0.448105</td>
      <td>-0.309486</td>
      <td>-0.363641</td>
      <td>-0.373975</td>
      <td>-0.184969</td>
      <td>0.090032</td>
      <td>Low</td>
      <td>Low</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>31283</td>
      <td>26488</td>
      <td>AAACAACGAATAGTTC-1</td>
      <td>VB1_0</td>
      <td>VB1_804</td>
      <td>VB1_1719</td>
      <td>VB1_654</td>
      <td>VB1_1936</td>
      <td>VB1_492</td>
      <td>VB1_2478</td>
      <td>4</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VB1_AAACAAGTATCTCCCA-1</td>
      <td>VB1</td>
      <td>VB1</td>
      <td>3817</td>
      <td>2574</td>
      <td>round1</td>
      <td>3824</td>
      <td>2572</td>
      <td>1</td>
      <td>1</td>
      <td>VB1</td>
      <td>3.117370</td>
      <td>vasc_high</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>4</td>
      <td>MHCI_high</td>
      <td>MHCI_high</td>
      <td>Tumor_low_MHCI_high</td>
      <td>Tumor_low</td>
      <td>vasc_high</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC</td>
      <td>1_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>0.124843</td>
      <td>0.508282</td>
      <td>-0.370904</td>
      <td>-0.208191</td>
      <td>-0.135008</td>
      <td>0.254889</td>
      <td>0.354699</td>
      <td>0.618041</td>
      <td>0.174323</td>
      <td>0.644896</td>
      <td>0.223238</td>
      <td>0.065979</td>
      <td>0.010680</td>
      <td>0.645669</td>
      <td>0.955596</td>
      <td>-0.325863</td>
      <td>0.259486</td>
      <td>-0.095793</td>
      <td>-0.280748</td>
      <td>-0.611245</td>
      <td>0.417938</td>
      <td>-0.344595</td>
      <td>-0.302494</td>
      <td>-0.264353</td>
      <td>0.223238</td>
      <td>0.694163</td>
      <td>Low</td>
      <td>Low</td>
      <td>vasc_high_MN4_low</td>
      <td>vasc_low</td>
      <td>13243</td>
      <td>8788</td>
      <td>AAACAAGTATCTCCCA-1</td>
      <td>VB1_1</td>
      <td>VB1_724</td>
      <td>VB1_270</td>
      <td>VB1_1481</td>
      <td>VB1_2524</td>
      <td>VB1_336</td>
      <td>VB1_1646</td>
      <td>6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VB1_AAACACCAATAACTGC-1</td>
      <td>VB1</td>
      <td>VB1</td>
      <td>6137</td>
      <td>3622</td>
      <td>round1</td>
      <td>4990</td>
      <td>3585</td>
      <td>0</td>
      <td>0</td>
      <td>VB1</td>
      <td>0.966779</td>
      <td>vasc_low</td>
      <td>17</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>4</td>
      <td>MHCI_high</td>
      <td>MHCI_high</td>
      <td>Tumor_high_MHCI_high</td>
      <td>Tumor_high_MHCI_high</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Lymphocytes B Cells (and Plasma Cells)</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>-0.223284</td>
      <td>0.124065</td>
      <td>0.148279</td>
      <td>-0.482537</td>
      <td>-0.294707</td>
      <td>-0.344041</td>
      <td>-0.235293</td>
      <td>0.083866</td>
      <td>-0.567738</td>
      <td>-0.239274</td>
      <td>-0.014502</td>
      <td>0.044275</td>
      <td>0.105805</td>
      <td>0.062049</td>
      <td>-0.386939</td>
      <td>-0.488036</td>
      <td>-0.447651</td>
      <td>0.072677</td>
      <td>0.316386</td>
      <td>0.408941</td>
      <td>0.290778</td>
      <td>-0.211729</td>
      <td>-0.364546</td>
      <td>-0.379476</td>
      <td>-0.014502</td>
      <td>0.143517</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>10072</td>
      <td>25947</td>
      <td>AAACACCAATAACTGC-1</td>
      <td>VB1_2</td>
      <td>VB1_410</td>
      <td>VB1_1610</td>
      <td>VB1_2358</td>
      <td>VB1_2145</td>
      <td>VB1_1193</td>
      <td>VB1_1279</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
</div>




```python
# turn nn_ids into a list in a new column
vb1_nn_list_temp = ['nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id']
vb1_nn['nn_list'] = vb1_nn[vb1_nn_list_temp].values.tolist()
vb1_nn.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>num_adjacent_nns</th>
      <th>nn_list</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VB1_AAACAACGAATAGTTC-1</td>
      <td>VB1</td>
      <td>VB1</td>
      <td>7174</td>
      <td>3783</td>
      <td>round1</td>
      <td>5068</td>
      <td>3601</td>
      <td>6</td>
      <td>6</td>
      <td>VB1</td>
      <td>0.000000</td>
      <td>vasc_low</td>
      <td>20</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_neg</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.409190</td>
      <td>-0.293434</td>
      <td>0.015793</td>
      <td>-0.623296</td>
      <td>-0.294338</td>
      <td>-0.538970</td>
      <td>0.160997</td>
      <td>-0.440288</td>
      <td>-0.049946</td>
      <td>-0.219445</td>
      <td>-0.184969</td>
      <td>0.010123</td>
      <td>0.001746</td>
      <td>0.523236</td>
      <td>-0.383538</td>
      <td>0.263978</td>
      <td>-0.317542</td>
      <td>-0.189507</td>
      <td>0.171521</td>
      <td>-0.565283</td>
      <td>-0.448105</td>
      <td>-0.309486</td>
      <td>-0.363641</td>
      <td>-0.373975</td>
      <td>-0.184969</td>
      <td>0.090032</td>
      <td>Low</td>
      <td>Low</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>31283</td>
      <td>26488</td>
      <td>AAACAACGAATAGTTC-1</td>
      <td>VB1_0</td>
      <td>VB1_804</td>
      <td>VB1_1719</td>
      <td>VB1_654</td>
      <td>VB1_1936</td>
      <td>VB1_492</td>
      <td>VB1_2478</td>
      <td>4</td>
      <td>[VB1_804, VB1_1719, VB1_654, VB1_1936, VB1_492...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VB1_AAACAAGTATCTCCCA-1</td>
      <td>VB1</td>
      <td>VB1</td>
      <td>3817</td>
      <td>2574</td>
      <td>round1</td>
      <td>3824</td>
      <td>2572</td>
      <td>1</td>
      <td>1</td>
      <td>VB1</td>
      <td>3.117370</td>
      <td>vasc_high</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>4</td>
      <td>MHCI_high</td>
      <td>MHCI_high</td>
      <td>Tumor_low_MHCI_high</td>
      <td>Tumor_low</td>
      <td>vasc_high</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC</td>
      <td>1_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>0.124843</td>
      <td>0.508282</td>
      <td>-0.370904</td>
      <td>-0.208191</td>
      <td>-0.135008</td>
      <td>0.254889</td>
      <td>0.354699</td>
      <td>0.618041</td>
      <td>0.174323</td>
      <td>0.644896</td>
      <td>0.223238</td>
      <td>0.065979</td>
      <td>0.010680</td>
      <td>0.645669</td>
      <td>0.955596</td>
      <td>-0.325863</td>
      <td>0.259486</td>
      <td>-0.095793</td>
      <td>-0.280748</td>
      <td>-0.611245</td>
      <td>0.417938</td>
      <td>-0.344595</td>
      <td>-0.302494</td>
      <td>-0.264353</td>
      <td>0.223238</td>
      <td>0.694163</td>
      <td>Low</td>
      <td>Low</td>
      <td>vasc_high_MN4_low</td>
      <td>vasc_low</td>
      <td>13243</td>
      <td>8788</td>
      <td>AAACAAGTATCTCCCA-1</td>
      <td>VB1_1</td>
      <td>VB1_724</td>
      <td>VB1_270</td>
      <td>VB1_1481</td>
      <td>VB1_2524</td>
      <td>VB1_336</td>
      <td>VB1_1646</td>
      <td>6</td>
      <td>[VB1_724, VB1_270, VB1_1481, VB1_2524, VB1_336...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VB1_AAACACCAATAACTGC-1</td>
      <td>VB1</td>
      <td>VB1</td>
      <td>6137</td>
      <td>3622</td>
      <td>round1</td>
      <td>4990</td>
      <td>3585</td>
      <td>0</td>
      <td>0</td>
      <td>VB1</td>
      <td>0.966779</td>
      <td>vasc_low</td>
      <td>17</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>4</td>
      <td>MHCI_high</td>
      <td>MHCI_high</td>
      <td>Tumor_high_MHCI_high</td>
      <td>Tumor_high_MHCI_high</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Lymphocytes B Cells (and Plasma Cells)</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>-0.223284</td>
      <td>0.124065</td>
      <td>0.148279</td>
      <td>-0.482537</td>
      <td>-0.294707</td>
      <td>-0.344041</td>
      <td>-0.235293</td>
      <td>0.083866</td>
      <td>-0.567738</td>
      <td>-0.239274</td>
      <td>-0.014502</td>
      <td>0.044275</td>
      <td>0.105805</td>
      <td>0.062049</td>
      <td>-0.386939</td>
      <td>-0.488036</td>
      <td>-0.447651</td>
      <td>0.072677</td>
      <td>0.316386</td>
      <td>0.408941</td>
      <td>0.290778</td>
      <td>-0.211729</td>
      <td>-0.364546</td>
      <td>-0.379476</td>
      <td>-0.014502</td>
      <td>0.143517</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>10072</td>
      <td>25947</td>
      <td>AAACACCAATAACTGC-1</td>
      <td>VB1_2</td>
      <td>VB1_410</td>
      <td>VB1_1610</td>
      <td>VB1_2358</td>
      <td>VB1_2145</td>
      <td>VB1_1193</td>
      <td>VB1_1279</td>
      <td>6</td>
      <td>[VB1_410, VB1_1610, VB1_2358, VB1_2145, VB1_11...</td>
    </tr>
  </tbody>
</table>
</div>



# - - - - - - - - - -
# VA2 VA2 VA2 VA2 VA2 VA2
# - - - - - - - - - -

# Some neighborhoods include spots that are farther than the immediate ring... ( distances > 8 ) these must be spots near the edge. I should remove these from the analysis. 


```python
# make into dfs
va2_id_df = pd.DataFrame(va2_nearest_ind, columns=['id','nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id'])
va2_dist_df = pd.DataFrame(va2_nearest_dist[:, 1:], columns=['nn1_dist','nn2_dist','nn3_dist','nn4_dist','nn5_dist','nn6_dist'])
print(va2_id_df.shape)
print(va2_dist_df.shape)
va2_id_df = va2_id_df.merge(va2_dist_df, left_index=True, right_index=True)
print(va2_id_df.shape)
va2_id_df.head(3)
```

    (4316, 7)
    (4316, 6)
    (4316, 13)





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>nn1_dist</th>
      <th>nn2_dist</th>
      <th>nn3_dist</th>
      <th>nn4_dist</th>
      <th>nn5_dist</th>
      <th>nn6_dist</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>1217</td>
      <td>467</td>
      <td>2483</td>
      <td>4204</td>
      <td>565</td>
      <td>2751</td>
      <td>407.977115</td>
      <td>407.977115</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>681</td>
      <td>2687</td>
      <td>3589</td>
      <td>2004</td>
      <td>3930</td>
      <td>2163</td>
      <td>407.977115</td>
      <td>407.977115</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>2444</td>
      <td>4259</td>
      <td>3524</td>
      <td>1222</td>
      <td>2630</td>
      <td>3521</td>
      <td>407.977115</td>
      <td>407.977115</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
    </tr>
  </tbody>
</table>
</div>




```python
# count the number of nearest neighbors directly adjacent
va2_id_df['num_adjacent_nns'] = va2_id_df[['nn1_dist','nn2_dist','nn3_dist','nn4_dist','nn5_dist','nn6_dist']].apply(lambda x: x[x < 420].count(), axis=1)
va2_id_df.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>nn1_dist</th>
      <th>nn2_dist</th>
      <th>nn3_dist</th>
      <th>nn4_dist</th>
      <th>nn5_dist</th>
      <th>nn6_dist</th>
      <th>num_adjacent_nns</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>1217</td>
      <td>467</td>
      <td>2483</td>
      <td>4204</td>
      <td>565</td>
      <td>2751</td>
      <td>407.977115</td>
      <td>407.977115</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>6</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>681</td>
      <td>2687</td>
      <td>3589</td>
      <td>2004</td>
      <td>3930</td>
      <td>2163</td>
      <td>407.977115</td>
      <td>407.977115</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>2444</td>
      <td>4259</td>
      <td>3524</td>
      <td>1222</td>
      <td>2630</td>
      <td>3521</td>
      <td>407.977115</td>
      <td>407.977115</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
</div>




```python
# function to add tissue slice identifier to an id column
def add_tissue_slice_id(df, tissue_slice_id, col_list):
    for col in col_list:
        df[col] = tissue_slice_id + '_' + df[col].astype(str)
    return(df)
```


```python
# Add tissue slice identifier VA to all ids
va2_temp = va2_id_df.copy()

col_list = ['id','nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id']
va2_temp = add_tissue_slice_id(va2_temp, 'VA2', col_list)
display(va2_temp.head(3))

va2_id_df = va2_temp.copy()
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>nn1_dist</th>
      <th>nn2_dist</th>
      <th>nn3_dist</th>
      <th>nn4_dist</th>
      <th>nn5_dist</th>
      <th>nn6_dist</th>
      <th>num_adjacent_nns</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VA2_0</td>
      <td>VA2_1217</td>
      <td>VA2_467</td>
      <td>VA2_2483</td>
      <td>VA2_4204</td>
      <td>VA2_565</td>
      <td>VA2_2751</td>
      <td>407.977115</td>
      <td>407.977115</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>6</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VA2_1</td>
      <td>VA2_681</td>
      <td>VA2_2687</td>
      <td>VA2_3589</td>
      <td>VA2_2004</td>
      <td>VA2_3930</td>
      <td>VA2_2163</td>
      <td>407.977115</td>
      <td>407.977115</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VA2_2</td>
      <td>VA2_2444</td>
      <td>VA2_4259</td>
      <td>VA2_3524</td>
      <td>VA2_1222</td>
      <td>VA2_2630</td>
      <td>VA2_3521</td>
      <td>407.977115</td>
      <td>407.977115</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>409.382217</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
</div>



```python
# drop neighborhoods with spots farther than 420 in nn6_dist
va2_id_df_sub1 = va2_id_df[va2_id_df['nn6_dist']<420]
print(va2_id_df_sub1.shape)

# allow one farther (so 5 neighbors immediately next to center spot)
va2_id_df_sub2 = va2_id_df[va2_id_df['nn5_dist']<420]
print(va2_id_df_sub2.shape)

# allow two farther (so 4 neighbors immediately next to center spot)
va2_id_df_sub3 = va2_id_df[va2_id_df['nn4_dist']<420]
print(va2_id_df_sub3.shape)
```

    (4021, 14)
    (4107, 14)
    (4237, 14)



```python
# let's require all 6 neighbors
va2_nn = va2_nn_pre.merge(va2_id_df[['id','nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id','num_adjacent_nns']], left_index=True, right_index=True)
print(va2_nn.shape)
va2_nn.head(3)
```

    (4316, 68)





<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>num_adjacent_nns</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VA2_AAACAAGTATCTCCCA-1</td>
      <td>VA2</td>
      <td>VA2</td>
      <td>45867</td>
      <td>8515</td>
      <td>round2</td>
      <td>24897</td>
      <td>7041</td>
      <td>5</td>
      <td>5</td>
      <td>VA2</td>
      <td>1.296336</td>
      <td>vasc_low</td>
      <td>23</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC</td>
      <td>5_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>-0.450968</td>
      <td>-0.448330</td>
      <td>-0.691477</td>
      <td>-0.374355</td>
      <td>-0.693455</td>
      <td>-0.806181</td>
      <td>-0.424467</td>
      <td>-0.261943</td>
      <td>-0.684682</td>
      <td>-0.765453</td>
      <td>-0.199923</td>
      <td>-0.133973</td>
      <td>-0.213757</td>
      <td>-0.290269</td>
      <td>0.029304</td>
      <td>-0.261245</td>
      <td>-0.741362</td>
      <td>-0.395326</td>
      <td>-0.297578</td>
      <td>0.182372</td>
      <td>-0.583172</td>
      <td>0.484041</td>
      <td>-0.332599</td>
      <td>-0.704982</td>
      <td>-0.199923</td>
      <td>-0.156565</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>7764.085688</td>
      <td>20534.804310</td>
      <td>AAACAAGTATCTCCCA-1</td>
      <td>VA2_0</td>
      <td>VA2_1217</td>
      <td>VA2_467</td>
      <td>VA2_2483</td>
      <td>VA2_4204</td>
      <td>VA2_565</td>
      <td>VA2_2751</td>
      <td>6</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VA2_AAACACCAATAACTGC-1</td>
      <td>VA2</td>
      <td>VA2</td>
      <td>11849</td>
      <td>5082</td>
      <td>round2</td>
      <td>22823</td>
      <td>5881</td>
      <td>0</td>
      <td>0</td>
      <td>VA2</td>
      <td>1.223823</td>
      <td>vasc_low</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>6</td>
      <td>MHCI_high</td>
      <td>MHCI_high</td>
      <td>Tumor_low_MHCI_high</td>
      <td>Tumor_low</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Lymphocytes B Cells (and Plasma Cells)</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>0.465755</td>
      <td>0.476869</td>
      <td>0.801786</td>
      <td>0.685619</td>
      <td>0.728101</td>
      <td>-0.440744</td>
      <td>-0.355418</td>
      <td>0.451091</td>
      <td>0.679239</td>
      <td>0.565557</td>
      <td>0.301495</td>
      <td>0.099595</td>
      <td>0.168763</td>
      <td>0.394440</td>
      <td>0.279309</td>
      <td>0.396732</td>
      <td>0.534097</td>
      <td>0.368466</td>
      <td>0.618810</td>
      <td>0.647113</td>
      <td>0.167148</td>
      <td>-0.005985</td>
      <td>0.455600</td>
      <td>0.641533</td>
      <td>0.301495</td>
      <td>0.153607</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>24672.634008</td>
      <td>23846.305076</td>
      <td>AAACACCAATAACTGC-1</td>
      <td>VA2_1</td>
      <td>VA2_681</td>
      <td>VA2_2687</td>
      <td>VA2_3589</td>
      <td>VA2_2004</td>
      <td>VA2_3930</td>
      <td>VA2_2163</td>
      <td>6</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VA2_AAACAGAGCGACTCCT-1</td>
      <td>VA2</td>
      <td>VA2</td>
      <td>43834</td>
      <td>8662</td>
      <td>round2</td>
      <td>24885</td>
      <td>7403</td>
      <td>5</td>
      <td>5</td>
      <td>VA2</td>
      <td>0.205495</td>
      <td>vasc_low</td>
      <td>24</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC</td>
      <td>5_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>-0.587803</td>
      <td>-0.647125</td>
      <td>-0.203480</td>
      <td>-0.608944</td>
      <td>-0.603551</td>
      <td>-0.840184</td>
      <td>-0.328385</td>
      <td>-0.589335</td>
      <td>-0.659993</td>
      <td>-0.799260</td>
      <td>-0.441499</td>
      <td>-0.277116</td>
      <td>-0.349001</td>
      <td>-0.654693</td>
      <td>0.151233</td>
      <td>-0.594247</td>
      <td>-0.666882</td>
      <td>-0.478262</td>
      <td>-0.676907</td>
      <td>-0.817864</td>
      <td>-0.660876</td>
      <td>0.140619</td>
      <td>-0.768692</td>
      <td>-0.746999</td>
      <td>-0.441499</td>
      <td>-0.126745</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>9484.342696</td>
      <td>7768.555088</td>
      <td>AAACAGAGCGACTCCT-1</td>
      <td>VA2_2</td>
      <td>VA2_2444</td>
      <td>VA2_4259</td>
      <td>VA2_3524</td>
      <td>VA2_1222</td>
      <td>VA2_2630</td>
      <td>VA2_3521</td>
      <td>6</td>
    </tr>
  </tbody>
</table>
</div>




```python
# turn nn_ids into a list in a new column
va2_nn_list_temp = ['nn1_id','nn2_id','nn3_id','nn4_id','nn5_id','nn6_id']
va2_nn['nn_list'] = va2_nn[va2_nn_list_temp].values.tolist()
va2_nn.head(3)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>num_adjacent_nns</th>
      <th>nn_list</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VA2_AAACAAGTATCTCCCA-1</td>
      <td>VA2</td>
      <td>VA2</td>
      <td>45867</td>
      <td>8515</td>
      <td>round2</td>
      <td>24897</td>
      <td>7041</td>
      <td>5</td>
      <td>5</td>
      <td>VA2</td>
      <td>1.296336</td>
      <td>vasc_low</td>
      <td>23</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC</td>
      <td>5_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>-0.450968</td>
      <td>-0.448330</td>
      <td>-0.691477</td>
      <td>-0.374355</td>
      <td>-0.693455</td>
      <td>-0.806181</td>
      <td>-0.424467</td>
      <td>-0.261943</td>
      <td>-0.684682</td>
      <td>-0.765453</td>
      <td>-0.199923</td>
      <td>-0.133973</td>
      <td>-0.213757</td>
      <td>-0.290269</td>
      <td>0.029304</td>
      <td>-0.261245</td>
      <td>-0.741362</td>
      <td>-0.395326</td>
      <td>-0.297578</td>
      <td>0.182372</td>
      <td>-0.583172</td>
      <td>0.484041</td>
      <td>-0.332599</td>
      <td>-0.704982</td>
      <td>-0.199923</td>
      <td>-0.156565</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>7764.085688</td>
      <td>20534.804310</td>
      <td>AAACAAGTATCTCCCA-1</td>
      <td>VA2_0</td>
      <td>VA2_1217</td>
      <td>VA2_467</td>
      <td>VA2_2483</td>
      <td>VA2_4204</td>
      <td>VA2_565</td>
      <td>VA2_2751</td>
      <td>6</td>
      <td>[VA2_1217, VA2_467, VA2_2483, VA2_4204, VA2_56...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VA2_AAACACCAATAACTGC-1</td>
      <td>VA2</td>
      <td>VA2</td>
      <td>11849</td>
      <td>5082</td>
      <td>round2</td>
      <td>22823</td>
      <td>5881</td>
      <td>0</td>
      <td>0</td>
      <td>VA2</td>
      <td>1.223823</td>
      <td>vasc_low</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>6</td>
      <td>MHCI_high</td>
      <td>MHCI_high</td>
      <td>Tumor_low_MHCI_high</td>
      <td>Tumor_low</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>Lymphocytes B Cells (and Plasma Cells)</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>0_Lymphocytes_B_Cells</td>
      <td>0.465755</td>
      <td>0.476869</td>
      <td>0.801786</td>
      <td>0.685619</td>
      <td>0.728101</td>
      <td>-0.440744</td>
      <td>-0.355418</td>
      <td>0.451091</td>
      <td>0.679239</td>
      <td>0.565557</td>
      <td>0.301495</td>
      <td>0.099595</td>
      <td>0.168763</td>
      <td>0.394440</td>
      <td>0.279309</td>
      <td>0.396732</td>
      <td>0.534097</td>
      <td>0.368466</td>
      <td>0.618810</td>
      <td>0.647113</td>
      <td>0.167148</td>
      <td>-0.005985</td>
      <td>0.455600</td>
      <td>0.641533</td>
      <td>0.301495</td>
      <td>0.153607</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>24672.634008</td>
      <td>23846.305076</td>
      <td>AAACACCAATAACTGC-1</td>
      <td>VA2_1</td>
      <td>VA2_681</td>
      <td>VA2_2687</td>
      <td>VA2_3589</td>
      <td>VA2_2004</td>
      <td>VA2_3930</td>
      <td>VA2_2163</td>
      <td>6</td>
      <td>[VA2_681, VA2_2687, VA2_3589, VA2_2004, VA2_39...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VA2_AAACAGAGCGACTCCT-1</td>
      <td>VA2</td>
      <td>VA2</td>
      <td>43834</td>
      <td>8662</td>
      <td>round2</td>
      <td>24885</td>
      <td>7403</td>
      <td>5</td>
      <td>5</td>
      <td>VA2</td>
      <td>0.205495</td>
      <td>vasc_low</td>
      <td>24</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_low</td>
      <td>NaN</td>
      <td>NaN</td>
      <td>SCLC</td>
      <td>5_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>-0.587803</td>
      <td>-0.647125</td>
      <td>-0.203480</td>
      <td>-0.608944</td>
      <td>-0.603551</td>
      <td>-0.840184</td>
      <td>-0.328385</td>
      <td>-0.589335</td>
      <td>-0.659993</td>
      <td>-0.799260</td>
      <td>-0.441499</td>
      <td>-0.277116</td>
      <td>-0.349001</td>
      <td>-0.654693</td>
      <td>0.151233</td>
      <td>-0.594247</td>
      <td>-0.666882</td>
      <td>-0.478262</td>
      <td>-0.676907</td>
      <td>-0.817864</td>
      <td>-0.660876</td>
      <td>0.140619</td>
      <td>-0.768692</td>
      <td>-0.746999</td>
      <td>-0.441499</td>
      <td>-0.126745</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>9484.342696</td>
      <td>7768.555088</td>
      <td>AAACAGAGCGACTCCT-1</td>
      <td>VA2_2</td>
      <td>VA2_2444</td>
      <td>VA2_4259</td>
      <td>VA2_3524</td>
      <td>VA2_1222</td>
      <td>VA2_2630</td>
      <td>VA2_3521</td>
      <td>6</td>
      <td>[VA2_2444, VA2_4259, VA2_3524, VA2_1222, VA2_2...</td>
    </tr>
  </tbody>
</table>
</div>



# Combine VA1, VB1, and VA2 by concatenating them


```python
ab1_a2_nn = pd.concat([va1_nn, vb1_nn, va2_nn])
print(va1_nn.shape)
print(vb1_nn.shape)
print(va2_nn.shape)
print(ab1_a2_nn.shape)
```

    (2806, 69)
    (2589, 69)
    (4316, 69)
    (9711, 69)



```python
display(ab1_a2_nn.head(2))
```


<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>num_adjacent_nns</th>
      <th>nn_list</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VA1_AAACAACGAATAGTTC-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>3625</td>
      <td>2571</td>
      <td>round1</td>
      <td>3649</td>
      <td>2571</td>
      <td>10</td>
      <td>10</td>
      <td>VA1</td>
      <td>2.648104</td>
      <td>vasc_high</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>2</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_high</td>
      <td>Respiratory epithelium</td>
      <td>Respiratory epithelium</td>
      <td>Respiratory Epithelium (Normal Lung)</td>
      <td>10_Respiratory_Epithelium_Lung</td>
      <td>7_10_Normal_Lung</td>
      <td>-0.056807</td>
      <td>-0.259847</td>
      <td>0.296748</td>
      <td>-0.038431</td>
      <td>0.513527</td>
      <td>-0.007136</td>
      <td>0.412196</td>
      <td>0.404657</td>
      <td>-0.158465</td>
      <td>0.690663</td>
      <td>-0.046498</td>
      <td>0.110800</td>
      <td>0.102126</td>
      <td>0.654483</td>
      <td>-0.176035</td>
      <td>-0.321725</td>
      <td>-0.469550</td>
      <td>0.189499</td>
      <td>0.139214</td>
      <td>-0.165127</td>
      <td>0.015559</td>
      <td>-0.755217</td>
      <td>0.633032</td>
      <td>-0.326533</td>
      <td>-0.046498</td>
      <td>0.031719</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_high_MN4_high</td>
      <td>vasc_high_MN4_upper80</td>
      <td>30730.0</td>
      <td>26443.0</td>
      <td>AAACAACGAATAGTTC-1</td>
      <td>VA1_0</td>
      <td>VA1_854</td>
      <td>VA1_1863</td>
      <td>VA1_2089</td>
      <td>VA1_517</td>
      <td>VA1_2692</td>
      <td>VA1_2472</td>
      <td>3</td>
      <td>[VA1_854, VA1_1863, VA1_2089, VA1_517, VA1_269...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VA1_AAACAAGTATCTCCCA-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>5804</td>
      <td>3329</td>
      <td>round1</td>
      <td>4956</td>
      <td>3308</td>
      <td>6</td>
      <td>6</td>
      <td>VA1</td>
      <td>1.001716</td>
      <td>vasc_low</td>
      <td>13</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>2</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_low</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.432040</td>
      <td>-0.363696</td>
      <td>-0.665610</td>
      <td>-0.380249</td>
      <td>-0.272964</td>
      <td>-0.616570</td>
      <td>-0.133588</td>
      <td>-0.277191</td>
      <td>-0.532077</td>
      <td>-0.284190</td>
      <td>-0.442197</td>
      <td>-0.206272</td>
      <td>-0.292946</td>
      <td>-0.220110</td>
      <td>0.816411</td>
      <td>-0.422218</td>
      <td>-0.542976</td>
      <td>-0.323644</td>
      <td>-0.431187</td>
      <td>-0.713442</td>
      <td>-0.122593</td>
      <td>0.079799</td>
      <td>-0.348075</td>
      <td>-0.420942</td>
      <td>-0.442197</td>
      <td>-0.484106</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>12698.0</td>
      <td>8743.0</td>
      <td>AAACAAGTATCTCCCA-1</td>
      <td>VA1_1</td>
      <td>VA1_284</td>
      <td>VA1_771</td>
      <td>VA1_352</td>
      <td>VA1_1783</td>
      <td>VA1_2738</td>
      <td>VA1_1604</td>
      <td>6</td>
      <td>[VA1_284, VA1_771, VA1_352, VA1_1783, VA1_2738...</td>
    </tr>
  </tbody>
</table>
</div>


# We already have the GSVA enrichments in the dataset. Let's continue. 


```python
cols_to_extract = ['Macrophage', 'T_Cell', 'CD8_T_Cell', 'NK_Cell', 'B_Cell', 'DC', 'Endothelial_Cell', 'LEC', 'MN4_EC_Phenotype']
```


```python
# function to extract enrichment values for each geneset of interest, for each nearest neighbor
def extract_enrichments_per_nn(df, cols_to_extract):
    # initialize results dataframe
    res_df = pd.DataFrame()
    # cycle through each col_to_extract, i.e. the current enrichment
    # cycle through each id, and extract values for each nearest neighbor. 
    for curr_geneset in cols_to_extract:
        print('- - - - - - - -')
        print(curr_geneset)
        print(' ')
        
        # initialize results columns for the current geneset as a column for each nn
        for i in range(1,7):
            df[curr_geneset + '_nn_' + str(i)] = np.nan

        # initialize results column with averages of the list
        df[curr_geneset + '_nn_mean'] = np.nan
        
        for curr_id in df['id'].unique(): 
            
            # get values for each neighbor
            # get current neighbor list to cycle through...
            # need to do .values[0] because it comes out as a numpy array which contains the list we want...
            curr_neighbor_list = df[df['id'] == curr_id]['nn_list'].values[0]
            for i, curr_neighbor in enumerate(curr_neighbor_list):
                # get values for each neighbor and place them into a list in a new column
                # add it to the row of the current id, in the current geneset column. 
                temp_enrichment_value = df[df['id']==curr_neighbor][curr_geneset].values[0]
                df.loc[df['id']==curr_id,curr_geneset + '_nn_' + str(i+1)] = temp_enrichment_value
                #res_list_temp.append(df[df['id']==curr_neighbor][curr_geneset].values[0])
            
        # combine the results columns into a list 
        df[curr_geneset + '_nn_vals'] = df[[curr_geneset + '_nn_1', curr_geneset + '_nn_2', curr_geneset + '_nn_3', curr_geneset + '_nn_4', curr_geneset + '_nn_5', curr_geneset + '_nn_6']].values.tolist()
        # find averages of the values and put in the mean results column for this geneset
        df[curr_geneset + '_nn_mean'] = df[[curr_geneset + '_nn_1', curr_geneset + '_nn_2', curr_geneset + '_nn_3', curr_geneset + '_nn_4', curr_geneset + '_nn_5', curr_geneset + '_nn_6']].values.mean(axis=1)
    
    return df
```


```python
test = extract_enrichments_per_nn(ab1_a2_nn, cols_to_extract)
display(test.head(6))
```

    - - - - - - - -
    Macrophage
     
    - - - - - - - -
    T_Cell
     
    - - - - - - - -
    CD8_T_Cell
     
    - - - - - - - -
    NK_Cell
     
    - - - - - - - -
    B_Cell
     
    - - - - - - - -
    DC
     
    - - - - - - - -
    Endothelial_Cell
     
    - - - - - - - -
    LEC
     
    - - - - - - - -
    MN4_EC_Phenotype
     



<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Barcode</th>
      <th>Tissue_Slice</th>
      <th>orig.ident</th>
      <th>nCount_Spatial</th>
      <th>nFeature_Spatial</th>
      <th>visium_round</th>
      <th>nCount_SCT</th>
      <th>nFeature_SCT</th>
      <th>SCT_snn_res.0.8</th>
      <th>seurat_clusters</th>
      <th>slice_ident</th>
      <th>vasc_summed_expr</th>
      <th>vasc_label</th>
      <th>Tumor_score</th>
      <th>Tumor_label1</th>
      <th>Tumor_label</th>
      <th>MHCI_score</th>
      <th>MHCI_label1</th>
      <th>MHCI_label</th>
      <th>Tumor_MHCI_label</th>
      <th>Tumor_MHCI_label2</th>
      <th>vasc_label1</th>
      <th>Navin_annotations</th>
      <th>Navin_annotations_simplified</th>
      <th>Navin_Clusters1</th>
      <th>Navin_Clusters</th>
      <th>Navin_Clusters2</th>
      <th>All_Immune</th>
      <th>Angiogenesis</th>
      <th>AyersIFNG</th>
      <th>B_Cell</th>
      <th>CD8_T_Cell</th>
      <th>DC</th>
      <th>Endothelial_Activation</th>
      <th>Endothelial_Cell</th>
      <th>Endothelial_Chemokines</th>
      <th>Exhaustion</th>
      <th>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th>
      <th>GOBP_LEUKOCYTE_MIGRATION</th>
      <th>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th>
      <th>ICB_Targets</th>
      <th>LEC</th>
      <th>M1_Macrophage</th>
      <th>M2_Macrophage</th>
      <th>MN4_EC_Phenotype</th>
      <th>MN4_EC_Phenotype_Top30</th>
      <th>Macrophage</th>
      <th>NK_Cell</th>
      <th>Proliferation</th>
      <th>T_Cell</th>
      <th>T_Reg</th>
      <th>Upregulated_by_2_3_CGAMP</th>
      <th>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th>
      <th>MN4_EC_Phenotype_Label</th>
      <th>MN4_EC_Percentile_Label</th>
      <th>vasc_MN4_label</th>
      <th>vasc_MN4_percentile_label</th>
      <th>x</th>
      <th>y</th>
      <th>cell</th>
      <th>id</th>
      <th>nn1_id</th>
      <th>nn2_id</th>
      <th>nn3_id</th>
      <th>nn4_id</th>
      <th>nn5_id</th>
      <th>nn6_id</th>
      <th>num_adjacent_nns</th>
      <th>nn_list</th>
      <th>Macrophage_nn_1</th>
      <th>Macrophage_nn_2</th>
      <th>Macrophage_nn_3</th>
      <th>Macrophage_nn_4</th>
      <th>Macrophage_nn_5</th>
      <th>Macrophage_nn_6</th>
      <th>Macrophage_nn_mean</th>
      <th>Macrophage_nn_vals</th>
      <th>T_Cell_nn_1</th>
      <th>T_Cell_nn_2</th>
      <th>T_Cell_nn_3</th>
      <th>T_Cell_nn_4</th>
      <th>T_Cell_nn_5</th>
      <th>T_Cell_nn_6</th>
      <th>T_Cell_nn_mean</th>
      <th>T_Cell_nn_vals</th>
      <th>CD8_T_Cell_nn_1</th>
      <th>CD8_T_Cell_nn_2</th>
      <th>CD8_T_Cell_nn_3</th>
      <th>CD8_T_Cell_nn_4</th>
      <th>CD8_T_Cell_nn_5</th>
      <th>CD8_T_Cell_nn_6</th>
      <th>CD8_T_Cell_nn_mean</th>
      <th>CD8_T_Cell_nn_vals</th>
      <th>NK_Cell_nn_1</th>
      <th>NK_Cell_nn_2</th>
      <th>NK_Cell_nn_3</th>
      <th>NK_Cell_nn_4</th>
      <th>NK_Cell_nn_5</th>
      <th>NK_Cell_nn_6</th>
      <th>NK_Cell_nn_mean</th>
      <th>NK_Cell_nn_vals</th>
      <th>B_Cell_nn_1</th>
      <th>B_Cell_nn_2</th>
      <th>B_Cell_nn_3</th>
      <th>B_Cell_nn_4</th>
      <th>B_Cell_nn_5</th>
      <th>B_Cell_nn_6</th>
      <th>B_Cell_nn_mean</th>
      <th>B_Cell_nn_vals</th>
      <th>DC_nn_1</th>
      <th>DC_nn_2</th>
      <th>DC_nn_3</th>
      <th>DC_nn_4</th>
      <th>DC_nn_5</th>
      <th>DC_nn_6</th>
      <th>DC_nn_mean</th>
      <th>DC_nn_vals</th>
      <th>Endothelial_Cell_nn_1</th>
      <th>Endothelial_Cell_nn_2</th>
      <th>Endothelial_Cell_nn_3</th>
      <th>Endothelial_Cell_nn_4</th>
      <th>Endothelial_Cell_nn_5</th>
      <th>Endothelial_Cell_nn_6</th>
      <th>Endothelial_Cell_nn_mean</th>
      <th>Endothelial_Cell_nn_vals</th>
      <th>LEC_nn_1</th>
      <th>LEC_nn_2</th>
      <th>LEC_nn_3</th>
      <th>LEC_nn_4</th>
      <th>LEC_nn_5</th>
      <th>LEC_nn_6</th>
      <th>LEC_nn_mean</th>
      <th>LEC_nn_vals</th>
      <th>MN4_EC_Phenotype_nn_1</th>
      <th>MN4_EC_Phenotype_nn_2</th>
      <th>MN4_EC_Phenotype_nn_3</th>
      <th>MN4_EC_Phenotype_nn_4</th>
      <th>MN4_EC_Phenotype_nn_5</th>
      <th>MN4_EC_Phenotype_nn_6</th>
      <th>MN4_EC_Phenotype_nn_mean</th>
      <th>MN4_EC_Phenotype_nn_vals</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>VA1_AAACAACGAATAGTTC-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>3625</td>
      <td>2571</td>
      <td>round1</td>
      <td>3649</td>
      <td>2571</td>
      <td>10</td>
      <td>10</td>
      <td>VA1</td>
      <td>2.648104</td>
      <td>vasc_high</td>
      <td>11</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>2</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_high</td>
      <td>Respiratory epithelium</td>
      <td>Respiratory epithelium</td>
      <td>Respiratory Epithelium (Normal Lung)</td>
      <td>10_Respiratory_Epithelium_Lung</td>
      <td>7_10_Normal_Lung</td>
      <td>-0.056807</td>
      <td>-0.259847</td>
      <td>0.296748</td>
      <td>-0.038431</td>
      <td>0.513527</td>
      <td>-0.007136</td>
      <td>0.412196</td>
      <td>0.404657</td>
      <td>-0.158465</td>
      <td>0.690663</td>
      <td>-0.046498</td>
      <td>0.110800</td>
      <td>0.102126</td>
      <td>0.654483</td>
      <td>-0.176035</td>
      <td>-0.321725</td>
      <td>-0.469550</td>
      <td>0.189499</td>
      <td>0.139214</td>
      <td>-0.165127</td>
      <td>0.015559</td>
      <td>-0.755217</td>
      <td>0.633032</td>
      <td>-0.326533</td>
      <td>-0.046498</td>
      <td>0.031719</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_high_MN4_high</td>
      <td>vasc_high_MN4_upper80</td>
      <td>30730.0</td>
      <td>26443.0</td>
      <td>AAACAACGAATAGTTC-1</td>
      <td>VA1_0</td>
      <td>VA1_854</td>
      <td>VA1_1863</td>
      <td>VA1_2089</td>
      <td>VA1_517</td>
      <td>VA1_2692</td>
      <td>VA1_2472</td>
      <td>3</td>
      <td>[VA1_854, VA1_1863, VA1_2089, VA1_517, VA1_269...</td>
      <td>-0.140374</td>
      <td>-0.683461</td>
      <td>-0.682427</td>
      <td>-0.025348</td>
      <td>-0.183020</td>
      <td>-0.689811</td>
      <td>-0.400740</td>
      <td>[-0.140373503285981, -0.683461428931692, -0.68...</td>
      <td>0.478743</td>
      <td>-0.114692</td>
      <td>-0.053042</td>
      <td>-0.077262</td>
      <td>0.058715</td>
      <td>0.163426</td>
      <td>0.075981</td>
      <td>[0.478743199137973, -0.114691753402848, -0.053...</td>
      <td>0.498593</td>
      <td>-0.155994</td>
      <td>-0.093356</td>
      <td>-0.117871</td>
      <td>-0.244188</td>
      <td>-0.240874</td>
      <td>-0.058948</td>
      <td>[0.498593160799533, -0.155993596157572, -0.093...</td>
      <td>0.049122</td>
      <td>0.468691</td>
      <td>0.145700</td>
      <td>0.489515</td>
      <td>-0.058071</td>
      <td>0.426193</td>
      <td>0.253525</td>
      <td>[0.0491215471122242, 0.468690624713802, 0.1456...</td>
      <td>0.180875</td>
      <td>0.452794</td>
      <td>0.212476</td>
      <td>0.420312</td>
      <td>0.267154</td>
      <td>0.500322</td>
      <td>0.338989</td>
      <td>[0.180875180364392, 0.452794442072525, 0.21247...</td>
      <td>-0.042891</td>
      <td>0.122026</td>
      <td>-0.513131</td>
      <td>-0.539717</td>
      <td>-0.384229</td>
      <td>0.003418</td>
      <td>-0.225754</td>
      <td>[-0.0428909507251363, 0.122025706475451, -0.51...</td>
      <td>-0.102464</td>
      <td>0.648017</td>
      <td>-0.245869</td>
      <td>0.397033</td>
      <td>-0.221405</td>
      <td>0.628478</td>
      <td>0.183965</td>
      <td>[-0.102463880569944, 0.648017492909985, -0.245...</td>
      <td>-0.168113</td>
      <td>-0.172881</td>
      <td>-0.078916</td>
      <td>0.954805</td>
      <td>-0.214143</td>
      <td>-0.202240</td>
      <td>0.019752</td>
      <td>[-0.168112618939413, -0.172881026928877, -0.07...</td>
      <td>0.087962</td>
      <td>0.045507</td>
      <td>-0.151581</td>
      <td>-0.066657</td>
      <td>0.088639</td>
      <td>0.121551</td>
      <td>0.020904</td>
      <td>[0.0879624949143438, 0.0455074885754431, -0.15...</td>
    </tr>
    <tr>
      <th>1</th>
      <td>VA1_AAACAAGTATCTCCCA-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>5804</td>
      <td>3329</td>
      <td>round1</td>
      <td>4956</td>
      <td>3308</td>
      <td>6</td>
      <td>6</td>
      <td>VA1</td>
      <td>1.001716</td>
      <td>vasc_low</td>
      <td>13</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>2</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_low</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.432040</td>
      <td>-0.363696</td>
      <td>-0.665610</td>
      <td>-0.380249</td>
      <td>-0.272964</td>
      <td>-0.616570</td>
      <td>-0.133588</td>
      <td>-0.277191</td>
      <td>-0.532077</td>
      <td>-0.284190</td>
      <td>-0.442197</td>
      <td>-0.206272</td>
      <td>-0.292946</td>
      <td>-0.220110</td>
      <td>0.816411</td>
      <td>-0.422218</td>
      <td>-0.542976</td>
      <td>-0.323644</td>
      <td>-0.431187</td>
      <td>-0.713442</td>
      <td>-0.122593</td>
      <td>0.079799</td>
      <td>-0.348075</td>
      <td>-0.420942</td>
      <td>-0.442197</td>
      <td>-0.484106</td>
      <td>Low</td>
      <td>Below20Percentile</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>12698.0</td>
      <td>8743.0</td>
      <td>AAACAAGTATCTCCCA-1</td>
      <td>VA1_1</td>
      <td>VA1_284</td>
      <td>VA1_771</td>
      <td>VA1_352</td>
      <td>VA1_1783</td>
      <td>VA1_2738</td>
      <td>VA1_1604</td>
      <td>6</td>
      <td>[VA1_284, VA1_771, VA1_352, VA1_1783, VA1_2738...</td>
      <td>-0.331374</td>
      <td>-0.253458</td>
      <td>-0.339693</td>
      <td>-0.704575</td>
      <td>-0.194675</td>
      <td>-0.698858</td>
      <td>-0.420439</td>
      <td>[-0.33137442702708, -0.253457547402507, -0.339...</td>
      <td>-0.365543</td>
      <td>0.087386</td>
      <td>-0.012532</td>
      <td>-0.199258</td>
      <td>0.075946</td>
      <td>-0.110388</td>
      <td>-0.087398</td>
      <td>[-0.365543343516556, 0.0873864569432376, -0.01...</td>
      <td>0.153244</td>
      <td>-0.266860</td>
      <td>-0.331099</td>
      <td>-0.202421</td>
      <td>-0.266560</td>
      <td>-0.151491</td>
      <td>-0.177531</td>
      <td>[0.153244482371191, -0.266860116069477, -0.331...</td>
      <td>0.546733</td>
      <td>-0.116703</td>
      <td>0.261950</td>
      <td>0.451097</td>
      <td>0.385720</td>
      <td>0.079316</td>
      <td>0.268019</td>
      <td>[0.546732833950551, -0.116702778523135, 0.2619...</td>
      <td>0.555634</td>
      <td>0.209008</td>
      <td>-0.146471</td>
      <td>0.079159</td>
      <td>0.620314</td>
      <td>-0.387013</td>
      <td>0.155105</td>
      <td>[0.555633961432688, 0.209008116797947, -0.1464...</td>
      <td>-0.141677</td>
      <td>-0.607476</td>
      <td>-0.188773</td>
      <td>-0.044288</td>
      <td>-0.523967</td>
      <td>-0.548042</td>
      <td>-0.342371</td>
      <td>[-0.141677294770587, -0.607476141948519, -0.18...</td>
      <td>0.165504</td>
      <td>-0.272119</td>
      <td>-0.291141</td>
      <td>0.141134</td>
      <td>-0.439560</td>
      <td>-0.067795</td>
      <td>-0.127330</td>
      <td>[0.165503689409276, -0.272119318571468, -0.291...</td>
      <td>-0.391778</td>
      <td>-0.259152</td>
      <td>-0.323065</td>
      <td>-0.190838</td>
      <td>-0.259652</td>
      <td>-0.210794</td>
      <td>-0.272546</td>
      <td>[-0.391778355670959, -0.259151830365927, -0.32...</td>
      <td>-0.052664</td>
      <td>-0.143471</td>
      <td>-0.302900</td>
      <td>-0.244650</td>
      <td>-0.359702</td>
      <td>-0.275860</td>
      <td>-0.229875</td>
      <td>[-0.0526644931955577, -0.143470613692345, -0.3...</td>
    </tr>
    <tr>
      <th>2</th>
      <td>VA1_AAACAATCTACTAGCA-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>4980</td>
      <td>3146</td>
      <td>round1</td>
      <td>4711</td>
      <td>3146</td>
      <td>6</td>
      <td>6</td>
      <td>VA1</td>
      <td>1.101286</td>
      <td>vasc_low</td>
      <td>19</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_low</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>SCLC + Effector Immune Cells (T+NK)</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>6_SCLC_Effector_Immune</td>
      <td>-0.312254</td>
      <td>-0.142246</td>
      <td>0.215974</td>
      <td>-0.522873</td>
      <td>-0.254499</td>
      <td>0.344170</td>
      <td>-0.096425</td>
      <td>0.015784</td>
      <td>-0.519466</td>
      <td>-0.273987</td>
      <td>-0.296916</td>
      <td>-0.112254</td>
      <td>-0.139708</td>
      <td>0.574189</td>
      <td>-0.242849</td>
      <td>-0.389510</td>
      <td>-0.529344</td>
      <td>-0.039327</td>
      <td>0.055921</td>
      <td>-0.707380</td>
      <td>-0.091447</td>
      <td>0.049511</td>
      <td>-0.304126</td>
      <td>0.924392</td>
      <td>-0.296916</td>
      <td>-0.466480</td>
      <td>High</td>
      <td>Low</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>29633.0</td>
      <td>20870.0</td>
      <td>AAACAATCTACTAGCA-1</td>
      <td>VA1_2</td>
      <td>VA1_2207</td>
      <td>VA1_1663</td>
      <td>VA1_1116</td>
      <td>VA1_1210</td>
      <td>VA1_920</td>
      <td>VA1_1124</td>
      <td>4</td>
      <td>[VA1_2207, VA1_1663, VA1_1116, VA1_1210, VA1_9...</td>
      <td>-0.102161</td>
      <td>-0.094041</td>
      <td>-0.714331</td>
      <td>-0.092432</td>
      <td>-0.111016</td>
      <td>-0.227783</td>
      <td>-0.223627</td>
      <td>[-0.102161065707172, -0.0940412510144283, -0.7...</td>
      <td>-0.175051</td>
      <td>0.162116</td>
      <td>0.372372</td>
      <td>0.515364</td>
      <td>0.144420</td>
      <td>-0.362806</td>
      <td>0.109403</td>
      <td>[-0.175050506937197, 0.162116360694768, 0.3723...</td>
      <td>0.555593</td>
      <td>0.546581</td>
      <td>0.324408</td>
      <td>-0.279168</td>
      <td>0.560360</td>
      <td>-0.298079</td>
      <td>0.234949</td>
      <td>[0.555593475471297, 0.546580604053197, 0.32440...</td>
      <td>0.024702</td>
      <td>0.028622</td>
      <td>0.330168</td>
      <td>-0.145092</td>
      <td>-0.052249</td>
      <td>0.390849</td>
      <td>0.096167</td>
      <td>[0.02470210181234, 0.0286224547429464, 0.33016...</td>
      <td>0.071004</td>
      <td>0.221824</td>
      <td>-0.069967</td>
      <td>-0.548883</td>
      <td>-0.525773</td>
      <td>-0.478993</td>
      <td>-0.221798</td>
      <td>[0.0710043326473108, 0.221823916630173, -0.069...</td>
      <td>0.000258</td>
      <td>-0.596759</td>
      <td>0.309931</td>
      <td>-0.372023</td>
      <td>-0.033602</td>
      <td>0.083418</td>
      <td>-0.101463</td>
      <td>[0.0002582966608133, -0.596759018177968, 0.309...</td>
      <td>-0.138867</td>
      <td>0.118531</td>
      <td>-0.449505</td>
      <td>-0.071069</td>
      <td>-0.414840</td>
      <td>0.317762</td>
      <td>-0.106331</td>
      <td>[-0.138866804624008, 0.118531428276318, -0.449...</td>
      <td>0.909088</td>
      <td>-0.175535</td>
      <td>0.780910</td>
      <td>-0.269854</td>
      <td>-0.224045</td>
      <td>-0.292559</td>
      <td>0.121334</td>
      <td>[0.90908809327848, -0.175535107021276, 0.78090...</td>
      <td>-0.205927</td>
      <td>-0.188143</td>
      <td>-0.105910</td>
      <td>-0.013501</td>
      <td>-0.098393</td>
      <td>0.137659</td>
      <td>-0.079036</td>
      <td>[-0.205927037792345, -0.188143473994638, -0.10...</td>
    </tr>
    <tr>
      <th>3</th>
      <td>VA1_AAACACCAATAACTGC-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>6761</td>
      <td>3762</td>
      <td>round1</td>
      <td>4936</td>
      <td>3641</td>
      <td>4</td>
      <td>4</td>
      <td>VA1</td>
      <td>0.907884</td>
      <td>vasc_low</td>
      <td>9</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>1</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_low_MHCI_low</td>
      <td>Tumor_low</td>
      <td>vasc_low</td>
      <td>SCLC+vasculature</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>4_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>-0.294997</td>
      <td>-0.159254</td>
      <td>-0.624650</td>
      <td>-0.709188</td>
      <td>-0.312688</td>
      <td>-0.128532</td>
      <td>0.547967</td>
      <td>-0.287383</td>
      <td>0.055645</td>
      <td>-0.309493</td>
      <td>-0.234070</td>
      <td>-0.035042</td>
      <td>-0.158471</td>
      <td>-0.261131</td>
      <td>-0.305161</td>
      <td>0.338145</td>
      <td>-0.330022</td>
      <td>-0.141812</td>
      <td>-0.190257</td>
      <td>-0.720670</td>
      <td>-0.219027</td>
      <td>-0.454912</td>
      <td>-0.372893</td>
      <td>-0.460346</td>
      <td>-0.234070</td>
      <td>0.077766</td>
      <td>Low</td>
      <td>Low</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>9524.0</td>
      <td>25898.0</td>
      <td>AAACACCAATAACTGC-1</td>
      <td>VA1_3</td>
      <td>VA1_1747</td>
      <td>VA1_425</td>
      <td>VA1_1286</td>
      <td>VA1_1385</td>
      <td>VA1_2332</td>
      <td>VA1_2567</td>
      <td>6</td>
      <td>[VA1_1747, VA1_425, VA1_1286, VA1_1385, VA1_23...</td>
      <td>-0.708070</td>
      <td>-0.172725</td>
      <td>-0.305942</td>
      <td>-0.710094</td>
      <td>-0.721270</td>
      <td>-0.723104</td>
      <td>-0.556867</td>
      <td>[-0.708069838814006, -0.172725483080612, -0.30...</td>
      <td>-0.233200</td>
      <td>0.078572</td>
      <td>0.475936</td>
      <td>-0.356947</td>
      <td>-0.373849</td>
      <td>-0.379486</td>
      <td>-0.131496</td>
      <td>[-0.23319977872425, 0.0785717073459754, 0.4759...</td>
      <td>-0.239341</td>
      <td>0.355265</td>
      <td>0.295795</td>
      <td>-0.281069</td>
      <td>0.314224</td>
      <td>-0.329398</td>
      <td>0.019246</td>
      <td>[-0.239340532575186, 0.355265387687536, 0.2957...</td>
      <td>-0.017243</td>
      <td>0.579839</td>
      <td>0.544887</td>
      <td>0.357208</td>
      <td>0.062353</td>
      <td>-0.258450</td>
      <td>0.211432</td>
      <td>[-0.0172430474396086, 0.579838912374754, 0.544...</td>
      <td>-0.691933</td>
      <td>-0.813996</td>
      <td>-0.805377</td>
      <td>-0.592536</td>
      <td>-0.680867</td>
      <td>-0.801487</td>
      <td>-0.731033</td>
      <td>[-0.69193292191416, -0.813995861950977, -0.805...</td>
      <td>-0.536209</td>
      <td>0.002456</td>
      <td>0.312658</td>
      <td>-0.152928</td>
      <td>-0.621705</td>
      <td>-0.625379</td>
      <td>-0.270185</td>
      <td>[-0.536209058824282, 0.0024558419296655, 0.312...</td>
      <td>-0.400199</td>
      <td>-0.277454</td>
      <td>0.170990</td>
      <td>-0.287088</td>
      <td>-0.455337</td>
      <td>-0.334244</td>
      <td>-0.263889</td>
      <td>[-0.400198776367015, -0.277454023705374, 0.170...</td>
      <td>-0.208142</td>
      <td>-0.257051</td>
      <td>-0.297359</td>
      <td>-0.273455</td>
      <td>-0.288158</td>
      <td>-0.324265</td>
      <td>-0.274738</td>
      <td>[-0.20814162832553, -0.25705141028191, -0.2973...</td>
      <td>-0.246047</td>
      <td>-0.084913</td>
      <td>-0.007454</td>
      <td>0.028760</td>
      <td>-0.320470</td>
      <td>-0.319042</td>
      <td>-0.158194</td>
      <td>[-0.246047435280206, -0.084913354050493, -0.00...</td>
    </tr>
    <tr>
      <th>4</th>
      <td>VA1_AAACAGCTTTCAGAAG-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>5132</td>
      <td>3219</td>
      <td>round1</td>
      <td>4797</td>
      <td>3218</td>
      <td>4</td>
      <td>4</td>
      <td>VA1</td>
      <td>3.751279</td>
      <td>vasc_high</td>
      <td>14</td>
      <td>Tumor_low</td>
      <td>Tumor_low</td>
      <td>4</td>
      <td>MHCI_high</td>
      <td>MHCI_high</td>
      <td>Tumor_low_MHCI_high</td>
      <td>Tumor_low</td>
      <td>vasc_high</td>
      <td>SCLC+vasculature</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>4_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>0.148657</td>
      <td>0.670108</td>
      <td>0.181029</td>
      <td>-0.682642</td>
      <td>0.511207</td>
      <td>0.351200</td>
      <td>0.298867</td>
      <td>0.536963</td>
      <td>-0.222596</td>
      <td>-0.272365</td>
      <td>-0.244035</td>
      <td>-0.024744</td>
      <td>-0.099322</td>
      <td>0.541333</td>
      <td>0.725377</td>
      <td>-0.416435</td>
      <td>0.129047</td>
      <td>0.159843</td>
      <td>0.433786</td>
      <td>0.630433</td>
      <td>0.598608</td>
      <td>-0.097602</td>
      <td>0.871973</td>
      <td>0.918392</td>
      <td>-0.244035</td>
      <td>-0.157300</td>
      <td>High</td>
      <td>Above80Percentile</td>
      <td>vasc_high_MN4_high</td>
      <td>vasc_high_MN4_upper80</td>
      <td>15282.0</td>
      <td>27943.0</td>
      <td>AAACAGCTTTCAGAAG-1</td>
      <td>VA1_4</td>
      <td>VA1_964</td>
      <td>VA1_1350</td>
      <td>VA1_2187</td>
      <td>VA1_320</td>
      <td>VA1_885</td>
      <td>VA1_1441</td>
      <td>6</td>
      <td>[VA1_964, VA1_1350, VA1_2187, VA1_320, VA1_885...</td>
      <td>-0.274701</td>
      <td>-0.359858</td>
      <td>-0.726771</td>
      <td>-0.261499</td>
      <td>-0.345378</td>
      <td>-0.100794</td>
      <td>-0.344833</td>
      <td>[-0.274700715316639, -0.35985792380019, -0.726...</td>
      <td>0.646957</td>
      <td>-0.377095</td>
      <td>-0.387680</td>
      <td>-0.367993</td>
      <td>-0.383631</td>
      <td>-0.173307</td>
      <td>-0.173791</td>
      <td>[0.646957089052535, -0.377094771708025, -0.387...</td>
      <td>0.739864</td>
      <td>0.286462</td>
      <td>-0.346208</td>
      <td>0.305705</td>
      <td>-0.337803</td>
      <td>0.555446</td>
      <td>0.200578</td>
      <td>[0.739864288138109, 0.286461826481311, -0.3462...</td>
      <td>0.568322</td>
      <td>-0.317762</td>
      <td>-0.291268</td>
      <td>0.557903</td>
      <td>-0.275162</td>
      <td>0.514957</td>
      <td>0.126165</td>
      <td>[0.568321679326659, -0.317762474933888, -0.291...</td>
      <td>-0.797920</td>
      <td>-0.702421</td>
      <td>-0.795172</td>
      <td>-0.494896</td>
      <td>-0.000282</td>
      <td>-0.156359</td>
      <td>-0.491175</td>
      <td>[-0.79792035031997, -0.702421272421103, -0.795...</td>
      <td>-0.509378</td>
      <td>-0.158298</td>
      <td>-0.631615</td>
      <td>0.537557</td>
      <td>-0.107571</td>
      <td>0.003694</td>
      <td>-0.144269</td>
      <td>[-0.50937767833692, -0.158298382085461, -0.631...</td>
      <td>-0.308972</td>
      <td>-0.134719</td>
      <td>-0.168973</td>
      <td>-0.050159</td>
      <td>-0.339376</td>
      <td>-0.145247</td>
      <td>-0.191241</td>
      <td>[-0.308972369596385, -0.134719262393709, -0.16...</td>
      <td>0.743919</td>
      <td>-0.346469</td>
      <td>-0.337367</td>
      <td>-0.289158</td>
      <td>0.746207</td>
      <td>-0.178436</td>
      <td>0.056449</td>
      <td>[0.743919132154805, -0.346469293858606, -0.337...</td>
      <td>-0.016872</td>
      <td>-0.129730</td>
      <td>-0.239732</td>
      <td>-0.094755</td>
      <td>-0.166601</td>
      <td>0.028572</td>
      <td>-0.103186</td>
      <td>[-0.0168724150989142, -0.129729812963318, -0.2...</td>
    </tr>
    <tr>
      <th>5</th>
      <td>VA1_AAACAGGGTCTATATT-1</td>
      <td>VA1</td>
      <td>VA1</td>
      <td>8196</td>
      <td>4380</td>
      <td>round1</td>
      <td>5147</td>
      <td>3829</td>
      <td>4</td>
      <td>4</td>
      <td>VA1</td>
      <td>0.000000</td>
      <td>vasc_low</td>
      <td>26</td>
      <td>Tumor_high</td>
      <td>Tumor_high</td>
      <td>3</td>
      <td>MHCI_low</td>
      <td>MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>Tumor_high_MHCI_low</td>
      <td>vasc_neg</td>
      <td>SCLC+vasculature</td>
      <td>SCLC</td>
      <td>SCLC</td>
      <td>4_SCLC</td>
      <td>1_4_5_11_SCLC</td>
      <td>-0.323709</td>
      <td>-0.471163</td>
      <td>-0.687675</td>
      <td>-0.694835</td>
      <td>-0.371723</td>
      <td>-0.246882</td>
      <td>-0.263471</td>
      <td>-0.369183</td>
      <td>-0.089673</td>
      <td>-0.370011</td>
      <td>-0.041501</td>
      <td>-0.092670</td>
      <td>-0.086812</td>
      <td>-0.322861</td>
      <td>-0.365973</td>
      <td>-0.506152</td>
      <td>0.126810</td>
      <td>-0.196939</td>
      <td>-0.169592</td>
      <td>-0.377017</td>
      <td>0.187969</td>
      <td>0.193422</td>
      <td>-0.372631</td>
      <td>-0.515952</td>
      <td>-0.041501</td>
      <td>-0.255055</td>
      <td>Low</td>
      <td>Low</td>
      <td>vasc_low</td>
      <td>vasc_low</td>
      <td>13841.0</td>
      <td>27122.0</td>
      <td>AAACAGGGTCTATATT-1</td>
      <td>VA1_5</td>
      <td>VA1_1045</td>
      <td>VA1_976</td>
      <td>VA1_1053</td>
      <td>VA1_2706</td>
      <td>VA1_1504</td>
      <td>VA1_190</td>
      <td>6</td>
      <td>[VA1_1045, VA1_976, VA1_1053, VA1_2706, VA1_15...</td>
      <td>-0.726844</td>
      <td>-0.714224</td>
      <td>-0.290575</td>
      <td>0.435964</td>
      <td>-0.728868</td>
      <td>-0.337161</td>
      <td>-0.393618</td>
      <td>[-0.726843721645307, -0.71422440541496, -0.290...</td>
      <td>-0.384312</td>
      <td>-0.357587</td>
      <td>-0.365868</td>
      <td>-0.071121</td>
      <td>-0.065675</td>
      <td>-0.002911</td>
      <td>-0.207912</td>
      <td>[-0.384311753424909, -0.357587207854525, -0.36...</td>
      <td>-0.341305</td>
      <td>-0.276666</td>
      <td>-0.286172</td>
      <td>0.242455</td>
      <td>-0.349810</td>
      <td>-0.322093</td>
      <td>-0.222265</td>
      <td>[-0.341304782869537, -0.276665999599592, -0.28...</td>
      <td>-0.284224</td>
      <td>-0.134279</td>
      <td>-0.158431</td>
      <td>0.248524</td>
      <td>-0.306371</td>
      <td>0.283083</td>
      <td>-0.058616</td>
      <td>[-0.284223962368656, -0.134278846170489, -0.15...</td>
      <td>-0.798838</td>
      <td>-0.815859</td>
      <td>-0.725968</td>
      <td>-0.608392</td>
      <td>-0.704284</td>
      <td>-0.568430</td>
      <td>-0.703629</td>
      <td>[-0.798838359330548, -0.815859282625153, -0.72...</td>
      <td>-0.629781</td>
      <td>-0.617194</td>
      <td>-0.619058</td>
      <td>0.491038</td>
      <td>-0.633611</td>
      <td>-0.161870</td>
      <td>-0.361746</td>
      <td>[-0.629781360981336, -0.617194119228131, -0.61...</td>
      <td>-0.346066</td>
      <td>-0.027284</td>
      <td>-0.065629</td>
      <td>-0.488348</td>
      <td>-0.327220</td>
      <td>0.128919</td>
      <td>-0.187605</td>
      <td>[-0.346066185525329, -0.0272838762524825, -0.0...</td>
      <td>-0.336067</td>
      <td>-0.269654</td>
      <td>-0.278356</td>
      <td>0.598862</td>
      <td>-0.342368</td>
      <td>-0.314463</td>
      <td>-0.157008</td>
      <td>[-0.336067213442525, -0.269653930786008, -0.27...</td>
      <td>-0.342761</td>
      <td>-0.202520</td>
      <td>-0.171931</td>
      <td>-0.023035</td>
      <td>-0.120509</td>
      <td>0.135843</td>
      <td>-0.120819</td>
      <td>[-0.342761482095657, -0.202520214305712, -0.17...</td>
    </tr>
  </tbody>
</table>
</div>


# Save neighborhoods


```python
test.to_csv('Processing/10_Neighborhood_Analysis/AB1_A2_Neighborhoods.csv')
```


```python
# subset to only neighborhoods with all 6 neighbors adjacent
test_sub_all_adjacent = test[test['num_adjacent_nns']==6]
print(test.shape)
print(test_sub_all_adjacent.shape)
```

    (9711, 141)
    (8522, 141)



```python
# save only all adjacent neighborhoods
test_sub_all_adjacent.to_csv('Processing/10_Neighborhood_Analysis/AB1_A2_Neighborhoods_All_Six_Adjacent.csv')
```


```python
# subset to only certain columns

```
