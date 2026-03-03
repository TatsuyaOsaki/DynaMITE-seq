# Try GSVA - save gmt and expression data and use old version of GSVA via jason's dockerhub. 


```R
R.home()
```


'/usr/local/lib/R'

installed.packages()[, c("Package", "LibPath")]

```R
.libPaths()
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'/usr/local/lib/R/site-library'</li><li>'/usr/local/lib/R/library'</li></ol>




```R
library(Seurat)
library(Matrix)
library(doParallel)
library(ggplot2)
library(org.Hs.eg.db)
library(msigdbr)
library(tidyverse)
library(dittoSeq)
library(pheatmap)
library(GSVA)
library(BaseSet)
library(readxl)
library(corrplot)
library(ggbeeswarm)
library(jsonlite)
```

    Loading required package: AnnotationDbi
    
    
    
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
    [31m✖[39m [34mpurrr[39m::[32mflatten_df()[39m      masks [34mhdf5r[39m::flatten_df()
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
    
    Attaching package: ‘BaseSet’
    
    
    The following object is masked from ‘package:lubridate’:
    
        union
    
    
    The following object is masked from ‘package:dplyr’:
    
        union
    
    
    The following object is masked from ‘package:tibble’:
    
        add_column
    
    
    The following object is masked from ‘package:AnnotationDbi’:
    
        select
    
    
    The following objects are masked from ‘package:IRanges’:
    
        active, union
    
    
    The following objects are masked from ‘package:S4Vectors’:
    
        active, union
    
    
    The following object is masked from ‘package:BiocGenerics’:
    
        union
    
    
    The following object is masked from ‘package:generics’:
    
        union
    
    
    The following object is masked from ‘package:stats’:
    
        filter
    
    
    The following object is masked from ‘package:base’:
    
        union
    
    
    corrplot 0.95 loaded
    
    
    Attaching package: ‘jsonlite’
    
    
    The following object is masked from ‘package:purrr’:
    
        flatten
    
    



```R
install.packages("hdf5r")
```

    Installing package into ‘/usr/local/lib/R/site-library’
    (as ‘lib’ is unspecified)
    



```R
library(hdf5r)
```

    
    Attaching package: ‘hdf5r’
    
    
    The following object is masked from ‘package:S4Vectors’:
    
        values
    
    

sessionInfo()

```R
# show more columns
options(repr.matrix.max.cols=100, repr.matrix.max.rows=100)
```


```R
# increase memory limit 
options(future.globals.maxSize = 1e9)
```

# load data
VA2 <- Load10X_Spatial("/immuno/ian/Projects/2024_09_Barbie_Visium_Round2/spaceranger_count/VisA_count/outs/")
VA1 <- Load10X_Spatial("/immuno/ian/Projects/2023_02_SCLC_10X_Visium/VA_count/outs/")
VB1 <- Load10X_Spatial("/immuno/ian/Projects/2023_02_SCLC_10X_Visium/VB_count/outs/")
VA2$orig.ident <- "VA2"
VA1$orig.ident <- "VA1"
VB1$orig.ident <- "VB1"
# set visium round, will use this for integration / batch correction
VA2$visium_round <- "round2"
VA1$visium_round <- "round1"
VB1$visium_round <- "round1"# Filter out spots that have zero UMIs, otherwise SCTransform won't run
VA2 <- subset(VA2, nCount_Spatial>0)
VA1 <- subset(VA1, nCount_Spatial>0)
VB1 <- subset(VB1, nCount_Spatial>0)#SCTransform normalization
VA2 <- SCTransform(VA2, assay = "Spatial", verbose = FALSE)
VA1 <- SCTransform(VA1, assay = "Spatial", verbose = FALSE)
VB1 <- SCTransform(VB1, assay = "Spatial", verbose = FALSE)VA2VA1VB1# Normalize the spatial assay
DefaultAssay(VA2) <- "Spatial"
DefaultAssay(VA1) <- "Spatial"
DefaultAssay(VB1) <- "Spatial"

VA2 = NormalizeData(VA2)
VA1 = NormalizeData(VA1)
VB1 = NormalizeData(VB1)# Get variable features
VA2 = FindVariableFeatures(VA2, assay='Spatial', nfeatures = 10000)
VA1 = FindVariableFeatures(VA1, assay='Spatial', nfeatures = 10000)
VB1 = FindVariableFeatures(VB1, assay='Spatial', nfeatures = 10000)#Scale data
VA2 <- ScaleData(VA2)
VA1 <- ScaleData(VA1)
VB1 <- ScaleData(VB1)VA2VA1VB1

```R
options(future.globals.maxSize = 3e+09)
```

# save these normalized datasets in case I need them again, won't have to redo SCTransform
saveRDS(VA2, 'Processing/norm_fromscratch_for_GSVA_20250625/VA2.rds')
saveRDS(VA1, 'Processing/norm_fromscratch_for_GSVA_20250625/VA1.rds')
saveRDS(VB1, 'Processing/norm_fromscratch_for_GSVA_20250625/VB1.rds')
# Load if needed


```R
VA2 <- readRDS('Processing/norm_fromscratch_for_GSVA_20250625/VA2.rds')
VA1 <- readRDS('Processing/norm_fromscratch_for_GSVA_20250625/VA1.rds')
VB1 <- readRDS('Processing/norm_fromscratch_for_GSVA_20250625/VB1.rds')
```

# prepare gene sets for GSVA


```R
# Gene lists added August 15, 2024
migration <- read.csv('/immuno/ian/Projects/2023_02_SCLC_10X_Visium/Seurat_Analysis/Processing_2024/GOBP_LEUKOCYTE_MIGRATION.csv')
reg_of_migration <- read.csv('/immuno/ian/Projects/2023_02_SCLC_10X_Visium/Seurat_Analysis/Processing_2024/GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION.csv')
migration %>% head(3)
reg_of_migration %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 1</caption>
<thead>
	<tr><th></th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>ABL1</td></tr>
	<tr><th scope=row>2</th><td>ABL2</td></tr>
	<tr><th scope=row>3</th><td>ADA </td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 3 × 1</caption>
<thead>
	<tr><th></th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>ABL1  </td></tr>
	<tr><th scope=row>2</th><td>ABL2  </td></tr>
	<tr><th scope=row>3</th><td>ADAM10</td></tr>
</tbody>
</table>




```R
# get them as named lists...
leuk_migration <- list('GOBP_LEUKOCYTE_MIGRATION' = migration$GOBP_LEUKOCYTE_MIGRATION,
                       'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION' = reg_of_migration$GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION)
```


```R
path <- "/immuno/ian/Projects/2023_02_SCLC_10X_Visium/Seurat_Analysis/Processing_2024/20240425_EC_Activation_Gene_Sets_lf.xlsx"
#geneset_df <- path %>% excel_sheets() %>% set_names() %>% map_dfr(read_excel, path = path, .id = "GeneSet")
geneset_df <- read_excel(path = path)
```


```R
geneset_df %>% head(2)
```


<table class="dataframe">
<caption>A tibble: 2 × 2</caption>
<thead>
	<tr><th scope=col>Gene</th><th scope=col>GeneSet</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><td>CCL5 </td><td>GOBP_LEUKOCYTE_ADHESION_TO_VASC</td></tr>
	<tr><td>FASLG</td><td>GOBP_LEUKOCYTE_ADHESION_TO_VASC</td></tr>
</tbody>
</table>




```R
# rename gene column to Symbol
geneset_df <- geneset_df %>% rename("Symbol" = "Gene")
# Split gene sets into lists of gene sets
geneset_df_list = split(x = geneset_df$Symbol, f = geneset_df$GeneSet)
```


```R
geneset_df_list$MN4_EC_Phenotype_Top30 <- geneset_df_list$MN4_EC_Phenotype[1:30]
```


```R
length(geneset_df_list$MN4_EC_Phenotype_Top30)
```


30



```R
# rename
navin_endo_geneset_lists <- geneset_df_list
```


```R
navin_endo_geneset_lists
```


<dl>
	<dt>$GOBP_LEUKOCYTE_ADHESION_TO_VASC</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CCL5'</li><li>'FASLG'</li><li>'SELE'</li><li>'TFP1'</li><li>'CX3CL1'</li><li>'B2M'</li><li>'CASP1'</li><li>'MMP9'</li><li>'TEK'</li><li>'VCAM1'</li><li>'TNFSF10'</li><li>'ICAM1'</li><li>'PTGIS'</li><li>'IL7'</li><li>'FGF2'</li><li>'TYMP'</li><li>'VWF'</li><li>'TGFB1'</li><li>'SOD1'</li><li>'SERPINE1'</li><li>'TJP1'</li><li>'OCLN'</li><li>'CLDN5'</li></ol>
</dd>
	<dt>$MN4_EC_Phenotype</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CXCL11'</li><li>'CXCL10'</li><li>'IFI44L'</li><li>'GBP4'</li><li>'IFI27'</li><li>'PECAM1'</li><li>'CD93'</li><li>'ECSCR'</li><li>'LGALS9'</li><li>'GBP2'</li><li>'MCAM'</li><li>'CLEC14A'</li><li>'MX1'</li><li>'CTSS'</li><li>'ADGRL4'</li><li>'MMRN1'</li><li>'HAPLN3'</li><li>'CCL2'</li><li>'HSPG2'</li><li>'PTGS2'</li><li>'THBS1'</li><li>'BST2'</li><li>'CDH5'</li><li>'CXCL9'</li><li>'EGFL7'</li><li>'PODXL'</li><li>'IFI6'</li><li>'IFIT3'</li><li>'IDO1'</li><li>'PRXL2A'</li><li>'GBP1'</li><li>'GBP5'</li><li>'MGP'</li><li>'EFEMP1'</li><li>'APOL1'</li><li>'TAP1'</li><li>'ISG20'</li><li>'SELE'</li><li>'ARHGAP29'</li><li>'FLT1'</li><li>'IL18BP'</li><li>'GIMAP4'</li><li>'IL32'</li><li>'PARP14'</li><li>'ENG'</li><li>'TGM2'</li><li>'TNFSF10'</li><li>'CALCRL'</li><li>'VWF'</li><li>'GOLM1'</li><li>'RNF213'</li><li>'SAMHD1'</li><li>'DDX58'</li><li>'BMPR2'</li><li>'ECE1'</li><li>'PLAAT4'</li><li>'IFIT1'</li><li>'PGF'</li><li>'PROCR'</li><li>'ICAM2'</li><li>'ICAM1'</li><li>'WARS'</li><li>'EDN1'</li><li>'RSAD2'</li><li>'WWTR1'</li><li>'OAS3'</li><li>'PRSS23'</li><li>'TRIM22'</li><li>'IRF1'</li><li>'LAMA4'</li><li>'PDGFB'</li><li>'LAP3'</li><li>'NCOA7'</li><li>'COL4A1'</li><li>'SERPINE1'</li><li>'EPSTI1'</li><li>'IL6ST'</li><li>'TCF4'</li><li>'XAF1'</li><li>'IFI35'</li><li>'CLDN5'</li><li>'PSMB9'</li><li>'IGFBP7'</li><li>'LGMN'</li><li>'IFIH1'</li><li>'USP18'</li><li>'TM4SF1'</li><li>'APP'</li><li>'RHOJ'</li><li>'PLPP3'</li><li>'VAMP5'</li><li>'APOL2'</li><li>'PNP'</li><li>'HHEX'</li><li>'OAS2'</li><li>'NNMT'</li><li>'UBE2L6'</li><li>'SAMD9L'</li><li>'SPARC'</li><li>'HYAL2'</li><li>'CD9'</li><li>'CRIM1'</li><li>'GNG11'</li><li>'SLC38A2'</li><li>'PXDN'</li><li>'SOD2'</li><li>'SPOCK1'</li><li>'FRMD4B'</li><li>'GALNT1'</li><li>'XIST'</li><li>'UTRN'</li><li>'ZC3HAV1'</li><li>'OPTN'</li><li>'TFPI2'</li><li>'ABL2'</li><li>'ITGA2'</li><li>'SOX18'</li><li>'APOL6'</li><li>'RAI14'</li><li>'COL18A1'</li><li>'TAPBP'</li><li>'MYCT1'</li><li>'STAT1'</li><li>'IFIT2'</li><li>'BMX'</li><li>'FN1'</li><li>'PSMB8'</li><li>'S100A16'</li><li>'IFI44'</li><li>'CAV1'</li><li>'ITGA5'</li><li>'PALMD'</li><li>'TGFBR2'</li><li>'MDK'</li><li>'SP100'</li><li>'OAS1'</li><li>'TXNDC5'</li><li>'DAAM1'</li><li>'ITGAV'</li><li>'CFLAR'</li><li>'COL4A2'</li><li>'ANGPT2'</li><li>'S100A13'</li><li>'PLS3'</li><li>'APLP2'</li><li>'JAG1'</li><li>'TMEM255B'</li><li>'HSPB1'</li><li>'PMAIP1'</li><li>'CAVIN1'</li><li>'SNTB2'</li><li>'NR2F2'</li><li>'DUSP6'</li><li>'NFIB'</li><li>'DDX60L'</li><li>'AKAP12'</li><li>'ROBO4'</li><li>'RGS5'</li><li>'ESAM'</li><li>'NRP2'</li><li>'ACTN1'</li><li>'ID3'</li><li>'CFH'</li><li>'HHIP'</li><li>'MACF1'</li><li>'CALD1'</li><li>'SPHK1'</li><li>'MMP2'</li><li>'RALA'</li><li>'MTUS1'</li><li>'LIPG'</li><li>'APOL3'</li><li>'ZNFX1'</li><li>'PRCP'</li><li>'TGFB1'</li><li>'CKAP4'</li><li>'NPDC1'</li><li>'SAT1'</li><li>'TAP2'</li><li>'OCIAD2'</li><li>'TFPI'</li><li>'MCUB'</li><li>'PDLIM5'</li><li>'FBLIM1'</li><li>'RASGRP3'</li><li>'ISG15'</li><li>'S1PR1'</li><li>'GSDMD'</li><li>'NMI'</li><li>'SOX4'</li><li>'KDR'</li><li>'TIE1'</li><li>'YBX3'</li><li>'LYL1'</li><li>'BACE2'</li><li>'SCARB2'</li><li>'CRIP2'</li><li>'GJA1'</li><li>'PLSCR1'</li><li>'MARCKS'</li><li>'RHOB'</li><li>'FNDC3B'</li><li>'IFITM1'</li><li>'LUZP1'</li><li>'EPAS1'</li><li>'IFITM3'</li><li>'IFI16'</li><li>'LAMB1'</li><li>'TNFSF13B'</li><li>'PDE4B'</li><li>'PALM2-AKAP2'</li><li>'PHACTR2'</li><li>'ITGB1'</li><li>'RAB13'</li><li>'FCGRT'</li><li>'SRPX'</li><li>'EHD4'</li><li>'PSMB10'</li><li>'TCIM'</li><li>'MYO10'</li><li>'ACVRL1'</li><li>'PKIG'</li><li>'CLDN14'</li><li>'FLNB'</li><li>'TUBB6'</li><li>'STAT2'</li><li>'ARHGAP21'</li><li>'NRCAM'</li><li>'SERPINH1'</li><li>'CD59'</li><li>'STOM'</li><li>'JAK1'</li><li>'PDGFA'</li><li>'CCND1'</li><li>'PML'</li><li>'GIMAP6'</li><li>'GRN'</li><li>'GIMAP7'</li><li>'KRT18'</li><li>'NECTIN2'</li><li>'PPFIBP1'</li><li>'GRASP'</li><li>'MX2'</li><li>'HLA-E'</li><li>'CTSB'</li><li>'SAMD9'</li><li>'MAST4'</li><li>'DDX60'</li><li>'RDX'</li><li>'GNB4'</li><li>'FZD4'</li><li>'RCN1'</li><li>'CTHRC1'</li><li>'SMAD6'</li><li>'NAMPT'</li><li>'TINAGL1'</li><li>'CCPG1'</li><li>'CCL14'</li><li>'WDFY1'</li><li>'HELZ2'</li><li>'VCAM1'</li><li>'VGLL4'</li><li>'LY6E'</li><li>'ADAR'</li><li>'PSME2'</li><li>'IGF2BP3'</li><li>'FKBP1A'</li><li>'BMP6'</li><li>'PCDH7'</li><li>'GNAI2'</li><li>'SPOCD1'</li><li>'CMPK2'</li><li>'SNCA'</li><li>'ODF2L'</li><li>'F2R'</li><li>'SPAG9'</li><li>'CBR3'</li><li>'EIF2AK2'</li><li>'IL1RL1'</li><li>'WSB1'</li><li>'FABP5'</li><li>'PCAT19'</li><li>'PTPRE'</li><li>'SASH1'</li><li>'SMURF2'</li><li>'CTNNA1'</li><li>'ESM1'</li><li>'MARCKSL1'</li><li>'FAM241A'</li><li>'RBMS1'</li><li>'DYRK4'</li><li>'PLOD2'</li><li>'GSTK1'</li><li>'SLC7A11'</li><li>'MEF2A'</li><li>'MARCH3'</li><li>'CLIC4'</li><li>'RTN4'</li><li>'MSN'</li><li>'NFE2L3'</li><li>'DIPK1B'</li><li>'PARP9'</li><li>'SP110'</li><li>'CARD16'</li><li>'PTTG1IP'</li><li>'CD99'</li><li>'CST3'</li><li>'RHOC'</li></ol>
</dd>
	<dt>$Upregulated_by_2_3_CGAMP</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CCL5'</li><li>'FASLG'</li><li>'SELE'</li><li>'TFP1'</li><li>'CX3CL1'</li><li>'B2M'</li><li>'CASP1'</li><li>'MMP9'</li><li>'TEK'</li><li>'VCAM1'</li><li>'TNFSF10'</li><li>'ICAM1'</li><li>'PTGIS'</li><li>'IL7'</li><li>'FGF2'</li><li>'TYMP'</li><li>'VWF'</li><li>'TGFB1'</li><li>'SOD1'</li><li>'SERPINE1'</li><li>'TJP1'</li><li>'OCLN'</li><li>'CLDN5'</li></ol>
</dd>
	<dt>$Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CCL5'</li><li>'TEK'</li><li>'VCAM1'</li><li>'SELE'</li><li>'IL7'</li><li>'FGF2'</li><li>'ICAM1'</li><li>'CXC3L1'</li><li>'TYMP'</li><li>'TJP1'</li><li>'OCLN'</li><li>'CLDN5'</li><li>'B2M'</li></ol>
</dd>
	<dt>$MN4_EC_Phenotype_Top30</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CXCL11'</li><li>'CXCL10'</li><li>'IFI44L'</li><li>'GBP4'</li><li>'IFI27'</li><li>'PECAM1'</li><li>'CD93'</li><li>'ECSCR'</li><li>'LGALS9'</li><li>'GBP2'</li><li>'MCAM'</li><li>'CLEC14A'</li><li>'MX1'</li><li>'CTSS'</li><li>'ADGRL4'</li><li>'MMRN1'</li><li>'HAPLN3'</li><li>'CCL2'</li><li>'HSPG2'</li><li>'PTGS2'</li><li>'THBS1'</li><li>'BST2'</li><li>'CDH5'</li><li>'CXCL9'</li><li>'EGFL7'</li><li>'PODXL'</li><li>'IFI6'</li><li>'IFIT3'</li><li>'IDO1'</li><li>'PRXL2A'</li></ol>
</dd>
</dl>



# let's do some other more updated genesets as well
# vasculature
vasculature <- c('PECAM1', 'VWF', 'CDH5', 'KDR', 'CLDN5', 'PLVAP')

# Ayers IFNG
AyersIFNG <- c('IDO1', 'CXCL10', 'CXCL9', 'HLA-DRA', 'STAT1', 'IFNG')

# Adhesion
Adhesion <- c('ICAM1', 'VCAM1', 'SELE')

vasc_list <- list("Vasculature"=vasculature, "AyersIFNG"=AyersIFNG, "Adhesion"=Adhesion)

```R
# Macrophages
Macrophage <- c("CD68", "CSF1R", "CD14", "ITGAM", "MARCO", "FCGR1A")
M1_Macrophage <- c("CD80", "CD86", "NOS2", "IRF5", "TNF", "IL1B", "IL6", "IL12A", "CXCL9", "CXCL10")
M2_Macrophage <- c("MRC1", "CD163", "MSR1", "ARG1", "IL10", "TGFBI", "RETNLA", "CHI3L3")
M1_Macrophage2 <- c("CSF1R", "CD14", "ITGAM", "CD80", "CD86", "NOS2", "IRF5", "TNF", "IL1B", "IL6", "IL12A", "CXCL9", "CXCL10")
M2_Macrophage2 <- c("CSF1R", "CD14", "ITGAM", "MRC1", "CD163", "MSR1", "ARG1", "IL10", "TGFBI", "RETNLA", "CHI3L3")
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

# Endothelial Cells
Endothelial_Cell <- c('CD34', 'MCAM', 'PECAM1', 'VWF', 'CDH5', 'KDR', 'VCAM1', 'FLT1', 'CLDN5', 'RAMP2', 'TEK', 'PLVAP', 'CLEC14A', 'ERG', 'ENG' ,'ACKR1', 'EGFL7')
LEC <- c('PROX1', 'LYVE1', 'PDPN', 'FLT4', 'PODXL')
Endothelial_Activation <- c('ICAM1', 'VCAM1', 'SELE', 'CXCL8', 'CCL2')
Endothelial_Chemokines <- c('CCL2', 'CCL3', 'CCL5', 'CCL19', 'CXCL8', 'CXCL9', 'CXCL10', 'CXCL11', 'CXCL12')
Angiogenesis <- c('VEGFA','KDR','ANGPT1','ANGPT2','TEK','FGF2','PDGFB','HIF1A','MMP2','MMP9','THBS1','PLGF','NOS3',
                  'CXCL12','DLL4','NOTCH1','EPAS1','SPHK1','TIMP1','TIMP2')

# Ayers IFNG
AyersIFNG <- c('IDO1', 'CXCL10', 'CXCL9', 'HLA-DRA', 'STAT1', 'IFNG')
# Exchaustion
Exhaustion <- c('LAG3', 'HAVCR2', 'TIGIT', 'TOX', 'EOMES', 'ENTPD1')

# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9578628/
# Table S3 gives genes correlated with proliferation...
# I'll take the top 14 plus the two popular ones MKI67 and TOP2A
Proliferation <- c('MKI67', 'TOP2A', 'MCM3', 'MCM2', 'MCM7', 'MCM4', 'MCM5', 'MCM6', 'RFC2', 'RFC4', 'SMC4', 'SMC2', 'CDK1', 'RFC3', 'KIF11', 'RRM1')

# ICB Targets
ICB_Targets <- c('CTLA4','PDCD1','CD274','LAG3','HAVCR2','TIGIT','TNFRSF4','TNFRSF18','IDO1','CD276','VTCN1','CD40','TNFRSF9')

# others
Epithelial_Cell <- c('EPCAM','KRT8','KRT18','KRT5','KRT19','KRT7','CDH1','CLDN4','TJP1','KRT14')
Fibroblast <- c('COL3A1','COL1A1','S100A4','ACTA2','COL1A2','DCN','THY1','FAP','PDGFRA','PDPN','LUM','TAGLN','VIM','PDGFRB','POSTN','FBLN1','FN1','DDR2','APOD')
```


```R
# make updated gene lists
updated_cell_type_lists <- list("All_Immune"=All_Immune, 
                                "Macrophage"=Macrophage, "M1_Macrophage"=M1_Macrophage2, "M2_Macrophage"=M2_Macrophage2,
                                "M1_Macrophage_old"=M1_Macrophage, "M2_Macrophage_old"=M2_Macrophage,
                                "T_Cell"=T_Cell, "CD8_T_Cell"=CD8_T_Cell, "T_Reg"=T_Reg,
                               "NK_Cell"=NK_Cell, "B_Cell"=B_Cell, "DC"=DC, "Epithelial_cell"=Epithelial_Cell, "Fibroblast"=Fibroblast)
Endothelial_lists <- list("Endothelial_Cell"=Endothelial_Cell, "LEC"=LEC, 
                          "Endothelial_Activation"=Endothelial_Activation, "Endothelial_Chemokines"=Endothelial_Chemokines,
                          "Angiogenesis"=Angiogenesis)

Activation_state_lists <- list("AyersIFNG"=AyersIFNG, "Exhaustion"=Exhaustion, "Proliferation"=Proliferation, "ICB_Targets"=ICB_Targets)
                                
updated_signature_lists <- c(updated_cell_type_lists, Endothelial_lists, navin_endo_geneset_lists, Activation_state_lists, leuk_migration)
```


```R
updated_signature_lists
```


<dl>
	<dt>$All_Immune</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'PTPRC'</li><li>'CD4'</li><li>'FOXP3'</li><li>'IFNG'</li><li>'CD68'</li><li>'CSF1R'</li><li>'CD14'</li><li>'ITGAM'</li><li>'MARCO'</li><li>'FCGR1A'</li><li>'CD80'</li><li>'CD86'</li><li>'NOS2'</li><li>'IRF5'</li><li>'TNF'</li><li>'IL1B'</li><li>'IL6'</li><li>'IL12A'</li><li>'CXCL9'</li><li>'CXCL10'</li><li>'MRC1'</li><li>'CD163'</li><li>'MSR1'</li><li>'ARG1'</li><li>'IL10'</li><li>'TGFBI'</li><li>'RETNLA'</li><li>'CHI3L3'</li><li>'FCGR1A'</li><li>'IL1B'</li><li>'AIF1'</li><li>'C1QC'</li><li>'C1QA'</li><li>'C1QB'</li><li>'CD3E'</li><li>'CD3D'</li><li>'CD3G'</li><li>'TRBC1'</li><li>'TRBC2'</li><li>'TRAC'</li><li>'CD28'</li><li>'CD247'</li><li>'CD5'</li><li>'CD27'</li><li>'CD8A'</li><li>'CD8B'</li><li>'GZMA'</li><li>'GZMB'</li><li>'GZMK'</li><li>'CD69'</li><li>'PRF1'</li><li>'NKG7'</li><li>'IL7R'</li><li>'CD3E'</li><li>'CD4'</li><li>'IL2RA'</li><li>'FOXP3'</li><li>'NCR1'</li><li>'FCGR3A'</li><li>'KLRK1'</li><li>'KLRD1'</li><li>'NCAM1'</li><li>'NKG7'</li><li>'GNLY'</li><li>'KLRC1'</li><li>'KLRF1'</li><li>'KLRB1'</li><li>'GZMB'</li><li>'CD160'</li><li>'XCL1'</li><li>'NCR3'</li><li>'KLRC2'</li><li>'XCL2'</li><li>'PRF1'</li><li>'KIR2DL3'</li><li>'MS4A1'</li><li>'CD19'</li><li>'CD79A'</li><li>'IGKC'</li><li>'IGHG3'</li><li>'MZB1'</li><li>'CR2'</li><li>'CD79B'</li><li>'BLNK'</li><li>'JCHAIN'</li><li>'IGHG2'</li><li>'IGHG1'</li><li>'CD22'</li><li>'BANK1'</li><li>'IGHD'</li><li>'IGHA1'</li><li>'IGHM'</li><li>'PAX5'</li><li>'BLK'</li><li>'LAMP3'</li><li>'ITGAX'</li><li>'ITGAM'</li><li>'HLA-DRA'</li><li>'CD80'</li><li>'CD86'</li><li>'BATF3'</li><li>'IRF4'</li><li>'IRF8'</li><li>'CD209'</li><li>'CEACAM8'</li><li>'ELANE'</li><li>'KIT'</li><li>'TPSAB1'</li><li>'CCR3'</li><li>'EPX'</li><li>'FCER1A'</li><li>'LYZ'</li></ol>
</dd>
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
<ol class=list-inline><li>'CSF1R'</li><li>'CD14'</li><li>'ITGAM'</li><li>'CD80'</li><li>'CD86'</li><li>'NOS2'</li><li>'IRF5'</li><li>'TNF'</li><li>'IL1B'</li><li>'IL6'</li><li>'IL12A'</li><li>'CXCL9'</li><li>'CXCL10'</li></ol>
</dd>
	<dt>$M2_Macrophage</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CSF1R'</li><li>'CD14'</li><li>'ITGAM'</li><li>'MRC1'</li><li>'CD163'</li><li>'MSR1'</li><li>'ARG1'</li><li>'IL10'</li><li>'TGFBI'</li><li>'RETNLA'</li><li>'CHI3L3'</li></ol>
</dd>
	<dt>$M1_Macrophage_old</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CD80'</li><li>'CD86'</li><li>'NOS2'</li><li>'IRF5'</li><li>'TNF'</li><li>'IL1B'</li><li>'IL6'</li><li>'IL12A'</li><li>'CXCL9'</li><li>'CXCL10'</li></ol>
</dd>
	<dt>$M2_Macrophage_old</dt>
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
	<dt>$Epithelial_cell</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'EPCAM'</li><li>'KRT8'</li><li>'KRT18'</li><li>'KRT5'</li><li>'KRT19'</li><li>'KRT7'</li><li>'CDH1'</li><li>'CLDN4'</li><li>'TJP1'</li><li>'KRT14'</li></ol>
</dd>
	<dt>$Fibroblast</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'COL3A1'</li><li>'COL1A1'</li><li>'S100A4'</li><li>'ACTA2'</li><li>'COL1A2'</li><li>'DCN'</li><li>'THY1'</li><li>'FAP'</li><li>'PDGFRA'</li><li>'PDPN'</li><li>'LUM'</li><li>'TAGLN'</li><li>'VIM'</li><li>'PDGFRB'</li><li>'POSTN'</li><li>'FBLN1'</li><li>'FN1'</li><li>'DDR2'</li><li>'APOD'</li></ol>
</dd>
	<dt>$Endothelial_Cell</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CD34'</li><li>'MCAM'</li><li>'PECAM1'</li><li>'VWF'</li><li>'CDH5'</li><li>'KDR'</li><li>'VCAM1'</li><li>'FLT1'</li><li>'CLDN5'</li><li>'RAMP2'</li><li>'TEK'</li><li>'PLVAP'</li><li>'CLEC14A'</li><li>'ERG'</li><li>'ENG'</li><li>'ACKR1'</li><li>'EGFL7'</li></ol>
</dd>
	<dt>$LEC</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'PROX1'</li><li>'LYVE1'</li><li>'PDPN'</li><li>'FLT4'</li><li>'PODXL'</li></ol>
</dd>
	<dt>$Endothelial_Activation</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'ICAM1'</li><li>'VCAM1'</li><li>'SELE'</li><li>'CXCL8'</li><li>'CCL2'</li></ol>
</dd>
	<dt>$Endothelial_Chemokines</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CCL2'</li><li>'CCL3'</li><li>'CCL5'</li><li>'CCL19'</li><li>'CXCL8'</li><li>'CXCL9'</li><li>'CXCL10'</li><li>'CXCL11'</li><li>'CXCL12'</li></ol>
</dd>
	<dt>$Angiogenesis</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'VEGFA'</li><li>'KDR'</li><li>'ANGPT1'</li><li>'ANGPT2'</li><li>'TEK'</li><li>'FGF2'</li><li>'PDGFB'</li><li>'HIF1A'</li><li>'MMP2'</li><li>'MMP9'</li><li>'THBS1'</li><li>'PLGF'</li><li>'NOS3'</li><li>'CXCL12'</li><li>'DLL4'</li><li>'NOTCH1'</li><li>'EPAS1'</li><li>'SPHK1'</li><li>'TIMP1'</li><li>'TIMP2'</li></ol>
</dd>
	<dt>$GOBP_LEUKOCYTE_ADHESION_TO_VASC</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CCL5'</li><li>'FASLG'</li><li>'SELE'</li><li>'TFP1'</li><li>'CX3CL1'</li><li>'B2M'</li><li>'CASP1'</li><li>'MMP9'</li><li>'TEK'</li><li>'VCAM1'</li><li>'TNFSF10'</li><li>'ICAM1'</li><li>'PTGIS'</li><li>'IL7'</li><li>'FGF2'</li><li>'TYMP'</li><li>'VWF'</li><li>'TGFB1'</li><li>'SOD1'</li><li>'SERPINE1'</li><li>'TJP1'</li><li>'OCLN'</li><li>'CLDN5'</li></ol>
</dd>
	<dt>$MN4_EC_Phenotype</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CXCL11'</li><li>'CXCL10'</li><li>'IFI44L'</li><li>'GBP4'</li><li>'IFI27'</li><li>'PECAM1'</li><li>'CD93'</li><li>'ECSCR'</li><li>'LGALS9'</li><li>'GBP2'</li><li>'MCAM'</li><li>'CLEC14A'</li><li>'MX1'</li><li>'CTSS'</li><li>'ADGRL4'</li><li>'MMRN1'</li><li>'HAPLN3'</li><li>'CCL2'</li><li>'HSPG2'</li><li>'PTGS2'</li><li>'THBS1'</li><li>'BST2'</li><li>'CDH5'</li><li>'CXCL9'</li><li>'EGFL7'</li><li>'PODXL'</li><li>'IFI6'</li><li>'IFIT3'</li><li>'IDO1'</li><li>'PRXL2A'</li><li>'GBP1'</li><li>'GBP5'</li><li>'MGP'</li><li>'EFEMP1'</li><li>'APOL1'</li><li>'TAP1'</li><li>'ISG20'</li><li>'SELE'</li><li>'ARHGAP29'</li><li>'FLT1'</li><li>'IL18BP'</li><li>'GIMAP4'</li><li>'IL32'</li><li>'PARP14'</li><li>'ENG'</li><li>'TGM2'</li><li>'TNFSF10'</li><li>'CALCRL'</li><li>'VWF'</li><li>'GOLM1'</li><li>'RNF213'</li><li>'SAMHD1'</li><li>'DDX58'</li><li>'BMPR2'</li><li>'ECE1'</li><li>'PLAAT4'</li><li>'IFIT1'</li><li>'PGF'</li><li>'PROCR'</li><li>'ICAM2'</li><li>'ICAM1'</li><li>'WARS'</li><li>'EDN1'</li><li>'RSAD2'</li><li>'WWTR1'</li><li>'OAS3'</li><li>'PRSS23'</li><li>'TRIM22'</li><li>'IRF1'</li><li>'LAMA4'</li><li>'PDGFB'</li><li>'LAP3'</li><li>'NCOA7'</li><li>'COL4A1'</li><li>'SERPINE1'</li><li>'EPSTI1'</li><li>'IL6ST'</li><li>'TCF4'</li><li>'XAF1'</li><li>'IFI35'</li><li>'CLDN5'</li><li>'PSMB9'</li><li>'IGFBP7'</li><li>'LGMN'</li><li>'IFIH1'</li><li>'USP18'</li><li>'TM4SF1'</li><li>'APP'</li><li>'RHOJ'</li><li>'PLPP3'</li><li>'VAMP5'</li><li>'APOL2'</li><li>'PNP'</li><li>'HHEX'</li><li>'OAS2'</li><li>'NNMT'</li><li>'UBE2L6'</li><li>'SAMD9L'</li><li>'SPARC'</li><li>'HYAL2'</li><li>'CD9'</li><li>'CRIM1'</li><li>'GNG11'</li><li>'SLC38A2'</li><li>'PXDN'</li><li>'SOD2'</li><li>'SPOCK1'</li><li>'FRMD4B'</li><li>'GALNT1'</li><li>'XIST'</li><li>'UTRN'</li><li>'ZC3HAV1'</li><li>'OPTN'</li><li>'TFPI2'</li><li>'ABL2'</li><li>'ITGA2'</li><li>'SOX18'</li><li>'APOL6'</li><li>'RAI14'</li><li>'COL18A1'</li><li>'TAPBP'</li><li>'MYCT1'</li><li>'STAT1'</li><li>'IFIT2'</li><li>'BMX'</li><li>'FN1'</li><li>'PSMB8'</li><li>'S100A16'</li><li>'IFI44'</li><li>'CAV1'</li><li>'ITGA5'</li><li>'PALMD'</li><li>'TGFBR2'</li><li>'MDK'</li><li>'SP100'</li><li>'OAS1'</li><li>'TXNDC5'</li><li>'DAAM1'</li><li>'ITGAV'</li><li>'CFLAR'</li><li>'COL4A2'</li><li>'ANGPT2'</li><li>'S100A13'</li><li>'PLS3'</li><li>'APLP2'</li><li>'JAG1'</li><li>'TMEM255B'</li><li>'HSPB1'</li><li>'PMAIP1'</li><li>'CAVIN1'</li><li>'SNTB2'</li><li>'NR2F2'</li><li>'DUSP6'</li><li>'NFIB'</li><li>'DDX60L'</li><li>'AKAP12'</li><li>'ROBO4'</li><li>'RGS5'</li><li>'ESAM'</li><li>'NRP2'</li><li>'ACTN1'</li><li>'ID3'</li><li>'CFH'</li><li>'HHIP'</li><li>'MACF1'</li><li>'CALD1'</li><li>'SPHK1'</li><li>'MMP2'</li><li>'RALA'</li><li>'MTUS1'</li><li>'LIPG'</li><li>'APOL3'</li><li>'ZNFX1'</li><li>'PRCP'</li><li>'TGFB1'</li><li>'CKAP4'</li><li>'NPDC1'</li><li>'SAT1'</li><li>'TAP2'</li><li>'OCIAD2'</li><li>'TFPI'</li><li>'MCUB'</li><li>'PDLIM5'</li><li>'FBLIM1'</li><li>'RASGRP3'</li><li>'ISG15'</li><li>'S1PR1'</li><li>'GSDMD'</li><li>'NMI'</li><li>'SOX4'</li><li>'KDR'</li><li>'TIE1'</li><li>'YBX3'</li><li>'LYL1'</li><li>'BACE2'</li><li>'SCARB2'</li><li>'CRIP2'</li><li>'GJA1'</li><li>'PLSCR1'</li><li>'MARCKS'</li><li>'RHOB'</li><li>'FNDC3B'</li><li>'IFITM1'</li><li>'LUZP1'</li><li>'EPAS1'</li><li>'IFITM3'</li><li>'IFI16'</li><li>'LAMB1'</li><li>'TNFSF13B'</li><li>'PDE4B'</li><li>'PALM2-AKAP2'</li><li>'PHACTR2'</li><li>'ITGB1'</li><li>'RAB13'</li><li>'FCGRT'</li><li>'SRPX'</li><li>'EHD4'</li><li>'PSMB10'</li><li>'TCIM'</li><li>'MYO10'</li><li>'ACVRL1'</li><li>'PKIG'</li><li>'CLDN14'</li><li>'FLNB'</li><li>'TUBB6'</li><li>'STAT2'</li><li>'ARHGAP21'</li><li>'NRCAM'</li><li>'SERPINH1'</li><li>'CD59'</li><li>'STOM'</li><li>'JAK1'</li><li>'PDGFA'</li><li>'CCND1'</li><li>'PML'</li><li>'GIMAP6'</li><li>'GRN'</li><li>'GIMAP7'</li><li>'KRT18'</li><li>'NECTIN2'</li><li>'PPFIBP1'</li><li>'GRASP'</li><li>'MX2'</li><li>'HLA-E'</li><li>'CTSB'</li><li>'SAMD9'</li><li>'MAST4'</li><li>'DDX60'</li><li>'RDX'</li><li>'GNB4'</li><li>'FZD4'</li><li>'RCN1'</li><li>'CTHRC1'</li><li>'SMAD6'</li><li>'NAMPT'</li><li>'TINAGL1'</li><li>'CCPG1'</li><li>'CCL14'</li><li>'WDFY1'</li><li>'HELZ2'</li><li>'VCAM1'</li><li>'VGLL4'</li><li>'LY6E'</li><li>'ADAR'</li><li>'PSME2'</li><li>'IGF2BP3'</li><li>'FKBP1A'</li><li>'BMP6'</li><li>'PCDH7'</li><li>'GNAI2'</li><li>'SPOCD1'</li><li>'CMPK2'</li><li>'SNCA'</li><li>'ODF2L'</li><li>'F2R'</li><li>'SPAG9'</li><li>'CBR3'</li><li>'EIF2AK2'</li><li>'IL1RL1'</li><li>'WSB1'</li><li>'FABP5'</li><li>'PCAT19'</li><li>'PTPRE'</li><li>'SASH1'</li><li>'SMURF2'</li><li>'CTNNA1'</li><li>'ESM1'</li><li>'MARCKSL1'</li><li>'FAM241A'</li><li>'RBMS1'</li><li>'DYRK4'</li><li>'PLOD2'</li><li>'GSTK1'</li><li>'SLC7A11'</li><li>'MEF2A'</li><li>'MARCH3'</li><li>'CLIC4'</li><li>'RTN4'</li><li>'MSN'</li><li>'NFE2L3'</li><li>'DIPK1B'</li><li>'PARP9'</li><li>'SP110'</li><li>'CARD16'</li><li>'PTTG1IP'</li><li>'CD99'</li><li>'CST3'</li><li>'RHOC'</li></ol>
</dd>
	<dt>$Upregulated_by_2_3_CGAMP</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CCL5'</li><li>'FASLG'</li><li>'SELE'</li><li>'TFP1'</li><li>'CX3CL1'</li><li>'B2M'</li><li>'CASP1'</li><li>'MMP9'</li><li>'TEK'</li><li>'VCAM1'</li><li>'TNFSF10'</li><li>'ICAM1'</li><li>'PTGIS'</li><li>'IL7'</li><li>'FGF2'</li><li>'TYMP'</li><li>'VWF'</li><li>'TGFB1'</li><li>'SOD1'</li><li>'SERPINE1'</li><li>'TJP1'</li><li>'OCLN'</li><li>'CLDN5'</li></ol>
</dd>
	<dt>$Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CCL5'</li><li>'TEK'</li><li>'VCAM1'</li><li>'SELE'</li><li>'IL7'</li><li>'FGF2'</li><li>'ICAM1'</li><li>'CXC3L1'</li><li>'TYMP'</li><li>'TJP1'</li><li>'OCLN'</li><li>'CLDN5'</li><li>'B2M'</li></ol>
</dd>
	<dt>$MN4_EC_Phenotype_Top30</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CXCL11'</li><li>'CXCL10'</li><li>'IFI44L'</li><li>'GBP4'</li><li>'IFI27'</li><li>'PECAM1'</li><li>'CD93'</li><li>'ECSCR'</li><li>'LGALS9'</li><li>'GBP2'</li><li>'MCAM'</li><li>'CLEC14A'</li><li>'MX1'</li><li>'CTSS'</li><li>'ADGRL4'</li><li>'MMRN1'</li><li>'HAPLN3'</li><li>'CCL2'</li><li>'HSPG2'</li><li>'PTGS2'</li><li>'THBS1'</li><li>'BST2'</li><li>'CDH5'</li><li>'CXCL9'</li><li>'EGFL7'</li><li>'PODXL'</li><li>'IFI6'</li><li>'IFIT3'</li><li>'IDO1'</li><li>'PRXL2A'</li></ol>
</dd>
	<dt>$AyersIFNG</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'IDO1'</li><li>'CXCL10'</li><li>'CXCL9'</li><li>'HLA-DRA'</li><li>'STAT1'</li><li>'IFNG'</li></ol>
</dd>
	<dt>$Exhaustion</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'LAG3'</li><li>'HAVCR2'</li><li>'TIGIT'</li><li>'TOX'</li><li>'EOMES'</li><li>'ENTPD1'</li></ol>
</dd>
	<dt>$Proliferation</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'MKI67'</li><li>'TOP2A'</li><li>'MCM3'</li><li>'MCM2'</li><li>'MCM7'</li><li>'MCM4'</li><li>'MCM5'</li><li>'MCM6'</li><li>'RFC2'</li><li>'RFC4'</li><li>'SMC4'</li><li>'SMC2'</li><li>'CDK1'</li><li>'RFC3'</li><li>'KIF11'</li><li>'RRM1'</li></ol>
</dd>
	<dt>$ICB_Targets</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'CTLA4'</li><li>'PDCD1'</li><li>'CD274'</li><li>'LAG3'</li><li>'HAVCR2'</li><li>'TIGIT'</li><li>'TNFRSF4'</li><li>'TNFRSF18'</li><li>'IDO1'</li><li>'CD276'</li><li>'VTCN1'</li><li>'CD40'</li><li>'TNFRSF9'</li></ol>
</dd>
	<dt>$GOBP_LEUKOCYTE_MIGRATION</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'ABL1'</li><li>'ABL2'</li><li>'ADA'</li><li>'ADAM10'</li><li>'ADAM17'</li><li>'ADAM8'</li><li>'ADD2'</li><li>'ADGRE2'</li><li>'ADORA1'</li><li>'ADTRP'</li><li>'AGER'</li><li>'AIF1'</li><li>'AIMP1'</li><li>'AIRE'</li><li>'AKIRIN1'</li><li>'AKT1'</li><li>'ALOX5'</li><li>'ANO6'</li><li>'ANXA1'</li><li>'APOD'</li><li>'APP'</li><li>'ARHGEF5'</li><li>'ARTN'</li><li>'ASB2'</li><li>'ASCL2'</li><li>'AZU1'</li><li>'B4GALT1'</li><li>'BCR'</li><li>'BDKRB1'</li><li>'BMP5'</li><li>'BSG'</li><li>'BST1'</li><li>'C1QBP'</li><li>'C3AR1'</li><li>'C5'</li><li>'C5AR1'</li><li>'C5AR2'</li><li>'CALCA'</li><li>'CALR'</li><li>'CAMK1D'</li><li>'CCL1'</li><li>'CCL11'</li><li>'CCL13'</li><li>'CCL14'</li><li>'CCL15'</li><li>'CCL16'</li><li>'CCL18'</li><li>'CCL19'</li><li>'CCL2'</li><li>'CCL20'</li><li>'CCL21'</li><li>'CCL22'</li><li>'CCL23'</li><li>'CCL24'</li><li>'CCL25'</li><li>'CCL26'</li><li>'CCL28'</li><li>'CCL3'</li><li>'CCL3L3'</li><li>'CCL4'</li><li>'CCL4L2'</li><li>'CCL5'</li><li>'CCL7'</li><li>'CCL8'</li><li>'CCN3'</li><li>'CCR1'</li><li>'CCR2'</li><li>'CCR5'</li><li>'CCR6'</li><li>'CCR7'</li><li>'CD177'</li><li>'CD200'</li><li>'CD200R1'</li><li>'CD300A'</li><li>'CD300H'</li><li>'CD34'</li><li>'CD47'</li><li>'CD74'</li><li>'CD81'</li><li>'CD9'</li><li>'CD99'</li><li>'CD99L2'</li><li>'CDC42'</li><li>'CH25H'</li><li>'CHGA'</li><li>'CHST2'</li><li>'CHST4'</li><li>'CKLF'</li><li>'CMKLR1'</li><li>'CNN2'</li><li>'CNR2'</li><li>'CORO1A'</li><li>'CREB3'</li><li>'CRK'</li><li>'CRKL'</li><li>'CRTAM'</li><li>'CSF1'</li><li>'CSF1R'</li><li>'CSF3R'</li><li>'CTSG'</li><li>'CX3CL1'</li><li>'CX3CR1'</li><li>'CXADR'</li><li>'CXCL1'</li><li>'CXCL10'</li><li>'CXCL11'</li><li>'CXCL12'</li><li>'CXCL13'</li><li>'CXCL16'</li><li>'CXCL17'</li><li>'CXCL2'</li><li>'CXCL3'</li><li>'CXCL5'</li><li>'CXCL6'</li><li>'CXCL8'</li><li>'CXCL9'</li><li>'CXCR1'</li><li>'CXCR2'</li><li>'CXCR3'</li><li>'CXCR4'</li><li>'CXCR5'</li><li>'CYP19A1'</li><li>'CYP7B1'</li><li>'DAPK2'</li><li>'DBH'</li><li>'DDT'</li><li>'DEFA1'</li><li>'DEFA1B'</li><li>'DEFB104A'</li><li>'DEFB104B'</li><li>'DEFB124'</li><li>'DNM1L'</li><li>'DOCK8'</li><li>'DPEP1'</li><li>'DPP4'</li><li>'DUSP1'</li><li>'ECM1'</li><li>'EDN1'</li><li>'EDN2'</li><li>'EDN3'</li><li>'EDNRB'</li><li>'ELANE'</li><li>'EMILIN1'</li><li>'EMP2'</li><li>'EPS8'</li><li>'EPX'</li><li>'EXT1'</li><li>'F11R'</li><li>'F2RL1'</li><li>'F7'</li><li>'FADD'</li><li>'FCER1G'</li><li>'FER'</li><li>'FFAR2'</li><li>'FLT1'</li><li>'FOLR2'</li><li>'FOXJ1'</li><li>'FPR2'</li><li>'FUT4'</li><li>'FUT7'</li><li>'FUT9'</li><li>'FYN'</li><li>'GAS6'</li><li>'GATA3'</li><li>'GBA1'</li><li>'GBF1'</li><li>'GCNT1'</li><li>'GCSAM'</li><li>'GCSAML'</li><li>'GDF15'</li><li>'GOLPH3'</li><li>'GP1BA'</li><li>'GP2'</li><li>'GPR15'</li><li>'GPR15LG'</li><li>'GPR18'</li><li>'GPR183'</li><li>'GPSM3'</li><li>'GREM1'</li><li>'H2BC1'</li><li>'HCK'</li><li>'HMGB1'</li><li>'HMOX1'</li><li>'HOXA7'</li><li>'HSD3B7'</li><li>'ICAM1'</li><li>'IL10'</li><li>'IL12A'</li><li>'IL16'</li><li>'IL17A'</li><li>'IL17RA'</li><li>'IL17RC'</li><li>'IL1A'</li><li>'IL1B'</li><li>'IL1R1'</li><li>'IL23A'</li><li>'IL27RA'</li><li>'IL33'</li><li>'IL34'</li><li>'IL6'</li><li>'IL6R'</li><li>'IRAK4'</li><li>'ITGA1'</li><li>'ITGA2'</li><li>'ITGA2B'</li><li>'ITGA3'</li><li>'ITGA4'</li><li>'ITGA6'</li><li>'ITGA7'</li><li>'ITGA9'</li><li>'ITGAL'</li><li>'ITGB1'</li><li>'ITGB2'</li><li>'ITGB3'</li><li>'ITGB7'</li><li>'JAGN1'</li><li>'JAM2'</li><li>'JAM3'</li><li>'JAML'</li><li>'KIT'</li><li>'KITLG'</li><li>'KLRC4-KLRK1'</li><li>'KLRK1'</li><li>'LBP'</li><li>'LCK'</li><li>'LEP'</li><li>'LGALS3'</li><li>'LGALS9'</li><li>'LGMN'</li><li>'LRCH1'</li><li>'LYN'</li><li>'LYST'</li><li>'LYVE1'</li><li>'MADCAM1'</li><li>'MAPK1'</li><li>'MAPK3'</li><li>'MCOLN2'</li><li>'MCU'</li><li>'MDK'</li><li>'MED23'</li><li>'MIA3'</li><li>'MICOS10-NBL1'</li><li>'MIF'</li><li>'MIR128-1'</li><li>'MIR146A'</li><li>'MIR223'</li><li>'MIR24-1'</li><li>'MMP14'</li><li>'MMP2'</li><li>'MMP28'</li><li>'MOSPD2'</li><li>'MPP1'</li><li>'MSMP'</li><li>'MSN'</li><li>'MSTN'</li><li>'MTUS1'</li><li>'MYD88'</li><li>'MYH9'</li><li>'MYO1G'</li><li>'NBL1'</li><li>'NCKAP1L'</li><li>'NEDD9'</li><li>'NF1'</li><li>'NINJ1'</li><li>'NKX2-3'</li><li>'NLRP12'</li><li>'NOD2'</li><li>'NUP85'</li><li>'OXSR1'</li><li>'P2RX4'</li><li>'P2RY12'</li><li>'PADI2'</li><li>'PAFAH1B1'</li><li>'PDE4B'</li><li>'PDGFB'</li><li>'PDGFD'</li><li>'PECAM1'</li><li>'PERP'</li><li>'PF4'</li><li>'PF4V1'</li><li>'PGF'</li><li>'PIK3CD'</li><li>'PIK3CG'</li><li>'PIK3R1'</li><li>'PIKFYVE'</li><li>'PIP5K1C'</li><li>'PLA2G1B'</li><li>'PLA2G7'</li><li>'PLCB1'</li><li>'PLEC'</li><li>'PLG'</li><li>'PLVAP'</li><li>'PODXL2'</li><li>'PPBP'</li><li>'PPIA'</li><li>'PPIB'</li><li>'PREX1'</li><li>'PRSS56'</li><li>'PRTN3'</li><li>'PTAFR'</li><li>'PTGER4'</li><li>'PTK2'</li><li>'PTK2B'</li><li>'PTN'</li><li>'PTPRJ'</li><li>'PTPRO'</li><li>'PYCARD'</li><li>'RABGEF1'</li><li>'RAC1'</li><li>'RAC2'</li><li>'RAC3'</li><li>'RARRES2'</li><li>'RET'</li><li>'RHOA'</li><li>'RIN3'</li><li>'RIPK3'</li><li>'RIPOR2'</li><li>'ROCK1'</li><li>'ROR2'</li><li>'RPL13A'</li><li>'RPS19'</li><li>'S100A12'</li><li>'S100A14'</li><li>'S100A7'</li><li>'S100A8'</li><li>'S100A9'</li><li>'S1PR1'</li><li>'SAA1'</li><li>'SBDS'</li><li>'SCG2'</li><li>'SELE'</li><li>'SELENOK'</li><li>'SELL'</li><li>'SELP'</li><li>'SELPLG'</li><li>'SERPINE1'</li><li>'SFTPD'</li><li>'SIRPA'</li><li>'SLAMF1'</li><li>'SLAMF8'</li><li>'SLC12A2'</li><li>'SLC8B1'</li><li>'SLIT2'</li><li>'SMPD3'</li><li>'SOS1'</li><li>'SPN'</li><li>'SPNS2'</li><li>'SRC'</li><li>'SRP54'</li><li>'ST3GAL4'</li><li>'STAP1'</li><li>'STAT5B'</li><li>'STK10'</li><li>'STK39'</li><li>'SWAP70'</li><li>'SYK'</li><li>'TACR1'</li><li>'TAFA4'</li><li>'TBX21'</li><li>'TGFB1'</li><li>'TGFB2'</li><li>'THBS1'</li><li>'THBS4'</li><li>'THY1'</li><li>'TIRAP'</li><li>'TMEM102'</li><li>'TNF'</li><li>'TNFAIP6'</li><li>'TNFRSF11A'</li><li>'TNFRSF14'</li><li>'TNFRSF18'</li><li>'TNFSF11'</li><li>'TNFSF14'</li><li>'TNFSF18'</li><li>'TREM1'</li><li>'TREM2'</li><li>'TRIM55'</li><li>'TRPM2'</li><li>'TRPM4'</li><li>'TRPV4'</li><li>'UMOD'</li><li>'VAV1'</li><li>'VAV3'</li><li>'VCAM1'</li><li>'VEGFA'</li><li>'VEGFB'</li><li>'VEGFC'</li><li>'VEGFD'</li><li>'WASL'</li><li>'WDR1'</li><li>'WNK1'</li><li>'WNT5A'</li><li>'XCL1'</li><li>'XCL2'</li><li>'XG'</li><li>'YES1'</li><li>'ZAP70'</li><li>'ZNF580'</li><li>'ZP3'</li></ol>
</dd>
	<dt>$GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</dt>
		<dd><style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'ABL1'</li><li>'ABL2'</li><li>'ADAM10'</li><li>'ADAM17'</li><li>'ADAM8'</li><li>'AGER'</li><li>'AIF1'</li><li>'AKIRIN1'</li><li>'ANO6'</li><li>'APP'</li><li>'ASCL2'</li><li>'BDKRB1'</li><li>'C1QBP'</li><li>'C3AR1'</li><li>'C5AR1'</li><li>'CALR'</li><li>'CAMK1D'</li><li>'CCL1'</li><li>'CCL19'</li><li>'CCL20'</li><li>'CCL21'</li><li>'CCL24'</li><li>'CCL3'</li><li>'CCL4'</li><li>'CCL5'</li><li>'CCL7'</li><li>'CCL8'</li><li>'CCR1'</li><li>'CCR2'</li><li>'CCR6'</li><li>'CCR7'</li><li>'CD47'</li><li>'CD74'</li><li>'CD99'</li><li>'CD99L2'</li><li>'CMKLR1'</li><li>'CORO1A'</li><li>'CREB3'</li><li>'CSF1'</li><li>'CSF1R'</li><li>'CX3CL1'</li><li>'CX3CR1'</li><li>'CXCL10'</li><li>'CXCL12'</li><li>'CXCL13'</li><li>'CXCL17'</li><li>'CXCL8'</li><li>'DAPK2'</li><li>'DEFB124'</li><li>'DNM1L'</li><li>'DOCK8'</li><li>'EDN1'</li><li>'EDN2'</li><li>'EDN3'</li><li>'F2RL1'</li><li>'F7'</li><li>'FADD'</li><li>'FPR2'</li><li>'FUT4'</li><li>'FUT7'</li><li>'GAS6'</li><li>'GPSM3'</li><li>'HMGB1'</li><li>'ICAM1'</li><li>'IL12A'</li><li>'IL1A'</li><li>'IL1R1'</li><li>'IL23A'</li><li>'IL34'</li><li>'IL6'</li><li>'IL6R'</li><li>'ITGA2'</li><li>'ITGA2B'</li><li>'ITGA4'</li><li>'ITGB3'</li><li>'JAM2'</li><li>'JAM3'</li><li>'KITLG'</li><li>'LBP'</li><li>'LGALS3'</li><li>'LGALS9'</li><li>'LGMN'</li><li>'LYVE1'</li><li>'MADCAM1'</li><li>'MAPK1'</li><li>'MAPK3'</li><li>'MCU'</li><li>'MDK'</li><li>'MED23'</li><li>'MIA3'</li><li>'MMP14'</li><li>'MOSPD2'</li><li>'MSTN'</li><li>'NCKAP1L'</li><li>'NEDD9'</li><li>'OXSR1'</li><li>'P2RX4'</li><li>'P2RY12'</li><li>'PDGFD'</li><li>'PERP'</li><li>'PGF'</li><li>'PIK3R1'</li><li>'PLA2G7'</li><li>'PLVAP'</li><li>'PRSS56'</li><li>'PTAFR'</li><li>'PTK2'</li><li>'PTK2B'</li><li>'PTN'</li><li>'PTPRJ'</li><li>'PYCARD'</li><li>'RAC1'</li><li>'RAC2'</li><li>'RARRES2'</li><li>'RHOA'</li><li>'RIPOR2'</li><li>'S100A14'</li><li>'S100A7'</li><li>'SELE'</li><li>'SELENOK'</li><li>'SELP'</li><li>'SERPINE1'</li><li>'SLAMF1'</li><li>'SPN'</li><li>'STK39'</li><li>'SWAP70'</li><li>'TACR1'</li><li>'TGFB1'</li><li>'THBS1'</li><li>'THBS4'</li><li>'THY1'</li><li>'TIRAP'</li><li>'TMEM102'</li><li>'TNF'</li><li>'TNFRSF14'</li><li>'TNFRSF18'</li><li>'TNFSF14'</li><li>'TNFSF18'</li><li>'TREM2'</li><li>'TRPV4'</li><li>'VEGFA'</li><li>'VEGFB'</li><li>'VEGFC'</li><li>'VEGFD'</li><li>'WNK1'</li><li>'WNT5A'</li><li>'XCL1'</li><li>'XG'</li><li>'ZNF580'</li><li>'ZP3'</li></ol>
</dd>
</dl>



# Save the lists as a json so you can load them into a different script


```R
write(toJSON(updated_signature_lists, pretty=TRUE, auto_unbox=TRUE), file = "Processing/gene_signature_lists_20250502.json")
```

# save as GMT file (file must not exist, as it appends to it)


```R
out_file = 'Processing/old_GSVA_docker/updated_signature_lists.gmt'

# Check if the file exists and delete it if it does
if (file.exists(out_file)) {
    file.remove(out_file)
}

for (name_set in names(updated_signature_lists)){
    description = name_set
    genes = updated_signature_lists[[name_set]]
    genes[length(genes)] = paste(genes[length(genes)], "\n", sep="")
    cat(name_set, description, genes, sep="\t", file = out_file, append = TRUE)
}
```


TRUE


# save as a csv for marco


```R
# Find the maximum length of vectors in the list
max_length <- max(sapply(updated_signature_lists, length))

# Pad each list to the max length with NA values
padded_list <- lapply(updated_signature_lists, function(x) {
  length(x) <- max_length
  x
})

# Convert the padded list to a data frame
df <- as.data.frame(padded_list, stringsAsFactors = FALSE)

# Write to CSV file
write.csv(df, "Output/05_do_GSVA/GSVA_genesets.csv", row.names = FALSE)
```

# final prep for the data

# prepare variable genes data


```R
# get lists of variable genes: 
var_genes_VA2 <- VariableFeatures(VA2)
var_genes_VA1 <- VariableFeatures(VA1)
var_genes_VB1 <- VariableFeatures(VB1)
```


```R
# subset datasets to include only those variable genes: 
data_VA2_vargenes <- GetAssayData(VA2)[var_genes_VA2,]
data_VA1_vargenes <- GetAssayData(VA1)[var_genes_VA1,]
data_VB1_vargenes <- GetAssayData(VB1)[var_genes_VB1,]
```


```R
dim(data_VA2_vargenes)
dim(data_VA1_vargenes)
dim(data_VB1_vargenes)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>10000</li><li>4316</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>10000</li><li>2806</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>10000</li><li>2589</li></ol>




```R
# save the gene expression data for GSVA
write.csv(data_VA2_vargenes, 'Processing/old_GSVA_docker/data_VA2_vargenes.csv')
write.csv(data_VA1_vargenes, 'Processing/old_GSVA_docker/data_VA1_vargenes.csv')
write.csv(data_VB1_vargenes, 'Processing/old_GSVA_docker/data_VB1_vargenes.csv')
```


```R

```

# Then, run the GSVA docker separately. 

# ran in terminal for each
docker run --rm --user 63704:1048860 -v /immuno:/immuno vacation/gsva:latest GSVA --gmt /immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/updated_signature_lists.gmt /immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/data_VA1_vargenes.csv --output /immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/VA1_GSVA_output.csv

docker run --rm --user 63704:1048860 -v /immuno:/immuno vacation/gsva:latest GSVA --gmt /immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/updated_signature_lists.gmt /immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/data_VB1_vargenes.csv --output /immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/VB1_GSVA_output.csv

docker run --rm --user 63704:1048860 -v /immuno:/immuno vacation/gsva:latest GSVA --gmt /immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/updated_signature_lists.gmt /immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/data_VA2_vargenes.csv --output /immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/VA2_GSVA_output.csv

# Now, re-load in the GSVA results. 


```R
Endo_gsva_VA2 <- read_csv('/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/VA2_GSVA_output.csv')
Endo_gsva_VA1 <- read_csv('/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/VA1_GSVA_output.csv')
Endo_gsva_VB1 <- read_csv('/immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Processing/old_GSVA_docker/VB1_GSVA_output.csv')
```

    [1mRows: [22m[34m30[39m [1mColumns: [22m[34m4317[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m    (1): name
    [32mdbl[39m (4316): AAACAAGTATCTCCCA.1, AAACACCAATAACTGC.1, AAACAGAGCGACTCCT.1, AAAC...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.
    [1mRows: [22m[34m30[39m [1mColumns: [22m[34m2807[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m    (1): name
    [32mdbl[39m (2806): AAACAACGAATAGTTC.1, AAACAAGTATCTCCCA.1, AAACAATCTACTAGCA.1, AAAC...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.
    [1mRows: [22m[34m30[39m [1mColumns: [22m[34m2590[39m
    [36m──[39m [1mColumn specification[22m [36m────────────────────────────────────────────────────────[39m
    [1mDelimiter:[22m ","
    [31mchr[39m    (1): name
    [32mdbl[39m (2589): AAACAACGAATAGTTC.1, AAACAAGTATCTCCCA.1, AAACACCAATAACTGC.1, AAAC...
    
    [36mℹ[39m Use `spec()` to retrieve the full column specification for this data.
    [36mℹ[39m Specify the column types or set `show_col_types = FALSE` to quiet this message.



```R
# set name column to rownames
Endo_gsva_VA2 <- Endo_gsva_VA2 %>% column_to_rownames("name")
Endo_gsva_VA1 <- Endo_gsva_VA1 %>% column_to_rownames("name")
Endo_gsva_VB1 <- Endo_gsva_VB1 %>% column_to_rownames("name")
# replace . with - 
colnames(Endo_gsva_VA2) <- gsub("\\.", "-", colnames(Endo_gsva_VA2))
colnames(Endo_gsva_VA1) <- gsub("\\.", "-", colnames(Endo_gsva_VA1))
colnames(Endo_gsva_VB1) <- gsub("\\.", "-", colnames(Endo_gsva_VB1))
```


```R
dim(Endo_gsva_VA2)
dim(Endo_gsva_VA1)
dim(Endo_gsva_VB1)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>30</li><li>4316</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>30</li><li>2806</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>30</li><li>2589</li></ol>




```R
Endo_gsva_VA2 %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 4316</caption>
<thead>
	<tr><th></th><th scope=col>AAACAAGTATCTCCCA-1</th><th scope=col>AAACACCAATAACTGC-1</th><th scope=col>AAACAGAGCGACTCCT-1</th><th scope=col>AAACAGCTTTCAGAAG-1</th><th scope=col>AAACAGGGTCTATATT-1</th><th scope=col>AAACAGTGTTCCTGGG-1</th><th scope=col>AAACATGGTGAGAGGA-1</th><th scope=col>AAACATTTCCCGGATT-1</th><th scope=col>AAACCACTACACAGAT-1</th><th scope=col>AAACCGGGTAGGTACC-1</th><th scope=col>AAACCGTTCGTCCAGG-1</th><th scope=col>AAACCTAAGCAGCCGG-1</th><th scope=col>AAACCTCATGAAGTTG-1</th><th scope=col>AAACGAAGAACATACC-1</th><th scope=col>AAACGAAGATGGAGTA-1</th><th scope=col>AAACGACAGTCTTGCC-1</th><th scope=col>AAACGAGACGGTTGAT-1</th><th scope=col>AAACGCCCGAGATCGG-1</th><th scope=col>AAACGGGCGTACGGGT-1</th><th scope=col>AAACGGTTGCGAACTG-1</th><th scope=col>AAACGTGTTCGCCCTA-1</th><th scope=col>AAACTAACGTGGCGAC-1</th><th scope=col>AAACTCGGTTCGCAAT-1</th><th scope=col>AAACTCGTGATATAAG-1</th><th scope=col>AAACTGCTGGCTCCAA-1</th><th scope=col>AAACTTAATTGCACGC-1</th><th scope=col>AAACTTGCAAACGTAT-1</th><th scope=col>AAAGAATGACCTTAGA-1</th><th scope=col>AAAGAATGTGGACTAA-1</th><th scope=col>AAAGACATGAAGTTTA-1</th><th scope=col>AAAGACCCAAGTCGCG-1</th><th scope=col>AAAGACTGGGCGCTTT-1</th><th scope=col>AAAGGCCCTATAATAC-1</th><th scope=col>AAAGGCTACGGACCAT-1</th><th scope=col>AAAGGCTCTCGCGCCG-1</th><th scope=col>AAAGGGATGTAGCAAG-1</th><th scope=col>AAAGGGCAGCTTGAAT-1</th><th scope=col>AAAGGTAAGCTGTACC-1</th><th scope=col>AAAGGTCAACGACATG-1</th><th scope=col>AAAGTAGCATTGCTCA-1</th><th scope=col>AAAGTCACTGATGTAA-1</th><th scope=col>AAAGTCGACCCTCAGT-1</th><th scope=col>AAAGTGTGATTTATCT-1</th><th scope=col>AAAGTTGACTCCCGTA-1</th><th scope=col>AAATAACCATACGGGA-1</th><th scope=col>AAATACCTATAAGCAT-1</th><th scope=col>AAATAGCTTAGACTTT-1</th><th scope=col>AAATAGGGTGCTATTG-1</th><th scope=col>AAATCACTCCTAAACG-1</th><th scope=col>AAATCCGATACACGCC-1</th><th scope=col>⋯</th><th scope=col>TTGCTGCACCTATCCA-1</th><th scope=col>TTGCTGGCCGGGCTTC-1</th><th scope=col>TTGGAAGAATACAGTC-1</th><th scope=col>TTGGACATGTGGCTTA-1</th><th scope=col>TTGGACCATCTGGCAA-1</th><th scope=col>TTGGACCTATAACAGT-1</th><th scope=col>TTGGAGTCTCCCTTCT-1</th><th scope=col>TTGGATATCGTCTACG-1</th><th scope=col>TTGGATCGACTTCTGG-1</th><th scope=col>TTGGATTGGGTACCAC-1</th><th scope=col>TTGGCCAAATTGTATC-1</th><th scope=col>TTGGCCATCTTGCGCT-1</th><th scope=col>TTGGCCTAGAATTTCG-1</th><th scope=col>TTGGCTCGCATGAGAC-1</th><th scope=col>TTGGGAAGACGAGCCG-1</th><th scope=col>TTGGGACACTGCCCGC-1</th><th scope=col>TTGGGACGTAAGAGTT-1</th><th scope=col>TTGGGCGGCGGTTGCC-1</th><th scope=col>TTGGGTTTATTCAGCG-1</th><th scope=col>TTGGTCACACTCGTAA-1</th><th scope=col>TTGGTGCGGTGTTGAA-1</th><th scope=col>TTGGTTGCGGTGCGCG-1</th><th scope=col>TTGTAACTTCATAGCG-1</th><th scope=col>TTGTAAGGACCTAAGT-1</th><th scope=col>TTGTAAGGCCAGTTGG-1</th><th scope=col>TTGTAATCCGTACTCG-1</th><th scope=col>TTGTACACCTCGAACA-1</th><th scope=col>TTGTATCACACAGAAT-1</th><th scope=col>TTGTCACCGCGGTATC-1</th><th scope=col>TTGTCGTTCAGTTACC-1</th><th scope=col>TTGTCTCGGCAAGATG-1</th><th scope=col>TTGTGAACCTAATCCG-1</th><th scope=col>TTGTGAGGCATGACGC-1</th><th scope=col>TTGTGCAGCCACGTCA-1</th><th scope=col>TTGTGCGGAAGCGGAT-1</th><th scope=col>TTGTGGCCCTGACAGT-1</th><th scope=col>TTGTGGTAGGAGGGAT-1</th><th scope=col>TTGTGGTGGTACTAAG-1</th><th scope=col>TTGTGTATGCCACCAA-1</th><th scope=col>TTGTGTTTCCCGAAAG-1</th><th scope=col>TTGTTAGCAAATTCGA-1</th><th scope=col>TTGTTCAGTGTGCTAC-1</th><th scope=col>TTGTTCTAGATACGCT-1</th><th scope=col>TTGTTGGCAATGACTG-1</th><th scope=col>TTGTTGTGTGTCAAGA-1</th><th scope=col>TTGTTTCACATCCAGG-1</th><th scope=col>TTGTTTCATTAGTCTA-1</th><th scope=col>TTGTTTCCATACAACT-1</th><th scope=col>TTGTTTGTATTACACG-1</th><th scope=col>TTGTTTGTGTAAATTC-1</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>All_Immune</th><td>-0.4509680</td><td>0.4657553</td><td>-0.5878032</td><td> 0.1120339</td><td>-0.4662198</td><td> 0.4454450</td><td>0.4450079</td><td> 0.2332552</td><td>-0.4711797</td><td>-0.13139340</td><td>-0.0004884595</td><td>0.2887079</td><td> 0.1160454</td><td>-0.4613412</td><td>0.5752656</td><td>-0.5663370</td><td>-0.4149761</td><td>-0.4839030</td><td> 0.3810739</td><td> 0.2749017</td><td>-0.5003722</td><td>-0.3690133</td><td> 0.2596578</td><td>-0.5281809</td><td> 0.1924258</td><td> 0.3518484</td><td>-0.4881647</td><td> 0.23090858</td><td>0.3195120</td><td>-0.5477629</td><td> 0.2234073</td><td>0.290982477</td><td> 0.2795565</td><td>-0.0187372</td><td> 0.1327373</td><td>-0.3775131</td><td>-0.2859881</td><td>-0.4973546</td><td>-0.5472882</td><td>-0.3741172</td><td> 0.06713531</td><td>-0.3780488</td><td>-0.3230999</td><td>0.3831547</td><td>-0.6041380</td><td> 0.3987243</td><td>-0.32929234</td><td>0.4984399</td><td>-0.2577851</td><td> 0.3630783</td><td>⋯</td><td>-0.03241744</td><td>-0.5141867</td><td>-0.2047050</td><td>-0.32076961</td><td>-0.5492200</td><td> 0.1215178</td><td>0.2482258</td><td>0.4486909</td><td>-0.4522968</td><td>-0.2251460</td><td> 0.4115713</td><td> 0.09835899</td><td>-0.2779309</td><td> 0.1606963</td><td>0.3407962</td><td>-0.5565992</td><td>-0.4922029</td><td> 0.07429765</td><td>-0.15796540</td><td> 0.38195882</td><td>-0.4427441</td><td>-0.2903827</td><td>-0.2586180</td><td>-0.4693500</td><td>-0.20057738</td><td>0.5865052</td><td>0.6469476</td><td>-0.09809743</td><td>0.02163196</td><td>-0.3908079</td><td> 0.2516179</td><td>0.4337731</td><td>-0.4773460</td><td> 0.3015086</td><td>0.4055023</td><td>-0.3956882</td><td>-0.2568207</td><td> 0.3126212</td><td>0.32807644</td><td> 0.2824381</td><td> 0.17840732</td><td>-0.3698462</td><td>0.5452200</td><td> 0.2173142</td><td>-0.3798123</td><td>0.5156872</td><td>0.5086120</td><td>-0.5089129</td><td> 0.1652607</td><td>-0.09204677</td></tr>
	<tr><th scope=row>Angiogenesis</th><td>-0.4483304</td><td>0.4768693</td><td>-0.6471249</td><td>-0.4374814</td><td>-0.6181580</td><td> 0.6592338</td><td>0.3477385</td><td> 0.2867388</td><td>-0.5357210</td><td>-0.04825953</td><td> 0.3532120687</td><td>0.5603880</td><td>-0.3283320</td><td>-0.1280803</td><td>0.1252315</td><td>-0.3990534</td><td>-0.4263739</td><td>-0.4775946</td><td> 0.6005745</td><td> 0.5379371</td><td>-0.4358833</td><td>-0.1604814</td><td> 0.5940708</td><td>-0.5708542</td><td> 0.6853363</td><td> 0.1661448</td><td>-0.3859823</td><td>-0.03683954</td><td>0.3554227</td><td>-0.6782890</td><td> 0.3911964</td><td>0.003809179</td><td> 0.6652143</td><td> 0.3837041</td><td>-0.7425053</td><td>-0.4695797</td><td> 0.5820798</td><td>-0.5164538</td><td>-0.6684685</td><td>-0.5609936</td><td>-0.01050554</td><td>-0.2072883</td><td>-0.3519201</td><td>0.2639977</td><td>-0.5373630</td><td>-0.1151478</td><td> 0.08861987</td><td>0.1513845</td><td>-0.1887910</td><td> 0.4796114</td><td>⋯</td><td> 0.04292217</td><td>-0.4525517</td><td>-0.3394941</td><td>-0.54484527</td><td>-0.6033033</td><td> 0.6433510</td><td>0.3490590</td><td>0.2538351</td><td>-0.5528601</td><td>-0.2689861</td><td> 0.3797404</td><td> 0.45532202</td><td> 0.4880685</td><td> 0.2356250</td><td>0.4667740</td><td>-0.5805193</td><td>-0.4143644</td><td>-0.12517614</td><td> 0.02338835</td><td>-0.04241745</td><td>-0.6096635</td><td>-0.5337488</td><td>-0.4800794</td><td>-0.6857988</td><td>-0.19413319</td><td>0.1093516</td><td>0.2175783</td><td>-0.16962664</td><td>0.18272659</td><td>-0.3730728</td><td> 0.6632008</td><td>0.5887285</td><td>-0.3483416</td><td>-0.3885622</td><td>0.3017526</td><td>-0.6147769</td><td>-0.4970270</td><td> 0.1604724</td><td>0.01409044</td><td> 0.1502281</td><td> 0.04419888</td><td>-0.3011749</td><td>0.3065242</td><td>-0.4036376</td><td>-0.2779393</td><td>0.3471697</td><td>0.5494112</td><td>-0.5691454</td><td>-0.0697942</td><td> 0.03176488</td></tr>
	<tr><th scope=row>AyersIFNG</th><td>-0.6914766</td><td>0.8017860</td><td>-0.2034801</td><td> 0.7277362</td><td> 0.4140443</td><td>-0.3833363</td><td>0.4490267</td><td>-0.4455322</td><td>-0.6811725</td><td> 0.77145002</td><td> 0.4789586445</td><td>0.5370652</td><td> 0.7558619</td><td>-0.1791315</td><td>0.8499016</td><td>-0.6878752</td><td> 0.1736779</td><td>-0.3839386</td><td>-0.2348940</td><td>-0.3290578</td><td>-0.4511666</td><td>-0.6425570</td><td>-0.2214886</td><td>-0.7223890</td><td>-0.2150860</td><td>-0.3940930</td><td>-0.1717308</td><td>-0.22238896</td><td>0.4705512</td><td>-0.7312925</td><td>-0.5832333</td><td>0.870148059</td><td>-0.3360634</td><td>-0.4707224</td><td>-0.4164685</td><td>-0.1731308</td><td> 0.1729505</td><td>-0.6674670</td><td> 0.1341916</td><td>-0.3869600</td><td> 0.32316942</td><td>-0.5027292</td><td>-0.6406563</td><td>0.5950460</td><td>-0.5947391</td><td> 0.7910258</td><td>-0.61704682</td><td>0.4542967</td><td>-0.5174114</td><td>-0.4356620</td><td>⋯</td><td>-0.33192632</td><td>-0.7205882</td><td> 0.1973497</td><td>-0.07297787</td><td>-0.6074430</td><td>-0.3161954</td><td>0.2733838</td><td>0.2342038</td><td>-0.5965993</td><td>-0.2452820</td><td>-0.3548128</td><td>-0.61414566</td><td> 0.1185986</td><td>-0.1250549</td><td>0.8063616</td><td>-0.5027872</td><td>-0.4565646</td><td> 0.12944008</td><td> 0.35324561</td><td>-0.24457344</td><td>-0.6575630</td><td>-0.2655657</td><td> 0.1865323</td><td>-0.5345757</td><td>-0.04486833</td><td>0.2365367</td><td>0.6801954</td><td>-0.62635054</td><td>0.59521480</td><td>-0.2844245</td><td>-0.2751635</td><td>0.1300172</td><td> 0.6199615</td><td>-0.2892237</td><td>0.8095875</td><td>-0.6418567</td><td>-0.6674670</td><td>-0.4351008</td><td>0.52152081</td><td>-0.4978214</td><td>-0.07459969</td><td>-0.3704368</td><td>0.4990678</td><td>-0.2306923</td><td>-0.5593237</td><td>0.1629367</td><td>0.7891932</td><td>-0.4593479</td><td>-0.1998800</td><td> 0.03407006</td></tr>
</tbody>
</table>




```R
Endo_gsva_VA1 %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 2806</caption>
<thead>
	<tr><th></th><th scope=col>AAACAACGAATAGTTC-1</th><th scope=col>AAACAAGTATCTCCCA-1</th><th scope=col>AAACAATCTACTAGCA-1</th><th scope=col>AAACACCAATAACTGC-1</th><th scope=col>AAACAGCTTTCAGAAG-1</th><th scope=col>AAACAGGGTCTATATT-1</th><th scope=col>AAACCCGAACGAAATC-1</th><th scope=col>AAACCGGAAATGTTAA-1</th><th scope=col>AAACCGGGTAGGTACC-1</th><th scope=col>AAACCGTTCGTCCAGG-1</th><th scope=col>AAACCTCATGAAGTTG-1</th><th scope=col>AAACGAGACGGTTGAT-1</th><th scope=col>AAACGGGTTGGTATCC-1</th><th scope=col>AAACTTGCAAACGTAT-1</th><th scope=col>AAAGACTGGGCGCTTT-1</th><th scope=col>AAAGCTTGCCTACATA-1</th><th scope=col>AAAGGCTACGGACCAT-1</th><th scope=col>AAAGGCTCTCGCGCCG-1</th><th scope=col>AAAGGGATGTAGCAAG-1</th><th scope=col>AAAGGGCAGCTTGAAT-1</th><th scope=col>AAAGTAGCATTGCTCA-1</th><th scope=col>AAAGTCGACCCTCAGT-1</th><th scope=col>AAAGTGTGATTTATCT-1</th><th scope=col>AAAGTTGACTCCCGTA-1</th><th scope=col>AAATAAGGTAGTGCCC-1</th><th scope=col>AAATAGCTTAGACTTT-1</th><th scope=col>AAATAGGGTGCTATTG-1</th><th scope=col>AAATCGCGGAAGGAGT-1</th><th scope=col>AAATCGTGTACCACAA-1</th><th scope=col>AAATCTAGCCCTGCTA-1</th><th scope=col>AAATGCTCGTTACGTT-1</th><th scope=col>AAATGGCCCGTGCCCT-1</th><th scope=col>AAATGGTCAATGTGCC-1</th><th scope=col>AAATGTATCTTATCCC-1</th><th scope=col>AAATTAACGGGTAGCT-1</th><th scope=col>AAATTAATAAGCGCGA-1</th><th scope=col>AAATTACCTATCGATG-1</th><th scope=col>AAATTCCAGGTCCAAA-1</th><th scope=col>AAATTGATAGTCCTTT-1</th><th scope=col>AAATTTGCGGGTGTGG-1</th><th scope=col>AACAACTGGTAGTTGC-1</th><th scope=col>AACAATACATTGTCGA-1</th><th scope=col>AACAATTACTCTACGC-1</th><th scope=col>AACACACGCTCGCCGC-1</th><th scope=col>AACACGAGACGCGGCC-1</th><th scope=col>AACACGCGGCCGCGAA-1</th><th scope=col>AACAGCTGTGTGGCAA-1</th><th scope=col>AACAGGATGGGCCGCG-1</th><th scope=col>AACAGGTAGTATGGAT-1</th><th scope=col>AACATATCAACTGGTG-1</th><th scope=col>⋯</th><th scope=col>TTGAGAAGTTTAGCAT-1</th><th scope=col>TTGAGAGTACTGCTAA-1</th><th scope=col>TTGAGCAGCCCACGGT-1</th><th scope=col>TTGAGCGCCACGTGAT-1</th><th scope=col>TTGAGTCCCGCTGCTG-1</th><th scope=col>TTGATCTAACTTTGTC-1</th><th scope=col>TTGATTAGCTGTTTCT-1</th><th scope=col>TTGATTATGCAGATGA-1</th><th scope=col>TTGCACAATTCAGAAA-1</th><th scope=col>TTGCATGCTGATCACG-1</th><th scope=col>TTGCCAAGCAGAACCC-1</th><th scope=col>TTGCCATAGCCCGCTC-1</th><th scope=col>TTGCCCTGATCACGGG-1</th><th scope=col>TTGCGCTCTCTCGCTT-1</th><th scope=col>TTGCGGCATCAGAAAG-1</th><th scope=col>TTGCGTAGTTTGAGGA-1</th><th scope=col>TTGCGTCGGCCAACCG-1</th><th scope=col>TTGCTAGCTACCAATC-1</th><th scope=col>TTGCTGATCATGTTCG-1</th><th scope=col>TTGCTGCACCTATCCA-1</th><th scope=col>TTGGAAGAATACAGTC-1</th><th scope=col>TTGGACATGTGGCTTA-1</th><th scope=col>TTGGACCATCTGGCAA-1</th><th scope=col>TTGGACCTATAACAGT-1</th><th scope=col>TTGGATATCGTCTACG-1</th><th scope=col>TTGGATCGACTTCTGG-1</th><th scope=col>TTGGCCATCTTGCGCT-1</th><th scope=col>TTGGCGATCCGAATAT-1</th><th scope=col>TTGGCTCGCATGAGAC-1</th><th scope=col>TTGGGACGTAAGAGTT-1</th><th scope=col>TTGGTATGGCTTGTGT-1</th><th scope=col>TTGGTCACACTCGTAA-1</th><th scope=col>TTGGTTGCGGTGCGCG-1</th><th scope=col>TTGTAAGGACCTAAGT-1</th><th scope=col>TTGTAATCCGTACTCG-1</th><th scope=col>TTGTACACCTCGAACA-1</th><th scope=col>TTGTCACCGCGGTATC-1</th><th scope=col>TTGTGAACCTAATCCG-1</th><th scope=col>TTGTGAGGCATGACGC-1</th><th scope=col>TTGTGCAGCCACGTCA-1</th><th scope=col>TTGTGGTAGGAGGGAT-1</th><th scope=col>TTGTGGTATAGGTATG-1</th><th scope=col>TTGTGTATGCCACCAA-1</th><th scope=col>TTGTGTTTCCCGAAAG-1</th><th scope=col>TTGTTAGCAAATTCGA-1</th><th scope=col>TTGTTCAGTGTGCTAC-1</th><th scope=col>TTGTTCTAGATACGCT-1</th><th scope=col>TTGTTTCACATCCAGG-1</th><th scope=col>TTGTTTCATTAGTCTA-1</th><th scope=col>TTGTTTCCATACAACT-1</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>All_Immune</th><td>-0.05680737</td><td>-0.4320402</td><td>-0.3122538</td><td>-0.2949967</td><td>0.1486568</td><td>-0.3237087</td><td>-0.09053114</td><td>-0.2726288931</td><td>-0.02927195</td><td>0.2272304</td><td>-0.30533013</td><td>-0.09327672</td><td> 0.2370990181</td><td>-0.2045095</td><td>-0.4812392</td><td>-0.32047356</td><td>-0.35327980</td><td>0.4838026</td><td>-0.2073011</td><td>-0.03469408</td><td>-0.2743828</td><td>-0.3542925</td><td>-0.09232902</td><td> 0.43852279</td><td>0.05718737</td><td>-0.01091304</td><td> 0.238856534</td><td> 0.29074356</td><td> 0.1345291</td><td>-0.21412090</td><td>-0.1148071</td><td>-0.14994788</td><td> 0.04422197</td><td>-0.2702715</td><td>-0.1906549</td><td>-0.48872274</td><td>-0.005195924</td><td> 0.07597530</td><td> 0.46371797</td><td> 0.006989429</td><td>-0.13753233</td><td> 0.1763031</td><td>-0.1672269</td><td>0.4987293</td><td>-0.04657804</td><td>-0.1632526</td><td>-0.06586516</td><td>-0.00495179</td><td>-0.33884358</td><td>-0.4092121</td><td>⋯</td><td> 0.1260934</td><td>-0.2177819</td><td>-0.05243433</td><td>-0.18632366</td><td>-0.1864494</td><td> 0.2934724</td><td> 0.007656173</td><td> 0.21899751</td><td> 0.226152554</td><td>-0.02608061</td><td>0.1076984</td><td>-0.07136198</td><td>-0.1725218</td><td>-0.3075612</td><td>0.2161298</td><td>-0.3635338</td><td>-0.02181073</td><td>-0.1857171</td><td> 0.11989717</td><td>0.1846133</td><td>-0.2591182</td><td>-0.2463493</td><td> 0.1819760</td><td>0.16433917</td><td>0.02919652</td><td>-0.11172757</td><td>-0.2451591</td><td>-0.2913795</td><td> 0.04997410</td><td>-0.2995629</td><td>0.1113056</td><td>0.2536364</td><td>-0.4829210</td><td> 0.10798115</td><td>-0.05319885</td><td> 0.1184411</td><td> 0.45123804</td><td>-0.04169607</td><td>-0.1495193</td><td>0.4031398</td><td>-0.2522672</td><td>-0.3097344</td><td>0.2593603</td><td>0.3800144</td><td>-0.03631398</td><td>-0.1106731</td><td>-0.1716187</td><td>-0.2348298</td><td>0.2165386</td><td>-0.12441889</td></tr>
	<tr><th scope=row>Angiogenesis</th><td>-0.25984750</td><td>-0.3636956</td><td>-0.1422460</td><td>-0.1592542</td><td>0.6701079</td><td>-0.4711628</td><td>-0.17020046</td><td>-0.0003359958</td><td>-0.04696236</td><td>0.4667481</td><td>-0.04242067</td><td> 0.11057705</td><td>-0.0734384262</td><td>-0.1535245</td><td>-0.2085669</td><td>-0.06393342</td><td>-0.04797819</td><td>0.2628455</td><td>-0.4874150</td><td>-0.01397154</td><td>-0.1809975</td><td>-0.3770225</td><td> 0.01444596</td><td> 0.45600551</td><td>0.11433245</td><td> 0.61165446</td><td>-0.001063347</td><td>-0.03768425</td><td>-0.2768315</td><td> 0.01258583</td><td> 0.2140868</td><td>-0.07065291</td><td>-0.22365934</td><td>-0.4863134</td><td>-0.2287146</td><td> 0.02370807</td><td>-0.240516559</td><td>-0.04283355</td><td> 0.69073273</td><td> 0.396864672</td><td> 0.01099594</td><td> 0.3814791</td><td>-0.3625070</td><td>0.7135791</td><td> 0.27613940</td><td>-0.3672031</td><td> 0.14704333</td><td> 0.23477637</td><td> 0.09344184</td><td> 0.3074853</td><td>⋯</td><td> 0.4925855</td><td>-0.2999384</td><td>-0.15864633</td><td>-0.12685964</td><td>-0.1202967</td><td>-0.0272748</td><td>-0.137567244</td><td> 0.37578126</td><td>-0.111272160</td><td> 0.32920329</td><td>0.1995298</td><td> 0.38889971</td><td> 0.2968146</td><td>-0.3646225</td><td>0.4377732</td><td>-0.4556723</td><td> 0.33623210</td><td>-0.1197518</td><td>-0.42917881</td><td>0.5180128</td><td> 0.2405660</td><td> 0.3213159</td><td>-0.3760181</td><td>0.08873116</td><td>0.02769650</td><td> 0.22182691</td><td>-0.1239461</td><td> 0.2137728</td><td> 0.47897057</td><td>-0.2304918</td><td>0.1725143</td><td>0.3445495</td><td>-0.4385952</td><td>-0.04535517</td><td>-0.24804941</td><td>-0.3532265</td><td>-0.04401931</td><td> 0.41844767</td><td>-0.5120538</td><td>0.2153092</td><td> 0.1425695</td><td>-0.5727502</td><td>0.3149629</td><td>0.1899648</td><td>-0.19324511</td><td> 0.1090272</td><td>-0.3369768</td><td>-0.1661996</td><td>0.1390813</td><td>-0.04127355</td></tr>
	<tr><th scope=row>AyersIFNG</th><td> 0.29674831</td><td>-0.6656103</td><td> 0.2159743</td><td>-0.6246499</td><td>0.1810290</td><td>-0.6876751</td><td>-0.62730316</td><td>-0.7877663981</td><td> 0.74286043</td><td>0.8028823</td><td>-0.61464586</td><td>-0.55212085</td><td> 0.0002946493</td><td>-0.1950133</td><td>-0.6573680</td><td>-0.77626174</td><td> 0.46046916</td><td>0.3165642</td><td>-0.7180761</td><td>-0.11215213</td><td>-0.2624282</td><td>-0.6690676</td><td>-0.15508890</td><td>-0.01051954</td><td>0.13230134</td><td> 0.19610502</td><td> 0.145079087</td><td> 0.60961239</td><td>-0.0628416</td><td>-0.03534975</td><td> 0.2532034</td><td>-0.63565426</td><td> 0.05291134</td><td> 0.3533220</td><td>-0.6528295</td><td>-0.66446579</td><td>-0.011548229</td><td> 0.58643835</td><td>-0.08183722</td><td>-0.493297319</td><td>-0.69383627</td><td>-0.6529677</td><td>-0.6557881</td><td>0.7807610</td><td> 0.38407167</td><td> 0.4766335</td><td> 0.26448510</td><td> 0.01186199</td><td>-0.67936046</td><td>-0.6913253</td><td>⋯</td><td>-0.4287446</td><td>-0.6838169</td><td>-0.60686667</td><td>-0.09510901</td><td>-0.4422309</td><td> 0.5527008</td><td>-0.002440988</td><td>-0.00010004</td><td> 0.001859777</td><td>-0.69346534</td><td>0.1049262</td><td> 0.59756382</td><td>-0.2087985</td><td> 0.0104767</td><td>0.1177903</td><td>-0.6728691</td><td> 0.66883956</td><td>-0.1125032</td><td> 0.04726975</td><td>0.3238978</td><td>-0.5029012</td><td>-0.6826700</td><td>-0.6563273</td><td>0.13081036</td><td>0.02402746</td><td> 0.04090202</td><td>-0.2879813</td><td> 0.1530709</td><td>-0.04214743</td><td>-0.6214486</td><td>0.7017851</td><td>0.1144523</td><td>-0.3177399</td><td> 0.59137459</td><td> 0.60564618</td><td> 0.8377351</td><td> 0.63372685</td><td> 0.13757917</td><td>-0.3816552</td><td>0.6703940</td><td>-0.1803174</td><td>-0.7806919</td><td>0.1623660</td><td>0.7509149</td><td> 0.30125482</td><td>-0.6524792</td><td>-0.6524610</td><td>-0.5788132</td><td>0.6862709</td><td>-0.59841081</td></tr>
</tbody>
</table>




```R
Endo_gsva_VB1 %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 2589</caption>
<thead>
	<tr><th></th><th scope=col>AAACAACGAATAGTTC-1</th><th scope=col>AAACAAGTATCTCCCA-1</th><th scope=col>AAACACCAATAACTGC-1</th><th scope=col>AAACAGAGCGACTCCT-1</th><th scope=col>AAACAGCTTTCAGAAG-1</th><th scope=col>AAACAGGGTCTATATT-1</th><th scope=col>AAACATGGTGAGAGGA-1</th><th scope=col>AAACCCGAACGAAATC-1</th><th scope=col>AAACCGGGTAGGTACC-1</th><th scope=col>AAACCGTTCGTCCAGG-1</th><th scope=col>AAACCTCATGAAGTTG-1</th><th scope=col>AAACGAAGATGGAGTA-1</th><th scope=col>AAACGAGACGGTTGAT-1</th><th scope=col>AAACGGGTTGGTATCC-1</th><th scope=col>AAACTCGTGATATAAG-1</th><th scope=col>AAACTTAATTGCACGC-1</th><th scope=col>AAACTTGCAAACGTAT-1</th><th scope=col>AAAGAATGACCTTAGA-1</th><th scope=col>AAAGACTGGGCGCTTT-1</th><th scope=col>AAAGCTTGCCTACATA-1</th><th scope=col>AAAGGCCCTATAATAC-1</th><th scope=col>AAAGGCTACGGACCAT-1</th><th scope=col>AAAGGGCAGCTTGAAT-1</th><th scope=col>AAAGTAGCATTGCTCA-1</th><th scope=col>AAAGTCGACCCTCAGT-1</th><th scope=col>AAAGTGTGATTTATCT-1</th><th scope=col>AAAGTTGACTCCCGTA-1</th><th scope=col>AAATAACCATACGGGA-1</th><th scope=col>AAATACCTATAAGCAT-1</th><th scope=col>AAATAGGGTGCTATTG-1</th><th scope=col>AAATCGCGGAAGGAGT-1</th><th scope=col>AAATCTAGCCCTGCTA-1</th><th scope=col>AAATGGCCCGTGCCCT-1</th><th scope=col>AAATGTATCTTATCCC-1</th><th scope=col>AAATTAATAAGCGCGA-1</th><th scope=col>AAATTACACGACTCTG-1</th><th scope=col>AAATTACCTATCGATG-1</th><th scope=col>AAATTCCAGGTCCAAA-1</th><th scope=col>AAATTGATAGTCCTTT-1</th><th scope=col>AAATTGCGGCGGTTCT-1</th><th scope=col>AAATTGGTGAGAAGCA-1</th><th scope=col>AAATTTGCGGGTGTGG-1</th><th scope=col>AACAACTGGTAGTTGC-1</th><th scope=col>AACAATACATTGTCGA-1</th><th scope=col>AACAATTACTCTACGC-1</th><th scope=col>AACACACGCTCGCCGC-1</th><th scope=col>AACACGAGACGCGGCC-1</th><th scope=col>AACAGCTGTGTGGCAA-1</th><th scope=col>AACAGGATGGGCCGCG-1</th><th scope=col>AACAGGTAGTATGGAT-1</th><th scope=col>⋯</th><th scope=col>TTGACCATGTTCTCCG-1</th><th scope=col>TTGACCGTGTTAATGA-1</th><th scope=col>TTGACGCTCCATGAGC-1</th><th scope=col>TTGACTACCATATGGT-1</th><th scope=col>TTGACTATTGTCCGGC-1</th><th scope=col>TTGAGAAGTTTAGCAT-1</th><th scope=col>TTGAGAGTACTGCTAA-1</th><th scope=col>TTGAGCAGCCCACGGT-1</th><th scope=col>TTGAGCGCCACGTGAT-1</th><th scope=col>TTGAGTCCCGCTGCTG-1</th><th scope=col>TTGATCTAACTTTGTC-1</th><th scope=col>TTGATGTGTAGTCCCG-1</th><th scope=col>TTGATTAGCTGTTTCT-1</th><th scope=col>TTGATTATGCAGATGA-1</th><th scope=col>TTGCATGCTGATCACG-1</th><th scope=col>TTGCCAAGCAGAACCC-1</th><th scope=col>TTGCCATAGCCCGCTC-1</th><th scope=col>TTGCCCTGATCACGGG-1</th><th scope=col>TTGCCGCAGACCTACA-1</th><th scope=col>TTGCGCTCTCTCGCTT-1</th><th scope=col>TTGCGGCATCAGAAAG-1</th><th scope=col>TTGCGTAGTTTGAGGA-1</th><th scope=col>TTGCGTCGGCCAACCG-1</th><th scope=col>TTGCTGCACCTATCCA-1</th><th scope=col>TTGGAAGAATACAGTC-1</th><th scope=col>TTGGACATGTGGCTTA-1</th><th scope=col>TTGGACCATCTGGCAA-1</th><th scope=col>TTGGACCTATAACAGT-1</th><th scope=col>TTGGCCATCTTGCGCT-1</th><th scope=col>TTGGCCTAGAATTTCG-1</th><th scope=col>TTGGCGATCCGAATAT-1</th><th scope=col>TTGGGACGTAAGAGTT-1</th><th scope=col>TTGGGCGGCGGTTGCC-1</th><th scope=col>TTGGTCACACTCGTAA-1</th><th scope=col>TTGGTTCGCTCAAAGG-1</th><th scope=col>TTGGTTGCGGTGCGCG-1</th><th scope=col>TTGTAAGGACCTAAGT-1</th><th scope=col>TTGTAAGGCCAGTTGG-1</th><th scope=col>TTGTACACCTCGAACA-1</th><th scope=col>TTGTCACCGCGGTATC-1</th><th scope=col>TTGTGAACCTAATCCG-1</th><th scope=col>TTGTGAGGCATGACGC-1</th><th scope=col>TTGTGCGGAAGCGGAT-1</th><th scope=col>TTGTGGTAGGAGGGAT-1</th><th scope=col>TTGTGGTATAGGTATG-1</th><th scope=col>TTGTTCTAGATACGCT-1</th><th scope=col>TTGTTGTGTGTCAAGA-1</th><th scope=col>TTGTTTCACATCCAGG-1</th><th scope=col>TTGTTTCATTAGTCTA-1</th><th scope=col>TTGTTTCCATACAACT-1</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>⋯</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>All_Immune</th><td>-0.40919006</td><td> 0.1248427</td><td>-0.2232842</td><td> 0.05892970</td><td>-0.3376218</td><td>-0.2276048</td><td>0.150023760</td><td>-0.064940679</td><td>-0.02756923</td><td> 0.37183500</td><td>-0.2437399</td><td>0.13869075</td><td>-0.1418207</td><td>0.24986994</td><td>-0.09196317</td><td>-0.3382891</td><td> 0.2758066</td><td> 0.2307397</td><td>0.3744783</td><td>-0.13036900</td><td>-0.3850378</td><td>0.5381396</td><td> 0.004713382</td><td>0.3145202</td><td>-0.09439654</td><td> 0.1329558</td><td>-0.20991742</td><td> 0.27584561</td><td>-0.1097791</td><td> 0.08518754</td><td>-0.3941738</td><td>-0.12179428</td><td>-0.002110972</td><td> 0.16536854</td><td>-0.4184972</td><td> 0.2080019</td><td>-0.35397539</td><td>-0.09456265</td><td> 0.02705141</td><td>-0.03273352</td><td>-0.001946636</td><td>-0.03811018</td><td>-0.1322734</td><td>0.3881608</td><td> 0.04797777</td><td>-0.3907203</td><td>-0.2847535</td><td>0.2956885</td><td>-0.02232197</td><td> 0.1240692</td><td>⋯</td><td>-0.2221356</td><td>-0.06358689</td><td>-0.19617519</td><td> 0.1700201</td><td>-0.06619192</td><td>-0.1446155</td><td>-0.3209004</td><td>-0.2057163</td><td>-0.1182971</td><td>0.08331947</td><td>-0.2036227</td><td>-0.3809979</td><td> 0.32525563</td><td> 0.05790675</td><td>-0.03960102</td><td>-0.01721894</td><td>0.06612352</td><td>-0.2785050</td><td>-0.4159777</td><td>-0.05612973</td><td>0.3624167</td><td> 0.3665413</td><td>0.2060093</td><td>0.12799406</td><td>-0.36692701</td><td>0.3230871</td><td>0.007297871</td><td>-0.1040083</td><td>-0.04166428</td><td>-0.4272843</td><td> 0.1134232</td><td> 0.1339973</td><td>-0.1298473</td><td>0.3075708</td><td> 0.13946384</td><td>-0.3207854</td><td>-0.1742502</td><td> 0.1140986</td><td>0.1189916</td><td>-0.4418161</td><td> 0.04587000</td><td>-0.2617494</td><td>-0.2568759</td><td>-0.3723166</td><td>-0.2127830</td><td>-0.02391508</td><td>-0.2477704</td><td>-0.3752086</td><td>0.2921108</td><td>-0.09576479</td></tr>
	<tr><th scope=row>Angiogenesis</th><td>-0.29343431</td><td> 0.5082822</td><td> 0.1240651</td><td> 0.01065969</td><td>-0.3232331</td><td>-0.4420861</td><td>0.006480752</td><td>-0.289911508</td><td>-0.21551481</td><td> 0.03019649</td><td>-0.2751823</td><td>0.06400930</td><td>-0.4947367</td><td>0.02270643</td><td>-0.22439359</td><td>-0.3397088</td><td>-0.2587877</td><td>-0.2956168</td><td>0.6900810</td><td>-0.01675589</td><td>-0.5698917</td><td>0.2728353</td><td> 0.539241691</td><td>0.6375109</td><td> 0.07068666</td><td>-0.1024371</td><td>-0.08228291</td><td>-0.05503112</td><td>-0.3108076</td><td>-0.06966676</td><td>-0.6143806</td><td>-0.04969817</td><td> 0.121563709</td><td> 0.05065225</td><td>-0.2775049</td><td>-0.4287338</td><td> 0.04362570</td><td>-0.45617468</td><td>-0.22058037</td><td>-0.38629812</td><td> 0.055442844</td><td>-0.25790735</td><td>-0.3914054</td><td>0.5630306</td><td>-0.09846821</td><td>-0.2008856</td><td>-0.2997329</td><td>0.1840702</td><td> 0.09683285</td><td>-0.1847827</td><td>⋯</td><td>-0.3891204</td><td> 0.17052181</td><td>-0.09495074</td><td>-0.1744794</td><td> 0.58765200</td><td>-0.3751693</td><td> 0.1842594</td><td>-0.4827688</td><td>-0.1960242</td><td>0.23315471</td><td>-0.2259431</td><td>-0.3010001</td><td>-0.03315458</td><td>-0.34345912</td><td>-0.28109768</td><td>-0.33129724</td><td>0.23152291</td><td>-0.1428820</td><td> 0.1018331</td><td>-0.34126208</td><td>0.5326904</td><td>-0.6005726</td><td>0.3175016</td><td>0.04804074</td><td> 0.08302104</td><td>0.3083732</td><td>0.380427701</td><td>-0.5336194</td><td> 0.66093932</td><td>-0.2005885</td><td>-0.2213533</td><td>-0.3282462</td><td> 0.4557656</td><td>0.3209752</td><td> 0.26010423</td><td>-0.3384646</td><td> 0.2394631</td><td> 0.5359387</td><td>0.3672545</td><td>-0.5035941</td><td> 0.07933361</td><td>-0.3412286</td><td>-0.5237940</td><td>-0.3509303</td><td>-0.2196897</td><td>-0.53662302</td><td>-0.5375393</td><td>-0.3140421</td><td>0.3118755</td><td>-0.01556957</td></tr>
	<tr><th scope=row>AyersIFNG</th><td> 0.01579317</td><td>-0.3709040</td><td> 0.1482789</td><td>-0.26579646</td><td>-0.5166586</td><td>-0.4080378</td><td>0.800607953</td><td> 0.007074506</td><td> 0.29402565</td><td>-0.25074055</td><td> 0.5782802</td><td>0.08324654</td><td>-0.6375188</td><td>0.10064837</td><td> 0.18461169</td><td>-0.7427456</td><td> 0.5715500</td><td>-0.6164082</td><td>0.8485634</td><td>-0.14285954</td><td>-0.4746703</td><td>0.7923441</td><td>-0.409756897</td><td>0.6465311</td><td>-0.09175987</td><td>-0.4151426</td><td>-0.39498303</td><td> 0.76952471</td><td> 0.5267728</td><td> 0.25085320</td><td>-0.5229271</td><td> 0.05089711</td><td> 0.227026621</td><td>-0.26677553</td><td>-0.4956242</td><td>-0.2090104</td><td>-0.09980889</td><td>-0.26192290</td><td> 0.46049548</td><td> 0.02158198</td><td>-0.310288440</td><td> 0.15136781</td><td>-0.3309115</td><td>0.4191018</td><td> 0.13555177</td><td>-0.2766496</td><td>-0.4195038</td><td>0.4757149</td><td> 0.56814493</td><td>-0.1229787</td><td>⋯</td><td>-0.2773975</td><td>-0.40410594</td><td>-0.32564980</td><td> 0.2761680</td><td>-0.40345045</td><td> 0.5207036</td><td>-0.3577300</td><td>-0.6499889</td><td>-0.2282238</td><td>0.38590808</td><td> 0.1235393</td><td>-0.7914832</td><td> 0.29609071</td><td>-0.33189761</td><td> 0.22174787</td><td>-0.33851308</td><td>0.24292621</td><td>-0.7424374</td><td>-0.1705314</td><td> 0.30723262</td><td>0.7909079</td><td> 0.4155481</td><td>0.3731621</td><td>0.63091866</td><td>-0.32842575</td><td>0.7209236</td><td>0.684416915</td><td> 0.3247352</td><td>-0.44012555</td><td>-0.3828995</td><td>-0.3348807</td><td> 0.7438662</td><td>-0.3467305</td><td>0.8044889</td><td>-0.07524357</td><td> 0.1480547</td><td> 0.5779482</td><td>-0.3349577</td><td>0.5258457</td><td>-0.7448307</td><td>-0.34016890</td><td>-0.7453827</td><td> 0.2523891</td><td> 0.6028686</td><td>-0.6305304</td><td> 0.10041967</td><td>-0.3117132</td><td>-0.3209406</td><td>0.1411895</td><td> 0.20544283</td></tr>
</tbody>
</table>




```R
# transpose so spot names are rownames and signatures are column names
Endo_gsva_VA2 <- Endo_gsva_VA2 %>% t() %>% as.data.frame() %>% setNames(rownames(Endo_gsva_VA2))
dim(Endo_gsva_VA2)
Endo_gsva_VA2 %>% head(3)

Endo_gsva_VA1 <- Endo_gsva_VA1 %>% t() %>% as.data.frame() %>% setNames(rownames(Endo_gsva_VA1))
dim(Endo_gsva_VA1)
Endo_gsva_VA1 %>% head(3)

Endo_gsva_VB1 <- Endo_gsva_VB1 %>% t() %>% as.data.frame() %>% setNames(rownames(Endo_gsva_VB1))
dim(Endo_gsva_VB1)
Endo_gsva_VB1 %>% head(3)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>4316</li><li>30</li></ol>




<table class="dataframe">
<caption>A data.frame: 3 × 30</caption>
<thead>
	<tr><th></th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACAAGTATCTCCCA-1</th><td>-0.4509680</td><td>-0.4483304</td><td>-0.6914766</td><td>-0.3743552</td><td>-0.6934548</td><td>-0.8061806</td><td>-0.4244671</td><td>-0.2619432</td><td>-0.6846815</td><td>0.2752691</td><td>-0.7654531</td><td>-0.7305166</td><td>-0.1999233</td><td>-0.1339729</td><td>-0.2137573</td><td>-0.2902687</td><td>0.02930414</td><td>-0.1282588</td><td>-0.2612448</td><td>-0.1995247</td><td>-0.7413616</td><td>-0.3953259</td><td>-0.2975777</td><td> 0.1823718</td><td>-0.5831718</td><td> 0.484040526</td><td>-0.3325986</td><td>-0.7049820</td><td>-0.1999233</td><td>-0.1565646</td></tr>
	<tr><th scope=row>AAACACCAATAACTGC-1</th><td> 0.4657553</td><td> 0.4768693</td><td> 0.8017860</td><td> 0.6856195</td><td> 0.7281010</td><td>-0.4407441</td><td>-0.3554181</td><td> 0.4510906</td><td> 0.6792392</td><td>0.4830992</td><td> 0.5655572</td><td> 0.5727547</td><td> 0.3014952</td><td> 0.0995948</td><td> 0.1687630</td><td> 0.3944400</td><td>0.27930871</td><td> 0.5291796</td><td> 0.3967318</td><td> 0.6683280</td><td> 0.5340971</td><td> 0.3684659</td><td> 0.6188102</td><td> 0.6471133</td><td> 0.1671480</td><td>-0.005985246</td><td> 0.4556002</td><td> 0.6415331</td><td> 0.3014952</td><td> 0.1536071</td></tr>
	<tr><th scope=row>AAACAGAGCGACTCCT-1</th><td>-0.5878032</td><td>-0.6471249</td><td>-0.2034801</td><td>-0.6089440</td><td>-0.6035510</td><td>-0.8401840</td><td>-0.3283846</td><td>-0.5893354</td><td>-0.6599931</td><td>0.2542799</td><td>-0.7992599</td><td>-0.7901462</td><td>-0.4414986</td><td>-0.2771161</td><td>-0.3490007</td><td>-0.6546928</td><td>0.15123322</td><td>-0.6707561</td><td>-0.5942471</td><td>-0.7829926</td><td>-0.6668816</td><td>-0.4782615</td><td>-0.6769066</td><td>-0.8178636</td><td>-0.6608763</td><td> 0.140619235</td><td>-0.7686918</td><td>-0.7469988</td><td>-0.4414986</td><td>-0.1267453</td></tr>
</tbody>
</table>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>2806</li><li>30</li></ol>




<table class="dataframe">
<caption>A data.frame: 3 × 30</caption>
<thead>
	<tr><th></th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACAACGAATAGTTC-1</th><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td></tr>
	<tr><th scope=row>AAACAAGTATCTCCCA-1</th><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td></tr>
	<tr><th scope=row>AAACAATCTACTAGCA-1</th><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td></tr>
</tbody>
</table>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>2589</li><li>30</li></ol>




<table class="dataframe">
<caption>A data.frame: 3 × 30</caption>
<thead>
	<tr><th></th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACAACGAATAGTTC-1</th><td>-0.4091901</td><td>-0.2934343</td><td> 0.01579317</td><td>-0.6232956</td><td>-0.2943380</td><td>-0.5389701</td><td> 0.1609967</td><td>-0.44028788</td><td>-0.04994617</td><td>-0.8102586</td><td>-0.2194454</td><td>-0.6907327</td><td>-0.18496858</td><td>0.01012297</td><td>0.00174647</td><td>0.52323585</td><td>-0.3835384</td><td> 0.04831849</td><td> 0.2639777</td><td>-0.4365621</td><td>-0.3175423</td><td>-0.18950681</td><td> 0.1715214</td><td>-0.5652826</td><td>-0.4481054</td><td>-0.3094859</td><td>-0.3636408</td><td>-0.3739748</td><td>-0.18496858</td><td>0.09003211</td></tr>
	<tr><th scope=row>AAACAAGTATCTCCCA-1</th><td> 0.1248427</td><td> 0.5082822</td><td>-0.37090396</td><td>-0.2081911</td><td>-0.1350080</td><td> 0.2548888</td><td> 0.3546991</td><td> 0.61804117</td><td> 0.17432327</td><td>-0.7612695</td><td> 0.6448965</td><td> 0.1331139</td><td> 0.22323805</td><td>0.06597860</td><td>0.01068036</td><td>0.64566927</td><td> 0.9555956</td><td>-0.32369132</td><td>-0.3258627</td><td> 0.1627362</td><td> 0.2594864</td><td>-0.09579307</td><td>-0.2807481</td><td>-0.6112452</td><td> 0.4179383</td><td>-0.3445947</td><td>-0.3024936</td><td>-0.2643529</td><td> 0.22323805</td><td>0.69416330</td></tr>
	<tr><th scope=row>AAACACCAATAACTGC-1</th><td>-0.2232842</td><td> 0.1240651</td><td> 0.14827892</td><td>-0.4825369</td><td>-0.2947074</td><td>-0.3440408</td><td>-0.2352929</td><td> 0.08386593</td><td>-0.56773770</td><td> 0.1338223</td><td>-0.2392742</td><td> 0.1998769</td><td>-0.01450242</td><td>0.04427482</td><td>0.10580519</td><td>0.06204871</td><td>-0.3869387</td><td>-0.07191923</td><td>-0.4880361</td><td>-0.1204101</td><td>-0.4476507</td><td> 0.07267653</td><td> 0.3163856</td><td> 0.4089414</td><td> 0.2907777</td><td>-0.2117291</td><td>-0.3645462</td><td>-0.3794759</td><td>-0.01450242</td><td>0.14351748</td></tr>
</tbody>
</table>




```R
length(colnames(Endo_gsva_VA2))
length(colnames(Endo_gsva_VA1))
length(colnames(Endo_gsva_VB1))
colnames(Endo_gsva_VA2)
colnames(Endo_gsva_VA1)
colnames(Endo_gsva_VB1)
```


30



30



30



<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'All_Immune'</li><li>'Angiogenesis'</li><li>'AyersIFNG'</li><li>'B_Cell'</li><li>'CD8_T_Cell'</li><li>'DC'</li><li>'Endothelial_Activation'</li><li>'Endothelial_Cell'</li><li>'Endothelial_Chemokines'</li><li>'Epithelial_cell'</li><li>'Exhaustion'</li><li>'Fibroblast'</li><li>'GOBP_LEUKOCYTE_ADHESION_TO_VASC'</li><li>'GOBP_LEUKOCYTE_MIGRATION'</li><li>'GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION'</li><li>'ICB_Targets'</li><li>'LEC'</li><li>'M1_Macrophage'</li><li>'M1_Macrophage_old'</li><li>'M2_Macrophage'</li><li>'M2_Macrophage_old'</li><li>'MN4_EC_Phenotype'</li><li>'MN4_EC_Phenotype_Top30'</li><li>'Macrophage'</li><li>'NK_Cell'</li><li>'Proliferation'</li><li>'T_Cell'</li><li>'T_Reg'</li><li>'Upregulated_by_2_3_CGAMP'</li><li>'Upregulated_by_2_3_CGAMP_IFNb_OVERLAP'</li></ol>



# Combine GSVA results from different Tissue Slices back together


```R
# Need to add Tissue Slice back to the Visium Spot Barcodes
Endo_gsva_VA1 %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 30</caption>
<thead>
	<tr><th></th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACAACGAATAGTTC-1</th><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.3752852</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td></tr>
	<tr><th scope=row>AAACAAGTATCTCCCA-1</th><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.6186324</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td></tr>
</tbody>
</table>




```R
Endo_gsva_VA2$Barcode1 <- rownames(Endo_gsva_VA2)
Endo_gsva_VA1$Barcode1 <- rownames(Endo_gsva_VA1)
Endo_gsva_VB1$Barcode1 <- rownames(Endo_gsva_VB1)
```


```R
Endo_gsva_VA2$Barcode2 <- paste0('VA2_', Endo_gsva_VA2$Barcode1)
Endo_gsva_VA1$Barcode2 <- paste0('VA1_', Endo_gsva_VA1$Barcode1)
Endo_gsva_VB1$Barcode2 <- paste0('VB1_', Endo_gsva_VB1$Barcode1)
```


```R
rownames(Endo_gsva_VA2) <- Endo_gsva_VA2$Barcode2
rownames(Endo_gsva_VA1) <- Endo_gsva_VA1$Barcode2
rownames(Endo_gsva_VB1) <- Endo_gsva_VB1$Barcode2
```


```R
Endo_gsva_VA1 %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 32</caption>
<thead>
	<tr><th></th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>Barcode1</th><th scope=col>Barcode2</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.3752852</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>AAACAACGAATAGTTC-1</td><td>VA1_AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.6186324</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>AAACAAGTATCTCCCA-1</td><td>VA1_AAACAAGTATCTCCCA-1</td></tr>
</tbody>
</table>




```R
Endo_gsva <- rbind(Endo_gsva_VA1, Endo_gsva_VB1, Endo_gsva_VA2)
```


```R
# drop barcode columns
Endo_gsva <- Endo_gsva %>% select(-c(Barcode1,Barcode2))
```


```R
dim(Endo_gsva)
Endo_gsva %>% head(3)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>30</li></ol>




<table class="dataframe">
<caption>A data.frame: 3 × 30</caption>
<thead>
	<tr><th></th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td></tr>
</tbody>
</table>




```R
# save outputs
#write.csv(Endo_gsva, "Output/05_do_GSVA_by_Tissue_Slice/202504025_Updated_GeneSets_GSVA_by_Tissue_Slice.csv")
write.csv(Endo_gsva, "Output/05_do_GSVA_by_Tissue_Slice/202504025_Updated_GeneSets_GSVA_by_Tissue_Slice.csv")
```
colnames(Endo_gsva)

```R

```

# Load GSVA results if needed


```R
#Endo_gsva <- read.csv("Output/05_do_GSVA_by_Tissue_Slice/202504025_Updated_GeneSets_GSVA_by_Tissue_Slice.csv")
Endo_gsva <- read.csv("Output/05_do_GSVA_by_Tissue_Slice/202504025_Updated_GeneSets_GSVA_by_Tissue_Slice.csv")
```


```R
dim(Endo_gsva)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>31</li></ol>




```R
Endo_gsva %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 31</caption>
<thead>
	<tr><th></th><th scope=col>X</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td></tr>
	<tr><th scope=row>2</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td></tr>
	<tr><th scope=row>3</th><td>VA1_AAACAATCTACTAGCA-1</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td></tr>
</tbody>
</table>




```R
# save current colnames
genesets_colnames <- colnames(Endo_gsva)
# save visium spot as a column
Endo_gsva$Visium_spot <- Endo_gsva$X
# set rownames to visium spot
rownames(Endo_gsva) <- Endo_gsva$X
# remove the X column
Endo_gsva <- Endo_gsva %>% select(-c(X))
# check
Endo_gsva %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 31</caption>
<thead>
	<tr><th></th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>Visium_spot</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>VA1_AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>VA1_AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>VA1_AAACAATCTACTAGCA-1</td></tr>
</tbody>
</table>



# Melt data for plotting


```R
Endo_gsva_melt <- Endo_gsva %>% pivot_longer(!Visium_spot, names_to = "GeneSet", values_to = "GSVA_enrichment")
Endo_gsva_melt %>% head(10)
```


<table class="dataframe">
<caption>A tibble: 10 × 3</caption>
<thead>
	<tr><th scope=col>Visium_spot</th><th scope=col>GeneSet</th><th scope=col>GSVA_enrichment</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>All_Immune            </td><td>-0.056807367</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>Angiogenesis          </td><td>-0.259847496</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>AyersIFNG             </td><td> 0.296748306</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>B_Cell                </td><td>-0.038430718</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>CD8_T_Cell            </td><td> 0.513526640</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>DC                    </td><td>-0.007136273</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>Endothelial_Activation</td><td> 0.412195890</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>Endothelial_Cell      </td><td> 0.404656780</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>Endothelial_Chemokines</td><td>-0.158465265</td></tr>
	<tr><td>VA1_AAACAACGAATAGTTC-1</td><td>Epithelial_cell       </td><td>-0.268586268</td></tr>
</tbody>
</table>




```R
options(repr.plot.width=20, repr.plot.height=12)
```


```R
Endo_gsva_plot1 <- ggplot(Endo_gsva_melt, aes(x=GeneSet, y=GSVA_enrichment)) + 
  geom_violin(width=1.5) + 
  theme(text=element_text(size=20)) + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1, size = 15)) + 
  labs(title = "Endothelial GeneSet Enrichment, All Visium Spots") + 
  theme(plot.title = element_text(size=22))
Endo_gsva_plot1
```

    Warning message:
    “[1m[22m`position_dodge()` requires non-overlapping [32mx[39m intervals.”



    
![png](output_87_1.png)
    


# Check correlation between gene sets


```R
# prep
Endo_gsva2 <- Endo_gsva %>% select(-c(Visium_spot))
#rownames(Endo_gsva2) <- Endo_gsva2$GeneSet

# computer correlation matrix
Endo_gsva_corr <- cor(Endo_gsva2)
```


```R
Endo_gsva_corr[1:5,1:5]
```


<table class="dataframe">
<caption>A matrix: 5 × 5 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th></tr>
</thead>
<tbody>
	<tr><th scope=row>All_Immune</th><td>1.0000000</td><td>0.6078798</td><td>0.5094556</td><td>0.8208330</td><td>0.6535781</td></tr>
	<tr><th scope=row>Angiogenesis</th><td>0.6078798</td><td>1.0000000</td><td>0.2960592</td><td>0.4713806</td><td>0.4033234</td></tr>
	<tr><th scope=row>AyersIFNG</th><td>0.5094556</td><td>0.2960592</td><td>1.0000000</td><td>0.3336670</td><td>0.2835691</td></tr>
	<tr><th scope=row>B_Cell</th><td>0.8208330</td><td>0.4713806</td><td>0.3336670</td><td>1.0000000</td><td>0.4283263</td></tr>
	<tr><th scope=row>CD8_T_Cell</th><td>0.6535781</td><td>0.4033234</td><td>0.2835691</td><td>0.4283263</td><td>1.0000000</td></tr>
</tbody>
</table>




```R
corrplot(Endo_gsva_corr, type = "upper", tl.col = "black", method="color")
```


    
![png](output_91_0.png)
    


# Load the integrated, annotated seurat object to merge these GSVA enrichments into


```R
data <- readRDS('Processing/AB1_A2_Annotated_step2Navin_20250326.rds')
```

# Annotate the spots with seurat clusters and Navin annotations


```R
data@meta.data[1:3,]
```


<table class="dataframe">
<caption>A data.frame: 3 × 27</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Barcode</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>VA1_AAACAACGAATAGTTC-1</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>VA1_AAACAAGTATCTCCCA-1</td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>VA1_AAACAATCTACTAGCA-1</td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td></tr>
</tbody>
</table>




```R
dim(data)
dim(data@meta.data)
dim(Endo_gsva)
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
<ol class=list-inline><li>9711</li><li>27</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>31</li></ol>




```R
Endo_gsva[1:3,]
```


<table class="dataframe">
<caption>A data.frame: 3 × 31</caption>
<thead>
	<tr><th></th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>Visium_spot</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td><td>VA1_AAACAACGAATAGTTC-1</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td><td>VA1_AAACAAGTATCTCCCA-1</td></tr>
	<tr><th scope=row>VA1_AAACAATCTACTAGCA-1</th><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td><td>VA1_AAACAATCTACTAGCA-1</td></tr>
</tbody>
</table>




```R
# save ordering in the metadata
barcode_ordering <- data@meta.data$Barcode
barcode_ordering[1:5]
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'VA1_AAACAACGAATAGTTC-1'</li><li>'VA1_AAACAAGTATCTCCCA-1'</li><li>'VA1_AAACAATCTACTAGCA-1'</li><li>'VA1_AAACACCAATAACTGC-1'</li><li>'VA1_AAACAGCTTTCAGAAG-1'</li></ol>




```R
# Merge new annotations with old
merged_annots <- merge(data@meta.data, Endo_gsva, by.x = 'Barcode', by.y = 'Visium_spot')
dim(merged_annots)
merged_annots %>% head(3)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>9711</li><li>57</li></ol>




<table class="dataframe">
<caption>A data.frame: 3 × 57</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.41219589</td><td> 0.40465678</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.37528524</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.18949894</td><td> 0.13921428</td><td>-0.1651271</td><td> 0.01555870</td><td>-0.75521724</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td></tr>
	<tr><th scope=row>2</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low </td><td>Tumor_low </td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low </td><td>Tumor_low          </td><td>vasc_low </td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.13358831</td><td>-0.27719093</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.61863236</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.32364399</td><td>-0.43118684</td><td>-0.7134420</td><td>-0.12259251</td><td> 0.07979950</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td></tr>
	<tr><th scope=row>3</th><td>VA1_AAACAATCTACTAGCA-1</td><td>VA1</td><td>4980</td><td>3146</td><td>round1</td><td>4711</td><td>3146</td><td>6 </td><td>6 </td><td>VA1</td><td>1.101286</td><td>vasc_low </td><td>19</td><td>Tumor_high</td><td>Tumor_high</td><td>3</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_high_MHCI_low</td><td>Tumor_high_MHCI_low</td><td>vasc_low </td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.31225376</td><td>-0.1422460</td><td> 0.2159743</td><td>-0.52287279</td><td>-0.2544987</td><td> 0.344170013</td><td>-0.09642489</td><td> 0.01578354</td><td>-0.5194661</td><td>-0.4129127</td><td>-0.2739867</td><td>-0.08923626</td><td>-0.29691635</td><td>-0.1122542</td><td>-0.1397081</td><td> 0.5741887</td><td>-0.2428486</td><td>-0.46012808</td><td>-0.3895103</td><td>-0.6018211</td><td>-0.5293438</td><td>-0.03932724</td><td> 0.05592052</td><td>-0.7073804</td><td>-0.09144675</td><td> 0.04951067</td><td>-0.3041263</td><td> 0.9243924</td><td>-0.29691635</td><td>-0.46647952</td></tr>
</tbody>
</table>




```R
# check if ordering is the same
table(data@meta.data$Barcode == merged_annots$Barcode)
```


    
    TRUE 
    9711 



```R
rownames(merged_annots) <- merged_annots$Barcode
```


```R
merged_annots %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 57</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.3752852</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.6186324</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td></tr>
</tbody>
</table>




```R
# set metadata to merged annots
data@meta.data <- merged_annots
```


```R
data@meta.data %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 57</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung      </td><td>-0.05680737</td><td>-0.2598475</td><td> 0.2967483</td><td>-0.03843072</td><td> 0.5135266</td><td>-0.007136273</td><td> 0.4121959</td><td> 0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td> 0.6906627</td><td> 0.3752852</td><td>-0.04649791</td><td> 0.1108001</td><td> 0.1021260</td><td> 0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.2823260</td><td>-0.4695504</td><td> 0.1894989</td><td> 0.1392143</td><td>-0.1651271</td><td> 0.0155587</td><td>-0.7552172</td><td> 0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td> 0.03171949</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1_AAACAAGTATCTCCCA-1</td><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td><td>SCLC                  </td><td>VA1</td><td>SCLC                  </td><td>SCLC + Effector Immune Cells (T+NK) </td><td>6_SCLC_Effector_Immune        </td><td>6_SCLC_Effector_Immune</td><td>-0.43204017</td><td>-0.3636956</td><td>-0.6656103</td><td>-0.38024920</td><td>-0.2729638</td><td>-0.616569655</td><td>-0.1335883</td><td>-0.2771909</td><td>-0.5320766</td><td>-0.6453954</td><td>-0.2841903</td><td>-0.6186324</td><td>-0.44219743</td><td>-0.2062721</td><td>-0.2929464</td><td>-0.2201101</td><td> 0.8164109</td><td>-0.48611256</td><td>-0.4222185</td><td>-0.6161525</td><td>-0.5429759</td><td>-0.3236440</td><td>-0.4311868</td><td>-0.7134420</td><td>-0.1225925</td><td> 0.0797995</td><td>-0.3480752</td><td>-0.4209421</td><td>-0.44219743</td><td>-0.48410628</td></tr>
</tbody>
</table>



# Save updated seurat object and metadata annotations


```R
write.csv(merged_annots, 'Processing/Annotations_05_do_GSVA_by_Tissue_Slice_20250425.csv')
saveRDS(data, 'Processing/AB1_A2_Annotated_05_do_GSVA_by_Tissue_Slice_20250425.rds')
```

# Now make a melted version for easy plotting 


```R
merged_annots_melt <- merge(data@meta.data, Endo_gsva_melt, by.x = 'Barcode', by.y = 'Visium_spot')
merged_annots_melt %>% head(3)
```


<table class="dataframe">
<caption>A data.frame: 3 × 59</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>GeneSet</th><th scope=col>GSVA_enrichment</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>-0.03843072</td><td>0.5135266</td><td>-0.007136273</td><td>0.4121959</td><td>0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.282326</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.1651271</td><td>0.0155587</td><td>-0.7552172</td><td>0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>All_Immune  </td><td>-0.05680737</td></tr>
	<tr><th scope=row>2</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>-0.03843072</td><td>0.5135266</td><td>-0.007136273</td><td>0.4121959</td><td>0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.282326</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.1651271</td><td>0.0155587</td><td>-0.7552172</td><td>0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>Angiogenesis</td><td>-0.25984750</td></tr>
	<tr><th scope=row>3</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>-0.03843072</td><td>0.5135266</td><td>-0.007136273</td><td>0.4121959</td><td>0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.282326</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.1651271</td><td>0.0155587</td><td>-0.7552172</td><td>0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>AyersIFNG   </td><td> 0.29674831</td></tr>
</tbody>
</table>



# save melted annot version


```R
write.csv(merged_annots_melt, 'Processing/Annotations_melted_05_do_GSVA_by_Tissue_Slice_20250425.csv')
```


```R
merged_annots_melt %>% head(4)
```


<table class="dataframe">
<caption>A data.frame: 4 × 59</caption>
<thead>
	<tr><th></th><th scope=col>Barcode</th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th><th scope=col>Navin_annotations</th><th scope=col>Tissue_Slice</th><th scope=col>Navin_annotations_simplified</th><th scope=col>Navin_Clusters1</th><th scope=col>Navin_Clusters</th><th scope=col>Navin_Clusters2</th><th scope=col>All_Immune</th><th scope=col>Angiogenesis</th><th scope=col>AyersIFNG</th><th scope=col>B_Cell</th><th scope=col>CD8_T_Cell</th><th scope=col>DC</th><th scope=col>Endothelial_Activation</th><th scope=col>Endothelial_Cell</th><th scope=col>Endothelial_Chemokines</th><th scope=col>Epithelial_cell</th><th scope=col>Exhaustion</th><th scope=col>Fibroblast</th><th scope=col>GOBP_LEUKOCYTE_ADHESION_TO_VASC</th><th scope=col>GOBP_LEUKOCYTE_MIGRATION</th><th scope=col>GOBP_POSITIVE_REGULATION_OF_LEUKOCYTE_MIGRATION</th><th scope=col>ICB_Targets</th><th scope=col>LEC</th><th scope=col>M1_Macrophage</th><th scope=col>M1_Macrophage_old</th><th scope=col>M2_Macrophage</th><th scope=col>M2_Macrophage_old</th><th scope=col>MN4_EC_Phenotype</th><th scope=col>MN4_EC_Phenotype_Top30</th><th scope=col>Macrophage</th><th scope=col>NK_Cell</th><th scope=col>Proliferation</th><th scope=col>T_Cell</th><th scope=col>T_Reg</th><th scope=col>Upregulated_by_2_3_CGAMP</th><th scope=col>Upregulated_by_2_3_CGAMP_IFNb_OVERLAP</th><th scope=col>GeneSet</th><th scope=col>GSVA_enrichment</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>-0.03843072</td><td>0.5135266</td><td>-0.007136273</td><td>0.4121959</td><td>0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.282326</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.1651271</td><td>0.0155587</td><td>-0.7552172</td><td>0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>All_Immune  </td><td>-0.05680737</td></tr>
	<tr><th scope=row>2</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>-0.03843072</td><td>0.5135266</td><td>-0.007136273</td><td>0.4121959</td><td>0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.282326</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.1651271</td><td>0.0155587</td><td>-0.7552172</td><td>0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>Angiogenesis</td><td>-0.25984750</td></tr>
	<tr><th scope=row>3</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>-0.03843072</td><td>0.5135266</td><td>-0.007136273</td><td>0.4121959</td><td>0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.282326</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.1651271</td><td>0.0155587</td><td>-0.7552172</td><td>0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>AyersIFNG   </td><td> 0.29674831</td></tr>
	<tr><th scope=row>4</th><td>VA1_AAACAACGAATAGTTC-1</td><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td><td>Respiratory epithelium</td><td>VA1</td><td>Respiratory epithelium</td><td>Respiratory Epithelium (Normal Lung)</td><td>10_Respiratory_Epithelium_Lung</td><td>7_10_Normal_Lung</td><td>-0.05680737</td><td>-0.2598475</td><td>0.2967483</td><td>-0.03843072</td><td>0.5135266</td><td>-0.007136273</td><td>0.4121959</td><td>0.4046568</td><td>-0.1584653</td><td>-0.2685863</td><td>0.6906627</td><td>0.3752852</td><td>-0.04649791</td><td>0.1108001</td><td>0.102126</td><td>0.6544833</td><td>-0.1760352</td><td>-0.06733032</td><td>-0.3217252</td><td>-0.282326</td><td>-0.4695504</td><td>0.1894989</td><td>0.1392143</td><td>-0.1651271</td><td>0.0155587</td><td>-0.7552172</td><td>0.6330324</td><td>-0.3265327</td><td>-0.04649791</td><td>0.03171949</td><td>B_Cell      </td><td>-0.03843072</td></tr>
</tbody>
</table>




```R

```

# Plot GSVA enrichments vs some labels


```R
options(repr.plot.width=24, repr.plot.height=16)

Endo_gsva_plot1 <- ggplot(merged_annots_melt, aes(x=Navin_Clusters, y=GSVA_enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5)) + 
  #geom_beeswarm(cex = 0.5, size=0.5) +
  facet_wrap(vars(GeneSet)) +
  labs(title = "GSVA Enrichments, by Navin Labeled Clusters") + 
  theme_bw() + 
  theme(text=element_text(size=20),
        plot.title = element_text(size=22),
       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Endo_gsva_plot1

# save plot
#png("/immuno/ian/Projects/2023_02_SCLC_10X_Visium/Seurat_Analysis/Output_CD31_Focus_20240307/VA_Endo_GSVA_bySeuratCluster.png",
#    width = 1600,
#    height = 800)
#print(Endo_gsva_plot1)
#dev.off()
```


    
![png](output_114_0.png)
    


# Plot GSVA results grouped by Navin Annotations


```R
options(repr.plot.width=24, repr.plot.height=16)

Endo_gsva_plot2 <- ggplot(merged_annots_melt, aes(x=Navin_annotations, y=GSVA_enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5)) + 
  #geom_beeswarm(cex = 0.5, size=0.5) +
  facet_wrap(vars(GeneSet)) +
  labs(title = "GSVA Enrichments, by Navin Labeled Annotations") + 
  theme_bw() + 
  theme(text=element_text(size=20),
        plot.title = element_text(size=22),
       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Endo_gsva_plot2

# save plot
#png("/immuno/ian/Projects/2023_02_SCLC_10X_Visium/Seurat_Analysis/Output_CD31_Focus_20240307/VA_Endo_GSVA_byNavinAnnot.png",
#    width = 1600,
#    height = 800)
#print(VA_Endo_gsva_plot2)
#dev.off()
```


    
![png](output_116_0.png)
    


# Plot GSVA results grouped by Tissue Slice


```R
options(repr.plot.width=24, repr.plot.height=16)

Endo_gsva_plot3 <- ggplot(merged_annots_melt, aes(x=Tissue_Slice, y=GSVA_enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5)) + 
  #geom_beeswarm(cex = 0.5, size=0.5) +
  facet_wrap(vars(GeneSet)) +
  labs(title = "GSVA Enrichments, by Tissue Slice") + 
  theme_bw() + 
  theme(text=element_text(size=20),
        plot.title = element_text(size=22),
       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Endo_gsva_plot3

# save plot
#png("/immuno/ian/Projects/2023_02_SCLC_10X_Visium/Seurat_Analysis/Output_CD31_Focus_20240307/VA_Endo_GSVA_byNavinAnnot.png",
#    width = 1600,
#    height = 800)
#print(VA_Endo_gsva_plot2)
#dev.off()
```


    
![png](output_118_0.png)
    


# Plot GSVA results grouped by vasc_label


```R
options(repr.plot.width=24, repr.plot.height=16)

Endo_gsva_plot4 <- ggplot(merged_annots_melt, aes(x=vasc_label, y=GSVA_enrichment)) + 
  geom_hline(yintercept=0) + 
  geom_violin(draw_quantiles = c(0.5)) + 
  #geom_beeswarm(cex = 0.5, size=0.5) +
  facet_wrap(vars(GeneSet)) +
  labs(title = "GSVA Enrichments, by vasc_label") + 
  theme_bw() + 
  theme(text=element_text(size=20),
        plot.title = element_text(size=22),
       axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

Endo_gsva_plot4

# save plot
#png("/immuno/ian/Projects/2023_02_SCLC_10X_Visium/Seurat_Analysis/Output_CD31_Focus_20240307/VA_Endo_GSVA_byNavinAnnot.png",
#    width = 1600,
#    height = 800)
#print(VA_Endo_gsva_plot2)
#dev.off()
```


    
![png](output_120_0.png)
    



```R

```


```R

```


```R

```


```R

```
