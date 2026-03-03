## February 2025
## Round 1 REVISIONS for SCLC 10X Visium
## Combine Round 1 analysis with Sample A of Round 2
## Annotate spots - Neg, Low, or High for certain markers
    - use marker combinations to annotate spots for things like Tumor
    - Incorporate Navin's annotations to VA1 to see how they correlate with the new clusters. 


```R
R.home()
```


'/usr/local/lib/R'



```R
library(Seurat)
#library(SeuratExtend)
library(Matrix)
library(doParallel)
library(ggplot2)
library(org.Hs.eg.db)
library(msigdbr)
library(tidyverse)
library(dittoSeq)
library(pheatmap)
library(RColorBrewer)
library(scales)
library(GGally)
library(rlang)
library(data.table)

options(repr.matrix.max.cols=50, repr.matrix.max.rows=100)

```

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
    
    
    Registered S3 method overwritten by 'GGally':
      method from   
      +.gg   ggplot2
    
    
    Attaching package: ‘rlang’
    
    
    The following objects are masked from ‘package:purrr’:
    
        %@%, flatten, flatten_chr, flatten_dbl, flatten_int, flatten_lgl,
        flatten_raw, invoke, splice
    
    
    The following object is masked from ‘package:Biobase’:
    
        exprs
    
    
    
    Attaching package: ‘data.table’
    
    
    The following object is masked from ‘package:rlang’:
    
        :=
    
    
    The following objects are masked from ‘package:lubridate’:
    
        hour, isoweek, mday, minute, month, quarter, second, wday, week,
        yday, year
    
    
    The following objects are masked from ‘package:dplyr’:
    
        between, first, last
    
    
    The following object is masked from ‘package:purrr’:
    
        transpose
    
    
    The following object is masked from ‘package:IRanges’:
    
        shift
    
    
    The following objects are masked from ‘package:S4Vectors’:
    
        first, second
    
    


# read data


```R
AB1_A2 <- readRDS('Processing/AB1_2_integrated_20250207.rds')
```

# check data


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
GetAssayData(AB1_A2, layer='data')[1:5,1:5]
```


    5 x 5 sparse Matrix of class "dgCMatrix"
            VA2_AAACAAGTATCTCCCA-1 VA2_AAACACCAATAACTGC-1 VA2_AAACAGAGCGACTCCT-1
    SAMD11               .                      .                      .        
    NOC2L                0.5032359              0.6119117              0.8623837
    KLHL17               .                      .                      .        
    PLEKHN1              .                      .                      .        
    PERM1                .                      .                      .        
            VA2_AAACAGCTTTCAGAAG-1 VA2_AAACAGGGTCTATATT-1
    SAMD11               .                      .        
    NOC2L                0.4945215              0.9601302
    KLHL17               .                      .        
    PLEKHN1              .                      .        
    PERM1                .                      .        



```R
test_var <- 'testing123'
test_var
sym(test_var)
#!!sym(test_var)

test_var2 <- paste0(test_var,'_var2')
test_var2
sym(test_var2)
```


'testing123'



    testing123



'testing123_var2'



    testing123_var2



```R
AB1_A2@meta.data %>% head(5)
```


<table class="dataframe">
<caption>A data.frame: 5 × 8</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA2_AAACAAGTATCTCCCA-1</th><td>VA2</td><td>45867</td><td>8515</td><td>round2</td><td>24897</td><td>7041</td><td>5</td><td>5</td></tr>
	<tr><th scope=row>VA2_AAACACCAATAACTGC-1</th><td>VA2</td><td>11849</td><td>5082</td><td>round2</td><td>22823</td><td>5881</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>VA2_AAACAGAGCGACTCCT-1</th><td>VA2</td><td>43834</td><td>8662</td><td>round2</td><td>24885</td><td>7403</td><td>5</td><td>5</td></tr>
	<tr><th scope=row>VA2_AAACAGCTTTCAGAAG-1</th><td>VA2</td><td>46896</td><td>8897</td><td>round2</td><td>24782</td><td>7399</td><td>0</td><td>0</td></tr>
	<tr><th scope=row>VA2_AAACAGGGTCTATATT-1</th><td>VA2</td><td>18610</td><td>6201</td><td>round2</td><td>22740</td><td>6158</td><td>4</td><td>4</td></tr>
</tbody>
</table>




```R
table(AB1_A2$orig.ident)
```


    
     VA1  VA2  VB1 
    2806 4316 2589 


# Let's add orig.ident as a factor column instead, so we can control the ordering in plots. 


```R
AB1_A2@meta.data$slice_ident <- factor(AB1_A2@meta.data$orig.ident, levels = c("VA2", "VA1", "VB1"))
```


```R
AB1_A2@meta.data %>% head(5)
```


<table class="dataframe">
<caption>A data.frame: 5 × 9</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA2_AAACAAGTATCTCCCA-1</th><td>VA2</td><td>45867</td><td>8515</td><td>round2</td><td>24897</td><td>7041</td><td>5</td><td>5</td><td>VA2</td></tr>
	<tr><th scope=row>VA2_AAACACCAATAACTGC-1</th><td>VA2</td><td>11849</td><td>5082</td><td>round2</td><td>22823</td><td>5881</td><td>0</td><td>0</td><td>VA2</td></tr>
	<tr><th scope=row>VA2_AAACAGAGCGACTCCT-1</th><td>VA2</td><td>43834</td><td>8662</td><td>round2</td><td>24885</td><td>7403</td><td>5</td><td>5</td><td>VA2</td></tr>
	<tr><th scope=row>VA2_AAACAGCTTTCAGAAG-1</th><td>VA2</td><td>46896</td><td>8897</td><td>round2</td><td>24782</td><td>7399</td><td>0</td><td>0</td><td>VA2</td></tr>
	<tr><th scope=row>VA2_AAACAGGGTCTATATT-1</th><td>VA2</td><td>18610</td><td>6201</td><td>round2</td><td>22740</td><td>6158</td><td>4</td><td>4</td><td>VA2</td></tr>
</tbody>
</table>



# Function to annotate neg, low, high for a list of markers

# Based on the previous checks (02__pre_annotation_function_development.ipynb), VA2 normalized data has a lower peak that is moving its median down and contibuting to a lower median threshold which results in more spots being "high". 
# Because of this, I think we should set a global threshold that we apply to all samples. The data is normalized, after all. 
## Let's add functionality to have different score thresholds for low and high depending on the tissue slice. This will be defined by a named list input, basically a dict. This will be a fraction... default will be 0.5, or the half of the max possible score. 


```R
# by slice version, will calculate scores within slice only to try to account for batch effects. 
# annotate_spots_scored_markers will take a marker list, and create a score for each marker. 
# 0 is negative expression, 1 is below nonzero median, and 2 is above nonzero median. 
# seu_obj is the seurat object we're working with
# marker_list is the list of genes we want to use to annotate. 
# We'll add the scores of all these genes together to get a combined score for the gene list. 
# Then we'll split the combined score at the median to create a final low or high score. 
# annot_name is the name of the annotation label we want for the combined marker labels. 
# cutoff is 50% of max possible score, unless otherwise specified. use a decimal... 0.5 = 50%
# score_per_slice determines whether the gene scores will be calculated on a per-tissue-slice basis. 
# by default, score_per_slice will be false, so it'll calculate based on the global mean. 
# However, for some genes that are low expressing, it might be better to do it by tissue slice. 
annotate_spots_scored_markers <- function(seu_obj, marker_list, annot_name, check_markers = TRUE, check_score=TRUE, assay='Spatial', 
                                          custom_score_cutoff_per_slice=c("VA1"=0.5, "VB1"=0.5, "VA2"=0.5), score_per_slice=FALSE) {

    # Error handling: Ensure all tissue slices have a cutoff value
    tissue_slices <- unique(seu_obj@meta.data$orig.ident)
    missing_slices <- setdiff(tissue_slices, names(custom_score_cutoff_per_slice))
    
    if (length(missing_slices) > 0) {
        warning("The following tissue slices are missing in custom_score_cutoff_per_slice and will be assigned a default cutoff of 0.5: ", 
                paste(missing_slices, collapse = ", "))
        # Add missing tissue slices with default cutoff value
        custom_score_cutoff_per_slice[missing_slices] <- 0.5
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # STEP 1: Get Normalized Expression Values
    
    # get expression for markers in marker_list, using the normalized data
    expr_df = FetchData(seu_obj, vars = marker_list, assay=assay_use, layer = 'data')
    # Merge orig.ident and slice_ident into the expr_df
    expr_df <- merge(expr_df, select(seu_obj@meta.data, c(orig.ident, slice_ident)), by='row.names', all=TRUE)
    
    # check the markers if desired
    if (check_markers == TRUE) {
        options(repr.plot.width=22, repr.plot.height=8)
        # copy expr_df
        expr_df_check <- expr_df
        for (curr_marker in marker_list) {
            # plot histogram split by slice
            print(
            ggplot(expr_df_check,aes(x = .data[[curr_marker]]))+geom_histogram()+facet_grid(~slice_ident)+theme_bw()+ggtitle(curr_marker)
                )
            # plot spatial feature plot
            print(SpatialFeaturePlot(seu_obj, features = c(curr_marker), pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=7))))
            # plot combined histogram
            print(
            ggplot(expr_df_check[expr_df[[curr_marker]] > 0,],aes(x = .data[[curr_marker]]))+geom_histogram()+theme_bw()+ggtitle(curr_marker)
                )
        }
    }
    
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # STEP 2: Split at Nonzero Median and Score 0, 1, or 2
    
    # Score each marker 0 (neg), 1 (below nonzero median), or 2 (above nonzero median)
    # initialize a list to store the score columns
    # (and check values per tissue slice (orig.ident))
    #scored_marker_list <- c()
    scored_marker_list <- vector(mode="character", length=length(marker_list))
    i=1
    for (curr_marker in marker_list) {

        # We'll use global median to do the thresholding. 
        # split at median
        curr_median <- median(expr_df[[curr_marker]][expr_df[[curr_marker]] > 0])
        # create a new column to score this marker
        curr_marker_score <- paste0(curr_marker,'_score')
        # add curr_marker_score to the scored_marker_list to easily access these columns later
        scored_marker_list[[i]] <- curr_marker_score
        # iterate the counter
        i = i + 1
        # default negative
        expr_df[[curr_marker_score]] <- 0

        if (score_per_slice) {
            # Use a tissue-slice specific median to score
            # calculate median per tissue slice and store it in a temporary column
            expr_df <- expr_df %>%
            group_by(orig.ident) %>%
            mutate(curr_median = median(.data[[curr_marker]][.data[[curr_marker]] > 0])) %>%
            ungroup()
            # use the slice-specific median to threshold score
            expr_df[[curr_marker_score]][expr_df[[curr_marker]] > 0.1] <- 1
            expr_df[[curr_marker_score]][expr_df[[curr_marker]] > expr_df$curr_median] <- 2
        } else {
            # Use global median to score
            # threshold at 0.1 to remove the zero population
            expr_df[[curr_marker_score]][expr_df[[curr_marker]] > 0.1] <- 1
            # threshold at at median to label high group
            expr_df[[curr_marker_score]][expr_df[[curr_marker]] > curr_median] <- 2
        }
        
        # just for QA, we can check per slice medians if desired
        if (check_markers == TRUE) {
            # Split, get median per slice, and concatenate (zeros exlcuded)
            curr_marker_median_label <- paste0(curr_marker,'_per_slice_median')
            expr_df <- expr_df %>%
              group_by(orig.ident) %>%
              mutate(!!curr_marker_median_label := median(.data[[curr_marker]][.data[[curr_marker]] > 0])) %>%
              ungroup()
    
            # Get median including zeros as well for comparison...
            curr_marker_median_zeros_label <- paste0(curr_marker,'_per_slice_median_with_zeros')
            expr_df <- expr_df %>%
              group_by(orig.ident) %>%
              mutate(!!curr_marker_median_zeros_label := median(.data[[curr_marker]])) %>%
              ungroup()

            # Calculate the fraction of zero values per group
            curr_marker_zero_fraction_label <- paste0(curr_marker, '_per_slice_zero_fraction')
            expr_df <- expr_df %>%
              group_by(orig.ident) %>%
              mutate(!!curr_marker_zero_fraction_label := sum(.data[[curr_marker]] == 0) / n()) %>%
              ungroup()
        }
        
        # check median value and label...
        print('---------------------')
        print('- - - - - - - - - - -')
        print(curr_marker)
        print(' - - global nonzero median - - -')
        print(curr_median)
        if (check_markers == TRUE) {
            print(' - - median per slice excluding zeros - - ')
            print(table(expr_df[[curr_marker_median_label]], expr_df$orig.ident))
            print(' - - median per slice including zeros - - ')
            print(table(expr_df[[curr_marker_median_zeros_label]], expr_df$orig.ident))
            print(' - - fraction of spots equal to zero per slice - - ')
            print(table(expr_df[[curr_marker_zero_fraction_label]], expr_df$orig.ident))
        }
        print(' - - resulting scores - - ')
        print('overall: ')
        print(table(expr_df[[curr_marker_score]]))
        print('by slice: ')
        print(table(expr_df[[curr_marker_score]], expr_df$orig.ident))
    }

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # STEP 3: Label signatures as Low or High based on scores
    # Use half the max possible score as the cutoff between Low and High. 
    
    # now sum together all scored markers
    annot_name_scored = paste0(annot_name,'_score')
    expr_df <- expr_df %>% mutate(!!sym(annot_name_scored) := select(., all_of(scored_marker_list)) %>% rowSums(na.rm = TRUE))
    # since the max score is 2 * length(scored_marker_list), the halfway point will be the cutoff for low vs high
    # Apply custom cutoff per tissue slice
    expr_df <- expr_df %>%
      rowwise() %>%
      mutate(cutoff = custom_score_cutoff_per_slice[as.character(orig.ident)] * 2 * length(scored_marker_list)) %>%
      ungroup()

    # fix the merge from above
    expr_df <- as.data.frame(expr_df)
    # set Row.names to rownames
    rownames(expr_df) <- expr_df$Row.names
    # remove Row.names column
    expr_df <- select(expr_df, -c('Row.names'))
    
    # annotate neg, low, and high as label1
    annot_name_label1 <- paste0(annot_name,'_label1')
    # default neg
    expr_df[[annot_name_label1]] <- paste0(annot_name,'_neg')
    # if above zero, label as low
    expr_df[[annot_name_label1]][expr_df[[annot_name_scored]] > 0] <- paste0(annot_name,'_low')
    # threshold at custom cutoff to label high group
    expr_df[[annot_name_label1]][expr_df[[annot_name_scored]] > expr_df$cutoff] <- paste0(annot_name,'_high')

    # annotate neg, low, and high as label
    annot_name_label <- paste0(annot_name,'_label')
    # default low
    expr_df[[annot_name_label]] <- paste0(annot_name,'_low')
    # threshold at custom cutoff to label high group
    expr_df[[annot_name_label]][expr_df[[annot_name_scored]] > expr_df$cutoff] <- paste0(annot_name,'_high')

    # Merge with the seurat object metadata
    seu_obj@meta.data <- merge(seu_obj@meta.data, select(expr_df, all_of(c(annot_name_scored, annot_name_label1, annot_name_label))), by = 'row.names', all = TRUE)
    
    # set Row.names to rownames
    rownames(seu_obj@meta.data) <- seu_obj@meta.data$Row.names
    # remove Row.names column
    seu_obj@meta.data <- select(seu_obj@meta.data, -c('Row.names'))
    
    
    # check score distribution
    if (check_score == TRUE) {
        options(repr.plot.width=22, repr.plot.height=8)
        # copy expr_df
        expr_df_check <- expr_df
        options(repr.plot.width=22, repr.plot.height=8)
        # plot histogram
        print(
        ggplot(expr_df_check,aes(x = .data[[annot_name_scored]]))+geom_histogram()+facet_grid(~slice_ident)+theme_bw()+ggtitle(annot_name_scored)
            )
        # plot scores on spatial feature plot
        options(repr.plot.width=22, repr.plot.height=8)
        print(SpatialFeaturePlot(seu_obj, features = c(annot_name_scored), pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=7))))
        # plot labels on spatial dim plot
        options(repr.plot.width=22, repr.plot.height=8)
        clrs_use0 <- setNames(c('#9c9c9eFF','#F39600FF','#A60E00FF'),c(paste0(annot_name,'_neg'),paste0(annot_name,'_low'),paste0(annot_name,'_high')))
        clrs_use <- setNames(c('#F39600FF','#A60E00FF'),c(paste0(annot_name,'_low'),paste0(annot_name,'_high')))
        old_idents <- Idents(seu_obj)
        Idents(seu_obj) <- annot_name_label1
        print(SpatialDimPlot(seu_obj, cols = clrs_use0, pt.size.factor=2) & guides(fill = guide_legend(override.aes = list(size=7))))
        Idents(seu_obj) <- old_idents
        # print scores and labels
        options(repr.plot.width=22, repr.plot.height=8)
        print(' - - resulting final scores - - ')
        print('overall: ')
        print(table(expr_df_check[[annot_name_scored]]))
        print('by slice: ')
        scored_table <- table(expr_df_check[[annot_name_scored]], expr_df_check$orig.ident)
        label_table <- table(expr_df_check[[annot_name_label1]], expr_df_check$orig.ident)
        print(scored_table)
        print(label_table)
        print('proportions by slice: ')
        print(prop.table(scored_table, margin = 2))  # Use margin = 2 for column proportions
        print(prop.table(label_table, margin = 2))  # Use margin = 2 for column proportions
        
        options(repr.plot.width=24, repr.plot.height=24)
        options(repr.plot.width=24, repr.plot.height=24)
        # use GGally package to do a pair plot of marker values
        plot1 <- ggpairs(expr_df_check %>% select(all_of(marker_list)),
                        upper = list(continuous = wrap("cor", size = 9))) 
        
        plot1 <- plot1 + 
                    theme_bw() +
                    theme(strip.text.x = element_text(size = 20),
                          strip.text.y = element_text(size = 20))
        print(plot1)
        options(repr.plot.width=24, repr.plot.height=24)
    }
    
    return(seu_obj)
}

#test <- annotate_spots_scored_markers(AB1_A2, c('PECAM1', 'VWF', 'CDH5', 'KDR', 'CLDN5', 'PLVAP'), 'vasc', check_markers=FALSE, assay='Spatial')
```

# Update June 9, 2025
# Use a simpler labeling method for the vasculature, we're losing too many spots. 
# Just take all the genes in the signature, sum their expression together, and threshold that for low vs high. 
# This will be less selective and capture more spots than with the scoring method. 
# Still use the scoring method for tumor annotation. 


```R
# simpler version that just summs gene expression and thresholds instead of using scores. 
annotate_spots_summed_expr_markers_simple <- function(
    seu_obj, marker_list, annot_name, assay='Spatial', per_slice=FALSE, check_markers=TRUE, make_plots=TRUE
) {
    # STEP 1: Fetch expression values for markers
    expr_df <- FetchData(seu_obj, vars = marker_list, assay = assay, layer = 'data')
    expr_df$orig.ident <- seu_obj@meta.data$orig.ident
    expr_df$slice_ident <- seu_obj@meta.data$slice_ident

    # STEP 2: Sum expression values across markers for each spot
    expr_df$summed_expr <- rowSums(expr_df[, marker_list], na.rm = TRUE)

    # STEP 3: Determine threshold (nonzero median), optionally per slice
    if (per_slice) {
        expr_df <- expr_df %>%
            group_by(orig.ident) %>%
            mutate(summed_expr_median = median(summed_expr[summed_expr > 0])) %>%
            ungroup()
        expr_df$label <- paste0(annot_name, '_low')
        expr_df$label[expr_df$summed_expr == 0] <- paste0(annot_name, '_neg')
        expr_df$label[expr_df$summed_expr > expr_df$summed_expr_median] <- paste0(annot_name, '_high')
    } else {
        nonzero_median <- median(expr_df$summed_expr[expr_df$summed_expr > 0])
        expr_df$label <- paste0(annot_name, '_low')
        expr_df$label[expr_df$summed_expr == 0] <- paste0(annot_name, '_neg')
        expr_df$label[expr_df$summed_expr > nonzero_median] <- paste0(annot_name, '_high')
    }

    # STEP 4: Merge annotation back into Seurat object metadata
    seu_obj@meta.data[[paste0(annot_name, "_summed_expr")]] <- expr_df$summed_expr
    seu_obj@meta.data[[paste0(annot_name, "_label")]] <- expr_df$label

    # STEP 5: Optional QC and spatial plots
    if (check_markers) {
        print(ggplot(expr_df, aes(x = summed_expr)) +
            geom_histogram() +
            facet_grid(~slice_ident) +
            theme_bw() +
            ggtitle(paste0("Summed expression for ", annot_name)))
        print(table(expr_df$label, expr_df$orig.ident))
    }
    
    if (make_plots) {
        options(repr.plot.width=22, repr.plot.height=8)
        # Plot summed expression spatially
        print(SpatialFeaturePlot(
            seu_obj, 
            features = paste0(annot_name, "_summed_expr"),
            pt.size.factor = 2
        ) + ggtitle(paste0("Spatial plot: summed expression for ", annot_name)))
        
        # Plot annotation labels spatially
        # Set colors for neg, low, high
        clrs_use <- setNames(
            c('#9c9c9eFF','#F39600FF','#A60E00FF'),
            c(paste0(annot_name,'_neg'), paste0(annot_name,'_low'), paste0(annot_name,'_high'))
        )
        old_idents <- Idents(seu_obj)
        Idents(seu_obj) <- paste0(annot_name, "_label")
        print(SpatialDimPlot(
            seu_obj, 
            cols = clrs_use,
            pt.size.factor = 2, 
            label = TRUE, 
            label.size = 3
        ) + ggtitle(paste0("Spatial plot: ", annot_name, " labels")))
        Idents(seu_obj) <- old_idents
    }

    return(seu_obj)
}

#test <- annotate_spots_summed_expr_markers_simple(
#  AB1_A2,
#  marker_list = c('PECAM1', 'VWF', 'CDH5', 'KDR', 'CLDN5', 'PLVAP'),
#  annot_name = 'vasc',
#  assay = 'Spatial',
#  per_slice = FALSE,
#  check_markers = TRUE
#)

```


```R

AB1_A2 <- annotate_spots_summed_expr_markers_simple(
  AB1_A2,
  marker_list = c('PECAM1', 'VWF', 'CDH5', 'KDR', 'CLDN5', 'PLVAP'),
  annot_name = 'vasc',
  assay = 'Spatial',
  per_slice = FALSE,
  check_markers = TRUE
)
```

    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.


               
                 VA1  VA2  VB1
      vasc_high  771 1733  665
      vasc_low   532 2080  558
      vasc_neg  1503  503 1366



    
![png](output_19_2.png)
    


    [1m[22mScale for [32mfill[39m is already present.
    Adding another scale for [32mfill[39m, which will replace the existing scale.
    [1m[22mScale for [32mfill[39m is already present.
    Adding another scale for [32mfill[39m, which will replace the existing scale.
    [1m[22mScale for [32mfill[39m is already present.
    Adding another scale for [32mfill[39m, which will replace the existing scale.



    
![png](output_19_4.png)
    



    
![png](output_19_5.png)
    


# Annotate Tumor Markers...


```R
# calculate cutoffs for custom_score_cutoff_per_slice... should be desired score divied by max possible score. 
c("VA1"=7.2, "VB1"=7.2, "VA2"=9) / 18
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>VA1</dt><dd>0.4</dd><dt>VB1</dt><dd>0.4</dd><dt>VA2</dt><dd>0.5</dd></dl>




```R
Tumor_and_E2F_list <- c('EPCAM','ASCL1','NKX2-1','CDK1','CHEK1','DLGAP5','HELLS','LMNB1','MCM2','MKI67','MYBL2','NASP','NUDT21','PLK1','RPA1','RFC2','UBE2T','UBR7')
```


```R

AB1_A2 <- annotate_spots_scored_markers(AB1_A2, Tumor_and_E2F_list, 'Tumor', check_markers=TRUE, assay='Spatial', custom_score_cutoff_per_slice=c("VA1"=0.4, "VB1"=0.4, "VA2"=0.44))
```

    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_1.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_3.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_5.png)
    



    
![png](output_23_6.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_8.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_10.png)
    



    
![png](output_23_11.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_13.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_15.png)
    



    
![png](output_23_16.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_18.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_20.png)
    



    
![png](output_23_21.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_23.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_25.png)
    



    
![png](output_23_26.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_28.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_30.png)
    



    
![png](output_23_31.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_33.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_35.png)
    



    
![png](output_23_36.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_38.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_40.png)
    



    
![png](output_23_41.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_43.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_45.png)
    



    
![png](output_23_46.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_48.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_50.png)
    



    
![png](output_23_51.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_53.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_55.png)
    



    
![png](output_23_56.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_58.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_60.png)
    



    
![png](output_23_61.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_63.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_65.png)
    



    
![png](output_23_66.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_68.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_70.png)
    



    
![png](output_23_71.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_73.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_75.png)
    



    
![png](output_23_76.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_78.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_80.png)
    



    
![png](output_23_81.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_83.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_85.png)
    



    
![png](output_23_86.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_88.png)
    


    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "EPCAM"
    [1] " - - global nonzero median - - -"
    [1] 1.780555
    [1] " - - median per slice excluding zeros - - "
                      
                        VA1  VA2  VB1
      1.58973313326251    0    0 2589
      1.63056832072081 2806    0    0
      1.89447113241465    0 4316    0
    [1] " - - median per slice including zeros - - "
                      
                        VA1  VA2  VB1
      1.3015178561061  2806    0    0
      1.35162227922205    0    0 2589
      1.8430047601044     0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.108897126969416    0 4316    0
      0.300888373889533    0    0 2589
      0.353884533143264 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    2242 3735 3734 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0  993  470  779
      1 1114 1424 1197
      2  699 2422  613
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "ASCL1"
    [1] " - - global nonzero median - - -"
    [1] 1.93906
    [1] " - - median per slice excluding zeros - - "
                      
                        VA1  VA2  VB1
      1.44939008287489 2806    0    0
      1.47252138040221    0    0 2589
      2.19168887724791    0 4316    0
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0.902187749586346 2806    0    0
      1.03386839523212     0    0 2589
      2.15640705385094     0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                        
                          VA1  VA2  VB1
      0.0764596848934198    0 4316    0
      0.421784472769409     0    0 2589
      0.466500356379187  2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    2731 3490 3490 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1309  330 1092
      1 1244 1015 1231
      2  253 2971  266
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "NKX2-1"
    [1] " - - global nonzero median - - -"
    [1] 1.858969
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.379513914776521    0 4316    0
      2.01271565717147  2806    0    0
      2.09933724744713     0    0 2589
    [1] " - - median per slice including zeros - - "
                      
                        VA1  VA2  VB1
      0                   0 4316    0
      1.87118900869168 2806    0    0
      2.02025167330389    0    0 2589
    [1] " - - fraction of spots equal to zero per slice - - "
                        
                          VA1  VA2  VB1
      0.0961761297798378    0    0 2589
      0.16535994297933   2806    0    0
      0.701807228915663     0 4316    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    3743 2984 2984 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0  464 3030  249
      1  919 1256  809
      2 1423   30 1531
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "CDK1"
    [1] " - - global nonzero median - - -"
    [1] 0.9006222
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.596254055152501    0 4316    0
      1.32791174714174  2806    0    0
      1.33441936329812     0    0 2589
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.405340147331745    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.360055607043559    0 4316    0
      0.548860563924295    0    0 2589
      0.591945830363507 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    4636 2538 2537 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1661 1554 1421
      1   95 2346   97
      2 1050  416 1071
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "CHEK1"
    [1] " - - global nonzero median - - -"
    [1] 0.8164323
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.622354128568846    0 4316    0
      1.22882458640562  2806    0    0
      1.23622668766573     0    0 2589
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.416349674493078    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.361214087117702    0 4316    0
      0.702974121282348    0    0 2589
      0.727369921596579 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    5420 2146 2145 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 2041 1559 1820
      1   39 2073   34
      2  726  684  735
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "DLGAP5"
    [1] " - - global nonzero median - - -"
    [1] 0.7615843
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.467935437951272    0 4316    0
      1.23952942733293     0    0 2589
      1.24241375479555  2806    0    0
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.215109837685384    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.465940685820204    0 4316    0
      0.682116647354191    0    0 2589
      0.696364932287954 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    5731 1990 1990 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1954 2011 1766
      1   17 1960   13
      2  835  345  810
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "HELLS"
    [1] " - - global nonzero median - - -"
    [1] 1.134671
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.793354594326016    0 4316    0
      1.54747591874701     0    0 2589
      1.55590374710844  2806    0    0
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0.640126340388917    0 4316    0
      1.13941973554525  2806    0    0
      1.2237189639569      0    0 2589
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.274096385542169    0 4316    0
      0.361143298570877    0    0 2589
      0.40128296507484  2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    3244 3234 3233 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1126 1183  935
      1  270 2715  249
      2 1410  418 1405
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "LMNB1"
    [1] " - - global nonzero median - - -"
    [1] 0.9927918
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.681552278132147    0 4316    0
      1.39566711098461     0    0 2589
      1.39749811192432  2806    0    0
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0    0
      0.531132127510214    0 4316    0
      0.864975727538744    0    0 2589
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.308387395736793    0 4316    0
      0.483198146002317    0    0 2589
      0.528510334996436 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    4065 2823 2823 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1483 1331 1251
      1  179 2490  154
      2 1144  495 1184
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "MCM2"
    [1] " - - global nonzero median - - -"
    [1] 0.9861554
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.766408941666347    0 4316    0
      1.32934023152993  2806    0    0
      1.34173020814965     0    0 2589
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.595264075430689    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.300741427247451    0 4316    0
      0.544611819235226    0    0 2589
      0.564148253741982 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    4291 2710 2710 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1583 1298 1410
      1  196 2323  191
      2 1027  695  988
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "MKI67"
    [1] " - - global nonzero median - - -"
    [1] 0.7684883
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.578105047538749    0 4316    0
      1.1943757147706   2806    0    0
      1.23122362036216     0    0 2589
    [1] " - - median per slice including zeros - - "
                      
                        VA1  VA2  VB1
      0                2806    0 2589
      0.36114534110023    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.379981464318814    0 4316    0
      0.736191579760525    0    0 2589
      0.742694226657163 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    5630 2041 2040 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 2084 1640 1906
      1   23 2002   16
      2  699  674  667
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "MYBL2"
    [1] " - - global nonzero median - - -"
    [1] 1.075136
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.804491487360202    0 4316    0
      1.42926076798589     0    0 2589
      1.44851798147451  2806    0    0
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0.635933003407319    0 4316    0
      0.873747042545319 2806    0    0
      0.942988681240424    0    0 2589
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.287534754402224    0 4316    0
      0.455774430281962    0    0 2589
      0.476122594440485 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    3757 2977 2977 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1336 1241 1180
      1  268 2480  229
      2 1202  595 1180
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "NASP"
    [1] " - - global nonzero median - - -"
    [1] 0.8561179
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.467130381442009    0 4316    0
      1.30064217465042  2806    0    0
      1.33632889224464     0    0 2589
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.233749284230873    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.44578313253012     0 4316    0
      0.581691772885284    0    0 2589
      0.616892373485389 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    5161 2275 2275 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1731 1924 1506
      1   65 2148   62
      2 1010  244 1021
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "NUDT21"
    [1] " - - global nonzero median - - -"
    [1] 1.189743
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.883705346479826    0 4316    0
      1.53165274644865     0    0 2589
      1.56289873402099  2806    0    0
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0.746534781172756    0 4316    0
      1.23553400256524  2806    0    0
      1.24566570662461     0    0 2589
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.258341056533828    0 4316    0
      0.346852066434917    0    0 2589
      0.353171774768354 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    3004 3354 3353 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0  991 1115  898
      1  361 2668  325
      2 1454  533 1366
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "PLK1"
    [1] " - - global nonzero median - - -"
    [1] 0.9046752
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.709401038355717    0 4316    0
      1.23835383021788  2806    0    0
      1.25778558076166     0    0 2589
    [1] " - - median per slice including zeros - - "
                      
                        VA1  VA2  VB1
      0                2806    0 2589
      0.52583676493121    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.305607043558851    0 4316    0
      0.682109764789736 2806    0    0
      0.682889146388567    0    0 2589
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    5001 2355 2355 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1914 1319 1768
      1  106 2169   80
      2  786  828  741
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "RPA1"
    [1] " - - global nonzero median - - -"
    [1] 0.7851232
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.546686647379621    0 4316    0
      1.26094786591725     0    0 2589
      1.27359517159765  2806    0    0
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.357976891279811    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.367701575532901    0 4316    0
      0.66928011404134  2806    0    0
      0.675936655079181    0    0 2589
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    5215 2248 2248 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1878 1587 1750
      1   34 2193   21
      2  894  536  818
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "RFC2"
    [1] " - - global nonzero median - - -"
    [1] 0.9235899
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.741621113913497    0 4316    0
      1.2727599991359   2806    0    0
      1.27461380144348     0    0 2589
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.588387308101852    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.295644114921223    0 4316    0
      0.630745461568173    0    0 2589
      0.671062009978617 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    4792 2460 2459 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1883 1276 1633
      1  120 2229  111
      2  803  811  845
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "UBE2T"
    [1] " - - global nonzero median - - -"
    [1] 0.9143303
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.664629550779886    0 4316    0
      1.29287082177518     0    0 2589
      1.30454346239692  2806    0    0
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.482541631710196    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.333178869323448    0 4316    0
      0.593665507918115    0    0 2589
      0.621525302922309 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    4719 2496 2496 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1744 1438 1537
      1  113 2278  105
      2  949  600  947
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "UBR7"
    [1] " - - global nonzero median - - -"
    [1] 1.073725
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.926210019227167    0 4316    0
      1.32903391276481  2806    0    0
      1.35673848119704     0    0 2589
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.798854262680656    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.241890639481001    0 4316    0
      0.558903051371186    0    0 2589
      0.563079116179615 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    4071 2820 2820 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1580 1044 1447
      1  304 2269  247
      2  922 1003  895


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_23_91.png)
    



    
![png](output_23_92.png)
    



    
![png](output_23_93.png)
    


    [1] " - - resulting final scores - - "
    [1] "overall: "
    
      0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19 
    148  25 169  89 258 178 302 216 379 309 352 359 392 396 509 449 484 539 652 767 
     20  21  22  23  24  25  26  27  28  29  30  31 
    783 671 488 316 185 134  75  38  25  16   6   2 
    [1] "by slice: "
        
         VA1 VA2 VB1
      0   85  54   9
      1    2  22   1
      2   41 104  24
      3   12  65  12
      4   89 137  32
      5   51  94  33
      6  112 136  54
      7   58 103  55
      8  130 135 114
      9  104 107  98
      10 135 100 117
      11 133  95 131
      12 127 113 152
      13 138 106 152
      14 202 136 171
      15 138 147 164
      16 150 177 157
      17 142 223 174
      18 175 318 159
      19 147 456 164
      20 135 507 141
      21 111 423 137
      22  97 292  99
      23  78 156  82
      24  52  77  56
      25  66  25  43
      26  40   8  27
      27  21   0  17
      28  19   0   6
      29  10   0   6
      30   4   0   2
      31   2   0   0
                
                  VA1  VA2  VB1
      Tumor_high 1387 2662 1434
      Tumor_low  1334 1600 1146
      Tumor_neg    85   54    9
    [1] "proportions by slice: "
        
                  VA1          VA2          VB1
      0  0.0302922309 0.0125115848 0.0034762457
      1  0.0007127584 0.0050973123 0.0003862495
      2  0.0146115467 0.0240963855 0.0092699884
      3  0.0042765502 0.0150602410 0.0046349942
      4  0.0317177477 0.0317423540 0.0123599846
      5  0.0181753386 0.0217794254 0.0127462341
      6  0.0399144690 0.0315106580 0.0208574739
      7  0.0206699929 0.0238646895 0.0212437234
      8  0.0463292944 0.0312789620 0.0440324450
      9  0.0370634355 0.0247914736 0.0378524527
      10 0.0481111903 0.0231696015 0.0451911935
      11 0.0473984319 0.0220111214 0.0505986868
      12 0.0452601568 0.0261816497 0.0587099266
      13 0.0491803279 0.0245597776 0.0587099266
      14 0.0719885959 0.0315106580 0.0660486674
      15 0.0491803279 0.0340593142 0.0633449208
      16 0.0534568781 0.0410101946 0.0606411742
      17 0.0506058446 0.0516682113 0.0672074160
      18 0.0623663578 0.0736793327 0.0614136732
      19 0.0523877406 0.1056533828 0.0633449208
      20 0.0481111903 0.1174698795 0.0544611819
      21 0.0395580898 0.0980074143 0.0529161839
      22 0.0345687812 0.0676552363 0.0382387022
      23 0.0277975766 0.0361445783 0.0316724604
      24 0.0185317177 0.0178405931 0.0216299730
      25 0.0235210264 0.0057924004 0.0166087292
      26 0.0142551675 0.0018535681 0.0104287370
      27 0.0074839629 0.0000000000 0.0065662418
      28 0.0067712046 0.0000000000 0.0023174971
      29 0.0035637919 0.0000000000 0.0023174971
      30 0.0014255167 0.0000000000 0.0007724990
      31 0.0007127584 0.0000000000 0.0000000000
                
                         VA1         VA2         VB1
      Tumor_high 0.494297933 0.616774791 0.553881808
      Tumor_low  0.475409836 0.370713624 0.442641947
      Tumor_neg  0.030292231 0.012511585 0.003476246



    
![png](output_23_95.png)
    



    
![png](output_23_96.png)
    


# Annotate MHCI Markers...


```R
# calculate cutoffs for custom_score_cutoff_per_slice... should be desired score divied by max possible score. 
c("VA1"=3, "VB1"=3, "VA2"=3) / 6
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>VA1</dt><dd>0.5</dd><dt>VB1</dt><dd>0.5</dd><dt>VA2</dt><dd>0.5</dd></dl>




```R

AB1_A2 <- annotate_spots_scored_markers(AB1_A2, c('TAP1','TAP2','B2M'), 'MHCI', check_markers=TRUE, assay='Spatial')
```

    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_26_1.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_26_3.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_26_5.png)
    



    
![png](output_26_6.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_26_8.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_26_10.png)
    



    
![png](output_26_11.png)
    


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_26_13.png)
    


    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "TAP1"
    [1] " - - global nonzero median - - -"
    [1] 1.051045
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.711593247622797    0 4316    0
      1.39996814823546     0    0 2589
      1.40673748343183  2806    0    0
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.564966007320299    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.239110287303058    0 4316    0
      0.518531717747684 2806    0    0
      0.518733101583623    0    0 2589
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    3830 2941 2940 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 1455 1032 1343
      1  250 2462  229
      2 1101  822 1017
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "TAP2"
    [1] " - - global nonzero median - - -"
    [1] 0.7990337
    [1] " - - median per slice excluding zeros - - "
                       
                         VA1  VA2  VB1
      0.589078619332677    0 4316    0
      1.22304180762793  2806    0    0
      1.22805569582272     0    0 2589
    [1] " - - median per slice including zeros - - "
                       
                         VA1  VA2  VB1
      0                 2806    0 2589
      0.412461850952901    0 4316    0
    [1] " - - fraction of spots equal to zero per slice - - "
                       
                         VA1  VA2  VB1
      0.315106580166821    0 4316    0
      0.729239088451139    0    0 2589
      0.733784746970777 2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
    5307 2202 2202 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0 2059 1360 1888
      1   40 2145   17
      2  707  811  684
    [1] "---------------------"
    [1] "- - - - - - - - - - -"
    [1] "B2M"
    [1] " - - global nonzero median - - -"
    [1] 3.290594
    [1] " - - median per slice excluding zeros - - "
                      
                        VA1  VA2  VB1
      3.16863012269802    0 4316    0
      3.38579700415636    0    0 2589
      3.40094185554013 2806    0    0
    [1] " - - median per slice including zeros - - "
                      
                        VA1  VA2  VB1
      3.16328400352503    0 4316    0
      3.38248419562585 2806    0    0
      3.3848389182514     0    0 2589
    [1] " - - fraction of spots equal to zero per slice - - "
                         
                           VA1  VA2  VB1
      0.00347624565469293    0    0 2589
      0.00880444856348471    0 4316    0
      0.022808267997149   2806    0    0
    [1] " - - resulting scores - - "
    [1] "overall: "
    
       0    1    2 
     111 4800 4800 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0   64   38    9
      1 1123 2588 1089
      2 1619 1690 1491


    [1m[22m`stat_bin()` using `bins = 30`. Pick better value with `binwidth`.



    
![png](output_26_16.png)
    



    
![png](output_26_17.png)
    



    
![png](output_26_18.png)
    


    [1] " - - resulting final scores - - "
    [1] "overall: "
    
       0    1    2    3    4    5    6 
     108 1250 2068 2663 2302  676  644 
    [1] "by slice: "
       
         VA1  VA2  VB1
      0   63   38    7
      1  471  309  470
      2  723  698  647
      3  452 1757  454
      4  720  918  664
      5  148  383  145
      6  229  213  202
               
                 VA1  VA2  VB1
      MHCI_high 1097 1514 1011
      MHCI_low  1646 2764 1571
      MHCI_neg    63   38    7
    [1] "proportions by slice: "
       
                VA1         VA2         VB1
      0 0.022451889 0.008804449 0.002703747
      1 0.167854597 0.071594069 0.181537273
      2 0.257662153 0.161723818 0.249903438
      3 0.161083393 0.407089898 0.175357281
      4 0.256593015 0.212696942 0.256469679
      5 0.052744120 0.088739574 0.056006180
      6 0.081610834 0.049351251 0.078022402
               
                        VA1         VA2         VB1
      MHCI_high 0.390947969 0.350787766 0.390498262
      MHCI_low  0.586600143 0.640407785 0.606797992
      MHCI_neg  0.022451889 0.008804449 0.002703747



    
![png](output_26_20.png)
    



    
![png](output_26_21.png)
    


# Function to combine labels, e.g. Tumor_High_MHCI_High


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
combine_labels <- function(seu_obj, label_col1, label_col2) {
    # Remove "_label" from end of the column names to combine them
    label_col1_sub <- gsub("_label","",label_col1)
    print(label_col1_sub)
    label_col2_sub <- gsub("_label","",label_col2)
    print(label_col2_sub)
    # combine labels
    combined_label_colname <- paste0(label_col1_sub,'_',label_col2_sub,'_label')
    print(combined_label_colname)

    # Logic to combine the labels
    seu_obj@meta.data[[combined_label_colname]] <- paste0(seu_obj@meta.data[[label_col1]], '_', seu_obj@meta.data[[label_col2]])

    return(seu_obj)
}
```

# check original tumor marker list


```R
# check original tumor marker list
test <- combine_labels(AB1_A2, 'Tumor_label', 'MHCI_label')

print(test@meta.data %>% head(3))
print(table(test@meta.data$Tumor_label))
print(table(test@meta.data$MHCI_label))
print(table(test@meta.data$Tumor_MHCI_label))
# plot labels on spatial dim plot
options(repr.plot.width=28, repr.plot.height=8)
old_idents <- Idents(test)
Idents(test) <- 'Tumor_MHCI_label'
print(SpatialDimPlot(test))
Idents(test) <- old_idents
```

    [1] "Tumor"
    [1] "MHCI"
    [1] "Tumor_MHCI_label"
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
                           MHCI_label    Tumor_MHCI_label
    VA1_AAACAACGAATAGTTC-1   MHCI_low  Tumor_low_MHCI_low
    VA1_AAACAAGTATCTCCCA-1   MHCI_low  Tumor_low_MHCI_low
    VA1_AAACAATCTACTAGCA-1   MHCI_low Tumor_high_MHCI_low
    
    Tumor_high  Tumor_low 
          5483       4228 
    
    MHCI_high  MHCI_low 
         3622      6089 
    
    Tumor_high_MHCI_high  Tumor_high_MHCI_low  Tumor_low_MHCI_high 
                    1903                 3580                 1719 
      Tumor_low_MHCI_low 
                    2509 



    
![png](output_31_1.png)
    



```R
# do with actual dataset
AB1_A2 <- combine_labels(AB1_A2, 'Tumor_label', 'MHCI_label')
```

    [1] "Tumor"
    [1] "MHCI"
    [1] "Tumor_MHCI_label"


# Make all Tumor Low in grey, then show Tumor High MHCI low, and Tumor High MHCI High as different colors. Similar to before. 


```R
#paletteer::paletteer_d("beyonce::X58")
#clrs <- paletteer::paletteer_d("RColorBrewer::OrRd")
clrs <- paletteer::paletteer_d("beyonce::X58")
show_col(clrs)
clrs
```


    <colors>
    [37m[48;5;232m#01000EFF[49m[39m [37m[48;5;166m#D42D24FF[49m[39m [30m[48;5;208m#EE7E24FF[49m[39m [30m[48;5;214m#FBA724FF[49m[39m [30m[48;5;223m#DABD8EFF[49m[39m [30m[48;5;250m#C9C8C6FF[49m[39m 



    
![png](output_34_1.png)
    



```R
clrs <- paletteer::paletteer_d("trekcolors::breen2")
show_col(clrs)
clrs
```


    <colors>
    [37m[48;5;88m#7E0500FF[49m[39m [37m[48;5;124m#A60E00FF[49m[39m [37m[48;5;160m#CE1800FF[49m[39m [30m[48;5;202m#DE4300FF[49m[39m [30m[48;5;208m#EF6F00FF[49m[39m [30m[48;5;214m#F39600FF[49m[39m [30m[48;5;220m#F7BD00FF[49m[39m [30m[48;5;222m#EDD072FF[49m[39m [30m[48;5;253m#E4E4E4FF[49m[39m 



    
![png](output_35_1.png)
    



```R
clrs_use <- c('Tumor_low_MHCI_low' = '#9c9c9eFF', 
              'Tumor_low_MHCI_high' = '#9c9c9eFF', 
              'Tumor_high_MHCI_low' = '#F39600FF',
              'Tumor_high_MHCI_high' = '#A60E00FF')
show_col(clrs_use)
```


    
![png](output_36_0.png)
    



```R
# Make another Tumor_MHCI_label2 which combines low groups. 
AB1_A2@meta.data$Tumor_MHCI_label2 <- "Tumor_low"
AB1_A2@meta.data$Tumor_MHCI_label2[AB1_A2@meta.data$Tumor_MHCI_label == "Tumor_high_MHCI_low"] <- "Tumor_high_MHCI_low"
AB1_A2@meta.data$Tumor_MHCI_label2[AB1_A2@meta.data$Tumor_MHCI_label == "Tumor_high_MHCI_high"] <- "Tumor_high_MHCI_high"
```


```R
# spatial dimplot colored by navin annotations

clrs_use2 <- c('Tumor_low' = '#9c9c9eFF', 
              'Tumor_high_MHCI_low' = '#F39600FF',
              'Tumor_high_MHCI_high' = '#A60E00FF')

Idents(AB1_A2) <- "Tumor_MHCI_label2"

#options(repr.plot.width=11, repr.plot.height=9)

plot <- SpatialDimPlot(AB1_A2, cols = clrs_use2, pt.size.factor=1.7, alpha=0.7) & guides(fill = guide_legend(override.aes = list(size=7)))
#+ labs(title="Tumor_MHCI_label")
plot

# save
png(paste0('Output/02_Annotate_spots/Tumor_MHCI_label2_SpatialDimPlot.png'),
    width = 2200,
    height = 800)
print(plot)
dev.off()

pdf(paste0('Output/02_Annotate_spots/Tumor_MHCI_label2_SpatialDimPlot.pdf'),
width = 22,
height = 8)
print(plot)
dev.off()

Idents(AB1_A2) <- "orig.ident"
```


<strong>pdf:</strong> 2



<strong>pdf:</strong> 2



    
![png](output_38_2.png)
    


# Let's relabel vasc_label to vasc_label1 and then use vasc_label to combine neg and low


```R
unique(AB1_A2@meta.data$vasc_label)
unique(AB1_A2@meta.data$Tumor_label1)
unique(AB1_A2@meta.data$Tumor_label)
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
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high'</li><li>'Tumor_neg'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high'</li></ol>




```R
AB1_A2@meta.data$vasc_label1 <- AB1_A2@meta.data$vasc_label
AB1_A2@meta.data$vasc_label[AB1_A2@meta.data$vasc_label1 == 'vasc_high'] <- 'vasc_high'
AB1_A2@meta.data$vasc_label[AB1_A2@meta.data$vasc_label1 == 'vasc_low'] <- 'vasc_low'
AB1_A2@meta.data$vasc_label[AB1_A2@meta.data$vasc_label1 == 'vasc_neg'] <- 'vasc_low'
```


```R
unique(AB1_A2@meta.data$vasc_label1)
unique(AB1_A2@meta.data$vasc_label)
unique(AB1_A2@meta.data$Tumor_label1)
unique(AB1_A2@meta.data$Tumor_label)
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




```R
dplyr::count(AB1_A2@meta.data, vasc_label1)
dplyr::count(AB1_A2@meta.data, vasc_label)
```


<table class="dataframe">
<caption>A data.frame: 3 × 2</caption>
<thead>
	<tr><th scope=col>vasc_label1</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>vasc_high</td><td>3169</td></tr>
	<tr><td>vasc_low </td><td>3170</td></tr>
	<tr><td>vasc_neg </td><td>3372</td></tr>
</tbody>
</table>




<table class="dataframe">
<caption>A data.frame: 2 × 2</caption>
<thead>
	<tr><th scope=col>vasc_label</th><th scope=col>n</th></tr>
	<tr><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;int&gt;</th></tr>
</thead>
<tbody>
	<tr><td>vasc_high</td><td>3169</td></tr>
	<tr><td>vasc_low </td><td>6542</td></tr>
</tbody>
</table>



# save seurat object


```R
# save seurat object for further analysis
saveRDS(AB1_A2, 'Processing/AB1_A2_Annotated_step1_20250326.rds')
```

# Save some plots


```R
colnames(AB1_A2@meta.data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'orig.ident'</li><li>'nCount_Spatial'</li><li>'nFeature_Spatial'</li><li>'visium_round'</li><li>'nCount_SCT'</li><li>'nFeature_SCT'</li><li>'SCT_snn_res.0.8'</li><li>'seurat_clusters'</li><li>'slice_ident'</li><li>'vasc_summed_expr'</li><li>'vasc_label'</li><li>'Tumor_score'</li><li>'Tumor_label1'</li><li>'Tumor_label'</li><li>'MHCI_score'</li><li>'MHCI_label1'</li><li>'MHCI_label'</li><li>'Tumor_MHCI_label'</li><li>'Tumor_MHCI_label2'</li><li>'vasc_label1'</li></ol>




```R
# /immuno/ian/Projects/2025_02_Barbie_Visium_Round1_Revisions/Output/02_Annotate_spots
```


```R
strsplit("MHCI_label", "_")[[1]][1]
strsplit("MHCI_label", "_")[[1]][2]
```


'MHCI'



'label'



```R
plot_list <- c('vasc_score', 'vasc_label1', 'vasc_label', 'Tumor_score', 'Tumor_label1', 'Tumor_label', 'MHCI_score', 'MHCI_label1', 'MHCI_label')
test_strings <- sapply(plot_list, function(x) {
  split_parts <- strsplit(x, "_")[[1]]
  remaining_parts <- split_parts[-length(split_parts)]
  paste(remaining_parts, collapse = "_")
})
test_strings
```


<style>
.dl-inline {width: auto; margin:0; padding: 0}
.dl-inline>dt, .dl-inline>dd {float: none; width: auto; display: inline-block}
.dl-inline>dt::after {content: ":\0020"; padding-right: .5ex}
.dl-inline>dt:not(:first-of-type) {padding-left: .5ex}
</style><dl class=dl-inline><dt>vasc_score</dt><dd>'vasc'</dd><dt>vasc_label1</dt><dd>'vasc'</dd><dt>vasc_label</dt><dd>'vasc'</dd><dt>Tumor_score</dt><dd>'Tumor'</dd><dt>Tumor_label1</dt><dd>'Tumor'</dd><dt>Tumor_label</dt><dd>'Tumor'</dd><dt>MHCI_score</dt><dd>'MHCI'</dd><dt>MHCI_label1</dt><dd>'MHCI'</dd><dt>MHCI_label</dt><dd>'MHCI'</dd></dl>




```R
'1' %in% 'vasc_label1'
```


FALSE



```R
grepl('1', 'vasc_label1', fixed = TRUE)
```


TRUE



```R
AB1_A2@meta.data %>% head(2)
```


<table class="dataframe">
<caption>A data.frame: 2 × 20</caption>
<thead>
	<tr><th></th><th scope=col>orig.ident</th><th scope=col>nCount_Spatial</th><th scope=col>nFeature_Spatial</th><th scope=col>visium_round</th><th scope=col>nCount_SCT</th><th scope=col>nFeature_SCT</th><th scope=col>SCT_snn_res.0.8</th><th scope=col>seurat_clusters</th><th scope=col>slice_ident</th><th scope=col>vasc_summed_expr</th><th scope=col>vasc_label</th><th scope=col>Tumor_score</th><th scope=col>Tumor_label1</th><th scope=col>Tumor_label</th><th scope=col>MHCI_score</th><th scope=col>MHCI_label1</th><th scope=col>MHCI_label</th><th scope=col>Tumor_MHCI_label</th><th scope=col>Tumor_MHCI_label2</th><th scope=col>vasc_label1</th></tr>
	<tr><th></th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>VA1_AAACAACGAATAGTTC-1</th><td>VA1</td><td>3625</td><td>2571</td><td>round1</td><td>3649</td><td>2571</td><td>10</td><td>10</td><td>VA1</td><td>2.648104</td><td>vasc_high</td><td>11</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_high</td></tr>
	<tr><th scope=row>VA1_AAACAAGTATCTCCCA-1</th><td>VA1</td><td>5804</td><td>3329</td><td>round1</td><td>4956</td><td>3308</td><td>6 </td><td>6 </td><td>VA1</td><td>1.001716</td><td>vasc_low </td><td>13</td><td>Tumor_low</td><td>Tumor_low</td><td>2</td><td>MHCI_low</td><td>MHCI_low</td><td>Tumor_low_MHCI_low</td><td>Tumor_low</td><td>vasc_low </td></tr>
</tbody>
</table>




```R
unique(AB1_A2@meta.data$vasc_label)
unique(AB1_A2@meta.data$Tumor_label1)
unique(AB1_A2@meta.data$Tumor_label)
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
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high'</li><li>'Tumor_neg'</li></ol>




<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'Tumor_low'</li><li>'Tumor_high'</li></ol>




```R
plot_list <- c('vasc_label1', 'vasc_label', 'Tumor_score', 'Tumor_label1', 'Tumor_label', 'MHCI_score', 'MHCI_label1', 'MHCI_label')

options(repr.plot.width=22, repr.plot.height=8)

for (plot_name in plot_list) {
    print('- - - - - - - - - - - - - - - - - - - - -')
    print(plot_name)
    
    clrs_use2 <- c('Tumor_low' = '#9c9c9eFF', 
                  'Tumor_high_MHCI_low' = '#F39600FF',
                  'Tumor_high_MHCI_high' = '#A60E00FF')
    
    # remove _score and _label from plot_name to make annot_name
    annot_name <- gsub('_score','',plot_name)
    annot_name <- gsub('_label1','',annot_name)
    annot_name <- gsub('_label','',annot_name)
    print(annot_name)
    
    Idents(AB1_A2) <- plot_name
    
    if (grepl('score', plot_name, fixed = TRUE)) {
        print('score')
        plot <- SpatialFeaturePlot(AB1_A2, features = c(plot_name), pt.size.factor=1.7, alpha=0.7) & guides(fill = guide_legend(override.aes = list(size=7)))
    } else if (grepl('label1', plot_name, fixed = TRUE)) {
        print('label1')
        clrs_use2 <- setNames(c('#9c9c9eFF','#F39600FF','#A60E00FF'),c(paste0(annot_name,'_neg'),paste0(annot_name,'_low'),paste0(annot_name,'_high')))
        plot <- SpatialDimPlot(AB1_A2, cols = clrs_use2, pt.size.factor=1.7, alpha=0.7) & guides(fill = guide_legend(override.aes = list(size=7)))
    } else {
        print('label')
        clrs_use2 <- setNames(c('#F39600FF','#A60E00FF'),c(paste0(annot_name,'_low'),paste0(annot_name,'_high')))
        plot <- SpatialDimPlot(AB1_A2, cols = clrs_use2, pt.size.factor=1.7, alpha=0.7) & guides(fill = guide_legend(override.aes = list(size=7)))
    }
    
    # save
    png(paste0('Output/02_Annotate_spots/',plot_name,'_SpatialDimPlot.png'),
        width = 2200,
        height = 800)
    print(plot)
    dev.off()
        
    pdf(paste0('Output/02_Annotate_spots/',plot_name,'_SpatialDimPlot.pdf'),
    width = 22,
    height = 8)
    print(plot)
    dev.off()
    
    Idents(AB1_A2) <- "orig.ident"
    
    }
```

    [1] "- - - - - - - - - - - - - - - - - - - - -"
    [1] "vasc_label1"
    [1] "vasc"
    [1] "label1"
    [1] "- - - - - - - - - - - - - - - - - - - - -"
    [1] "vasc_label"
    [1] "vasc"
    [1] "label"
    [1] "- - - - - - - - - - - - - - - - - - - - -"
    [1] "Tumor_score"
    [1] "Tumor"
    [1] "score"
    [1] "- - - - - - - - - - - - - - - - - - - - -"
    [1] "Tumor_label1"
    [1] "Tumor"
    [1] "label1"
    [1] "- - - - - - - - - - - - - - - - - - - - -"
    [1] "Tumor_label"
    [1] "Tumor"
    [1] "label"
    [1] "- - - - - - - - - - - - - - - - - - - - -"
    [1] "MHCI_score"
    [1] "MHCI"
    [1] "score"
    [1] "- - - - - - - - - - - - - - - - - - - - -"
    [1] "MHCI_label1"
    [1] "MHCI"
    [1] "label1"
    [1] "- - - - - - - - - - - - - - - - - - - - -"
    [1] "MHCI_label"
    [1] "MHCI"
    [1] "label"



```R

```


```R

```


```R

```


```R

```
