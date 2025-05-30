Step 2: Establishing a workflow for ATAC-seq-QC–Part A (step 2.1): What
seurat normalization method should be used?
================
Xenia Sirodzha
2025-02-06

- [Libraries](#libraries)
- [1. Loading and inspecting](#1-loading-and-inspecting)
- [2. Pseudo bulking and
  normalization](#2-pseudo-bulking-and-normalization)
- [3 Assessement of transformation by visualization by
  celltype](#3-assessement-of-transformation-by-visualization-by-celltype)
  - [3.1 Data wrangling for
    visualization](#31-data-wrangling-for-visualization)
  - [3.2 Visualization](#32-visualization)
    - [3.2.1 Cell counts histogram before
      transformation](#321-cell-counts-histogram-before-transformation)
    - [3.2.2 Visualization of “CLR”](#322-visualization-of-clr)
      - [3.2.2.1 Distribution after pseudo-bulking by
        celltype](#3221-distribution-after-pseudo-bulking-by-celltype)
      - [3.2.2.2 Sum of all peaks](#3222-sum-of-all-peaks)
      - [3.2.2.3 Coefficient of
        variation](#3223-coefficient-of-variation)
      - [3.2.2.4 Count distribution by cell
        type](#3224-count-distribution-by-cell-type)
    - [3.2.3 Visualization of
      “LogNormalize”](#323-visualization-of-lognormalize)
      - [3.2.3.1 Distribution after pseudo-bulking by
        celltype](#3231-distribution-after-pseudo-bulking-by-celltype)
      - [3.2.3.2 Sum of all peaks](#3232-sum-of-all-peaks)
      - [3.2.3.3 Coefficient of
        variation](#3233-coefficient-of-variation)
      - [3.2.3.4 Count distribution by cell
        type](#3234-count-distribution-by-cell-type)
- [4. Summary](#4-summary)
- [Session info](#session-info)

Here, I aimed to establish a workflow for ATAC-Seq-QC by performing
typical ATAC-seq-QC steps on a reduced set of data. All steps were
approved by Alexander Sasse.

To be precise, in this step, the data was normalized and pseudo-bulked
using Seurat’s `AggregateExpression` function. To assess the effects of
these transformations the data was partially visualized. Overall, the
goal was to determine what normalization method, which is integrated
into the `AggregateExpression` function from the Seurat package to use.

# Libraries

These are the packages required for the code below.

``` r
library(tidyverse) 
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ dplyr     1.1.4     ✔ readr     2.1.5
    ## ✔ forcats   1.0.0     ✔ stringr   1.5.1
    ## ✔ ggplot2   3.5.1     ✔ tibble    3.2.1
    ## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
    ## ✔ purrr     1.0.4     
    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(binom)
library(lubridate)
library(broom)
```

    ## Warning: Paket 'broom' wurde unter R Version 4.4.3 erstellt

``` r
library(scales) #for scale transformations
```

    ## 
    ## Attache Paket: 'scales'
    ## 
    ## Das folgende Objekt ist maskiert 'package:purrr':
    ## 
    ##     discard
    ## 
    ## Das folgende Objekt ist maskiert 'package:readr':
    ## 
    ##     col_factor

``` r
library(cowplot) #arranging plots 
```

    ## 
    ## Attache Paket: 'cowplot'
    ## 
    ## Das folgende Objekt ist maskiert 'package:lubridate':
    ## 
    ##     stamp

``` r
library(ggpubr)
```

    ## 
    ## Attache Paket: 'ggpubr'
    ## 
    ## Das folgende Objekt ist maskiert 'package:cowplot':
    ## 
    ##     get_legend

``` r
library(ggridges)
library(ggthemes) #stylising plots
```

    ## 
    ## Attache Paket: 'ggthemes'
    ## 
    ## Das folgende Objekt ist maskiert 'package:cowplot':
    ## 
    ##     theme_map

``` r
library(ggplot2)
library(ggsci) 
library(wesanderson)
library(readxl) #for loading excel files
```

    ## Warning: Paket 'readxl' wurde unter R Version 4.4.3 erstellt

``` r
library(Seurat) #for rna-seq data, atac-seq data
```

    ## Lade nötiges Paket: SeuratObject

    ## Warning: Paket 'SeuratObject' wurde unter R Version 4.4.3 erstellt

    ## Lade nötiges Paket: sp
    ## 
    ## Attache Paket: 'SeuratObject'
    ## 
    ## Die folgenden Objekte sind maskiert von 'package:base':
    ## 
    ##     intersect, t

``` r
library(Signac) #for atac seq data
```

# 1. Loading and inspecting

First, I loaded the sn-ATAC-seq Seurat Object, which contains the count
matrix of 5 cell types and 1000 randomly selected peaks. More details
about this Seurat Object can be found in
“02-06-ATAC-Seq-Analysis-Step-1.1”.

``` r
#loading file
sub_zfish_atac <- readRDS("Data/Lin_2023_zfish_snATAC-seq/sub_zfish_atac.rds")

#overview
glimpse(sub_zfish_atac) #overview over file
```

    ## Formal class 'Seurat' [package "SeuratObject"] with 13 slots
    ##   ..@ assays      :List of 1
    ##   .. ..$ atac:Formal class 'Assay5' [package "SeuratObject"] with 8 slots
    ##   ..@ meta.data   :'data.frame': 5514 obs. of  27 variables:
    ##   .. ..$ orig.ident       : Factor w/ 6 levels "10hpf","12hpf",..: 5 5 5 5 5 5 5 5 5 5 ...
    ##   .. ..$ nCount_atac      : num [1:5514] 168 181 194 137 93 120 136 103 106 113 ...
    ##   .. ..$ nFeature_atac    : int [1:5514] 83 95 91 58 50 54 61 47 50 59 ...
    ##   .. ..$ ...1             : chr [1:5514] "5hpf_1#5.25hpf_1_BC0021_N02" "5hpf_1#5.25hpf_1_BC0012_N01" "5hpf_1#5.25hpf_1_BC0018_N01" "5hpf_1#5.25hpf_1_BC0078_N02" ...
    ##   .. ..$ Sample           : chr [1:5514] "5hpf_1" "5hpf_1" "5hpf_1" "5hpf_1" ...
    ##   .. ..$ TSSEnrichment    : num [1:5514] 8.7 10.07 9.01 8.62 8.5 ...
    ##   .. ..$ ReadsInTSS       : num [1:5514] 3175 3524 2745 2755 1652 ...
    ##   .. ..$ ReadsInPromoter  : num [1:5514] 14232 14938 11951 11721 7048 ...
    ##   .. ..$ ReadsInBlacklist : num [1:5514] 0 0 0 0 0 0 0 0 0 0 ...
    ##   .. ..$ PromoterRatio    : num [1:5514] 0.0787 0.0861 0.0835 0.0988 0.0615 ...
    ##   .. ..$ PassQC           : num [1:5514] 1 1 1 1 1 1 1 1 1 1 ...
    ##   .. ..$ NucleosomeRatio  : num [1:5514] 1.16 1.23 1.08 1.22 1.2 ...
    ##   .. ..$ nMultiFrags      : num [1:5514] 17542 13440 12844 11903 11281 ...
    ##   .. ..$ nMonoFrags       : num [1:5514] 41913 38847 34432 26658 26007 ...
    ##   .. ..$ nFrags           : num [1:5514] 90396 86794 71572 59314 57328 ...
    ##   .. ..$ nDiFrags         : num [1:5514] 30941 34507 24296 20753 20040 ...
    ##   .. ..$ BlacklistRatio   : num [1:5514] 0 0 0 0 0 0 0 0 0 0 ...
    ##   .. ..$ stage            : chr [1:5514] "5hpf" "5hpf" "5hpf" "5hpf" ...
    ##   .. ..$ Clusters         : chr [1:5514] "C1" "C1" "C1" "C1" ...
    ##   .. ..$ ReadsInPeaks     : num [1:5514] 78617 77133 68699 61051 46438 ...
    ##   .. ..$ FRIP             : num [1:5514] 0.435 0.444 0.48 0.515 0.405 ...
    ##   .. ..$ celltype         : chr [1:5514] "YSL/presumptive endoderm" "YSL/presumptive endoderm" "YSL/presumptive endoderm" "YSL/presumptive endoderm" ...
    ##   .. ..$ predictedCell    : chr [1:5514] "6h_3 CELL4337_N1 _" "6h_2_CELL3455_N1" "6h_2_CELL3657_N1" "LC-5.3h-2_CELL2003_N1_LC-5.3h-2" ...
    ##   .. ..$ predictedGroup   : chr [1:5514] "YSL/endoderm" "YSL/endoderm" "YSL/endoderm" "YSL/endoderm" ...
    ##   .. ..$ predictedScore   : num [1:5514] 0.919 0.949 0.806 0.631 0.92 ...
    ##   .. ..$ DoubletScore     : num [1:5514] 22.2 44.5 323.3 293.3 0 ...
    ##   .. ..$ DoubletEnrichment: num [1:5514] 3.3 4.4 51.5 12.4 0.3 1.4 2.1 0.2 24.1 11.6 ...
    ##   ..@ active.assay: chr "atac"
    ##   ..@ active.ident: Factor w/ 6 levels "10hpf","12hpf",..: 5 5 5 5 5 5 5 5 5 5 ...
    ##   .. ..- attr(*, "names")= chr [1:5514] "5hpf_1#5.25hpf_1_BC0021_N02" "5hpf_1#5.25hpf_1_BC0012_N01" "5hpf_1#5.25hpf_1_BC0018_N01" "5hpf_1#5.25hpf_1_BC0078_N02" ...
    ##   ..@ graphs      : list()
    ##   ..@ neighbors   : list()
    ##   ..@ reductions  : list()
    ##   ..@ images      : list()
    ##   ..@ project.name: chr "SeuratProject"
    ##   ..@ misc        : list()
    ##   ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
    ##   .. ..$ : int [1:3] 5 0 2
    ##   ..@ commands    : list()
    ##   ..@ tools       : list()

``` r
sub_zfish_atac@assays$atac$counts[1:10,1:2] #overview count matrix
```

    ## 10 x 2 sparse Matrix of class "dgCMatrix"
    ##                         5hpf_1#5.25hpf_1_BC0021_N02 5hpf_1#5.25hpf_1_BC0012_N01
    ## chr7:20868664-20869164                            1                           .
    ## chr25:7596198-7596698                             .                           .
    ## chr24:28306284-28306784                           2                           .
    ## chr14:5775997-5776497                             .                           .
    ## chr14:21287634-21288134                           .                           .
    ## chr3:2152506-2153006                              .                           1
    ## chr19:3481940-3482440                             .                           .
    ## chr18:918337-918837                               .                           .
    ## chr7:71301133-71301633                            .                           .
    ## chr11:14025086-14025586                           .                           .

Overall, this Seurat object contains both the count matrix and the meta
data.

# 2. Pseudo bulking and normalization

At the start, I performed pseudo-bulking using Seurat’s
`AggregateExpression` function. As, I was not sure which Seurat provided
normalization method to use, I performed both logarithmic and
centralized logarithmic normalization and looked whether and how the
output differed.

To ensure that the randomized normalizations from Seurat are
reproducible, I set a seed.

``` r
## ensuring always the same random result
set.seed(59762489)
 
#testing if setting the seed worked
sample(1000,1) #running multiple times, same result
```

    ## [1] 684

Next, I pseudo-bulked and normalized the Seurat Object `sub_zfish_atac`.
The values were aggregated by cell type.

``` r
# Normalizing with Seurat
psd.bulk.sub_atc.clr= AggregateExpression(sub_zfish_atac, group.by = "celltype", normalization.method = "CLR", return.seurat = T)
```

    ## Normalizing across features

    ## Centering and scaling data matrix

``` r
psd.bulk.sub_atc.lgn= AggregateExpression(sub_zfish_atac, group.by = "celltype", normalization.method = "LogNormalize", return.seurat = T)
```

    ## Centering and scaling data matrix

Then, I checked whether the pseudo-bulking worked as intended and
whether the output differed. I expect that there will be only 5 columns
left for the 5 cell types and that the output will differ

``` r
#before pseudo bulking
dim(sub_zfish_atac@assays$atac) 
```

    ## [1] 1000 5514

``` r
#after pseudo bulking
dim(psd.bulk.sub_atc.clr$atac$counts) #values are stored in counts layer
```

    ## [1] 1000    5

``` r
dim(psd.bulk.sub_atc.lgn$atac$counts) 
```

    ## [1] 1000    5

``` r
#testing if the output is identical
identical(psd.bulk.sub_atc.clr@assays$atac$data, psd.bulk.sub_atc.lgn@assays$atac$data)
```

    ## [1] FALSE

``` r
#overview output
glimpse(psd.bulk.sub_atc.clr)
```

    ## Formal class 'Seurat' [package "SeuratObject"] with 13 slots
    ##   ..@ assays      :List of 1
    ##   .. ..$ atac:Formal class 'Assay5' [package "SeuratObject"] with 8 slots
    ##   ..@ meta.data   :'data.frame': 5 obs. of  2 variables:
    ##   .. ..$ orig.ident: chr [1:5] "lateral plate mesoderm" "neural crest" "periderm/epidermis" "primary neuron" ...
    ##   .. ..$ celltype  : chr [1:5] "lateral plate mesoderm" "neural crest" "periderm/epidermis" "primary neuron" ...
    ##   ..@ active.assay: chr "atac"
    ##   ..@ active.ident: Factor w/ 5 levels "lateral plate mesoderm",..: 1 2 3 4 5
    ##   .. ..- attr(*, "names")= chr [1:5] "lateral plate mesoderm" "neural crest" "periderm/epidermis" "primary neuron" ...
    ##   ..@ graphs      : list()
    ##   ..@ neighbors   : list()
    ##   ..@ reductions  : list()
    ##   ..@ images      : list()
    ##   ..@ project.name: chr "Aggregate"
    ##   ..@ misc        : list()
    ##   ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
    ##   .. ..$ : int [1:3] 5 0 2
    ##   ..@ commands    :List of 1
    ##   .. ..$ ScaleData.atac:Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
    ##   ..@ tools       : list()

``` r
glimpse(psd.bulk.sub_atc.lgn)
```

    ## Formal class 'Seurat' [package "SeuratObject"] with 13 slots
    ##   ..@ assays      :List of 1
    ##   .. ..$ atac:Formal class 'Assay5' [package "SeuratObject"] with 8 slots
    ##   ..@ meta.data   :'data.frame': 5 obs. of  2 variables:
    ##   .. ..$ orig.ident: chr [1:5] "lateral plate mesoderm" "neural crest" "periderm/epidermis" "primary neuron" ...
    ##   .. ..$ celltype  : chr [1:5] "lateral plate mesoderm" "neural crest" "periderm/epidermis" "primary neuron" ...
    ##   ..@ active.assay: chr "atac"
    ##   ..@ active.ident: Factor w/ 5 levels "lateral plate mesoderm",..: 1 2 3 4 5
    ##   .. ..- attr(*, "names")= chr [1:5] "lateral plate mesoderm" "neural crest" "periderm/epidermis" "primary neuron" ...
    ##   ..@ graphs      : list()
    ##   ..@ neighbors   : list()
    ##   ..@ reductions  : list()
    ##   ..@ images      : list()
    ##   ..@ project.name: chr "Aggregate"
    ##   ..@ misc        : list()
    ##   ..@ version     :Classes 'package_version', 'numeric_version'  hidden list of 1
    ##   .. ..$ : int [1:3] 5 0 2
    ##   ..@ commands    :List of 1
    ##   .. ..$ ScaleData.atac:Formal class 'SeuratCommand' [package "SeuratObject"] with 5 slots
    ##   ..@ tools       : list()

The outputs meet the expectations.

# 3 Assessement of transformation by visualization by celltype

## 3.1 Data wrangling for visualization

Next, to assess if the transformation is suitable for quality control, I
analysed and visualised the data.

To achieve this, the normalised counts needed to be extracted and
converted into a ggplot tidyverse compatible format: a data frame. I
performed analysis and visualization for both counts in an analogical
manner.

First, I calculated the sum of all peaks and the coefficient of variance
(cv).

``` r
## extracting tibbles 

# CLR
df_psd.bulk_atc.clr=as_tibble(psd.bulk.sub_atc.clr@assays$atac$data, rownames="rownames")

# LogNormalize

df_psd.bulk_atc.lgn=as_tibble(psd.bulk.sub_atc.lgn@assays$atac$data, rownames="rownames")


## calculating sum over all peaks and coefficient of variance of all peaks

# CLR:
df_psd.bulk_atc.clr %>% 
  pivot_longer(cols=!contains("rownames"), ,values_to = "pk_counts", names_to = "celltype") %>% 
  mutate(pk_count_sum= sum(pk_counts)) %>% 
  mutate(cv=sd(pk_counts) / mean(pk_counts)) %>% 
  group_by(celltype) %>% 
  mutate(pk_count_sum_celltype= sum(pk_counts)) %>% 
  mutate(cv_celltype=sd(pk_counts) / mean(pk_counts)) %>% 
  mutate(all="all celltypes")-> df_psd.bulk_atc.clr

# LogNormalize

df_psd.bulk_atc.lgn %>% 
  pivot_longer(cols=!contains("rownames"), ,values_to = "pk_counts", names_to = "celltype") %>% 
  mutate(pk_count_sum= sum(pk_counts)) %>% 
  mutate(cv=sd(pk_counts) / mean(pk_counts)) %>% 
  group_by(celltype) %>% 
  mutate(pk_count_sum_celltype= sum(pk_counts)) %>% 
  mutate(cv_celltype=sd(pk_counts) / mean(pk_counts)) %>% 
  mutate(all="all celltypes")-> df_psd.bulk_atc.lgn
```

This code results in two data sets which have identical column names.

``` r
head(df_psd.bulk_atc.clr)
```

    ## # A tibble: 6 × 8
    ## # Groups:   celltype [5]
    ##   rownames           celltype pk_counts pk_count_sum    cv pk_count_sum_celltype
    ##   <chr>              <chr>        <dbl>        <dbl> <dbl>                 <dbl>
    ## 1 chr7:20868664-208… lateral…     0.653        3511. 0.586                  563.
    ## 2 chr7:20868664-208… neural …     0.737        3511. 0.586                  703.
    ## 3 chr7:20868664-208… perider…     0.612        3511. 0.586                  619.
    ## 4 chr7:20868664-208… primary…     0.542        3511. 0.586                  737.
    ## 5 chr7:20868664-208… YSL/pre…     0.959        3511. 0.586                  889.
    ## 6 chr25:7596198-759… lateral…     0.769        3511. 0.586                  563.
    ## # ℹ 2 more variables: cv_celltype <dbl>, all <chr>

``` r
head(df_psd.bulk_atc.lgn)
```

    ## # A tibble: 6 × 8
    ## # Groups:   celltype [5]
    ##   rownames           celltype pk_counts pk_count_sum    cv pk_count_sum_celltype
    ##   <chr>              <chr>        <dbl>        <dbl> <dbl>                 <dbl>
    ## 1 chr7:20868664-208… lateral…      4.65        8934. 0.542                 1730.
    ## 2 chr7:20868664-208… neural …      4.63        8934. 0.542                 1833.
    ## 3 chr7:20868664-208… perider…      4.53        8934. 0.542                 1767.
    ## 4 chr7:20868664-208… primary…      4.23        8934. 0.542                 1847.
    ## 5 chr7:20868664-208… YSL/pre…      4.62        8934. 0.542                 1757.
    ## 6 chr25:7596198-759… lateral…      1.76        8934. 0.542                 1730.
    ## # ℹ 2 more variables: cv_celltype <dbl>, all <chr>

To plot a histogram, I also require the counts of all cells in a long
format. This is achieved with the code below.

``` r
df_pk_mtrx= as_tibble(sub_zfish_atac$atac$counts, rownames="rownames")
df_pk_mtrx %>% pivot_longer(cols=!contains("rownames"), ,values_to = "pk_counts", names_to = "id") -> df_pk_mtrx.l
head(df_pk_mtrx.l)
```

    ## # A tibble: 6 × 3
    ##   rownames               id                          pk_counts
    ##   <chr>                  <chr>                           <dbl>
    ## 1 chr7:20868664-20869164 5hpf_1#5.25hpf_1_BC0021_N02         1
    ## 2 chr7:20868664-20869164 5hpf_1#5.25hpf_1_BC0012_N01         0
    ## 3 chr7:20868664-20869164 5hpf_1#5.25hpf_1_BC0018_N01         2
    ## 4 chr7:20868664-20869164 5hpf_1#5.25hpf_1_BC0078_N02         1
    ## 5 chr7:20868664-20869164 5hpf_1#5.25hpf_1_BC0044_N01         0
    ## 6 chr7:20868664-20869164 5hpf_1#5.25hpf_1_BC0065_N02         0

## 3.2 Visualization

Here, the data was visualized to assess what normalization method to use
for ATAC-Seq-QC. All plots were generated using the tidyverse ggplot
package.

### 3.2.1 Cell counts histogram before transformation

I plotted counts before transformation as a histogram.

``` r
df_pk_mtrx.l %>% 
  ggplot(aes(x=pk_counts))+geom_histogram(binwidth = 0.7, color="gray15", alpha=0.8)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none")
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/histogram%20all%20cells%20before%20transformation-1.png)<!-- -->

Most values were 0 therefore, I excluded them before plotting again.

``` r
df_pk_mtrx.l %>% 
  filter(pk_counts!=0) %>% 
  ggplot(aes(x=pk_counts))+geom_histogram(binwidth = 0.5, color="gray15", alpha=0.8)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none")
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/histogram%20all%20cells%20before%20transformation%20without%200s-1.png)<!-- -->

Overall, most peaks have a count of 2 and the distribution is close to a
normal distribution.

### 3.2.2 Visualization of “CLR”

In this section, the CLR normalized data was visualized.

#### 3.2.2.1 Distribution after pseudo-bulking by celltype

Counts were plotted in a histogram overall and then by cell type.

``` r
# Overall

df_psd.bulk_atc.clr %>% 
 ggplot(aes(x=pk_counts))+geom_histogram(binwidth = 0.2, color="gray15", alpha=0.8)+
 scale_y_continuous(breaks = scales::breaks_pretty(n=10))+
  coord_fixed(ratio = 0.01)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") ->dis.counts.clr


# by cell type

 ggplot(df_psd.bulk_atc.clr, aes(x=pk_counts, fill = celltype))+
  geom_histogram(binwidth = 0.2, color="gray15", alpha=0.8)+
  scale_y_continuous(breaks = scales::breaks_pretty(n=10))+
  facet_wrap(facets = "celltype")+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") ->dis.counts.cell.clr
```

They were arranged below.

``` r
dis.counts.clr|dis.counts.cell.clr
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/vis%20clr:%20displaying%20histograms%20after%20transformation-1.png)<!-- -->

Generally, counts are somewhat normally distributed. There are no cell
type specific trends. However, most counts are between the values 0 and
1.5 which is not typical for well normalized ATAC-seq-data.

#### 3.2.2.2 Sum of all peaks

I counted and visualized the sum of all peaks by cell type and in
general.

In addition, I extracted the value for the sum of all peaks. It is
3527.945.

``` r
unique(df_psd.bulk_atc.clr$pk_count_sum)
```

    ## [1] 3510.855

The sums of all peaks were visualized in a barplot, analogically to the
visualization by histogram

``` r
# Overall

df_psd.bulk_atc.clr %>% 
  ungroup() %>% 
  distinct(all, pk_count_sum) %>% 
  ggplot(aes(x=all, y=pk_count_sum))+
  geom_col(alpha=0.8, color="gray10")+
  coord_fixed(ratio = 0.0005)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> sum.peaks.all.clr


# by cell type

df_psd.bulk_atc.clr %>% 
  group_by(celltype) %>% 
  distinct(pk_count_sum_celltype, pk_count_sum) %>% 
  ggplot(aes(x=celltype, y=pk_count_sum_celltype, fill=celltype))+
  geom_col(alpha=0.8, color="gray10")+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> sum.peaks.celltype.clr

sum.peaks.all.clr|sum.peaks.celltype.clr
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/vis%20clr:%20bar%20plot%20sum%20of%20all%20peaks-1.png)<!-- -->

The YSL/presumptive mesoderm cells have the highest sum of peaks with
around 900. But overall they do not differ much.

#### 3.2.2.3 Coefficient of variation

The coefficient of variation was visualized as the sum of all peaks
previously.

``` r
# Overall

df_psd.bulk_atc.clr %>% 
  ungroup() %>% 
  distinct(all, cv) %>% 
  ggplot(aes(x=all, y=cv ))+
  geom_col(alpha=0.8, color="gray10")+
  scale_y_continuous(n.breaks=10)+
  coord_fixed(ratio = 5)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> cv.all.clr


# by cell type

df_psd.bulk_atc.clr %>% 
  group_by(celltype) %>% 
  distinct(cv_celltype, cv) %>% 
  ggplot(aes(x=celltype, y=cv_celltype, fill=celltype))+
  geom_col(alpha=0.8, color="gray10")+
  scale_y_continuous(n.breaks=10)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") ->cv.celltype.clr

cv.all.clr|cv.celltype.clr
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/vis%20clr:%20bar%20plot%20coefficient%20of%20variation-1.png)<!-- -->

Overall, coefficients of variation of each cell type are very similar to
one another and to the coefficient of variation of all cells. However,
periderm/epidermis cells have the highest coefficient of variation.

#### 3.2.2.4 Count distribution by cell type

Then, just as previously count distribution was visualized with a violin
plot.

``` r
df_psd.bulk_atc.clr %>% 
  mutate(all="all") %>% 
 ggplot( aes( x=all, y=pk_counts))+
  geom_violin(color="grey20", alpha=0.3, fill="gray30", trim=F)+
  coord_fixed(ratio = 0.5)+
   scale_y_continuous(breaks = seq(0,10,1))+ #y-axis more easily readable
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> v.plot.all.clr


 ggplot(df_psd.bulk_atc.clr, aes(x=celltype, y=pk_counts, fill = celltype))+
  geom_violin(color="grey20", alpha=0.3, trim=F)+
    scale_y_continuous(breaks = seq(0,10,1))+ #y-axis more easily readable
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> v.plot.cells.clr
 
 v.plot.all.clr|v.plot.cells.clr
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/vis%20clr:%20violin%20plot%20counts-1.png)<!-- -->

To visualise median and quantiles, a boxplot was added.

``` r
df_psd.bulk_atc.clr %>% 
  mutate(all="all") %>% 
 ggplot( aes( x=all, y=pk_counts))+
  geom_violin(color="grey20", alpha=0.3, fill="gray30", trim=F)+
  geom_boxplot(width = 0.07)+
  stat_boxplot(geom ='errorbar', width=0.2, position=position_dodge(0.5))+
   scale_y_continuous(breaks = seq(0,10,1))+ #y-axis more easily readable
  coord_fixed(ratio = 0.5)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> v.plot.all.box.clr


 ggplot(df_psd.bulk_atc.clr, aes(x=celltype, y=pk_counts, fill = celltype))+
  geom_violin(color="grey20", alpha=0.3, trim=F)+
  geom_boxplot(width = 0.07)+
  stat_boxplot(geom ='errorbar', width=0.2, position=position_dodge(0.5))+
    scale_y_continuous(breaks = seq(0,10,1))+ #y-axis more easily readable
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> v.plot.cells.box.clr
 
 v.plot.all.box.clr|v.plot.cells.box.clr
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/vis%20clr:%20violin%20plot%20counts%20plus%20boxplot-1.png)<!-- -->

The distribution of counts of all celltypes is very similar to each
other and the distribution of counts overall. The YSL/presumptive
endoderm cells have slightly higher counts overall. Most counts range
from 0.5 to 1 and are centered around approximately 0.6.

### 3.2.3 Visualization of “LogNormalize”

In this section, the LogNormalize normalized data was visualized
analogically to CLR normalized data (chapter 3.2.2).

#### 3.2.3.1 Distribution after pseudo-bulking by celltype

Counts were plotted in a histogram overall and then by cell type.

``` r
# Overall 

df_psd.bulk_atc.lgn %>% 
 ggplot(aes(x=pk_counts))+geom_histogram(binwidth = 0.2, color="gray15", alpha=0.8)+
 scale_y_continuous(breaks = scales::breaks_pretty(n=10))+
  coord_fixed(ratio = 0.025)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") ->dis.counts.lgn


# by cell type

 ggplot(df_psd.bulk_atc.lgn, aes(x=pk_counts, fill = celltype))+
  geom_histogram(binwidth = 0.2, color="gray15", alpha=0.8)+
  scale_y_continuous(breaks = scales::breaks_pretty(n=10))+
  facet_wrap(facets = "celltype")+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") ->dis.counts.cell.lgn

dis.counts.lgn|dis.counts.cell.lgn
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/vis%20lgn:%20histograms%20after%20transformation-1.png)<!-- -->

#### 3.2.3.2 Sum of all peaks

The sum of all peaks is 9091.372, which is much higher than when the
data is normalized with CLR.

``` r
unique(df_psd.bulk_atc.lgn$pk_count_sum)
```

    ## [1] 8934.439

Bar plot sum of all peaks.

``` r
# Overall

df_psd.bulk_atc.lgn %>% 
  ungroup() %>% 
  distinct(all, pk_count_sum) %>% 
  ggplot(aes(x=all, y=pk_count_sum))+
  geom_col(alpha=0.8, color="gray10")+
  coord_fixed(ratio = 0.0005)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> sum.peaks.all.lgn


#by cell type

df_psd.bulk_atc.lgn %>% 
  group_by(celltype) %>% 
  distinct(pk_count_sum_celltype, pk_count_sum) %>% 
  ggplot(aes(x=celltype, y=pk_count_sum_celltype, fill=celltype))+
  geom_col(alpha=0.8, color="gray10")+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> sum.peaks.celltype.lgn

sum.peaks.all.lgn|sum.peaks.celltype.lgn
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/vis%20lgn:%20bar%20plot%20sum%20of%20all%20peaks-1.png)<!-- -->

Overall, the sum of all peaks is much higher than when the data is
normalized with CLR. The sum of all peaks by cell type differs less
between each cell type than in CLR-normalized data.

#### 3.2.3.3 Coefficient of variation

Bar plot coefficient of variation.

``` r
# Overall

df_psd.bulk_atc.lgn %>% 
  ungroup() %>% 
  distinct(all, cv) %>% 
  ggplot(aes(x=all, y=cv ))+
  geom_col(alpha=0.8, color="gray10")+
  scale_y_continuous(n.breaks=10)+
  coord_fixed(ratio = 5)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> cv.all.lgn


# by cell type

df_psd.bulk_atc.lgn %>% 
  group_by(celltype) %>% 
  distinct(cv_celltype, cv) %>% 
  ggplot(aes(x=celltype, y=cv_celltype, fill=celltype))+
  geom_col(alpha=0.8, color="gray10")+
  scale_y_continuous(n.breaks=10)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") ->cv.celltype.lgn

cv.all.lgn|cv.celltype.lgn
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/vis%20lgn:%20bar%20plot%20coefficient%20of%20variation-1.png)<!-- -->

The coefficients of variation are very similar to CLR-normalized
coefficients of variation although they are slightly reduced overall.

#### 3.2.3.4 Count distribution by cell type

Violin plots

``` r
# Overall 

df_psd.bulk_atc.lgn %>% 
  mutate(all="all") %>% 
 ggplot( aes( x=all, y=pk_counts))+
  geom_violin(color="grey20", alpha=0.3, fill="gray30", trim=F)+
  coord_fixed(ratio = 0.5)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> v.plot.all.lgn


# by cell type

 ggplot(df_psd.bulk_atc.lgn, aes(x=celltype, y=pk_counts, fill = celltype))+
  geom_violin(color="grey20", alpha=0.3, trim=F)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> v.plot.cells.lgn
 
 v.plot.all.lgn|v.plot.cells.lgn
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/vis%20lgn:%20violin%20plot%20counts-1.png)<!-- -->

Violin plots with a boxplot inside.

``` r
# Overall

df_psd.bulk_atc.lgn %>% 
  mutate(all="all") %>% 
 ggplot( aes( x=all, y=pk_counts))+
  geom_violin(color="grey20", alpha=0.3, fill="gray30", trim=F)+
  geom_boxplot(width = 0.07)+
  stat_boxplot(geom ='errorbar', width=0.2, position=position_dodge(0.5))+
   scale_y_continuous(breaks = seq(0,10,1))+ #y-axis more easily readable
  coord_fixed(ratio = 0.5)+
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> v.plot.all.box.lgn


# by cell type

 ggplot(df_psd.bulk_atc.lgn, aes(x=celltype, y=pk_counts, fill = celltype))+
  geom_violin(color="grey20", alpha=0.3, trim=F)+
  geom_boxplot(width = 0.07)+
  stat_boxplot(geom ='errorbar', width=0.2, position=position_dodge(0.5))+
    scale_y_continuous(breaks = seq(0,10,1))+ #y-axis more easily readable
  theme_bw()+ theme(text = element_text(size = 9, family="sans", color="black"), #stylisation
        axis.text = element_text(size=9, color="black"),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title=element_text(size=9, face = "bold"),
        plot.subtitle=element_text(size=10.5),
        plot.title=element_text(size=11, face = "bold"), 
        legend.position ="none") -> v.plot.cells.box.lgn
 
 v.plot.all.box.lgn|v.plot.cells.box.lgn
```

![](02-06-ATAC-Seq-Analysis-Step-2.1-seurat-norm_files/figure-gfm/vis%20lgn:%20violin%20plot%20counts%20plus%20boxplot-1.png)<!-- -->

The counts are similarily distributed, however values are more centered
around count value 1.5. Most counts range from 1 to 3. In contrast to
CLR-normalized data the YSL/presumptive endoderm cells do not have
slightly higher counts here.

# 4. Summary

For the most part visualization itself works well. Pseudobulking using
the Seurat package is successful as well. However, histograms and violin
plots of both CLR and LogNormalize from the seurat package show that the
data was not sufficiently normalized as most counts are in a range of
0.5 to 2 or 1 to 3. Whereas if the counts were well normalized one would
expect most of them to range from 1 to 4. Moreover, out of both Seurat
provided options “LogNormalize” yields the results closest to the
expectation of well normalized ATAC-seq data.

# Session info

``` r
sessionInfo()
```

    ## R version 4.4.2 (2024-10-31 ucrt)
    ## Platform: x86_64-w64-mingw32/x64
    ## Running under: Windows 11 x64 (build 26100)
    ## 
    ## Matrix products: default
    ## 
    ## 
    ## locale:
    ## [1] LC_COLLATE=German_Germany.utf8  LC_CTYPE=German_Germany.utf8   
    ## [3] LC_MONETARY=German_Germany.utf8 LC_NUMERIC=C                   
    ## [5] LC_TIME=German_Germany.utf8    
    ## 
    ## time zone: Europe/Berlin
    ## tzcode source: internal
    ## 
    ## attached base packages:
    ## [1] stats     graphics  grDevices utils     datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] Signac_1.14.0      Seurat_5.2.1       SeuratObject_5.0.2 sp_2.2-0          
    ##  [5] readxl_1.4.5       wesanderson_0.3.7  ggsci_3.2.0        ggthemes_5.1.0    
    ##  [9] ggridges_0.5.6     ggpubr_0.6.0       cowplot_1.1.3      scales_1.3.0      
    ## [13] broom_1.0.8        binom_1.1-1.1      lubridate_1.9.4    forcats_1.0.0     
    ## [17] stringr_1.5.1      dplyr_1.1.4        purrr_1.0.4        readr_2.1.5       
    ## [21] tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1      tidyverse_2.0.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RcppAnnoy_0.0.22        splines_4.4.2           later_1.4.1            
    ##   [4] bitops_1.0-9            cellranger_1.1.0        polyclip_1.10-7        
    ##   [7] fastDummies_1.7.5       lifecycle_1.0.4         rstatix_0.7.2          
    ##  [10] globals_0.17.0          lattice_0.22-6          MASS_7.3-64            
    ##  [13] backports_1.5.0         magrittr_2.0.3          plotly_4.10.4          
    ##  [16] rmarkdown_2.29          yaml_2.3.10             httpuv_1.6.15          
    ##  [19] sctransform_0.4.1       spam_2.11-1             spatstat.sparse_3.1-0  
    ##  [22] reticulate_1.40.0       pbapply_1.7-2           RColorBrewer_1.1-3     
    ##  [25] abind_1.4-8             zlibbioc_1.52.0         Rtsne_0.17             
    ##  [28] GenomicRanges_1.58.0    BiocGenerics_0.52.0     GenomeInfoDbData_1.2.13
    ##  [31] IRanges_2.40.1          S4Vectors_0.44.0        ggrepel_0.9.6          
    ##  [34] irlba_2.3.5.1           listenv_0.9.1           spatstat.utils_3.1-2   
    ##  [37] goftest_1.2-3           RSpectra_0.16-2         spatstat.random_3.3-2  
    ##  [40] fitdistrplus_1.2-2      parallelly_1.43.0       codetools_0.2-20       
    ##  [43] RcppRoll_0.3.1          tidyselect_1.2.1        UCSC.utils_1.2.0       
    ##  [46] farver_2.1.2            matrixStats_1.5.0       stats4_4.4.2           
    ##  [49] spatstat.explore_3.3-4  jsonlite_1.8.9          progressr_0.15.1       
    ##  [52] Formula_1.2-5           survival_3.8-3          tools_4.4.2            
    ##  [55] ica_1.0-3               Rcpp_1.0.14             glue_1.8.0             
    ##  [58] gridExtra_2.3           xfun_0.52               GenomeInfoDb_1.42.3    
    ##  [61] withr_3.0.2             fastmap_1.2.0           digest_0.6.37          
    ##  [64] timechange_0.3.0        R6_2.6.1                mime_0.12              
    ##  [67] colorspace_2.1-1        scattermore_1.2         tensor_1.5             
    ##  [70] spatstat.data_3.1-6     utf8_1.2.4              generics_0.1.3         
    ##  [73] data.table_1.16.4       httr_1.4.7              htmlwidgets_1.6.4      
    ##  [76] uwot_0.2.2              pkgconfig_2.0.3         gtable_0.3.6           
    ##  [79] lmtest_0.9-40           XVector_0.46.0          htmltools_0.5.8.1      
    ##  [82] carData_3.0-5           dotCall64_1.2           png_0.1-8              
    ##  [85] spatstat.univar_3.1-1   knitr_1.50              rstudioapi_0.17.1      
    ##  [88] tzdb_0.4.0              reshape2_1.4.4          nlme_3.1-167           
    ##  [91] zoo_1.8-12              KernSmooth_2.23-26      parallel_4.4.2         
    ##  [94] miniUI_0.1.1.1          pillar_1.10.2           grid_4.4.2             
    ##  [97] vctrs_0.6.5             RANN_2.6.2              promises_1.3.2         
    ## [100] car_3.1-3               xtable_1.8-4            cluster_2.1.8          
    ## [103] evaluate_1.0.3          cli_3.6.1               compiler_4.4.2         
    ## [106] Rsamtools_2.22.0        rlang_1.1.4             crayon_1.5.3           
    ## [109] future.apply_1.11.3     ggsignif_0.6.4          labeling_0.4.3         
    ## [112] plyr_1.8.9              stringi_1.8.4           viridisLite_0.4.2      
    ## [115] deldir_2.0-4            BiocParallel_1.40.2     munsell_0.5.1          
    ## [118] Biostrings_2.74.1       lazyeval_0.2.2          spatstat.geom_3.3-5    
    ## [121] Matrix_1.7-1            RcppHNSW_0.6.0          hms_1.1.3              
    ## [124] patchwork_1.3.0         future_1.40.0           shiny_1.10.0           
    ## [127] ROCR_1.0-11             igraph_2.1.4            fastmatch_1.1-6
