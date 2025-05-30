Step 1
================
Xenia Sirodzha
2025-02-06

- [Libraries](#libraries)
- [1. ATAC-seq data loading and
  inspection](#1-atac-seq-data-loading-and-inspection)
  - [1.1 Loading the data set into R](#11-loading-the-data-set-into-r)
- [2. Data transformation](#2-data-transformation)
  - [2.1 Reducing the data sets to run
    locally](#21-reducing-the-data-sets-to-run-locally)
- [Side-note: How to install Signac](#side-note-how-to-install-signac)
- [Session info](#session-info)

Overall, the goal here was it to explore the data set and reduce it to
run the subsequent ATAC-Seq Quality Control steps (ATAC-Seq-QC) locally.

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
library(readxl) #for loading excel files
```

    ## Warning: Paket 'readxl' wurde unter R Version 4.4.3 erstellt

``` r
library(ggsci)
library(Seurat) #for RNA-seq and ATAC-seq data
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

# 1. ATAC-seq data loading and inspection

## 1.1 Loading the data set into R

All data was sourced from: Lin X, Yang X, Chen C, et al. Single-nucleus
chromatin landscapes during zebrafish early embryogenesis. Sci Data.
2023;10(1):464. <doi:10.1038/s41597-023-02373-y> . It was kindly
compiled and pre-processed by Lauren Saunders.

First, the data was loaded into R using the code below.

``` r
readRDS("Data/Lin_2023_zfish_snATAC-seq/Cell-type_Peak_Matrix.rds") -> zfish_snATAC_seq_pk_mtrx #count matrix 
read_csv("Data/Lin_2023_zfish_snATAC-seq/atac_all.metaData.csv")-> zfish_mta_data # meta data
```

    ## New names:
    ## Rows: 50637 Columns: 24
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (7): ...1, Sample, stage, Clusters, celltype, predictedCell, predictedG... dbl
    ## (17): TSSEnrichment, ReadsInTSS, ReadsInPromoter, ReadsInBlacklist, Prom...
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

``` r
read_csv("Data/Lin_2023_zfish_snATAC-seq/SupplT2.csv")-> zfish_snATAC_seq_cell_type_mrkrs #markers for cell types
```

    ## New names:
    ## Rows: 8684 Columns: 12
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (3): group_name, seqnames, name dbl (9): ...1, group, start, end, strand, idx,
    ## Log2FC, FDR, MeanDiff
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

Then, the count matrix was explored:

``` r
class(zfish_snATAC_seq_pk_mtrx) #what class does the object have
```

    ## [1] "dgCMatrix"
    ## attr(,"package")
    ## [1] "Matrix"

``` r
dim(zfish_snATAC_seq_pk_mtrx) #gives dimensions
```

    ## [1] 370058  50637

``` r
zfish_snATAC_seq_pk_mtrx[1:10, 1:2]
```

    ## 10 x 2 sparse Matrix of class "dgCMatrix"
    ##                  24hpf_1#24hpf_1_BC00224_N06 24hpf_1#24hpf_1_BC00014_N03
    ## chr1:5232-5732                             .                           .
    ## chr1:5787-6287                             .                           .
    ## chr1:10088-10588                           .                           .
    ## chr1:10991-11491                           .                           .
    ## chr1:11895-12395                           .                           2
    ## chr1:12474-12974                           1                           .
    ## chr1:14016-14516                           .                           .
    ## chr1:14704-15204                           1                           .
    ## chr1:16672-17172                           3                           .
    ## chr1:18404-18904                           3                           .

The object `zfish_snATAC_seq_pk_mtrx` is a compressed matrix consisting
of 370058 rows, whose names represent the locations on the chromosome
(“peak”), and 50637 columns, whose names represent all single cell
samples, which were collected. The sample names include information
about how many hours post fertilisation (hpf) the embryos were
collected.

The data frame `zfish_snATAC_seq_pk_mtrx` contains the meta data
information for every sample of the count matrix, such as cell type,
reads in peaks and many more.

``` r
class(zfish_mta_data)
```

    ## [1] "spec_tbl_df" "tbl_df"      "tbl"         "data.frame"

``` r
dim(zfish_mta_data)
```

    ## [1] 50637    24

``` r
head(zfish_mta_data)
```

    ## # A tibble: 6 × 24
    ##   ...1          Sample TSSEnrichment ReadsInTSS ReadsInPromoter ReadsInBlacklist
    ##   <chr>         <chr>          <dbl>      <dbl>           <dbl>            <dbl>
    ## 1 3hpf_1#3hpf_… 3hpf_1          6.50       1452            7341                0
    ## 2 3hpf_1#3hpf_… 3hpf_1          4.21        829            4887                0
    ## 3 3hpf_1#3hpf_… 3hpf_1          8.49       1171            5702                0
    ## 4 3hpf_1#3hpf_… 3hpf_1          5.02        684            3907                0
    ## 5 3hpf_1#3hpf_… 3hpf_1          6.57       1337            6664                0
    ## 6 3hpf_1#3hpf_… 3hpf_1          4.14        868            5580                0
    ## # ℹ 18 more variables: PromoterRatio <dbl>, PassQC <dbl>,
    ## #   NucleosomeRatio <dbl>, nMultiFrags <dbl>, nMonoFrags <dbl>, nFrags <dbl>,
    ## #   nDiFrags <dbl>, BlacklistRatio <dbl>, stage <chr>, Clusters <chr>,
    ## #   ReadsInPeaks <dbl>, FRIP <dbl>, celltype <chr>, predictedCell <chr>,
    ## #   predictedGroup <chr>, predictedScore <dbl>, DoubletScore <dbl>,
    ## #   DoubletEnrichment <dbl>

The data frame `zfish_snATAC_seq_cell_type_mrkrs` contains further
information on all peaks, which could have been used to assign cell
types.

``` r
class(zfish_snATAC_seq_cell_type_mrkrs)
```

    ## [1] "spec_tbl_df" "tbl_df"      "tbl"         "data.frame"

``` r
dim(zfish_snATAC_seq_cell_type_mrkrs)
```

    ## [1] 8684   12

``` r
head(zfish_snATAC_seq_cell_type_mrkrs)
```

    ## # A tibble: 6 × 12
    ##    ...1 group group_name seqnames    start      end strand name      idx Log2FC
    ##   <dbl> <dbl> <chr>      <chr>       <dbl>    <dbl>  <dbl> <chr>   <dbl>  <dbl>
    ## 1  3258     1 C1         chr15    43625549 43604031      2 ctsc      435   3.64
    ## 2   146     1 C1         chr1     12126648 12107682      2 mttp      146   4.45
    ## 3  9559     1 C1         chr4      4804032  4832149      1 slc13a4    87   3.94
    ## 4  1927     1 C1         chr13     5004934  5029013      1 psap       69   3.70
    ## 5 11538     1 C1         chr7     24236282 24210086      2 slc7a8a   230   3.33
    ## 6  3610     1 C1         chr16    29509127 29513434      1 ctss2.1   305   4.73
    ## # ℹ 2 more variables: FDR <dbl>, MeanDiff <dbl>

Additionally the two tables below were loaded and explored.

Supplementary table 1 consists mostly of descriptions of the isolation
of each sample.

``` r
read_excel("Data/Lin_2023_zfish_snATAC-seq/SupplT11 summary table of data deposited.xlsx") %>% head() #further experiment data Supplementary Table 1
```

    ## New names:
    ## • `` -> `...10`

    ## # A tibble: 6 × 12
    ##   project_ID sample_name(*the same number…¹ fq_name library_strategy description
    ##   <chr>      <chr>                          <chr>   <chr>            <chr>      
    ## 1 CNP0002827 snATAC_zf3.3hpf_1              3.3hpf… snATAC-seq       single nuc…
    ## 2 CNP0002827 snATAC_zf3.3hpf_1              3.3hpf… snATAC-seq       single nuc…
    ## 3 CNP0002827 snATAC_zf3.3hpf_2              3.3hpf… snATAC-seq       single nuc…
    ## 4 CNP0002827 snATAC_zf3.3hpf_2              3.3hpf… snATAC-seq       single nuc…
    ## 5 CNP0002827 snATAC_zf5.25hpf_1             5.25hp… snATAC-seq       single nuc…
    ## 6 CNP0002827 snATAC_zf5.25hpf_1             5.25hp… snATAC-seq       single nuc…
    ## # ℹ abbreviated name:
    ## #   ¹​`sample_name(*the same number represents the same library)`
    ## # ℹ 7 more variables: fq1_name <chr>, fq2_name <chr>,
    ## #   experiment_accession <chr>, run_accession <chr>, ...10 <lgl>,
    ## #   `accession link` <chr>, `NCBI SRA Accesion Number` <chr>

Supplementary table 3 contains further information on large fraction of
the peaks.

``` r
read_csv("Data/Lin_2023_zfish_snATAC-seq/SupplT13-all developmental stages intergrated cell type marker of scRNA-seq.csv") %>% head() #further experiment data Supplementary Table 3
```

    ## New names:
    ## Rows: 28907 Columns: 8
    ## ── Column specification
    ## ──────────────────────────────────────────────────────── Delimiter: "," chr
    ## (2): ...1, gene dbl (6): p_val, avg_log2FC, pct.1, pct.2, p_val_adj, cluster
    ## ℹ Use `spec()` to retrieve the full column specification for this data. ℹ
    ## Specify the column types or set `show_col_types = FALSE` to quiet this message.
    ## • `` -> `...1`

    ## # A tibble: 6 × 8
    ##   ...1    p_val avg_log2FC pct.1 pct.2 p_val_adj cluster gene   
    ##   <chr>   <dbl>      <dbl> <dbl> <dbl>     <dbl>   <dbl> <chr>  
    ## 1 h2af1al     0       4.58 0.991 0.053         0       0 h2af1al
    ## 2 cldnd       0       4.37 0.988 0.1           0       0 cldnd  
    ## 3 h1m         0       4.23 0.999 0.218         0       0 h1m    
    ## 4 ccna1       0       3.79 0.989 0.158         0       0 ccna1  
    ## 5 cldng       0       3.73 0.984 0.093         0       0 cldng  
    ## 6 acp5a       0       3.71 0.991 0.146         0       0 acp5a

# 2. Data transformation

## 2.1 Reducing the data sets to run locally

Next, I reduced the count matrix `zfish_snATAC_seq_pk_mtrx` and meta
data `zfish_mta_data` to perform QC steps locally before applying them
to the data sets containing all samples and all peaks. To achieve this,
I selected cell types Lauren Saunders deemed interesting and sampled
1000 peaks randomly.

``` r
## subset cells for 5 celltypes
unique(zfish_mta_data$celltype)
```

    ##  [1] "blastomere"               "YSL/presumptive endoderm"
    ##  [3] "EVL"                      "epiblast"                
    ##  [5] "hypoblast"                "YSL"                     
    ##  [7] "lateral plate mesoderm"   "anterior/posterior axis" 
    ##  [9] "neural keel"              "neural crest"            
    ## [11] "periderm/epidermis"       "segmental plate"         
    ## [13] "UND"                      "mesenchyme cell"         
    ## [15] "central nervous system"   "integument"              
    ## [17] "musculature system"       "primary neuron"          
    ## [19] "immature eye"             "forebrain"               
    ## [21] "erythroid lineage cell"   "digestive system"        
    ## [23] "neural stem cell"

``` r
## counting cells
zfish_mta_data %>% 
  group_by(celltype) %>% 
  tally() %>% 
  arrange(-n) %>% 
  filter(celltype %in% c("neural crest", "periderm/epidermis", "lateral plate mesoderm", "primary neuron","YSL/presumptive endoderm")) %>% #one from each germ layer plus neural crest and neurons
  mutate(sum(n)) %>% head() #testing if the cell types add up to ~ 5000
```

    ## # A tibble: 5 × 3
    ##   celltype                     n `sum(n)`
    ##   <chr>                    <int>    <int>
    ## 1 primary neuron            1598     5514
    ## 2 neural crest              1177     5514
    ## 3 periderm/epidermis         992     5514
    ## 4 lateral plate mesoderm     944     5514
    ## 5 YSL/presumptive endoderm   803     5514

``` r
## filter selected celltypes
sub_mta_data.df = zfish_mta_data %>% 
  filter(celltype %in% c("neural crest", "periderm/epidermis", "lateral plate mesoderm", "primary neuron","YSL/presumptive endoderm"))

##selecting 1000 peaks randomly
sub_pk_mtrx = zfish_snATAC_seq_pk_mtrx[sample(rownames(zfish_snATAC_seq_pk_mtrx), 1000), c(unique(sub_mta_data.df$...1))]

##testing if reduction was successful
dim(sub_pk_mtrx)
```

    ## [1] 1000 5514

``` r
sub_pk_mtrx[1:10,1:10] 
```

    ## 10 x 10 sparse Matrix of class "dgCMatrix"

    ##   [[ suppressing 10 column names '5hpf_1#5.25hpf_1_BC0021_N02', '5hpf_1#5.25hpf_1_BC0012_N01', '5hpf_1#5.25hpf_1_BC0018_N01' ... ]]

    ##                                            
    ## chr17:21082296-21082796 . . . . . . . . . .
    ## chr2:30749626-30750126  1 . . . . . 2 . . .
    ## chr18:11118728-11119228 . . 2 . . . . 1 . .
    ## chr16:31922878-31923378 1 2 1 1 . . . 2 . .
    ## chr12:48094992-48095492 . . . . . . . 1 . .
    ## chr9:14538913-14539413  . . . . . . . . . .
    ## chr7:22773984-22774484  1 2 . 2 . . . . 2 .
    ## chr2:1643095-1643595    . . . . . . . . . .
    ## chr20:32402182-32402682 . . . 2 . . . . . .
    ## chr20:34137531-34138031 . . . . . . . . . .

There are around 5000 cells after subsetting cells of 5 cell types from
the data set.

Lastly, as preparation for ATAC-Seq-QC on the reduced data set, I
created a Seurat Object containing the reduced count matrix and meta
data and saved it as an RDS file.

``` r
# creating Seurat Object 
sub_zfish_atac=CreateSeuratObject(counts = sub_pk_mtrx, assay = "atac", meta.data = sub_mta_data.df)
```

    ## Warning: Setting row names on a tibble is deprecated.

``` r
#save as RDS file
saveRDS(sub_zfish_atac, file = "Data/Lin_2023_zfish_snATAC-seq/sub_zfish_atac.rds")
```

# Side-note: How to install Signac

I also installed the signac package. It was created by the [Stuart
Lab](https://www.nature.com/articles/s41592-021-01282-5). Signac is a
package, which was specially developed for ATAC-seq data sets.

``` r
setRepositories(ind=1:3) # needed to automatically install Bioconductor dependencies
install.packages("Signac")
```

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
    ##  [1] Seurat_5.2.1       SeuratObject_5.0.2 sp_2.2-0           ggsci_3.2.0       
    ##  [5] readxl_1.4.5       broom_1.0.8        binom_1.1-1.1      lubridate_1.9.4   
    ##  [9] forcats_1.0.0      stringr_1.5.1      dplyr_1.1.4        purrr_1.0.4       
    ## [13] readr_2.1.5        tidyr_1.3.1        tibble_3.2.1       ggplot2_3.5.1     
    ## [17] tidyverse_2.0.0   
    ## 
    ## loaded via a namespace (and not attached):
    ##   [1] RColorBrewer_1.1-3     rstudioapi_0.17.1      jsonlite_1.8.9        
    ##   [4] magrittr_2.0.3         spatstat.utils_3.1-2   farver_2.1.2          
    ##   [7] rmarkdown_2.29         vctrs_0.6.5            ROCR_1.0-11           
    ##  [10] spatstat.explore_3.3-4 htmltools_0.5.8.1      cellranger_1.1.0      
    ##  [13] sctransform_0.4.1      parallelly_1.43.0      KernSmooth_2.23-26    
    ##  [16] htmlwidgets_1.6.4      ica_1.0-3              plyr_1.8.9            
    ##  [19] plotly_4.10.4          zoo_1.8-12             igraph_2.1.4          
    ##  [22] mime_0.12              lifecycle_1.0.4        pkgconfig_2.0.3       
    ##  [25] Matrix_1.7-1           R6_2.6.1               fastmap_1.2.0         
    ##  [28] fitdistrplus_1.2-2     future_1.40.0          shiny_1.10.0          
    ##  [31] digest_0.6.37          colorspace_2.1-1       patchwork_1.3.0       
    ##  [34] tensor_1.5             RSpectra_0.16-2        irlba_2.3.5.1         
    ##  [37] progressr_0.15.1       spatstat.sparse_3.1-0  timechange_0.3.0      
    ##  [40] httr_1.4.7             polyclip_1.10-7        abind_1.4-8           
    ##  [43] compiler_4.4.2         bit64_4.6.0-1          withr_3.0.2           
    ##  [46] backports_1.5.0        fastDummies_1.7.5      MASS_7.3-64           
    ##  [49] tools_4.4.2            lmtest_0.9-40          httpuv_1.6.15         
    ##  [52] future.apply_1.11.3    goftest_1.2-3          glue_1.8.0            
    ##  [55] nlme_3.1-167           promises_1.3.2         grid_4.4.2            
    ##  [58] Rtsne_0.17             cluster_2.1.8          reshape2_1.4.4        
    ##  [61] generics_0.1.3         gtable_0.3.6           spatstat.data_3.1-6   
    ##  [64] tzdb_0.4.0             data.table_1.16.4      hms_1.1.3             
    ##  [67] utf8_1.2.4             spatstat.geom_3.3-5    RcppAnnoy_0.0.22      
    ##  [70] ggrepel_0.9.6          RANN_2.6.2             pillar_1.10.2         
    ##  [73] spam_2.11-1            RcppHNSW_0.6.0         vroom_1.6.5           
    ##  [76] later_1.4.1            splines_4.4.2          lattice_0.22-6        
    ##  [79] survival_3.8-3         bit_4.6.0              deldir_2.0-4          
    ##  [82] tidyselect_1.2.1       miniUI_0.1.1.1         pbapply_1.7-2         
    ##  [85] knitr_1.50             gridExtra_2.3          scattermore_1.2       
    ##  [88] xfun_0.52              matrixStats_1.5.0      stringi_1.8.4         
    ##  [91] lazyeval_0.2.2         yaml_2.3.10            evaluate_1.0.3        
    ##  [94] codetools_0.2-20       cli_3.6.1              uwot_0.2.2            
    ##  [97] xtable_1.8-4           reticulate_1.40.0      munsell_0.5.1         
    ## [100] Rcpp_1.0.14            globals_0.17.0         spatstat.random_3.3-2 
    ## [103] png_0.1-8              spatstat.univar_3.1-1  parallel_4.4.2        
    ## [106] dotCall64_1.2          listenv_0.9.1          viridisLite_0.4.2     
    ## [109] scales_1.3.0           ggridges_0.5.6         crayon_1.5.3          
    ## [112] rlang_1.1.4            cowplot_1.1.3
