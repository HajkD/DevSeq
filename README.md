# DevSeq Project

## Reproducible Scripts for the Publication:

Schuster C, Gabel A, Drost H-G, Grosse I, Meyerowitz E. DevSeq: 


### Reference Genome and Proteome Retrieval

The CDS and proteome files of _Arabidopsis thaliana_ have been retrieved from ENSEMBL
via the [biomartr](https://github.com/HajkD/biomartr) package. First, users need to
[install the biomartr](https://github.com/HajkD/biomartr#installation) package.

```r
# install.packages("biomartr")

# download CDS for Arabidopsis thaliana
biomartr::getCDS(db = "ensemblgenomes", organism = "Arabidopsis thaliana", path = "data/CDS")

# download CDS for Tarenaya hassleriana
biomartr::getCDS(db = "refseq", organism = "Tarenaya hassleriana", path = "data/CDS/subject_species")
```

The CDS and Proteome files for `` have been downloaded from [Phytozome V11](https://phytozome.jgi.doe.gov/pz/portal.html) on 17 Nov 2016.

### Perform Orthology Inference and Generate dN/dS tables

### Installing the orthologr package

For dNdS computations first users need to [install the orthologr](https://github.com/HajkD/orthologr#installation-guide) package. 
As a next step, please also [install BLAST+](https://github.com/HajkD/orthologr/blob/master/vignettes/Install.Rmd#install-blast) into the `/usr/local/bin` directory on your computer. 

```r
# compute dN/dS table of A. thaliana vs. all other species
orthologr::map.generator(
               query_file      = "data/CDS/Arabidopsis_thaliana.TAIR10.cds.all.fa.gz",
               subjects.folder = "data/CDS/subject_species",
               eval            = "1E-5", 
               ortho_detection = "RBH",
               aa_aln_type      = "pairwise",
               aa_aln_tool      = "NW", 
               codon_aln_tool   = "pal2nal", 
               dnds_est.method  = "Comeron",
               output.folder    = "data/dNdS_maps",
               comp_cores       = 12
               )

```

Import all dNdS maps:

```r
# Import all dNdS maps and store each pairwise comparison as list element
# of map.list
map.list <- lapply(list.files("data/dNdS_maps/"), function(map) {
    
    readr::read_delim(
        file.path("data/dNdS_maps/",map),
        col_names = TRUE,
        delim = ";",
        col_types = readr::cols("query_id" = readr::col_character(),
                                "subject_id"= readr::col_character(),
                                "dN" = readr::col_double(),
                                "dS"  = readr::col_double(),
                                "dNdS" = readr::col_double())
    )
})

# rename list elements 
names(map.list) <- paste0("Ath_vs_", c("Alyr", "Crub", "Esals", "Mtrunc", "Thassl"))

# look at import
map.list
```

```
$Ath_vs_Alyr
 A tibble: 21,654 × 5
        query_id subject_id       dN     dS    dNdS
           <chr>      <chr>    <dbl>  <dbl>   <dbl>
1  AT1G01010.1.1     333554 0.106400 0.2537 0.41950
2  AT1G01020.1.1     470181 0.040230 0.1037 0.38790
3  AT1G01030.1.1     470180 0.014990 0.1265 0.11850
4  AT1G01040.1.1     333551 0.013470 0.1165 0.11560
5  AT1G01050.2.1     909874 0.000000 0.1750 0.00000
6  AT1G01060.8.1     470177 0.044950 0.1133 0.39670
7  AT1G01070.1.1     918864 0.018300 0.1059 0.17280
8  AT1G01080.1.1     909871 0.033980 0.1056 0.32170
9  AT1G01090.1.1     470171 0.009104 0.2181 0.04174
10 AT1G01110.2.1     333544 0.032480 0.1220 0.26620
 ... with 21,644 more rows

$Ath_vs_Crub
 A tibble: 21,246 × 5
        query_id      subject_id      dN     dS    dNdS
           <chr>           <chr>   <dbl>  <dbl>   <dbl>
1  AT1G01010.1.1 Carubv10009049m 0.11180 0.2760 0.40520
2  AT1G01020.1.1 Carubv10011984m 0.08264 0.1810 0.45650
3  AT1G01030.1.1 Carubv10009540m 0.02824 0.2248 0.12560
4  AT1G01040.1.1 Carubv10008073m 0.02317 0.2349 0.09863
5  AT1G01050.2.1 Carubv10010288m 0.01063 0.2386 0.04453
6  AT1G01060.4.1 Carubv10008551m 0.05169 0.2002 0.25820
7  AT1G01070.1.1 Carubv10009303m 0.08801 0.3165 0.27800
8  AT1G01080.1.1 Carubv10009913m 0.07032 0.2134 0.32940
9  AT1G01090.1.1 Carubv10009201m 0.01581 0.3457 0.04573
10 AT1G01110.2.1 Carubv10008796m 0.05539 0.2266 0.24440
 ... with 21,236 more rows

$Ath_vs_Esals
 A tibble: 19,982 × 5
        query_id      subject_id       dN     dS    dNdS
           <chr>           <chr>    <dbl>  <dbl>   <dbl>
1  AT1G01010.1.1 Thhalv10007665m 0.181900 0.3618 0.50280
2  AT1G01020.1.1 Thhalv10008618m 0.101300 0.2806 0.36100
3  AT1G01030.1.1 Thhalv10008068m 0.025890 0.2531 0.10230
4  AT1G01040.1.1 Thhalv10006531m 0.033420 0.3129 0.10680
5  AT1G01050.2.1 Thhalv10008767m 0.007183 0.3848 0.01867
6  AT1G01060.8.1 Thhalv10007029m 0.096500 0.2749 0.35100
7  AT1G01070.1.1 Thhalv10009501m 0.085800 0.3232 0.26550
8  AT1G01080.1.1 Thhalv10008355m 0.063880 0.2909 0.21960
9  AT1G01090.1.1 Thhalv10007695m 0.022560 0.4576 0.04930
10 AT1G01110.2.1 Thhalv10007299m 0.060020 0.3247 0.18480
 ... with 19,972 more rows

$Ath_vs_Mtrunc
 A tibble: 12,863 × 5
        query_id      subject_id      dN    dS    dNdS
           <chr>           <chr>   <dbl> <dbl>   <dbl>
1  AT1G01040.1.1 Medtr7g118350.1 0.14510 1.675 0.08662
2  AT1G01050.2.1 Medtr8g024050.1 0.08582 1.318 0.06513
3  AT1G01060.4.1 Medtr7g118330.2 0.33950 1.605 0.21150
4  AT1G01080.1.1 Medtr7g118230.1 0.49390    NA      NA
5  AT1G01090.1.1 Medtr8g024310.1 0.12610 1.752 0.07198
6  AT1G01110.2.1 Medtr8g024540.1 0.26390    NA      NA
7  AT1G01120.1.1 Medtr7g118170.1 0.17790 2.074 0.08577
8  AT1G01140.1.1 Medtr8g024600.1 0.15170 1.523 0.09955
9  AT1G01170.2.1 Medtr7g091770.2 0.21350 1.760 0.12130
10 AT1G01180.1.1 Medtr8g024680.1 0.25220 1.820 0.13860
 ... with 12,853 more rows

$Ath_vs_Thassl
 A tibble: 15,791 × 5
        query_id                                  subject_id      dN     dS    dNdS
           <chr>                                       <chr>   <dbl>  <dbl>   <dbl>
1  AT1G01020.1.1 lcl|NW_010965696.1_cds_XP_010542814.1_17939 0.23950 0.9072 0.26400
2  AT1G01030.1.1 lcl|NW_010967707.1_cds_XP_010556651.1_30297 0.14160 1.1620 0.12190
3  AT1G01040.1.1 lcl|NW_010967707.1_cds_XP_010556639.1_30293 0.08732 0.7440 0.11740
4  AT1G01040.2.1 lcl|NW_010967707.1_cds_XP_010556638.1_30296 0.08784 0.7449 0.11790
5  AT1G01050.2.1 lcl|NW_010965696.1_cds_XP_010542812.1_17937 0.03153 0.7757 0.04065
6  AT1G01060.6.1 lcl|NW_010967707.1_cds_XP_010556647.1_30291 0.23920 0.8518 0.28080
7  AT1G01060.8.1 lcl|NW_010967707.1_cds_XP_010556646.1_30289 0.20380 0.8670 0.23500
8  AT1G01070.1.1 lcl|NW_010969548.1_cds_XP_010522883.1_36546 0.22380 1.3040 0.17160
9  AT1G01080.1.1 lcl|NW_010967707.1_cds_XP_010556629.1_30271 0.15890 0.8563 0.18550
10 AT1G01090.1.1 lcl|NW_010967707.1_cds_XP_010556630.1_30270 0.07587 1.0110 0.07507
 ... with 15,781 more rows
```

Individual pairwise dNdS files can then be selected by typing:

```r
# example: how to select individual dNdS maps
# here: dNdS map for A thaliana vs E salsugineum
map.list$Ath_vs_Esals
```

```
 A tibble: 19,982 × 5
        query_id      subject_id       dN     dS    dNdS
           <chr>           <chr>    <dbl>  <dbl>   <dbl>
1  AT1G01010.1.1 Thhalv10007665m 0.181900 0.3618 0.50280
2  AT1G01020.1.1 Thhalv10008618m 0.101300 0.2806 0.36100
3  AT1G01030.1.1 Thhalv10008068m 0.025890 0.2531 0.10230
4  AT1G01040.1.1 Thhalv10006531m 0.033420 0.3129 0.10680
5  AT1G01050.2.1 Thhalv10008767m 0.007183 0.3848 0.01867
6  AT1G01060.8.1 Thhalv10007029m 0.096500 0.2749 0.35100
7  AT1G01070.1.1 Thhalv10009501m 0.085800 0.3232 0.26550
8  AT1G01080.1.1 Thhalv10008355m 0.063880 0.2909 0.21960
9  AT1G01090.1.1 Thhalv10007695m 0.022560 0.4576 0.04930
10 AT1G01110.2.1 Thhalv10007299m 0.060020 0.3247 0.18480
 ... with 19,972 more rows
```

Analogous all other individual dNdS maps can be selected by typing `map.list$Ath_vs_Alyr`, `map.list$Ath_vs_Crub`, `map.list$Ath_vs_Thassl`, etc.

### Detection of all `A. thaliana` genes that have intersecting orthologs with all other species

```r
# first rename colnames of individual dNdS maps
for (i in seq_along(map.list)) {
   colnames(map.list[[i]])[2:5] <- paste0(names(map.list)[i], c("_subject_id","_dN", "_dS", "_dNdS"))
}

# combine all geneids into one file
all.maps <- dplyr::bind_rows(lapply(map.list, function(x) dplyr::select(x, query_id)))

# detect genes that have orthologs in all other species
length(all.maps$query_id[which(table(all.maps$query_id) == length(map.list))])
```

```
9280
```

Thus, `9280 A. thaliana` genes have orthologs in all other species.
Now, we combine all dN, dS, and dNdS information for these `9280` genes.

```r
# store all intersecting orthologs in tibble
all.orthologs <- tibble::as_tibble(all.maps$query_id[which(table(all.maps$query_id) == length(map.list))])
colnames(all.orthologs) <- "query_id"

lapply(map.list, function(x) dplyr::inner_join(all.orthologs, x, by = "query_id"))
```



The [DevSeqR package](https://github.com/HajkD/DevSeqR) allows users to reproduce all analyses and to perform
additional exploratory data analysis using the DevSeq dataset.

## Install `DevSeqR` package

```r
# install.packages("devtools")

# install the current version of DevSeqR on your system
library(devtools)
install_github("HajkD/DevSeqR", build_vignettes = TRUE, dependencies = TRUE)
```

## Getting Started

```r
# load DevSeqR package
library(DevSeqR)
# import experimental design information for the DevSeq dataset
DevSeqSample <- load.sample.info()

DevSeqSample
```

```
   Sample_ID DevSeq_ID    Species Tissue                                        Description
       <dbl>     <chr>      <chr>  <chr>                                              <chr>
1          1     10-12 A.thaliana   root                       root, root tip (top 0.5-1mm)
2          2     13-15 A.thaliana   root root, maturation zone (2-3mm piece above root tip)
3          3   106-108 A.thaliana   root    root, whole root (3mm piece including root tip)
4          4       1-3 A.thaliana   root                                   root, whole root
5          5       7-9 A.thaliana   root                                   root, whole root
6          6       4-6 A.thaliana   root                                   root, whole root
7          7     16-18 A.thaliana   stem                                          hypocotyl
8          8     22-24 A.thaliana   stem                        3rd internode (0.5cm piece)
9          9     19-21 A.thaliana   stem                          2nd internode (1cm piece)
10        10     25-27 A.thaliana   stem                          1st internode (1cm piece)
 ... with 130 more rows, and 5 more variables: Age <chr>, Substrate <chr>,
   Photoperiod <chr>, Replicates <dbl>, Sample_size <chr>
```
