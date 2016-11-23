# DevSeq Project

## Reproducible Scripts for the Publication:

Schuster C, Gabel A, Drost H-G, Grosse I, Meyerowitz E. DevSeq: 


### Reference Genome and Proteome Retrieval

The CDS and proteome files of _Arabidopsis thaliana_ have been retrieved from ENSEMBL
via the [biomartr](https://github.com/HajkD/biomartr) package. First, users need to
install the [biomartr](https://github.com/HajkD/biomartr#installation) package.

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

For dNdS computations first users need to install the [orthologr](https://github.com/HajkD/orthologr#installation-guide) package. 

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