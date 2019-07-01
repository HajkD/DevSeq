# DevSeq Project

## Reproducible Scripts for the Publication:

> Schuster C, Gabel A, Drost H-G, Leyser O, Grosse I, Meyerowitz E. Comparative developmental transcriptome atlases across nine plant body plan organs.


* [1. Resource Retrieval](#Resource Retrieval)
    - [1.1 Reference Genome and CDS Retrieval](#reference-genome-and-cds-retrieval)
    - [1.2 Software Retrieval](#reference-genome-and-cds-retrieval)
* [2. Orthology Inference and dNdS Estimation](#orthology-inference-and-dnds-estimation)
    - [2.1 BLAST best reciprocal hit](#blast-best-reciprocal-hit)
    - [2.2 Orthogroup Inference with Orthofinder2](#orthogroup-inference-with-orthofinder2)

## Resource Retrieval

### Reference Genome and CDS Retrieval

The reference genomes used for all analyses were obtained from diverse databases:



```r
# install.packages("biomartr")

# download CDS for Arabidopsis thaliana
biomartr::getCDS(db = "ensemblgenomes",
                 organism = "Arabidopsis thaliana",
                 path = "data/CDS")
                 
# download Proteome for Arabidopsis thaliana
biomartr::getProteome(db = "ensemblgenomes",
                      organism = "Arabidopsis thaliana",
                      path = "data/Proteome")
                 
# download CDS for Tarenaya hassleriana
biomartr::getCDS(db = "refseq",
                 organism = "Tarenaya hassleriana",
                 path = "data/CDS/subject_species")
# download Proteome for Tarenaya hassleriana
biomartr::getProteome(db = "refseq",
                      organism = "Tarenaya hassleriana",
                      path = "data/Proteome")
```

The CDS and Proteome files for `A. lyrata`, `B. distachyon`, `C. rubella`, `E. salsugineum`, `M. truncatula` have been downloaded from [Phytozome V11](https://phytozome.jgi.doe.gov/pz/portal.html) on 17 Nov 2016. The CDS sequences for `Picea abies` were downloaded from ftp://plantgenie.org/Data/ConGenIE/Picea_abies/v1.0/FASTA/GenePrediction/Pabies1.0-all-cds.fna.gz on 26 Mar 2018.

### Software Retrieval

#### Installing the orthologr package

For dNdS computations first users need to [install the orthologr](https://github.com/HajkD/orthologr#installation-guide) package. 
As a next step, please also [install BLAST+](https://github.com/HajkD/orthologr/blob/master/vignettes/Install.Rmd#install-blast) into the `/usr/local/bin` directory on your computer. 
Please be aware that the following dN/dS computations might run several hours with 12 - 20 cores or even up to 1 - 1.5 days with 4 - 8 cores.

The coding sequence files of _Arabidopsis thaliana_ have been retrieved from ENSEMBL
via the [biomartr](https://github.com/HajkD/biomartr) package. First, users need to
[install the biomartr](https://github.com/HajkD/biomartr#installation) package.


### Orthology Inference and dNdS Estimation


### BLAST best reciprocal hit

## Generate 1:1 orthologs tables

### Generate 1:1 orthologs tables for A. thaliana

```r
# compute dN/dS table of A. thaliana vs. all other species
orthologr::map_generator(
               query_file      = "data/DevSeq_CDS/Query_files/Athaliana.fa",
               subjects.folder = "data/CDS/Athaliana_subject_files",
               eval            = "1E-5", # e value threshold for ortholog detection
               ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
               aa_aln_type      = "pairwise",
               aa_aln_tool      = "NW", # use Needleman-Wunsch Algorithm for global codon alignment
               codon_aln_tool   = "pal2nal", 
               dnds_est.method  = "Comeron", # use robust dN/dS estimation (Comeron's method)
               output.folder    = "data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/",
               comp_cores       = 12
               )

```

Import all dNdS maps:

```r
# Import all A thaliana dNdS maps and store each pairwise comparison as list element
Ath_pairwise_dNdS_maps <- lapply(list.files("data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/"), function(map) {
    
    readr::read_delim(
        file.path("data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/", map),
        col_names = TRUE,
        delim = ";"
    )
})

# rename list elements 
names(Ath_pairwise_dNdS_maps) <- paste0("Athaliana_vs_", c("Alyrata", "Bdistachyon","Crubella", "Esalsugineum", "Mtruncatula", "Thassleriana"))

# look at import
Ath_pairwise_dNdS_maps
```

```
$Athaliana_vs_Alyrata
# A tibble: 23,844 x 15
   query_id subject_id      dN    dS   dNdS perc_identity alig_length
   <chr>    <chr>        <dbl> <dbl>  <dbl>         <dbl>       <dbl>
 1 AT1G010… AL1G11530… 0.106   0.254 0.420           74.0         469
 2 AT1G010… AL1G11510… 0.0402  0.104 0.388           91.1         246
 3 AT1G010… AL1G11500… 0.0565  0.210 0.270           96.3         298
 4 AT1G010… AL1G11480… 0.0134  0.117 0.115           96.2        1915
 5 AT1G010… AL1G11470… 0       0.175 0              100           213
 6 AT1G010… AL1G11460… 0.0444  0.117 0.379           88.5         654
 7 AT1G010… AL1G11450… 0.0183  0.106 0.173           95.1         366
 8 AT1G010… AL1G11440… 0.0340  0.110 0.310           82.9         327
 9 AT1G010… AL1G11430… 0.00910 0.218 0.0417          96.8         434
10 AT1G011… AL1G11410… 0.0325  0.122 0.266           93.6         528
# … with 23,834 more rows, and 8 more variables: mismatches <dbl>,
#   gap_openings <dbl>, q_start <dbl>, q_end <dbl>, s_start <dbl>,
#   s_end <dbl>, evalue <dbl>, bit_score <dbl>

$Athaliana_vs_Bdistachyon
# A tibble: 12,148 x 15
   query_id subject_id    dN    dS    dNdS perc_identity alig_length mismatches
   <chr>    <chr>      <dbl> <dbl>   <dbl>         <dbl>       <dbl>      <dbl>
 1 AT1G010… Bradi3g16… 0.549 NA    NA               38.0         397        198
 2 AT1G010… Bradi3g13… 0.457 NA    NA               53.4         189         87
 3 AT1G010… Bradi5g01… 0.156  1.98  0.0790          86.7         354         44
 4 AT1G011… Bradi1g68… 0.297 NA    NA               60           550        186
 5 AT1G011… Bradi1g76… 0.198  1.92  0.103           72.6         445        115
 6 AT1G011… Bradi1g50… 0.197  1.20  0.165           72.6          73         20
 7 AT1G011… Bradi4g35… 0.474 NA    NA               56.5         248        101
 8 AT1G012… Bradi1g35… 0.235  1.52  0.155           65.4         214         63
 9 AT1G012… Bradi1g77… 0.255  2.77  0.0919          63.2        1060        365
10 AT1G012… Bradi4g35… 0.263 NA    NA               66.3         249         77
# … with 12,138 more rows, and 7 more variables: gap_openings <dbl>,
#   q_start <dbl>, q_end <dbl>, s_start <dbl>, s_end <dbl>, evalue <dbl>,
#   bit_score <dbl>

$Athaliana_vs_Crubella
# A tibble: 22,721 x 15
   query_id subject_id     dN    dS   dNdS perc_identity alig_length mismatches
   <chr>    <chr>       <dbl> <dbl>  <dbl>         <dbl>       <dbl>      <dbl>
 1 AT1G010… Carubv100… 0.112  0.276 0.405           68.5         482         88
 2 AT1G010… Carubv100… 0.0826 0.181 0.456           84.0         244         38
 3 AT1G010… Carubv100… 0.0282 0.225 0.126           93.9         361         19
 4 AT1G010… Carubv100… 0.0232 0.236 0.0980          94.5        1915         94
 5 AT1G010… Carubv100… 0.0106 0.239 0.0445          97.7         213          5
 6 AT1G010… Carubv100… 0.0540 0.202 0.267           87.9         651         68
 7 AT1G010… Carubv100… 0.0880 0.316 0.278           82.2         371         60
 8 AT1G010… Carubv100… 0.0703 0.213 0.329           85.1         295         35
 9 AT1G010… Carubv100… 0.0158 0.346 0.0457          95.2         433         15
10 AT1G011… Carubv100… 0.0554 0.227 0.244           84.9         543         63
# … with 22,711 more rows, and 7 more variables: gap_openings <dbl>,
#   q_start <dbl>, q_end <dbl>, s_start <dbl>, s_end <dbl>, evalue <dbl>,
#   bit_score <dbl>

$Athaliana_vs_Esalsugineum
# A tibble: 22,831 x 15
   query_id subject_id      dN    dS   dNdS perc_identity alig_length
   <chr>    <chr>        <dbl> <dbl>  <dbl>         <dbl>       <dbl>
 1 AT1G010… Thhalv100… 0.182   0.362 0.503           64.9         445
 2 AT1G010… Thhalv100… 0.101   0.281 0.361           81.1         244
 3 AT1G010… Thhalv100… 0.0259  0.253 0.102           91.9         360
 4 AT1G010… Thhalv100… 0.0334  0.314 0.106           92.9        1923
 5 AT1G010… NewID_000… 0.00718 0.385 0.0187          98.6         213
 6 AT1G010… NewID_000… 0.0965  0.275 0.351           77.6         664
 7 AT1G010… NewID_000… 0.0855  0.323 0.265           84.2         367
 8 AT1G010… NewID_000… 0.0639  0.291 0.220           86.9         297
 9 AT1G010… NewID_000… 0.0226  0.458 0.0493          94.2         432
10 AT1G011… Thhalv100… 0.0600  0.325 0.185           84.3         547
# … with 22,821 more rows, and 8 more variables: mismatches <dbl>,
#   gap_openings <dbl>, q_start <dbl>, q_end <dbl>, s_start <dbl>,
#   s_end <dbl>, evalue <dbl>, bit_score <dbl>

$Athaliana_vs_Mtruncatula
# A tibble: 14,287 x 15
   query_id subject_id     dN    dS    dNdS perc_identity alig_length
   <chr>    <chr>       <dbl> <dbl>   <dbl>         <dbl>       <dbl>
 1 AT1G010… Medtr7g08… 0.378   1.93  0.196           51.2         240
 2 AT1G010… Medtr7g11… 0.145   1.67  0.0868          74.2        1909
 3 AT1G010… Medtr8g02… 0.0858  1.32  0.0651          88.3         213
 4 AT1G010… Medtr7g11… 0.340   1.60  0.212           43.8         493
 5 AT1G010… Medtr7g11… 0.494  NA    NA               44.7         197
 6 AT1G010… Medtr8g02… 0.126   1.75  0.0720          82           400
 7 AT1G011… Medtr8g02… 0.264  NA    NA               55.3         580
 8 AT1G011… Medtr7g11… 0.178   2.07  0.0858          73.0         523
 9 AT1G011… Medtr8g02… 0.143   2.39  0.0597          74.6         437
10 AT1G011… Medtr7g09… 0.214   1.76  0.121           68.4          76
# … with 14,277 more rows, and 8 more variables: mismatches <dbl>,
#   gap_openings <dbl>, q_start <dbl>, q_end <dbl>, s_start <dbl>,
#   s_end <dbl>, evalue <dbl>, bit_score <dbl>

$Athaliana_vs_Thassleriana
# A tibble: 17,249 x 15
   query_id subject_id     dN    dS   dNdS perc_identity alig_length mismatches
   <chr>    <chr>       <dbl> <dbl>  <dbl>         <dbl>       <dbl>      <dbl>
 1 AT1G010… XM_010544… 0.240  0.907 0.264           67.0         224         71
 2 AT1G010… XM_010558… 0.142  1.16  0.122           60.9         391         79
 3 AT1G010… XM_010558… 0.0873 0.744 0.117           81.6        1954        281
 4 AT1G010… XM_010544… 0.0315 0.776 0.0406          93.5         215         12
 5 AT1G010… XM_010558… 0.204  0.867 0.235           57.1         750        196
 6 AT1G010… XM_010558… 0.239  0.852 0.281           52.4         750        220
 7 AT1G010… XM_010524… 0.224  1.30  0.172           65.2         368        118
 8 AT1G010… XM_010558… 0.159  0.856 0.186           70.5         295         75
 9 AT1G010… XM_010524… 0.0673 0.831 0.0810          85.9         432         55
10 AT1G011… XM_010558… 0.153  1.38  0.111           68.8         382         82
# … with 17,239 more rows, and 7 more variables: gap_openings <dbl>,
#   q_start <dbl>, q_end <dbl>, s_start <dbl>, s_end <dbl>, evalue <dbl>,
#   bit_score <dbl>
```

Individual pairwise dNdS files can then be selected by typing:

```r
# example: how to select individual dNdS maps
# here: dNdS map for A thaliana vs E salsugineum
Ath_pairwise_dNdS_maps$Athaliana_vs_Alyrata
```

```
 # A tibble: 23,844 x 15
   query_id subject_id      dN    dS   dNdS perc_identity alig_length
   <chr>    <chr>        <dbl> <dbl>  <dbl>         <dbl>       <dbl>
 1 AT1G010… AL1G11530… 0.106   0.254 0.420           74.0         469
 2 AT1G010… AL1G11510… 0.0402  0.104 0.388           91.1         246
 3 AT1G010… AL1G11500… 0.0565  0.210 0.270           96.3         298
 4 AT1G010… AL1G11480… 0.0134  0.117 0.115           96.2        1915
 5 AT1G010… AL1G11470… 0       0.175 0              100           213
 6 AT1G010… AL1G11460… 0.0444  0.117 0.379           88.5         654
 7 AT1G010… AL1G11450… 0.0183  0.106 0.173           95.1         366
 8 AT1G010… AL1G11440… 0.0340  0.110 0.310           82.9         327
 9 AT1G010… AL1G11430… 0.00910 0.218 0.0417          96.8         434
10 AT1G011… AL1G11410… 0.0325  0.122 0.266           93.6         528
# … with 23,834 more rows, and 8 more variables: mismatches <dbl>,
#   gap_openings <dbl>, q_start <dbl>, q_end <dbl>, s_start <dbl>,
#   s_end <dbl>, evalue <dbl>, bit_score <dbl>
```

Analogous all other individual dNdS maps can be selected by typing `map.list$Ath_vs_Alyr`,`map.list$Ath_vs_Bdist`, `map.list$Ath_vs_Crub`, `map.list$Ath_vs_Thassl`, etc.

### Detection of all `A. thaliana` genes that have intersecting orthologs with all other species

```r
# first rename colnames of individual dNdS maps
for (i in seq_along(map.list)) {
   colnames(map.list[[i]])[2:5] <- paste0(names(map.list)[i], c("_subject_id","_dN", "_dS", "_dNdS"))
}

# combine all geneids into one file
all.maps <- dplyr::bind_rows(lapply(map.list, function(x) tibble::as_tibble(unique(x$query_id))))
colnames(all.maps) <- "query_id"

# detect genes that have orthologs in all other species
length(names(table(all.maps$query_id))[which(table(all.maps$query_id) == length(map.list))])
```

```
7276
```

Thus, `7276 A. thaliana` genes have orthologs in all other species.

Look at how many orthologs are shared between 1,2,3,.. species:

```r
table(table(all.maps$query_id))
```

```
   1    2    3    4    5    6 
4201 3235 4152 4035 3934 7276 
```

Now, we combine all dN, dS, and dNdS information for these `7276` genes.

```r
# store all intersecting orthologs in tibble
all.orthologs <- tibble::as_tibble(names(table(all.maps$query_id))[which(table(all.maps$query_id) == length(map.list))])
colnames(all.orthologs) <- "query_id"

# generate orthologs tables
orthologs <- lapply(map.list, function(x) dplyr::inner_join(all.orthologs, x, by = "query_id"))

# join orthologs tables to a final cross-species orthologs dNdS file
final.orthologs <- orthologs$Ath_vs_Alyr
for (i in (seq_along(map.list) - 1)) {
    final.orthologs <- dplyr::inner_join(final.orthologs, orthologs[[i + 1]], by = "query_id")
}

# filter table
final.orthologs <- dplyr::select(
    final.orthologs,
    -Ath_vs_Alyr_subject_id.x,
    -Ath_vs_Alyr_dN.x,
    -Ath_vs_Alyr_dS.x,
    -Ath_vs_Alyr_dNdS.x
)
colnames(final.orthologs)[2:5] <-
    c("Ath_vs_Alyr_subject_id",
      "Ath_vs_Alyr_dN",
      "Ath_vs_Alyr_dS",
      "Ath_vs_Alyr_dNdS")

# looking at the final table
final.orthologs
```

```
 A tibble: 7,276 × 25
        query_id Ath_vs_Alyr_subject_id Ath_vs_Alyr_dN Ath_vs_Alyr_dS Ath_vs_Alyr_dNdS
           <chr>                  <chr>          <dbl>          <dbl>            <dbl>
1  AT1G01080.1.1                 909871       0.033980        0.10560          0.32170
2  AT1G01090.1.1                 470171       0.009104        0.21810          0.04174
3  AT1G01120.1.1                 918858       0.003072        0.13260          0.02317
4  AT1G01170.2.1                 311317       0.000000        0.30640          0.00000
5  AT1G01180.1.1                 909860       0.038320        0.15400          0.24880
6  AT1G01200.1.1                 470156       0.019050        0.16750          0.11370
7  AT1G01225.1.1                 470154       0.013320        0.11220          0.11870
8  AT1G01230.1.1                 470153       0.002573        0.06798          0.03785
9  AT1G01260.3.1                 470147       0.016700        0.20880          0.07996
10 AT1G01280.1.1                 470146       0.026260        0.15670          0.16760
 ... with 7,266 more rows, and 20 more variables: Ath_vs_Bdist_subject_id <chr>,
   Ath_vs_Bdist_dN <dbl>, Ath_vs_Bdist_dS <dbl>, Ath_vs_Bdist_dNdS <dbl>,
   Ath_vs_Crub_subject_id <chr>, Ath_vs_Crub_dN <dbl>, Ath_vs_Crub_dS <dbl>,
   Ath_vs_Crub_dNdS <dbl>, Ath_vs_Esals_subject_id <chr>, Ath_vs_Esals_dN <dbl>,
   Ath_vs_Esals_dS <dbl>, Ath_vs_Esals_dNdS <dbl>, Ath_vs_Mtrunc_subject_id <chr>,
   Ath_vs_Mtrunc_dN <dbl>, Ath_vs_Mtrunc_dS <dbl>, Ath_vs_Mtrunc_dNdS <dbl>,
   Ath_vs_Thassl_subject_id <chr>, Ath_vs_Thassl_dN <dbl>, Ath_vs_Thassl_dS <dbl>,
   Ath_vs_Thassl_dNdS <dbl>
```

Store final table in `;` separated file:

```r
# create new folder "ortho_table"
dir.create("data/ortho_table")
# store final orthologs file in ortho_table folder
readr::write_delim(final.orthologs, "data/ortho_table/DevSeq_all_species_intersect_orthologs.csv", delim = ";")
```


## Generate large table with all genes

```r
# store all intersecting orthologs in tibble
all.orthologs <- tibble::as_tibble(names(table(all.maps$query_id))[which(table(all.maps$query_id) == length(map.list))])
colnames(all.orthologs) <- "query_id"

# generate full orthologs tables
orthologs_full <- lapply(map.list, function(x) dplyr::full_join(all.orthologs, x, by = "query_id"))

# join orthologs tables to a final cross-species orthologs dNdS file including all genes
final.orthologs_full <- orthologs_full$Ath_vs_Alyr
for (i in (seq_along(map.list) - 1)) {
    final.orthologs_full <- dplyr::full_join(final.orthologs_full, orthologs_full[[i + 1]], by = "query_id")
}

# filter table
final.orthologs_full <- dplyr::select(
    final.orthologs_full,
    -Ath_vs_Alyr_subject_id.x,
    -Ath_vs_Alyr_dN.x,
    -Ath_vs_Alyr_dS.x,
    -Ath_vs_Alyr_dNdS.x
)
colnames(final.orthologs_full)[2:5] <-
    c("Ath_vs_Alyr_subject_id",
      "Ath_vs_Alyr_dN",
      "Ath_vs_Alyr_dS",
      "Ath_vs_Alyr_dNdS")

# looking at the final table
final.orthologs_full
```

Store final table in `;` separated file:

```r
# create new folder "ortho_table" 
dir.create("data/ortho_table")
# store final orthologs file in ortho_table folder
readr::write_delim(final.orthologs_full, "data/ortho_table/DevSeq_all_species_fulljoin_orthologs.csv", delim = ";")
```

### Generate 1:1 orthologs tables for B. distachyon

```r
# compute dN/dS table of B. distachyon vs. M. truncatula and T. hassleriana
orthologr::map.generator(
               query_file      = "data/CDS/Bdistachyon_314_v3.1.cds.fa.gz",
               subjects.folder = "data/CDS/subject_species_Bdistachyon",
               eval            = "1E-5", # e value threshold for ortholog detection
               ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
               aa_aln_type      = "pairwise",
               aa_aln_tool      = "NW", # use Needleman-Wunsch Algorithm for global codon alignment
               codon_aln_tool   = "pal2nal", 
               dnds_est.method  = "Comeron", # use robust dN/dS estimation (Comeron's method)
               output.folder    = "data/dNdS_maps/Bdistachyon",
               comp_cores       = 6
               )

```


### Generate 1:1 orthologs tables for C. rubella

```r
# compute dN/dS table of C. rubella vs. E. salsugineum, M. truncatula,  T. hassleriana, and B. distachyon
orthologr::map.generator(
               query_file      = "data/CDS/Crubella_183_v1.0.cds.fa.gz",
               subjects.folder = "data/CDS/subject_species_Crubella/",
               eval            = "1E-5", # e value threshold for ortholog detection
               ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
               aa_aln_type      = "pairwise",
               aa_aln_tool      = "NW", # use Needleman-Wunsch Algorithm for global codon alignment
               codon_aln_tool   = "pal2nal", 
               dnds_est.method  = "Comeron", # use robust dN/dS estimation (Comeron's method)
               output.folder    = "data/dNdS_maps/",
               comp_cores       = 6
               )
```


### Orthogroup Inference with Orthofinder2


```r
orthologr::translate_cds_to_protein_all()

orthologr::retrieve_longest_isoforms_all(proteome_folder = "devseq_orthofinder_protein_2019_06_27", 
                           annotation_folder = "devseq_orthologr_gtf_2019_06_19", 
                           output_folder = "devseq_orthofinder_protein_longest_2019_06_27",
                           annotation_format = "gtf")


orthologr::orthofinder2(proteome_folder = "devseq_orthofinder_protein_longest_2019_06_27", comp_cores = 4)
```

## Install `DevSeqR` package

The [DevSeqR package](https://github.com/HajkD/DevSeqR) allows users to reproduce all analyses and to perform
additional exploratory data analysis using the DevSeq dataset.

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
