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

### Reference Genomes, CDS, lncRNAs, etc

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
orthologr::map_generator_dnds(
               query_file      = "data/DevSeq_CDS/Query_files/Athaliana.fa",
               subjects_folder = "data/DevSeq_CDS/Athaliana_subject_files",
               eval            = "1E-5", # e value threshold for ortholog detection
               ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
               aa_aln_type      = "pairwise",
               aa_aln_tool      = "NW", # use Needleman-Wunsch Algorithm for global codon alignment
               codon_aln_tool   = "pal2nal", 
               dnds_est_method  = "Comeron", # use robust dN/dS estimation (Comeron's method)
               min_qry_coverage_hsp = 50, # min query coverage of hsp >= 50% of initial query locus
               min_qry_perc_identity = 30, # min percent identity of hsp >= 30% to initial query locus
               output_folder    = "data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/",
               comp_cores       = 12
               )
               
               
# retrieve pairwise ortho tables and detect for each gene locus exactly one 
# representative splice variant that maximizes the sequence homology
# between the two orthologous loci
devseq_ortho_tables <- orthologr::generate_ortholog_tables_all(
    dNdS_folder = "data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana",
    annotation_file_query = "data/DevSeq_GTF/Query_files/Athaliana.gtf",
    annotation_folder_subject = "data/DevSeq_GTF/Athaliana_subject_files",
    output_folder = "data/DevSeq_orthologs/OrthologR/DevSeq_Plants/",
    output_type = "gene_locus",
    format = c("gtf", "gtf")
)


# retrieve the core set of orthologs in tidy data format
devseq_core_orthologs <-
    retrieve_core_orthologs(
        ortho_tables = devseq_ortho_tables,
        species_order = c(
            "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
        )
    )

# look at the results
devseq_core_orthologs
```



### Generate ortholog tables by gene locus 

The following procedures determine the best splice variant per gene locus
that maximizes the sequence homology between the query locus and the subject locus.
The respective splice variant in the query and subject will then be used to represent the 
homology relationship between the homologous gene loci.

```r
### QUERY: A thaliana

# Ath vs Aly
orthologr::generate_ortholog_tables(
  dNdS_file = "data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/map_q=Athaliana.fa_s=Alyrata.fa_orthodetec=RBH_eval=1E-5.csv",
  annotation_file_query = "data/DevSeq_GTF/Athaliana.gtf",
  annotation_file_subject = "data/DevSeq_GTF/Alyrata.gtf",
  output_folder = "data/DevSeq_orthologs/OrthologR/DevSeq_Plants/Query_Athaliana",
  format = c("gtf", "gtf")
)

# Ath vs Bdist
orthologr::generate_ortholog_tables(
  dNdS_file = "data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/map_q=Athaliana.fa_s=Bdistachyon.fa_orthodetec=RBH_eval=1E-5.csv,
  annotation_file_query = "data/DevSeq_GTF/Athaliana.gtf",
  annotation_file_subject = "data/DevSeq_GTF/Bdistachyon.gtf",
  output_folder = "data/DevSeq_orthologs/OrthologR/DevSeq_Plants/Query_Athaliana",
  format = c("gtf", "gtf")
)

# Ath vs Crubella
orthologr::generate_ortholog_tables(
  dNdS_file = "data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/map_q=Athaliana.fa_s=Crubella.fa_orthodetec=RBH_eval=1E-5.csv",
  annotation_file_query = "data/DevSeq_GTF/Athaliana.gtf",
  annotation_file_subject = "data/DevSeq_GTF/Crubella.gtf",
  output_folder = "data/DevSeq_orthologs/OrthologR/DevSeq_Plants/Query_Athaliana",
  format = c("gtf", "gtf")
)


# Ath vs Esalsugineum
orthologr::generate_ortholog_tables(
  dNdS_file = "data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/map_q=Athaliana.fa_s=Esalsugineum.fa_orthodetec=RBH_eval=1E-5.csv",
  annotation_file_query = "data/DevSeq_GTF/Athaliana.gtf",
  annotation_file_subject = "data/DevSeq_GTF/Esalsugineum.gtf",
  output_folder = "data/DevSeq_orthologs/OrthologR/DevSeq_Plants/Query_Athaliana",
  format = c("gtf", "gtf")
)


# Ath vs Mtruncatula
orthologr::generate_ortholog_tables(
  dNdS_file = "data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/map_q=Athaliana.fa_s=Mtruncatula.fa_orthodetec=RBH_eval=1E-5.csv",
  annotation_file_query = "data/DevSeq_GTF/Athaliana.gtf",
  annotation_file_subject = "data/DevSeq_GTF/Mtruncatula.gtf",
  output_folder = "data/DevSeq_orthologs/OrthologR/DevSeq_Plants/Query_Athaliana",
  format = c("gtf", "gtf")
)

# Ath vs Thassleriana
orthologr::generate_ortholog_tables(
  dNdS_file = "data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/map_q=Athaliana.fa_s=Thassleriana.fa_orthodetec=RBH_eval=1E-5.csv",
  annotation_file_query = "data/DevSeq_GTF/Athaliana.gtf",
  annotation_file_subject = "data/DevSeq_GTF/Thassleriana.gtf",
  output_folder = "data/DevSeq_orthologs/OrthologR/DevSeq_Plants/Query_Athaliana",
  format = c("gtf", "gtf")
)
```




```r

orthologs_by_gene_path <-
  file.path(
    "orthologs_by_gene_locus_qry_athaliana",
    list.files("orthologs_by_gene_locus_qry_athaliana/")
  )

# import and generate final ortholog table for A thaliana
ortholog_tbl_athaliana <-
  dplyr::bind_rows(sapply(orthologs_by_gene_path, function(x) {
    list(readr::read_delim(
      x,
      delim = ";",
      col_names = TRUE)
  }))


ortholog_tbl_athaliana <- dplyr::mutate(ortholog_tbl_athaliana, q_len = q_end - q_start + 1, scope = 1 - (abs(q_len - alig_length) / q_len))

# number of orthologs found in other species
ortholog_tbl_athaliana_pairwise <- dplyr::summarize(dplyr::group_by(ortholog_tbl_athaliana, subject_species),
                 n_orthologs = n())

ortholog_tbl_athaliana_pairwise$subject_species <- factor(
  ortholog_tbl_athaliana_pairwise$subject_species,
  levels = c(
    "Alyrata",
    "Crubella",
    "Esalsugineum",
    "Thassleriana",
    "Mtruncatula",
    "Bdistachyon"
  )
)
ortholog_tbl_athaliana_core <-
  dplyr::do(dplyr::group_by(ortholog_tbl_athaliana, query_id),
            retrieve_core_genes(.))

ortholog_tbl_athaliana_core <-
  dplyr::filter(ortholog_tbl_athaliana_core, !is.na(query_id))

```

### Detection of all `A. thaliana` genes that have intersecting orthologs with all other species

```r
# first rename colnames of individual dNdS maps
for (i in seq_along(Ath_pairwise_dNdS_maps)) {
   colnames(Ath_pairwise_dNdS_maps[[i]])[2:5] <- paste0(names(Ath_pairwise_dNdS_maps)[i], c("_subject_id","_dN", "_dS", "_dNdS"))
}

# combine all geneids into one file
all.maps <- dplyr::bind_rows(lapply(Ath_pairwise_dNdS_maps, function(x) tibble::as_tibble(unique(x$query_id))))
colnames(all.maps) <- "query_id"

# detect genes that have orthologs in all other species
length(names(table(all.maps$query_id))[which(table(all.maps$query_id) == length(Ath_pairwise_dNdS_maps))])
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

#### Running Orthofinder2 for DevSeq Species
```r
# translate all CDS sequences into protein sequences and check
# for all CDS sequences whether or not they are divisible by 3
orthologr::translate_cds_to_protein_all(input_folder = "data/DevSeq_CDS/Query_files", 
                                        output_folder = "data/DevSeq_Proteins",
                                        delete_corrupt_cds = FALSE)

# for each translated protein file retrieve the longest isoform per gene locus  
orthologr::retrieve_longest_isoforms_all(
                           proteome_folder = "data/DevSeq_Proteins", 
                           annotation_folder = "data/DevSeq_GTF/Query_files", 
                           output_folder = "data/DevSeq_Proteins_longest_isoforms",
                           annotation_format = "gtf")

# run orthofinder2 on protein sequences with longest splice variants 
orthologr::orthofinder2(proteome_folder = "data/DevSeq_Proteins_longest_isoforms", comp_cores = 4)
```

#### Running Orthofinder2 for Brawand Species



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
