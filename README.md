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
               query_file      = "../DevSeq_data/DevSeq_CDS/Query_files/Athaliana.fa",
               subjects_folder = "../DevSeq_data/DevSeq_CDS/Athaliana_subject_files",
               eval            = "1E-5", # e value threshold for ortholog detection
               ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
               aa_aln_type      = "pairwise",
               aa_aln_tool      = "NW", # use Needleman-Wunsch Algorithm for global codon alignment
               codon_aln_tool   = "pal2nal", 
               dnds_est_method  = "Comeron", # use robust dN/dS estimation (Comeron's method)
               min_qry_coverage_hsp = 70, # min query coverage of hsp >= 70% of initial query locus
               min_qry_perc_identity = 30, # min percent identity of hsp >= 30% to initial query locus
               output_folder    = "../DevSeq_data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana/",
               comp_cores       = 4
               )
               
               
# retrieve pairwise ortho tables and detect for each gene locus exactly one 
# representative splice variant that maximizes the sequence homology
# between the two orthologous loci
devseq_ortho_tables <- orthologr::generate_ortholog_tables_all(
    dNdS_folder = "../DevSeq_data/DevSeq_dNdS_maps/OrthologR/DevSeq_Plants/Query_Athaliana",
    annotation_file_query = "../DevSeq_data/DevSeq_GTF/Query_files/Athaliana.gtf",
    annotation_folder_subject = "../DevSeq_data/DevSeq_GTF/Athaliana_subject_files",
    output_folder = "../DevSeq_data/DevSeq_orthologs/OrthologR/DevSeq_Plants/",
    output_type = "gene_locus",
    format = c("gtf", "gtf")
)

readr::write_tsv(devseq_ortho_tables, "devseq_ortho_tables_all_splice_variants.tsv", col_names = TRUE)


# retrieve the core set of orthologs in tidy data format
devseq_core_orthologs <-
    orthologr::retrieve_core_orthologs(
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

readr::write_tsv(devseq_core_orthologs, "devseq_core_orthologs_representative_splice_variant.tsv", col_names = TRUE)

# look at the results
devseq_core_orthologs
```

Visualizing the inferred orthologs

```r
p_plot_pairwise_orthologs <- orthologr::plot_pairwise_orthologs(ortho_tables = devseq_ortho_tables,
        species_order = c(
            "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
        ), 
        n_core_orthologs = 7175)

cowplot::save_plot(
    "p_plot_pairwise_orthologs.pdf",
    p_plot_pairwise_orthologs,
    base_height = 9,
    base_width = 14
)
```


```r
# determining best set of thresholds for retaining orthologs
qcovhsp_70_perc_identity_30 <- dplyr::summarize(dplyr::group_by(dplyr::filter(devseq_ortho_tables,  qcovhsp >= 70, perc_identity >= 30), subject_species),
                                                n_orthologs = dplyr::n())
qcovhsp_70_perc_identity_50 <- dplyr::summarize(dplyr::group_by(dplyr::filter(devseq_ortho_tables,  qcovhsp >= 70, perc_identity >= 50), subject_species),
                                                n_orthologs = dplyr::n())
qcovhsp_70_perc_identity_70 <- dplyr::summarize(dplyr::group_by(dplyr::filter(devseq_ortho_tables,  qcovhsp >= 70, perc_identity >= 70), subject_species),
                                                n_orthologs = dplyr::n())
qcovhsp_90_perc_identity_70 <- dplyr::summarize(dplyr::group_by(dplyr::filter(devseq_ortho_tables,  qcovhsp >= 90, perc_identity >= 70), subject_species),
                                                n_orthologs = dplyr::n())
qcovhsp_90_perc_identity_90 <- dplyr::summarize(dplyr::group_by(dplyr::filter(devseq_ortho_tables,  qcovhsp >= 90, perc_identity >= 90), subject_species),
                                                n_orthologs = dplyr::n())


orthologr::plot_pairwise_orthologs(
        ortho_tables = devseq_ortho_tables,
        species_order = c(
                "Alyrata",
                "Crubella",
                "Esalsugineum",
                "Thassleriana",
                "Mtruncatula",
                "Bdistachyon"
        ),
        n_core_orthologs = NULL
) + ggplot2::geom_line(ggplot2::aes(x = subject_species,
                                    y = n_orthologs,
                                    group = 1), data = qcovhsp_70_perc_identity_30, size = 2)  + ggplot2::geom_point(size = 4, data = qcovhsp_70_perc_identity_30) +
        ggplot2::geom_text(
                ggplot2::aes(label = n_orthologs),
                data = qcovhsp_70_perc_identity_50,
                hjust = 0,
                vjust = 2,
                size = 3
        ) + ggplot2::geom_line(ggplot2::aes(x = subject_species,
                                            y = n_orthologs,
                                            group = 1), data = qcovhsp_70_perc_identity_50, size = 2)  + ggplot2::geom_point(size = 4, data = qcovhsp_70_perc_identity_50) +
        ggplot2::geom_text(
                ggplot2::aes(label = n_orthologs),
                data = qcovhsp_70_perc_identity_50,
                hjust = 0,
                vjust = 2,
                size = 3
        ) + ggplot2::geom_line(ggplot2::aes(x = subject_species,
                                            y = n_orthologs,
                                            group = 1), data = qcovhsp_70_perc_identity_70, size = 2)  + ggplot2::geom_point(size = 4, data = qcovhsp_70_perc_identity_70) +
        ggplot2::geom_text(
                ggplot2::aes(label = n_orthologs),
                data = qcovhsp_90_perc_identity_70,
                hjust = 0,
                vjust = 2,
                size = 3
        ) + ggplot2::geom_line(ggplot2::aes(x = subject_species,
                                            y = n_orthologs,
                                            group = 1), data = qcovhsp_90_perc_identity_70, size = 2)  + ggplot2::geom_point(size = 4, data = qcovhsp_90_perc_identity_70) +
        ggplot2::geom_text(
                ggplot2::aes(label = n_orthologs),
                data = qcovhsp_70_perc_identity_50,
                hjust = 0,
                vjust = 2,
                size = 3
        ) + ggplot2::geom_line(ggplot2::aes(x = subject_species,
                                            y = n_orthologs,
                                            group = 1), data = qcovhsp_90_perc_identity_90, size = 2)  + ggplot2::geom_point(size = 4, data = qcovhsp_90_perc_identity_90) +
        ggplot2::geom_text(
                ggplot2::aes(label = n_orthologs),
                data = qcovhsp_90_perc_identity_90,
                hjust = 0,
                vjust = 2,
                size = 3
        )


```


### Orthogroup Inference with Orthofinder2

#### Running Orthofinder2 for DevSeq Species

```r
# translate all CDS sequences into protein sequences and check
# for all CDS sequences whether or not they are divisible by 3
orthologr::translate_cds_to_protein_all(input_folder = "../DevSeq_data/DevSeq_CDS/Query_files", 
                                        output_folder = "../DevSeq_data/DevSeq_Proteins",
                                        delete_corrupt_cds = FALSE)

# for each translated protein file retrieve the longest isoform per gene locus  
orthologr::retrieve_longest_isoforms_all(
                           proteome_folder = "../DevSeq_data/DevSeq_Proteins", 
                           annotation_folder = "../DevSeq_data/DevSeq_GTF/Query_files", 
                           output_folder = "../DevSeq_data/DevSeq_Proteins_longest_isoforms",
                           annotation_format = "gtf")

# run orthofinder2 on protein sequences with longest splice variants 
orthologr::orthofinder2(proteome_folder = "../DevSeq_data/DevSeq_Proteins_longest_isoforms", comp_cores = 4)
```

#### Running Orthofinder2 for Brawand Species

```r
# translate all CDS sequences into protein sequences and check
# for all CDS sequences whether or not they are divisible by 3
orthologr::translate_cds_to_protein_all(input_folder = "../DevSeq_data/Brawand_CDS/Query_files", 
                                        output_folder = "../DevSeq_data/Brawand_Proteins",
                                        delete_corrupt_cds = FALSE)

# for each translated protein file retrieve the longest isoform per gene locus  
orthologr::retrieve_longest_isoforms_all(
                           proteome_folder = "../DevSeq_data/Brawand_Proteins", 
                           annotation_folder = "../DevSeq_data/Brawand_GTF/Query_files", 
                           output_folder = "../DevSeq_data/Brawand_Proteins_longest_isoforms",
                           annotation_format = "gtf")

# run orthofinder2 on protein sequences with longest splice variants 
orthologr::orthofinder2(proteome_folder = "../DevSeq_data/Brawand_Proteins_longest_isoforms", comp_cores = 4)
```

### Inference of orthologous lncRNAs

```r
# gene loci lncRNAs
orthologr::map_generator_lnc(
  query_file      = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Antisense_exonic/Athaliana.fa",
  subjects.folder = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Antisense_exonic/subject_files/",
  eval            = "1E-5", # e value threshold for ortholog detection
  ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
  output.folder    = "../DevSeq_data/DevSeq_lncRNA_orthologs/Antisense_exonic",
  comp_cores       = 4
)

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
