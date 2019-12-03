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


# Orthology Inference and dNdS Estimation

## DevSeq

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
               min_qry_coverage_hsp = 30, # min query coverage of hsp >= 30% of initial query locus
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

# add scope column
devseq_ortho_tables <- dplyr::mutate(devseq_ortho_tables, scope = 1 - (abs(q_len - alig_length) / q_len))
```

### Determining best homology threshold for orthologs

```r
p_all_ortho_thresholds_devseq <- orthologr::plot_diverse_homology_thresholds(devseq_ortho_tables, species_order = c(
                        "Alyrata",
                        "Crubella",
                        "Esalsugineum",
                        "Thassleriana",
                        "Mtruncatula",
                        "Bdistachyon"
                ))

cowplot::save_plot(
    "p_all_ortho_thresholds_devseq.pdf",
    p_all_ortho_thresholds_devseq,
    base_height = 14,
    base_width = 18
)
```


### Gene locus vs splice variant threshold assessment for orthologs

```r
p_core_sets_multi <- orthologr::plot_diverse_homology_thresholds_core_orthologs(devseq_ortho_tables, species_order = c(
                "Alyrata",
                "Crubella",
                "Esalsugineum",
                "Thassleriana",
                "Mtruncatula",
                "Bdistachyon"
            ), type = "gene_locus")

cowplot::save_plot(
    "p_core_sets_multi_devseq.pdf",
    p_core_sets_multi,
    base_height = 9,
    base_width = 14
)

p_core_sets_multsplice_variant_devseq <- orthologr::plot_diverse_homology_thresholds_core_orthologs(devseq_ortho_tables, species_order = c(
                "Alyrata",
                "Crubella",
                "Esalsugineum",
                "Thassleriana",
                "Mtruncatula",
                "Bdistachyon"
            ), type = "both")


cowplot::save_plot(
    "p_core_sets_multsplice_variant_devseq.pdf",
    p_core_sets_multsplice_variant_devseq,
    base_height = 9,
    base_width = 14
)
```


### Retrieve final DevSeq ortholog tables

```r
# use homology thresholds:
devseq_ortho_tables_final <- dplyr::filter(devseq_ortho_tables,  qcovhsp >= 70, perc_identity >= 30, scope >= 0.7)
readr::write_tsv(devseq_ortho_tables_final, "devseq_ortho_tables_all_splice_variants.tsv", col_names = TRUE)


# retrieve the core set of orthologs in tidy data format
devseq_core_orthologs <-
    orthologr::retrieve_core_orthologs(
        ortho_tables = devseq_ortho_tables_final,
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


# test that core set only includes unique query_IDs and subject_IDs
dplyr::summarize(
  dplyr::group_by(devseq_core_orthologs, subject_species),
  n_unique_query_gene_locus_id = length(unique(query_gene_locus_id)),
  n_non_unique_query_gene_locus_id = length((query_gene_locus_id)),
  n_unique_query_id = length(unique(query_id)),
  n_non_unique_query_id = length((query_id)),
  n_unique_subject_gene_locus_id = length(unique(subject_gene_locus_id)),
  n_non_unique_subject_gene_locus_id = length((subject_gene_locus_id)),
  n_unique_subject_id = length(unique(subject_id)),
  n_non_unique_subject_id = length((subject_id))
)
```

### Visualizing the inferred orthologs

```r
p_plot_pairwise_orthologs <- orthologr::plot_pairwise_orthologs(ortho_tables = devseq_ortho_tables_final,
        species_order = c(
            "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
        ), 
        n_core_orthologs = length(unique(devseq_core_orthologs$query_gene_locus_id)))

cowplot::save_plot(
    "p_plot_pairwise_orthologs_devseq.pdf",
    p_plot_pairwise_orthologs,
    base_height = 8,
    base_width = 12
)


p_dNdS_5_Athaliana <- metablastr::gg_species_dnds_blast_tbl(
  dplyr::filter(devseq_ortho_tables_final, !is.na(dNdS), dNdS <= 5),
  type = "dNdS",
  order = c(
    "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
  ), xlab = "dNdS"
)

p_alig_length <- metablastr::gg_species_feature_blast_tbl(devseq_ortho_tables_final,
  type = "alig_length",
  order = c(
    "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
  ), xlab = "Alignment length in amino acids (BLAST hits)"
)


p_perc_identity <- metablastr::gg_species_feature_blast_tbl(devseq_ortho_tables_final,
  type = "perc_identity",
  order = c(
    "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
  ), xlab = "Percent Identity (BLAST hits)"
)

p_ortholog_divergence_devseq <- gridExtra::grid.arrange(p_alig_length, p_dNdS_5_Athaliana, p_perc_identity, ncol = 1)

cowplot::save_plot(
    "p_ortholog_divergence_devseq.pdf",
    p_ortholog_divergence_devseq,
    base_height = 20,
    base_width = 12
)


p_dNdS_5_devseq <- metablastr::gg_species_dnds_blast_tbl(
  dplyr::filter(devseq_ortho_tables_final, !is.na(dNdS), dNdS <= 5),
  type = "dNdS",
  order = c(
    "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
  ), xlab = "dNdS"
)

p_dN_5_devseq<- metablastr::gg_species_dnds_blast_tbl(
  dplyr::filter(devseq_ortho_tables_final, !is.na(dNdS), dNdS <= 5),
  type = "dN",
  order = c(
    "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
  ), xlab = "dN"
)

p_dS_5_devseq <- metablastr::gg_species_dnds_blast_tbl(
  dplyr::filter(devseq_ortho_tables_final, !is.na(dNdS), dNdS <= 5),
  type = "dS",
  order = c(
    "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
  ), xlab = "dS"
)


p_dN_plus_dS_5_devseq <- metablastr::gg_species_dnds_blast_tbl(
  dplyr::filter(devseq_ortho_tables_final, !is.na(dNdS), dNdS <= 5),
  type = "dN+dS",
  order = c(
    "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
  ), xlab = "dN+dS"
)

p_all_dNdS_devseq <- gridExtra::grid.arrange(p_dN_5_devseq, p_dS_5_devseq, p_dNdS_5_devseq, p_dN_plus_dS_5_devseq, nrow = 2)

cowplot::save_plot(
    "p_all_dNdS_devseq.pdf",
    p_all_dNdS_devseq,
    base_height = 16,
    base_width = 20
)
```

## EvoSeq

```
# compute dN/dS table of A. thaliana vs. all other species
orthologr::map_generator_dnds(
               query_file      = "../DevSeq_data/EvoSeq_CDS/Query_files/Athaliana.fa",
               subjects_folder = "../DevSeq_data/EvoSeq_CDS/Athaliana_subject_files",
               eval            = "1E-5", # e value threshold for ortholog detection
               ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
               aa_aln_type      = "pairwise",
               aa_aln_tool      = "NW", # use Needleman-Wunsch Algorithm for global codon alignment
               codon_aln_tool   = "pal2nal", 
               dnds_est_method  = "Comeron", # use robust dN/dS estimation (Comeron's method)
               min_qry_coverage_hsp = 70, # min query coverage of hsp >= 70% of initial query locus
               min_qry_perc_identity = 30, # min percent identity of hsp >= 30% to initial query locus
               output_folder    = "../DevSeq_data/EvoSeq_dNdS_maps/OrthologR/EvoSeq_Plants/Query_Athaliana/",
               comp_cores       = 4
               )
               
               
# retrieve pairwise ortho tables and detect for each gene locus exactly one 
# representative splice variant that maximizes the sequence homology
# between the two orthologous loci
evoseq_ortho_tables <- orthologr::generate_ortholog_tables_all(
    dNdS_folder = "../DevSeq_data/EvoSeq_dNdS_maps/OrthologR/EvoSeq_Plants/Query_Athaliana/",
    annotation_file_query = "../DevSeq_data/EvoSeq_GTF/Query_files/Athaliana.gtf",
    annotation_folder_subject = "../DevSeq_data/EvoSeq_GTF/Athaliana_subject_files",
    output_folder = "../DevSeq_data/EvoSeq_orthologs/OrthologR/EvoSeq_Plants/",
    output_type = "gene_locus",
    format = c("gtf", "gtf")
)

readr::write_tsv(evoseq_ortho_tables, "evoseq_ortho_tables_all_splice_variants.tsv", col_names = TRUE)


# retrieve the core set of orthologs in tidy data format
evoseq_core_orthologs <-
    orthologr::retrieve_core_orthologs(
        ortho_tables = evoseq_ortho_tables,
        species_order = c(
            "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
        )
    )

readr::write_tsv(evoseq_core_orthologs, "evoseq_core_orthologs_representative_splice_variant.tsv", col_names = TRUE)

# look at the results
evoseq_core_orthologs
```

### Running Orthologr for Brawand Species

```r
# compute dN/dS table of Human vs. all other species
orthologr::map_generator_dnds(
               query_file      = "../DevSeq_data/Brawand_CDS/Query_files/Human.fa",
               subjects_folder = "../DevSeq_data/Brawand_CDS/Human_subject_files",
               eval            = "1E-5", # e value threshold for ortholog detection
               ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
               aa_aln_type      = "pairwise",
               aa_aln_tool      = "NW", # use Needleman-Wunsch Algorithm for global codon alignment
               codon_aln_tool   = "pal2nal", 
               dnds_est_method  = "Comeron", # use robust dN/dS estimation (Comeron's method)
               min_qry_coverage_hsp = 30, # min query coverage of hsp >= 30% of initial query locus
               min_qry_perc_identity = 30, # min percent identity of hsp >= 30% to initial query locus
               output_folder    = "../DevSeq_data/Brawand_dNdS_maps/OrthologR/Brawand_Animals/Query_Human/",
               comp_cores       = 4
               )
               
               
# retrieve pairwise ortho tables and detect for each gene locus exactly one 
# representative splice variant that maximizes the sequence homology
# between the two orthologous loci
brawand_ortho_tables <- orthologr::generate_ortholog_tables_all(
    dNdS_folder = "../DevSeq_data/Brawand_dNdS_maps/OrthologR/Brawand_Animals/Query_Human/",
    annotation_file_query = "../DevSeq_data/Brawand_GTF/Query_files/Human.gtf",
    annotation_folder_subject = "../DevSeq_data/Brawand_GTF/Human_subject_files",
    output_folder = "../DevSeq_data/Brawand_orthologs/OrthologR/Brawand_Animals/",
    output_type = "gene_locus",
    format = c("gtf", "gtf")
)

brawand_ortho_tables <- dplyr::mutate(brawand_ortho_tables, scope = 1 - (abs(q_len - alig_length) / q_len))
```

### Determining best homology threshold for orthologs

```r
p_all_ortho_thresholds_brawand <- orthologr::plot_diverse_homology_thresholds(brawand_ortho_tables, 
     species_order = c(
            "Chimp",
            "Bonobo",
            "Gorilla",
            "Orangutan",
            "Macaque",
            "Mouse",
            "Opossum",
            "Chicken"
        ))

cowplot::save_plot(
    "p_all_ortho_thresholds_brawand.pdf",
    p_all_ortho_thresholds_brawand,
    base_height = 14,
    base_width = 18
)
```

### Gene locus vs splice variant threshold assessment for orthologs

```r
p_core_sets_multi_brawand <- orthologr::plot_diverse_homology_thresholds_core_orthologs(brawand_ortho_tables, species_order = c(
            "Chimp",
            "Bonobo",
            "Gorilla",
            "Orangutan",
            "Macaque",
            "Mouse",
            "Opossum",
            "Chicken"
        ), type = "gene_locus")

cowplot::save_plot(
    "p_core_sets_multi_brawand.pdf",
    p_core_sets_multi_brawand,
    base_height = 9,
    base_width = 14
)

p_core_sets_multsplice_variant_brawand <- orthologr::plot_diverse_homology_thresholds_core_orthologs(brawand_ortho_tables, species_order = c(
            "Chimp",
            "Bonobo",
            "Gorilla",
            "Orangutan",
            "Macaque",
            "Mouse",
            "Opossum",
            "Chicken"
        ), type = "both")


cowplot::save_plot(
    "p_core_sets_multsplice_variant_brawand.pdf",
    p_core_sets_multsplice_variant_brawand,
    base_height = 9,
    base_width = 14
)
```
### Retrieve final Brawand ortholog tables

```r
# use homology thresholds:
brawand_ortho_tables_final <- dplyr::filter(brawand_ortho_tables,  qcovhsp >= 70, perc_identity >= 30, scope >= 0.7)
readr::write_tsv(brawand_ortho_tables, "brawand_ortho_tables_all_splice_variants.tsv", col_names = TRUE)

# retrieve the core set of orthologs in tidy data format
brawand_core_orthologs <-
    orthologr::retrieve_core_orthologs(
        ortho_tables = brawand_ortho_tables_final,
        species_order = c(
            "Chimp",
            "Bonobo",
            "Gorilla",
            "Orangutan",
            "Macaque",
            "Mouse",
            "Opossum",
            "Chicken"
        )
    )

readr::write_tsv(brawand_core_orthologs, "brawand_core_orthologs_representative_splice_variant.tsv", col_names = TRUE)

# look at the results
brawand_core_orthologs

# test that core set only includes unique query_IDs and subject_IDs
dplyr::summarize(
  dplyr::group_by(brawand_core_orthologs, subject_species),
  n_unique_query_gene_locus_id = length(unique(query_gene_locus_id)),
  n_non_unique_query_gene_locus_id = length((query_gene_locus_id)),
  n_unique_query_id = length(unique(query_id)),
  n_non_unique_query_id = length((query_id)),
  n_unique_subject_gene_locus_id = length(unique(subject_gene_locus_id)),
  n_non_unique_subject_gene_locus_id = length((subject_gene_locus_id)),
  n_unique_subject_id = length(unique(subject_id)),
  n_non_unique_subject_id = length((subject_id))
)
```

### Visualizing Orthologr results

```r
p_n_orthos_human <- orthologr::plot_pairwise_orthologs(ortho_tables = brawand_ortho_tables_final,
        species_order = c(
            "Chimp",
            "Bonobo",
            "Gorilla",
            "Orangutan",
            "Macaque",
            "Mouse",
            "Opossum",
            "Chicken"
        ), 
        n_core_orthologs = length(unique(brawand_core_orthologs$query_gene_locus_id)))


cowplot::save_plot(
    "p_n_orthos_human_brawand.pdf",
    p_n_orthos_human,
    base_height = 8,
    base_width = 12
)



p_dNdS_human <- metablastr::gg_species_dnds_blast_tbl(
  dplyr::filter(brawand_ortho_tables, !is.na(dNdS)),
  type = "dN/dS",
  order = c(
    "Chimp",
    "Bonobo",
    "Gorilla",
    "Orangutan",
    "Macaque",
    "Mouse",
    "Opossum",
    "Chicken"
  ), xlab = "dN/dS"
)

p_dNdS_5_Human_brawand <- metablastr::gg_species_dnds_blast_tbl(
  dplyr::filter(brawand_ortho_tables, !is.na(dNdS), dNdS <= 5),
  type = "dN/dS",
  order = c(
    "Chimp",
    "Bonobo",
    "Gorilla",
    "Orangutan",
    "Macaque",
    "Mouse",
    "Opossum",
    "Chicken"
  ), xlab = "dN/dS"
)
p_all_human <- gridExtra::grid.arrange(p_n_orthos_human, p_dNdS_human, p_dNdS_5_Human_brawand)



p_dN_5_human <- metablastr::gg_species_dnds_blast_tbl(
  dplyr::filter(brawand_ortho_tables, !is.na(dNdS), dNdS <= 5),
  type = "dN",
  order = c(
    "Chimp",
    "Bonobo",
    "Gorilla",
    "Orangutan",
    "Macaque",
    "Mouse",
    "Opossum",
    "Chicken"
  ), xlab = "dN"
)

p_dS_5_human <- metablastr::gg_species_dnds_blast_tbl(
  dplyr::filter(brawand_ortho_tables, !is.na(dNdS), dNdS <= 5),
  type = "dS",
  order = c(
    "Chimp",
    "Bonobo",
    "Gorilla",
    "Orangutan",
    "Macaque",
    "Mouse",
    "Opossum",
    "Chicken"
  ), xlab = "dS"
)


p_dN_plus_dS_5_human <- metablastr::gg_species_dnds_blast_tbl(
  dplyr::filter(brawand_ortho_tables, !is.na(dNdS), dNdS <= 5),
  type = "dN+dS",
  order = c(
    "Chimp",
    "Bonobo",
    "Gorilla",
    "Orangutan",
    "Macaque",
    "Mouse",
    "Opossum",
    "Chicken"
  ), xlab = "dN+dS"
)

p_all_dNdS_human <- gridExtra::grid.arrange(p_dN_5_human, p_dS_5_human, p_dNdS_5_human, p_dN_plus_dS_5_human, nrow = 2)

cowplot::save_plot(
    "p_all_dNdS_human.pdf",
    p_all_dNdS_human,
    base_height = 16,
    base_width = 20
)


p_alig_length_brawand <- metablastr::gg_species_feature_blast_tbl(brawand_ortho_tables,
  type = "alig_length",
  order = c(
    "Chimp",
    "Bonobo",
    "Gorilla",
    "Orangutan",
    "Macaque",
    "Mouse",
    "Opossum",
    "Chicken"
  ), xlab = "Alignment Length (BLAST hits)"
)


p_perc_identity_brawand <- metablastr::gg_species_feature_blast_tbl(brawand_ortho_tables,
  type = "perc_identity",
  order = c(
    "Chimp",
    "Bonobo",
    "Gorilla",
    "Orangutan",
    "Macaque",
    "Mouse",
    "Opossum",
    "Chicken"
  ), xlab = "Percent Identity (BLAST hits)"
)



p_ortholog_divergence_brawand <- gridExtra::grid.arrange(p_alig_length_brawand, p_dNdS_5_Human_brawand, p_perc_identity_brawand)

cowplot::save_plot(
    "p_ortholog_divergence_brawand.pdf",
    p_ortholog_divergence_brawand,
    base_height = 20,
    base_width = 12
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


### Visualize Pairwise Orthogroups
path_to_devseq_pairwise_orthos <- "../DevSeq_data/DevSeq_orthologs/Orthofinder2/DevSeq_Plants/Results_DevSeq_Proteins_longest_isoforms/Orthologues/Orthologues_Athaliana/"
devseq_pairwise_files <- list.files(path_to_devseq_pairwise_orthos)

devseq_pairwise_df <- lapply(file.path(path_to_devseq_pairwise_orthos, devseq_pairwise_files), readr::read_tsv)
names(devseq_pairwise_df) <- c("Alyrata", "Bdistachyon", "Crubella", "Esalsugineum", "Mtruncatula", "Thassleriana")

devseq_pairwise_n_orthogroups <- lapply(devseq_pairwise_df, nrow)

devseq_pairwise_df <- tibble::tibble(subject_species = names(devseq_pairwise_n_orthogroups), n_orthologs = unlist(devseq_pairwise_n_orthogroups))


devseq_pairwise_df$subject_species <- factor(devseq_pairwise_df$subject_species,
                                             levels = c("Alyrata",
                                                        "Crubella",
                                                        "Esalsugineum",
                                                        "Thassleriana",
                                                        "Mtruncatula",
                                                        "Bdistachyon"))
### retrieve core orthogroups
og_file_devseq <- "../DevSeq_data/DevSeq_orthologs/Orthofinder2/DevSeq_Plants/Results_DevSeq_Proteins_longest_isoforms/Orthogroups/Orthogroups.tsv"
sc_file_devseq <- "../DevSeq_data/DevSeq_orthologs/Orthofinder2/DevSeq_Plants/Results_DevSeq_Proteins_longest_isoforms/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"

DevSeq_OF2_CoreOrthologs <- orthologr::orthofinder2_retrieve_core_orthologs(orthogroups_file = og_file_devseq, single_copy_file = sc_file_devseq)

### visualize pairwise orthologs
devseq_plot <-
        metablastr::gg_pairwise_orthologs_line(devseq_pairwise_df,
                                               vline = nrow(DevSeq_OF2_CoreOrthologs),
                                               ymax = 26000,
                                               title = "Orthofinder2: Pairwise genome comparisons: A. thaliana vs. Subject Species")

# store DevSeq Orthofinder2 core orthogroups
readr::write_tsv(DevSeq_OF2_CoreOrthologs, "devseq_OF2_CoreOrthologs.tsv",col_names = TRUE)

cowplot::save_plot(
    "p_orthofinder2_orthologs_devseq.pdf",
    devseq_plot,
    base_height = 8,
    base_width = 12
)
```

#### Running Orthofinder2 for Brawand Species

Extracting coding sequences from genome and gtf files.

```r
# Bonobo
extractSeqsFromAnnotation(
  annotation = "../DevSeq_data/Brawand_GTF/Pan_paniscus.panpan1.1.92.gtf",
  genome_fasta = "../DevSeq_data/Brawand_Genomes/Pan_paniscus.panpan1.1.dna.toplevel.fa",
  outputFile = "../DevSeq_data/Brawand_CDS/Bonobo.fa",
  format = "gtf",
  metaFeature = "tx",
  feature = "CDS"
)

# Chicken
extractSeqsFromAnnotation(
  annotation = "../DevSeq_data/Brawand_GTF/Gallus_gallus.Gallus_gallus-5.0.92.gtf",
  genome_fasta = "../DevSeq_data/Brawand_Genomes/Gallus_gallus.Gallus_gallus-5.0.dna.toplevel.fa",
  outputFile = "../DevSeq_data/Brawand_CDS/Chicken.fa",
  format = "gtf",
  metaFeature = "tx",
  feature = "CDS"
)

# Chimp
extractSeqsFromAnnotation(
  annotation = "../DevSeq_data/Brawand_GTF/Pan_troglodytes.Pan_tro_3.0.92.gtf",
  genome_fasta = "../DevSeq_data/Brawand_Genomes/Pan_troglodytes.Pan_tro_3.0.dna.toplevel.fa",
  outputFile = "../DevSeq_data/Brawand_CDS/Chimp.fa",
  format = "gtf",
  metaFeature = "tx",
  feature = "CDS"
)

# Gorilla
extractSeqsFromAnnotation(
  annotation = "../DevSeq_data/Brawand_GTF/Gorilla_gorilla.gorGor4.92.gtf",
  genome_fasta = "../DevSeq_data/Brawand_Genomes/Gorilla_gorilla.gorGor4.dna.toplevel.fa",
  outputFile = "../DevSeq_data/Brawand_CDS/Gorilla.fa",
  format = "gtf",
  metaFeature = "tx",
  feature = "CDS"
)

# Human
extractSeqsFromAnnotation(
  annotation = "../DevSeq_data/Brawand_GTF/Homo_sapiens.GRCh38.92_ensembl.gtf",
  genome_fasta = "../DevSeq_data/Brawand_Genomes/Homo_sapiens.GRCh38.dna.toplevel.fa",
  outputFile = "../DevSeq_data/Brawand_CDS/Human.fa",
  format = "gtf",
  metaFeature = "tx",
  feature = "CDS"
)

# Macaque
extractSeqsFromAnnotation(
  annotation = "../DevSeq_data/Brawand_GTF/Macaca_mulatta.Mmul_8.0.1.92.gtf",
  genome_fasta = "../DevSeq_data/Brawand_Genomes/Macaca_mulatta.Mmul_8.0.1.dna.toplevel.fa",
  outputFile = "../DevSeq_data/Brawand_CDS/Macaque.fa",
  format = "gtf",
  metaFeature = "tx",
  feature = "CDS"
)

# Mouse
extractSeqsFromAnnotation(
  annotation = "../DevSeq_data/Brawand_GTF/Mus_musculus.GRCm38.92.gtf",
  genome_fasta = "../DevSeq_data/Brawand_Genomes/Mus_musculus.GRCm38.dna.toplevel.fa",
  outputFile = "../DevSeq_data/Brawand_CDS/Mouse.fa",
  format = "gtf",
  metaFeature = "tx",
  feature = "CDS"
)

# Orangutan
extractSeqsFromAnnotation(
  annotation = "../DevSeq_data/Brawand_GTF/Pongo_abelii.PPYG2.92.gtf",
  genome_fasta = "../DevSeq_data/Brawand_Genomes/Pongo_abelii.PPYG2.dna.toplevel.fa",
  outputFile = "../DevSeq_data/Brawand_CDS/Orangutan.fa",
  format = "gtf",
  metaFeature = "tx",
  feature = "CDS"
)

# Opossum
extractSeqsFromAnnotation(
  annotation = "../DevSeq_data/Brawand_GTF/Monodelphis_domestica.monDom5.92.gtf",
  genome_fasta = "../DevSeq_data/Brawand_Genomes/Monodelphis_domestica.monDom5.dna.toplevel.fa",
  outputFile = "../DevSeq_data/Brawand_CDS/Opossum.fa",
  format = "gtf",
  metaFeature = "tx",
  feature = "CDS"
)
```


```r
# translate all CDS sequences into protein sequences and check
# for all CDS sequences whether or not they are divisible by 3
orthologr::translate_cds_to_protein_all(input_folder = "../DevSeq_data/Brawand_CDS/Query_files", 
                                        output_folder = "../DevSeq_data/Brawand_Proteins",
                                        delete_corrupt_cds = FALSE)

file.rename(from = "../DevSeq_data/Brawand_GTF/Query_files/Gallus_gallus.Gallus_gallus-5.0.92.gtf", to = "../DevSeq_data/Brawand_GTF/Query_files/Chicken.gtf")
file.rename(from = "../DevSeq_data/Brawand_GTF/Query_files/Gorilla_gorilla.gorGor4.92.gtf", to = "../DevSeq_data/Brawand_GTF/Query_files/Gorilla.gtf")
file.rename(from = "../DevSeq_data/Brawand_GTF/Query_files/Homo_sapiens.GRCh38.92_ensembl.gtf", to = "../DevSeq_data/Brawand_GTF/Query_files/Human.gtf")
file.rename(from = "../DevSeq_data/Brawand_GTF/Query_files/Macaca_mulatta.Mmul_8.0.1.92.gtf", to = "../DevSeq_data/Brawand_GTF/Query_files/Macaque.gtf")
file.rename(from = "../DevSeq_data/Brawand_GTF/Query_files/Monodelphis_domestica.monDom5.92.gtf", to = "../DevSeq_data/Brawand_GTF/Query_files/Opossum.gtf")
file.rename(from = "../DevSeq_data/Brawand_GTF/Query_files/Mus_musculus.GRCm38.92.gtf", to = "../DevSeq_data/Brawand_GTF/Query_files/Mouse.gtf")
file.rename(from = "../DevSeq_data/Brawand_GTF/Query_files/Pan_paniscus.panpan1.1.92.gtf", to = "../DevSeq_data/Brawand_GTF/Query_files/Bonobo.gtf")
file.rename(from = "../DevSeq_data/Brawand_GTF/Query_files/Pan_troglodytes.Pan_tro_3.0.92.gtf", to = "../DevSeq_data/Brawand_GTF/Query_files/Chimp.gtf")
file.rename(from = "../DevSeq_data/Brawand_GTF/Query_files/Pongo_abelii.PPYG2.92.gtf", to = "../DevSeq_data/Brawand_GTF/Query_files/Orangutan.gtf")

# for each translated protein file retrieve the longest isoform per gene locus  
orthologr::retrieve_longest_isoforms_all(
                           proteome_folder = "../DevSeq_data/Brawand_Proteins", 
                           annotation_folder = "../DevSeq_data/Brawand_GTF/Query_files", 
                           output_folder = "../DevSeq_data/Brawand_Proteins_longest_isoforms",
                           annotation_format = "gtf")

# run orthofinder2 on protein sequences with longest splice variants 
orthologr::orthofinder2(proteome_folder = "../DevSeq_data/Brawand_Proteins_longest_isoforms", comp_cores = 4)

# Brawand Core Orthologs
og_file_brawand <- "../DevSeq_data/Brawand_orthologs/Orthofinder2/Results_Brawand_Proteins_longest_isoforms/Orthogroups/Orthogroups.tsv"
sc_file_brawand <- "../DevSeq_data/Brawand_orthologs/Orthofinder2/Results_Brawand_Proteins_longest_isoforms/Orthogroups/Orthogroups_SingleCopyOrthologues.txt"

Brawand_OF2_CoreOrthologs <- orthologr::orthofinder2_retrieve_core_orthologs(orthogroups_file = og_file_brawand, single_copy_file = sc_file_brawand)

readr::write_tsv(Brawand_OF2_CoreOrthologs, "brawand_OF2_CoreOrthologs.tsv",col_names = TRUE)

### Brawand Pairwise Orthologs
Brawand_Statistics_PerSpecies <- readr::read_tsv("~/Desktop/Projects/orthofinder/orthofinder_results/Animals_Result/Results_Jan30/Comparative_Genomics_Statistics/Statistics_PerSpecies.tsv", n_max = 10)
names(Brawand_Statistics_PerSpecies)[1] <- "Type"

path_to_brawand_pairwise_orthos <- "../DevSeq_data/Brawand_orthologs/Orthofinder2/Results_Brawand_Proteins_longest_isoforms/Orthologues/Orthologues_Human"
brawand_pairwise_files <- list.files(path_to_brawand_pairwise_orthos)

brawand_pairwise_df <- lapply(file.path(path_to_brawand_pairwise_orthos, brawand_pairwise_files), readr::read_tsv)
names(brawand_pairwise_df) <- c("Bonobo", "Chicken", "Chimp", "Gorilla", "Macaque", "Mouse", "Opossum", "Orangutan")

brawand_pairwise_n_orthogroups <- lapply(brawand_pairwise_df, nrow)

brawand_pairwise_df <- tibble::tibble(subject_species = names(brawand_pairwise_n_orthogroups), n_orthologs = unlist(brawand_pairwise_n_orthogroups))


brawand_pairwise_df$subject_species <- factor(brawand_pairwise_df$subject_species,
                                             levels = c("Chimp",
                                                        "Bonobo",
                                                        "Gorilla",
                                                        "Orangutan",
                                                        "Macaque",
                                                        "Mouse",
                                                        "Opossum",
                                                        "Chicken"))


brawand_plot <- metablastr::gg_pairwise_orthologs_line(brawand_pairwise_df, vline = nrow(Brawand_OF2_CoreOrthologs), title = "Orthofinder2: Pairwise genome comparisons: Human vs. Subject Species")

cowplot::save_plot(
    "p_orthofinder2_orthologs_brawand.pdf",
    brawand_plot,
    base_height = 8,
    base_width = 12
)
```



### Inference of Orthologous lncRNAs for DevSeq Data

We used `*_No_TE_genes_tpm_sample_names.csv` files.

```r
### Use E-value threshold 1E-2
# gene loci lncRNAs
lnc_map_1e2 <- orthologr::map_generator_lnc(
  query_file      = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Query_files/Athaliana.fa",
  subjects_folder = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Athaliana_subject_files/",
  eval            = "1E-2", # e value threshold for ortholog detection
  ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
  output_folder    = "../DevSeq_data/DevSeq_lncRNA_orthologs/lncRNA_gene_loci_1e2",
  min_qry_coverage_hsp = 30,
  min_qry_perc_identity = 30,
  comp_cores       = 4
)


devseq_core_orthologs_lnc_1e2 <-
    orthologr::lnc_map_core_orthologs(
        lnc_map = lnc_map_1e2,
        species_order = c(
            "Alyrata",
            "Crubella",
            "Esalsugineum",
            "Thassleriana",
            "Mtruncatula",
            "Bdistachyon"
        )
    )

all_lncRNAs_df_1e2 <- orthologr::lnc_map_counts(lnc_map_1e2, species_order = c(
    "Alyrata",
    "Crubella",
    "Esalsugineum",
    "Thassleriana",
    "Mtruncatula",
    "Bdistachyon"
  ))


### Use E-value threshold 1E-3
# gene loci lncRNAs
lnc_map_1e3 <- orthologr::map_generator_lnc(
  query_file      = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Query_files/Athaliana.fa",
  subjects_folder = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Athaliana_subject_files/",
  eval            = "1E-3", # e value threshold for ortholog detection
  ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
  output_folder    = "../DevSeq_data/DevSeq_lncRNA_orthologs/lncRNA_gene_loci_1e3",
  min_qry_coverage_hsp = 30,
  min_qry_perc_identity = 30,
  comp_cores       = 4
)

    
all_lncRNAs_df_1e3 <- orthologr::lnc_map_counts(lnc_map_1e3, species_order = c(
    "Alyrata",
    "Crubella",
    "Esalsugineum",
    "Thassleriana",
    "Mtruncatula",
    "Bdistachyon"
  ))

### Use E-value threshold 1E-5
# gene loci lncRNAs
lnc_map_1e5 <- orthologr::map_generator_lnc(
  query_file      = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Query_files/Athaliana.fa",
  subjects_folder = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Athaliana_subject_files/",
  eval            = "1E-5", # e value threshold for ortholog detection
  ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
  output_folder    = "../DevSeq_data/DevSeq_lncRNA_orthologs/lncRNA_gene_loci_1e5",
  min_qry_coverage_hsp = 30,
  min_qry_perc_identity = 30,
  comp_cores       = 4
)


all_lncRNAs_df_1e5 <- orthologr::lnc_map_counts(lnc_map_1e5, species_order = c(
    "Alyrata",
    "Crubella",
    "Esalsugineum",
    "Thassleriana",
    "Mtruncatula",
    "Bdistachyon"
  ))
  
### Use E-value threshold 1E-7
# gene loci lncRNAs
lnc_map_1e7 <- orthologr::map_generator_lnc(
  query_file      = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Query_files/Athaliana.fa",
  subjects_folder = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Athaliana_subject_files/",
  eval            = "1E-7", # e value threshold for ortholog detection
  ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
  output_folder    = "../DevSeq_data/DevSeq_lncRNA_orthologs/lncRNA_gene_loci_1e7",
  min_qry_coverage_hsp = 30,
  min_qry_perc_identity = 30,
  comp_cores       = 4
)


all_lncRNAs_df_1e7 <- orthologr::lnc_map_counts(lnc_map_1e7, species_order = c(
    "Alyrata",
    "Crubella",
    "Esalsugineum",
    "Thassleriana",
    "Mtruncatula",
    "Bdistachyon"
  ))
  

### Use E-value threshold 1E-10
# gene loci lncRNAs
lnc_map_1e10 <- orthologr::map_generator_lnc(
  query_file      = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Query_files/Athaliana.fa",
  subjects_folder = "../DevSeq_data/DevSeq_lncRNAs/lncRNA_gene_loci/Athaliana_subject_files/",
  eval            = "1E-10", # e value threshold for ortholog detection
  ortho_detection = "RBH", # use conservative method: BLAST best reciprocal hit
  output_folder    = "../DevSeq_data/DevSeq_lncRNA_orthologs/lncRNA_gene_loci_1e10",
  min_qry_coverage_hsp = 30,
  min_qry_perc_identity = 30,
  comp_cores       = 4
)


all_lncRNAs_df_1e10 <- orthologr::lnc_map_counts(lnc_map_1e10, species_order = c(
    "Alyrata",
    "Crubella",
    "Esalsugineum",
    "Thassleriana",
    "Mtruncatula",
    "Bdistachyon"
  ))
  
metablastr::gg_pairwise_orthologs_line(all_lncRNAs_df_1e2, title = "Number of orthologous lncRNAs: A. thaliana vs Subject Species") + 
ggplot2::geom_line(ggplot2::aes(x = subject_species,
                                    y = n_orthologs,
                                    group = 1), 
                                    data = all_lncRNAs_df_1e3, size = 2)  + 
                                    ggplot2::geom_point(size = 4, data = all_lncRNAs_df_1e3) +
        ggplot2::geom_text(
                ggplot2::aes(label = n_orthologs),
                data = all_lncRNAs_df_1e3,
                hjust = 0,
                vjust = 2,
                size = 2
        )  + 
ggplot2::geom_line(ggplot2::aes(x = subject_species,
                                    y = n_orthologs,
                                    group = 1), 
                                    data = all_lncRNAs_df_1e5, size = 2)  + 
                                    ggplot2::geom_point(size = 4, data = all_lncRNAs_df_1e5) +
        ggplot2::geom_text(
                ggplot2::aes(label = n_orthologs),
                data = all_lncRNAs_df_1e5,
                hjust = 0,
                vjust = 2,
                size = 2
        )  + 
ggplot2::geom_line(ggplot2::aes(x = subject_species,
                                    y = n_orthologs,
                                    group = 1), 
                                    data = all_lncRNAs_df_1e7, size = 2)  + 
                                    ggplot2::geom_point(size = 4, data = all_lncRNAs_df_1e7) +
        ggplot2::geom_text(
                ggplot2::aes(label = n_orthologs),
                data = all_lncRNAs_df_1e7,
                hjust = 0,
                vjust = 2,
                size = 2
        )  + 
ggplot2::geom_line(ggplot2::aes(x = subject_species,
                                    y = n_orthologs,
                                    group = 1), 
                                    data = all_lncRNAs_df_1e10, size = 2)  + 
                                    ggplot2::geom_point(size = 4, data = all_lncRNAs_df_1e10) +
        ggplot2::geom_text(
                ggplot2::aes(label = n_orthologs),
                data = all_lncRNAs_df_1e10,
                hjust = 0,
                vjust = 2,
                size = 2
        ) 

p_all_lncRNA_plot <- metablastr::gg_pairwise_orthologs_line(all_lncRNAs_df_1e2, title = "Number of orthologous lncRNAs (only NATs): A. thaliana vs Subject Species")

cowplot::save_plot(
  "p_all_lncRNA_plot.pdf",
  p_all_lncRNA_plot,
  base_height = 8,
  base_width = 14
)
```

## Promotor Evolution Analysis

### Extracting promotor sequences of all genes from all DevSeq species

```r
anno_files <- file.path("/Volumes/BKP_1/BKP_1/devseq_data/DevSeq_GTF/Query_files", list.files("/Volumes/BKP_1/BKP_1/devseq_data/DevSeq_GTF/Query_files"))
genome_files <- file.path("/Volumes/BKP_1/BKP_1/devseq_data/DevSeq_Genomes/Query_files", list.files("/Volumes/BKP_1/BKP_1/devseq_data/DevSeq_Genomes/Query_files"))
# c(250, 500, 1000, 1500, 2000)
promotor_length_vector <- c(250, 500, 1000)

working_dir <- getwd()

dir.create("Promotor_Sequences")
message("Starting promotor extraction for promotor lengths: ", paste0(promotor_length_vector, collapse = ", "))
message(" and species: ", paste0(genome_files, collapse = ", "), " ...")
for (j in seq_len(length(promotor_length_vector))) {
message("Processing promotor length ", promotor_length_vector[j], " ...")

if (!file.exists(file.path("Promotor_Sequences", paste0("Promotor_length_", promotor_length_vector[j]))))
   dir.create(file.path("Promotor_Sequences", paste0("Promotor_length_", promotor_length_vector[j])))
   
  for (i in seq_len(length(anno_files))) {
    species_name <-
      unlist(stringr::str_split(basename(genome_files[i]), "[.]"))[1]
      message("Processing species ", species_name, " ...")

   setwd(file.path("Promotor_Sequences", paste0("Promotor_length_", promotor_length_vector[j])))
   
    metablastr::extract_promotor_seqs_from_genome(
      annotation_file = anno_files[i],
      genome_file = genome_files[i],
      annotation_format = "gtf",
      promotor_length = promotor_length_vector[j]
      )
      
    setwd(working_dir)  
    cat("\n")
    cat("\n")
  }
  message("Promotor extraction terminated successfully!")
}
```

```r
promotor_homology_ath_500 <-
  orthologr::promotor_divergence_of_orthologous_genes(
        promotor_folder = "Promotor_Sequences/Promotor_length_500/",
        ortholog_tables_folder = "orthologs_by_gene_locus_qry_athaliana")
        
readr::write_delim(promotor_homology_ath_500, "ortholog_table_promotor_homology_ath_500.csv", delim = ";")
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
