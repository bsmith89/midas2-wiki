This page is still under development

***

# MIDAS DB

At its core, MIDAS DB is an interface between collections of genomes and a set of custom files needed for strain-level metagenomics analysis. The new implementation of MIDAS DB reads in a Table Of Contents (TOC) file, containing genome-to-species assignment and a choice of representative genome for each species. This new infrastructure dramatically simplifies the prior knowledge needed to build a custom MIDAS database, and enables the dynamic assignment of representative genomes. Presently, MIDAS DB contains two public databases readily available for users: UHGG v1 (4644 species with 286997 genomes) and GTDB v202 (47893 species with 258405 genomes) in a public AWS S3 bucket.  MIDAS is the first pipeline to allow users the flexibility to select representative genomes, a key component for the accuracy of SNP calling, as demonstrated in Supplementary Table S2. 


## UHGG

A collection of 286,997 genomes assembled from metagenomes, isolates and single cells from human stool samples has been clustered into 4,644 species in an effort similar to [IGGdb 1.0](https://github.com/snayfach/IGGdb).   We refer to this new collection as [IGGdb 2.0](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/), IGG 2.0, IGG+, or simply IGG.  Perhaps the most important difference with respect to the original IGGdb 1.0 is that this new collection contains only gut genomes.

This repository contains tools for building a [MIDAS](https://github.com/snayfach/MIDAS) database from IGG 2.0, geared to run on [AWS Batch](https://aws.amazon.com/batch/) as described in [PairANI](https://github.com/czbiohub/pairani/wiki).  The resulting database will be located at [s3://microbiome-igg/2.0](http://microbiome-igg.s3.amazonaws.com/2.0/README.TXT).

Our project to update [MIDAS for IGGdb 1.0](https://github.com/czbiohub/MIDAS-IGGdb/blob/master/README.md) is on hold.  We will first update MIDAS to run with this new database, directly from S3.

# Inputs

The IGG 2.0 genomes, genome-to-species assignments, and a choice of representative genome for each species, were provided by [Alexandre Almeida](https://www.ebi.ac.uk/about/people/alexandre-almeida) of EBI and mirrored in S3 by [Jason Shi](http://docpollard.org/people/jason-shi/), whose `s3://jason.shi-bucket/IGGdb2.0/clean_set/` serves as input to the tools in this repository. 

Six-digit numeric species ids were arbitrarily assigned by Jason Shi in `s3://jason.shi-bucket/IGGdb2.0/alt_species_ids.tsv`.


# Local Layout

```
Output                                          Producer            Meaning
------------------------------------------------------------------------------------------------------------
midas_iggdb                                     DB-related files    Mirror s3://miocriombe-igg/2.0/
```


# Target Layout in S3

A table of contents listing all genomes will be located at [s3://microbiome-igg/2.0/genomes.tsv](http://microbiome-igg.s3.amazonaws.com/2.0/genomes.tsv).  It will look like
```
genome              species      representative        genome_is_representative
GUT_GENOME138501    104351       GUT_GENOME269084      0
GUT_GENOME269084    104351       GUT_GENOME269084      1
...
```
with 286,997 rows below the headers (one per genome) and 4,644 distinct species.

Every species has a unique representative genome.

Column `genome_is_representative` is defined as `genome == representative`.

All inputs would be relocated under
```
s3://microbiome-igg/2.0/cleaned_genomes/<species>/<genome>.fa.lz4
```
for instance, to match the earlier example,
```
s3://microbiome-igg/2.0/cleaned_genomes/104351/GUT_GENOME138501.fa.lz4
s3://microbiome-igg/2.0/cleaned_genomes/104351/GUT_GENOME269084.fa.lz4
...
```
In the relocated genomes, all contig headers have been replaced by serially assigned ids of the form `<genome_id>_<contig_number>_<contig_length>_<sequence_hash>` where the contig_numbers increment from 1 in each genome, and sequence hash is for quick identification of duplicates.
```
>UHGG000001_C0_L209.5k_Hcfe300
>UHGG000001_C1_L163.2k_H1e6257
...
```

## Gene Annotations

Genes were identified by [Prokka](https://github.com/tseemann/prokka).

```
s3://microbiome-igg/2.0/gene_annotations/104351/GUT_GENOME138501/GUT_GENOME138501.{faa, ffn, fna, gff, tsv}.lz4
```

## Marker Genes Identification

Homologs of a specified marker gene set (typically universal single copy genes) are identified for each input genome, cross-referenced, collated, and indexed with hs-blastn to use in species abundance estimation for metagenomic samples.   This involves the following steps.

15 Phyeco single copy marker genes are available at:
```
s3://microbiome-igg/2.0/marker_gene_models/phyeco/marker_genes.hmm.lz4
s3://microbiome-igg/2.0/marker_gene_models/phyeco/marker_genes.mapping_cutoffs.lz4

```

For each genome, run hmmsearch against a selected marker gene set (typically the phyeco set mentioned above), yielding

```
s3://microbiome-igg/2.0/marker_genes/<marker_set>/temp/104351/GUT_GENOME138501/GUT_GENOME138501.{hmmsearch, markers.fa, markers.map}
```

where `GUT_GENOME{XXXXXX}.markers.fa` contains the prokka-annotated genes that map to marker genes, and `GUT_GENOME{XXXXXX}.markers.map` is a TSV map of prokka gene id to marker gene id.

These are concatenated, across all genomes, into monolithic `marker_genes.fa`, from which hs-blastn index is constructed.

```
s3://microbiome-igg/2.0/marker_genes/<marker_set>/marker_genes.fa
s3://microbiome-igg/2.0/marker_genes/<marker_set>/marker_genes.fa.{sa, bwt, sequence}
```

## Pan-Genomes

```
s3://microbiome-igg/2.0/pangenomes/104351/{genes.ffn, centroids.ffn, gene_info.txt}
s3://microbiome-igg/2.0/pangenomes/104351/temp/centroids.{99, 95, 90, 85, 80, 75}.ffn
s3://microbiome-igg/2.0/pangenomes/104351/temp/uclust.{99, 95, 90, 85, 80, 75}.txt
```
The species pangenome `{YYYYYY}/genes.ffn` is a concatenation of all `GUT_GENOME{XXXXXX}.ffn` produced by Prokka for the genomes in that `species_id`, with some minor cleanup: 

  * extra newlines have been removed from the gene sequences, so that each gene sequence occupies a single line below the gene header

  * DNA letters have been converted to uppercase, and 

  * degenerate genes (with empty sequences or headers containing "|") have been excluded

The temp files are produced by clustering `genes.ffn` with [vsearch](https://github.com/torognes/vsearch) to 99, 95, ... percent identity (PID), respectively.

The top level `centroids.ffn` file represents the 99 percent identity clusters -- with each cluster represented by the gene at its center.

For every gene from `genes.ffn`, the centroids of the clusters containing that gene are listed in `gene_info.txt`.  For example, in the following excerpt, we see 5 genes that belong to two different 99-PID clusters, but only to a single 95-PID cluster centered around the gene `UHGG000001_00001`.

```
gene_id                centroid_99            centroid_95
UHGG000001_00001       UHGG000001_00001       UHGG000001_00001
UHGG000001_00002       UHGG000001_00002       UHGG000001_00001
UHGG000001_00003       UHGG000001_00002       UHGG000001_00001
UHGG000001_00004       UHGG000001_00002       UHGG000001_00001
UHGG000001_00005       UHGG000001_00002       UHGG000001_00001
```


# Migrate MIDAS DB v1.2

MIDAS DB is driven by a TOC with four columns

```
genome     species                        representative   genome_is_representative
1190605.3  Enterovibrio_norvegicus_54866  1190605.3        1
349741.6   Akkermansia_muciniphila_55290  349741.6         1
1408424.3  Bacillus_bogoriensis_60417     1408424.3        1
```

## Phyeco Marker Genes Database

Need to create `marker_centroids` under the `marker_genes` folder, inside which is the ID mapping from the SGC marker genes of the representative genomes to the centroid_99 genes of the pan-genome.


## Representative Genomes Database

Rename the FASTA file into `{species_id}/{genome_id}/{genome_id}.fna`. Generate the `{genome_id}.ffn` for all the genes per rep genome.

## Pan Genomes Database

We need `centorids.ffn` and `gene_info.txt`.


# See Also

[Operations](https://github.com/czbiohub/iggtools/wiki/Operations)
