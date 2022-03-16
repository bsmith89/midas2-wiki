
# MIDAS DB

At its core, MIDAS DB is an interface between **collections of genomes** and a set of custom files needed for strain-level metagenomic analysis. The new implementation of MIDAS DB reads in a **Table Of Contents (TOC) file**, containing genome-to-species assignment and a choice of representative genome for each species. This new infrastructure simplifies the prior knowledge needed to build a custom MIDAS database, and enables the dynamic assignment of representative genomes. This repository contains tools for building a [MIDAS DB](https://github.com/snayfach/MIDAS) from any collection of genomes, geared to run on [AWS Batch](https://aws.amazon.com/batch/) as described in [PairANI](https://github.com/czbiohub/pairani/wiki).  


Presently, MIDAS DB contains two public databases readily available for users: UHGG v1 (4644 species with 286997 genomes) and GTDB v202 (47893 species with 258405 genomes) in a public AWS S3 bucket.  MIDAS is the first pipeline to allow users the flexibility to select representative genomes, a key component for the accuracy of SNP calling.


## UHGG

A collection of 286,997 genomes assembled from metagenomes, isolates and single cells from human stool samples has been clustered into 4,644 species in an effort similar to [IGGdb 1.0](https://github.com/snayfach/IGGdb).   We refer to this new collection as [UHGG](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/).  Perhaps the most important difference with respect to the IGGdb 1.0 is that UHGG contains only gut genomes.

The collection of UHGG genomes were provided by [Alexandre Almeida](https://www.ebi.ac.uk/about/people/alexandre-almeida) of EBI and mirrored in S3 by [Jason Shi](http://docpollard.org/people/jason-shi/), whose `s3://jason.shi-bucket/IGGdb2.0/clean_set/` serves as input to the tools in this repository. Six-digit numeric species ids were arbitrarily assigned by Jason Shi in `s3://jason.shi-bucket/IGGdb2.0/alt_species_ids.tsv`. 

## GTDB

ADD HERE.


# Target Layout in S3

Take UHGG as an example, the resulting MIDAS DB is located at: `uhgg`=[s3://microbiome-igg/2.0/](http://microbiome-igg.s3.amazonaws.com/2.0/README.TXT)


## Table of Contents

The `s3://microbiome-igg/2.0/genomes.tsv` listing all genomes has four columns, with per genome each row:

```
genome              species      representative        genome_is_representative
GUT_GENOME138501    104351       GUT_GENOME269084      0
GUT_GENOME269084    104351       GUT_GENOME269084      1
...
```

## Collections of Genomes

All input genome fasta files would be relocated under `cleaned_genomes/<species>/<genome>.fa.lz4`

For instance, to match the earlier example,

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

Genome annotations for each genome were identified by [Prokka](https://github.com/tseemann/prokka), and the annotated genes were kept under the directory of `genes_annotation/<species>/<genome>`

```
s3://microbiome-igg/2.0/gene_annotations/104351/GUT_GENOME138501/GUT_GENOME138501.{faa, ffn, fna, gff, tsv}.lz4
```

## Single Copy Marker Genes Identification

Homologs of a specified single copy marker gene set are identified for each input genome, cross-referenced, collated, and indexed with `hs-blastn` to use in species abundance estimation for metagenomic samples.   This involves the following steps.

The HMM model 15 Phyeco single copy marker genes are available at:
```
s3://microbiome-igg/2.0/marker_gene_models/phyeco/marker_genes.hmm.lz4
s3://microbiome-igg/2.0/marker_gene_models/phyeco/marker_genes.mapping_cutoffs.lz4

```

For each annotated genome, single copy marker genes (`{YYYYYY}.markers.fa`) were identified with `hmmsearch`, as well as the mapping of Prokka gene ID to corresponding marker gene ID (`{YYYYYY}.markers.map`), under `marker_genes/phyeco/temp/<species>/<genome>/<genome>`

```
s3://microbiome-igg/2.0/marker_genes/phyeco/temp/104351/GUT_GENOME138501/GUT_GENOME138501.{hmmsearch, markers.fa, markers.map}
```

The identified marker genes for all the representative genomes were concatenated into monolithic `marker_genes.fa`, from which hs-blastn index is constructed.

```
s3://microbiome-igg/2.0/marker_genes/phyeco/marker_genes.fa
s3://microbiome-igg/2.0/marker_genes/phyeco/marker_genes.fa.{sa, bwt, sequence}
```

## Pan-Genomes

Species-level pan-genes refer to the set of non-redundant set of genes that represent the genetic diversity within one species cluster. 

To start with, for each <species>, all the annotated genes from its all genome members were concatenated into `pangenomes/<species>/genes.ffn`, with some minor cleanup: 

  * extra newlines have been removed from the gene sequences, so that each gene sequence occupies a single line below the gene header

  * DNA letters have been converted to uppercase, and 

  * degenerate genes (with empty sequences or headers containing "|") have been excluded

```
s3://microbiome-igg/2.0/pangenomes/104351/{genes.ffn}
```


Second, the concatenated genes were clustered based on 99% percent identity (PID) using [vsearch](https://github.com/torognes/vsearch). A centroid gene (`centroids.99.ffn`) was selected to represent each cluster. The centroid.99 were further on clustered to 95, 90, ... PID, respectively.  The top level `centroids.ffn` file represents the 99 percent identity clusters -- with each cluster represented by the gene at its center.


```
s3://microbiome-igg/2.0/pangenomes/104351/{centroids.ffn}
s3://microbiome-igg/2.0/pangenomes/104351/temp/centroids.{99, 95, 90, 85, 80, 75}.ffn
s3://microbiome-igg/2.0/pangenomes/104351/temp/uclust.{99, 95, 90, 85, 80, 75}.txt
```

For every gene from `genes.ffn`, the centroids of the clusters containing that gene are listed in `gene_info.txt`.  

```
s3://microbiome-igg/2.0/pangenomes/104351/{gene_info.txt}
```

For example, in the following excerpt, we see 5 genes that belong to two different 99-PID clusters, but only to a single 95-PID cluster centered around the gene `UHGG000001_00001`.

```
gene_id                centroid_99            centroid_95
UHGG000001_00001       UHGG000001_00001       UHGG000001_00001
UHGG000001_00002       UHGG000001_00002       UHGG000001_00001
UHGG000001_00003       UHGG000001_00002       UHGG000001_00001
UHGG000001_00004       UHGG000001_00002       UHGG000001_00001
UHGG000001_00005       UHGG000001_00002       UHGG000001_00001
```


# Local MIDB DB Layout

MIDAS 2.0 users don't need to pre-download the entire MIDAS DB. Only species needed for the given analysis would be downloaded on the fly, and the local MIDAS DB mirror the S3 layout.

```
Output                      Producer              Meaning
---------------------------------------------------------------------------------
midasdb                     DB-related files      Mirror s3://miocriombe-igg/2.0/
```


# See Also

[Operations](https://github.com/czbiohub/iggtools/wiki/Operations)

