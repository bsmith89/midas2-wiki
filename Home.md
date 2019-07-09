# Microbiome - Integrated Gut Genome Tools

A collection of 286,997 metagenome-assembled genomes (MAGs) from human gut samples has been clustered into 4,644 species in an effort similar to [IGGdb 1.0](https://github.com/snayfach/IGGdb).   We refer to this new collection as IGGdb 2.0, IGG 2.0, IGG+, or simply IGG.

This repository contains tools for building a [MIDAS](https://github.com/snayfach/MIDAS) database from IGG 2.0, geared to run on [AWS Batch](https://aws.amazon.com/batch/) as described in [PairANI](https://github.com/czbiohub/pairani/wiki).  The resulting database will be located at [s3://microbiome-igg/2.0](http://microbiome-igg.s3.amazonaws.com/2.0/README.TXT).

Our project to update [MIDAS for IGGdb 1.0](https://github.com/czbiohub/MIDAS-IGGdb/blob/master/README.md) is on hold.  We will first update MIDAS to run with this new database, directly from S3.

# Data layout in S3

A table of contents listing all MAGs will be located at [s3://microbiome-igg/2.0/genomes.tsv](http://microbiome-igg.s3.amazonaws.com/2.0/genomes.tsv).  It will look like
```
genome              species      representative        genome_is_representative
GUT_GENOME138501    104351       GUT_GENOME269084      0
GUT_GENOME269084    104351       GUT_GENOME269084      1
...
```
with 4644 rows below the headers.


Column `genome_is_representative` is defined as `genome == representative`.

Every species has a single representative genome.  This representative genome could be used to identify the species.  While doing so would be elegant, for our purposes here, purely numeric species ids are more convenient.

The 6-digit species ids above were arbitrarily assigned by Jason Shi; see `s3://jason.shi-bucket/IGGdb2.0/alt_species_ids.tsv`.
