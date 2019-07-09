# Microbiome - Integrated Gut Genome Tools

A collection of 286,997 metagenome-assembled genomes (MAGs) from human gut samples has been clustered into 4,644 species in an effort similar to [IGGdb 1.0](https://github.com/snayfach/IGGdb).   We refer to this new collection as IGGdb 2.0, IGG 2.0, IGG+, or simply IGG.

This repository contains tools for building a [MIDAS](https://github.com/snayfach/MIDAS) database from IGG 2.0, geared to run on [AWS Batch](https://aws.amazon.com/batch/) as described in [PairANI](https://github.com/czbiohub/pairani/wiki).  The database will be located at [s3://microbiome-igg/2.0](http://microbiome-igg.s3.amazonaws.com/2.0/README.TXT)

Our project to update [MIDAS for IGGdb 1.0](https://github.com/czbiohub/MIDAS-IGGdb/blob/master/README.md) is on hold until this new database has been built.
