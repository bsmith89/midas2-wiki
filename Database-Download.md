## List MIDAS DB

MIDAS Reference Database (MIDAS DB) is comprised of species pangenomes, representative genomes, marker genes. We pre-constructed two MIDAS DBs from two public microbial genome collections and available for users to download:

```
$ midas2 database --list
uhgg 286997 genomes from 4644 species version 1.0
gtdb 258405 genomes from 47893 species version r202
```

## Initiate MIDAS DB
 
To download MIDAS DB, users need to specify two command-line options:
```
--midasdb_dir: /path/to/the/local/MIDASDB
--midasdb_name: the name prokaryotic genome database {uhgg, gtdb}
```

First time users need to initialize a local MIDAS DB with specified _midasdb_name_, which will downloaded the minimal information needed to run the database customization step.

```
midas2 database --init --midasdb_name gtdb --midasdb_dir ${midasdb_dir}
```

## Download MIDAS DB

MIDAS DB UHGG requires 93 GB of free space to download and decompress, and MIDAS DB GTDB requires 539 GB of free space to download and decompress. Users can download the entire MIDAS DB GTDB to local directory `${my_midasdb_gtdb}`:

```
midas2 database --download --midasdb_name gtdb --midasdb_dir ${my_midasdb_gtdb} -s all
```

## Slice Download MIDASDB 

Since only sufficiently abundant species can be accurately metagenotyped, the scale of genotype-able species in a deeply sequenced shotgun metagenomics sample would be around 100 to 150 species. Therefore, MIDAS 2.0 provide slice-download the reference database by species. 

Given a list of species that we want to perform strain-level analysis, we can save these species into a TSV file `my_species_list=/path/to/species/tsv`, as following:

```
102506
122337
```

Then before starting the SNV/CNV module analysis, we can download the MIDAS DB only for these species:

```
midas2 database --download --midasdb_name gtdb --midasdb_dir ${my_midasdb_gtdb} --species_list ${my_species_list}
```

For example, for users who have followed the across-samples abundant species detection step, users can generated the list of species to genotype from the `species_prevalence.tsv`:

```
# Generate the list of species that is present in more than one sample
awk '$6 > 1 {print $6}' ${midas_outdir}/species/`species_prevalence.tsv` > my_species_list.tsv

# Download MIDAS DB GTDB only for the species in my list
midas2 database --download --midasdb_name gtdb --midasdb_dir ${my_midasdb_gtdb} --species_list my_species_list.tsv
```
