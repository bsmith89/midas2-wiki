
## Custome Collection of Genomes

The new infrastructure significantly simplified the prior information needed to build a MIDAS Reference Database from any collections of genomes. Users needs to organize the genomes and produce the `genomes.tsv` according to the Inputs section.  

```
#TODO: refer to unit database testing
```

# Database Download

There are two ways to download above-mentioned, pre-computed MIDASDB-UHGG and MIDASDB-GTDB: centralized download and on-demand download. Locally downloaded MIDASDB mirror the target layout.

## Centralized Download

The default MIDASDBs can be downloaded as following:

- MIDASDB-UHGG: [download midasdb_uhgg.tar.gz](toadd link). The entire databases requires 93 GB of free space to download and decompress.
- MIDASDB-GTDB: [download midasdb_gtdb.tar.gz](toadd link). The entire databases requires 539 GB of free space to download and decompress.

Once download and unpack the tarball (`tar -zxvf midasdb_uhgg.tar.gz`), users need to specify two MIDASDB related input parameters to midas2 command:

- `--midasdb_dir`: the absolute path of the local MIDASDB.
- `--midasdb_name`: the name prokaryotic genome database.

For example: `midas2 run_snps outdir --midasdb_name uhgg --midasdb_dir midasdb_uhgg [options]`.

## On Demand Download

The above mentioned two MIDASDBs are also being hosted on a [public S3 bucket](s3://microbiome-pollardlab/). The default MIDAS reference databases contain the representative genomes (rep-genome) and pan-genome databases for **all** the species. Since only species with sufficiently high reads coverage would be genotyped, the scale of genotype-able species even in the deeply sequenced shotgun metagenomics samples would be around 100 to 150 species. Therefore, MIDAS 2.0 can **slice** the reference database efficiently based on the species needed for the current analysis, and furthermore only **download** the necessary files during the analysis. 

The requirements for the on-demand downloading is an AWS account with access to S3. In addition, the advanced features of MIDAS 2.0 that re-assign representative genomes for the single-sample SNPs module requests the users to on-demand download the rep-genome databases of newly assigned representative genomes.

Alternatively, users can slice-download the database files of `species_id=100001` from MIDASDB-UHGG to local path `single_midasdb_uhgg`:

```
$ midas2 database --download -s 100001 --midasdb_name uhgg --midasdb_dir $single_midasdb_uhgg
```

Another example is users can slice-download the database files of `species_id=100100` from MIDASDB-GTDB to local path `single_midasdb_gtdb`:

```
$ midas2 database --download -s 100001 --midasdb_name gtdb --midasdb_dir single_midasdb_gtdb
```