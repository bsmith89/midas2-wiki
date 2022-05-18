
MIDAS 2.0 users can locally build a new MIDAS DB for a custom collection of genomes. The target layout of MIDAS DB can refer to [this page]().  This page is focused specifically on the database construction commands. 

To start with, users need to organize the genomes in a specific format and produce the TOC `genomes.tsv` according to the [Inputs](https://github.com/czbiohub/MIDAS2.0/wiki/4.-MIDAS-Database#inputs)

We have prepared a toy collections of genomes.

```
git clone https://github.com/czbiohub/MIDAS2.0.git .
cd MIDAS2.0/tests/genomes
```

We will build the new MIDASDB under the directory of `MIDAS2.0/tests/genomes`. There are two command-line parameters that we need to pass:

```
--debug: keep the local file after successfully build the database
--force: re-build the database even if already locally exists
```

First, annotate all the genomes.

```
midas2 annotate_genome --species all --midasdb_name newdb --midasdb_dir my_new_midasdb --debug --force
midas2 build_midasdb --generate_gene_feature --genomes all --midasdb_name newdb --midasdb_dir my_new_midasdb --debug --force
```

Second, infer SCGs for all the genomes and build marker database
```
midas2 infer_markers --genomes all --midasdb_name newdb--midasdb_dir my_new_midasdb --debug --force
midas2 build_midasdb --build_markerdb --midasdb_name newdb--midasdb_dir my_new_midasdb --debug --force

```

Third, build species pangenomes
```
midas2 build_pangenome --species all --midasdb_name newdb --midasdb_dir my_new_midasdb --debug --force
midas2 build_midasdb --generate_cluster_info --species all --midasdb_name newdb --midasdb_dir my_new_midasdb --debug --force
```




