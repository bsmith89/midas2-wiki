# MIDAS2: Metagenomic Intra-Species Diversity Analysis

Metagenomic Intra-Species Diversity Analysis
([MIDAS](https://genome.cshlp.org/content/26/11/1612)) is an integrated set of
workflows for
**profiling strain-level genomic variations in shotgun metagenomic data**.
Specifically, MIDAS is designed to profile strain genotypes in two
complementary ways:

- single-nucleotide variants (SNVs) across polymorphic sites in a species reference
  genome
- gene copy number variation (CNVs) across the full species pangenome

<!--
TODO: Tried to rename MIDAS 2.0 as MIDAS2 everywhere.
If it's MIDAS 2.0, then what's the
next minor version? MIDAS 2.1? 2.0.2? MIDAS2 is also cleaner to my eye and
matches the program name
-->
MIDAS2 implements the same analyses as the original
(here referred to as [MIDAS1](https://github.com/snayfach/MIDAS)),
but re-engineered to (1) allow for multiple, alternative reference databases
(MIDASDBs), and (2) optimize scaling to collections of thousands of samples.

This documentation includes instructions for [installing software](#advanced-installation),
[downloading one or more reference databases](#module-database-management),
and running the three standard modules: [species selection](#module-species-selection),
[SNV profiling](#module-single-nucleotide-variant-analysis), and
[CNV profiling](#module-copy-number-variant-analysis).
In addition, the philosophy and details of how MIDAS operates are described,
along with instructions for more advanced usage
(e.g. [building custom database](#advanced-building-your-own-midasdb)).

# Quickstart

1. **Install MIDAS2.**

(NOTE: These instruction assume that [conda has already been installed](TODO-link-to-conda).)

<!--
# TODO: Where does midas2.yml come from?
# TODO: Is this really the best "quick-start" installation if it requires
# tricky stuff like cpanm?
# TODO: Why are we using pip, but just for midas2? Couldn't we have used conda?
-->
```
# Create a new conda environment with required software.
conda env create -n midas2 -f midas2.yml
# Temporary fix for Prokka:
cpanm Bio::SearchIO::hmmer --force
pip install midas2
```

2. **Download example data.**

<!--
TODO: wget URL
TODO: Download *two* example files so that we can run a full workflow.
-->
```
mkdir -p midas_example
cd midas_example
wget TODO
tar -zxf TODO
# Inflates into a reads/ directory with R1/R2 for each of two samples.
```

3. **Pre-download marker genes.**

```
# TODO: midas2 download_db  ... # Download the MIDASDB-UHGG to ./midas_db/uhgg
```

4. **Identify abundant species.**

(NOTE: This is designed only to select high-coverage species in each sample.
_It is not intended to quantify species abundance._)

<!--
TODO: It would be nice if there was a file in midasdb_dir that told it the
name of the database. Why does it even need to know the name? Couldn't I
just point it to an arbitrary path with the database?
TODO: I dropped num_cores to 2 from 8. A personal laptop is less likely to have
8 cores...
TODO: Any reason it's just R1 and not R2 as well?
-->
```
for sample_name in sample1 sample2
do
    midas2 run_species \
    --sample_name ${sample_name} \
    -1 reads/${sample_name}_R1.fastq.gz \
        --midasdb_name uhgg \
        --midasdb_dir ./midas_db/uhgg \
        --num_cores 2 \
        midas_output
done

```

5. **Run SNVs module.**

<!--
TODO: Can we drop the extra flags (e.g. chunksize, select_by, etc)?
Are there defaults that can be used for this quickstart?
-->
```
for sample_name in sample1 sample2
do
    midas2 run_snps \
        --sample_name ${sample_name} \
        -1 reads/${sample_name}_R1.fastq.gz \
        --midasdb_name uhgg \
        --midasdb_dir midas_db/uhgg \
        --chunk_size 500000 \
        --select_by median_marker_coverage,unique_fraction_covered \
        --select_threshold=5,0.5 \
        --num_cores 2 \
        midas_output
done
# TODO: midas2 merge_snps
```

The SNV results can be found in TODO,
summarizing the number of sites that have been identified
as polymorphic within and across samples for each species found to be at
sufficient coverage.

5. **Run CNVs module.**

<!--
TODO: Can we drop the extra flags (e.g. chunksize, select_by, etc)?
Are there defaults that can be used for this quickstart?
-->
```
for sample_name in sample1 sample2
do
    midas2 run_genes \
        --sample_name ${sample_name} \
        -1 reads/${sample_name}_R1.fastq.gz \
        --midasdb_name uhgg \
        --midasdb_dir midas_db/uhgg \
        --select_by median_marker_coverage,unique_fraction_covered \
        --select_threshold=5,0.5 \
        --num_cores 2 \
        midas_output
done
# TODO: midas2 merge_genes
```

The CNV results can be found in TODO,
summarizing TODO within and across samples for each species found to be at
sufficient coverage.

# Overview: The MIDAS Workflow

MIDAS contains two strain-level reads-to-table analysis modules: population
SNVs analysis (SNV module) and pan-genome CNVs analysis (CNV module).  Each
module includes two sequential steps: single-sample analysis and across-samples
analysis.

![Figure 1: MIDAS2 Analysis Modules](https://github.com/czbiohub/MIDAS2.0/blob/master/docs/images/Fig.Modules.png?raw=true "Title")

Before running these modules, however, the MIDAS workflow starts by identifying
species at high coverage in each sample (species module).
The results of this analysis are then be used to filter the databases used
in subsequent steps.
The goal is to optimize the Bowtie2 indexes used for each analysis
in order to maximize specificity and sensitivity.
This is important because fundamentally there is a trade-off between:

- "cannibalization" of reads by reference genomes with shared sequence
- incorrect mapping of reads to genomes of interest from different,
  closely related organisms found in the same

# Overview: The MIDASDBs

- TODO: What the difference between the two high-level DBs?
- TODO: DBs are downloaded when needed
- TODO: Users _can_ make custom DBs: link to that section
- TODO: The DB includes different parts for the different modules.
- TODO: DBs can be downloaded and shared across projects

# Overview: Maximizing Performance

- TODO: running samples in parallel
- TODO: when to use extra threads in the same process?
- TODO: pre-downloading databases
- TODO: how to choose a chunksize
- TODO: mention prebuilt bowtie2 indexes
- TODO: etc.

# Overview: Similar Tools

TODO

# Overview: The MIDAS Interface

## Common CLI options

TODO: output dir, sample name, database dir/name, species selection, etc.

## Input files

TODO

## Output files

TODO

# Module: Database Management

A MIDASDB is comprised of species pangenomes, representative genomes, and
marker genes.
We have already constructed and made available for download MIDASDBs for two
comprehensive public microbial genome collections.

TODO: Link to the websites for each of these.

```
$ midas2 database --list
uhgg 286997 genomes from 4644 species version 1.0
gtdb 258405 genomes from 47893 species version r202
```

Most MIDAS commands take as an argument a path to a local mirror of the MIDASDB,
as well as the name of the database.

```
--midasdb_dir: /path/to/the/local/MIDASDB
--midasdb_name: the name prokaryotic genome database (either 'uhgg' or 'gtdb')
```

For the purposes of this documentation we'll generally assume that we're working
with the prebuilt `gtdb` MIDASDB and that the local mirror is in a subdirectory
of the project root: `./midasdb`.
However, in practice, users may prefer to have a centralized
MIDASDB so that multiple projects and users can use previously cached files.

## Preloading species marker genes

We highly recommend that first time users initialize a local MIDASDB
This downloads the minimal information needed to run the species module
needed for species selection by the SNV and CNV modules.

```
midas2 database --init --midasdb_name gtdb --midasdb_dir ./midasdb
```

## Preloading all species

While it _is_ possible to download an entire MIDASDB using the following
command:

```
midas2 database --download --midasdb_name gtdb --midasdb_dir ./midasdb -s all
```

this requires a large amount of data transfer and storage:
93 GB for `uhgg` and 539 GB for `gtdb`.

## Preloading select species

Instead, we recommend that users take advantage of species-level database
sharding to download and decompress only the necessary portions of a
MIDASDB.

This can be done on-demand based on various species-filtering criteria (see TODO)
passed to the single-sample analyses of both the SNV and CNV modules.
Unfortunately, this does not play nicely with parallelizing computation
across samples, because multiple `run_snps` or `run_genes` processes
might simultaneously attempt to download the same species data.
Instead, we recommend that users consider pre-selecting a list of species
for download.

If we save the following list of species IDs (here an example with only two
species) to a plain text file named `./species.list`:

```
102506
122337
```

we can then run the following to preload all of the data needed for these species:

```
midas2 database \
    --download --midasdb_name gtdb --midasdb_dir ./midasdb \
    --species_list ./species.list
```

Advanced users may be interested in using the outputs of the `merge_species`
command used in the [species selection module](#module-species-selection)
to generate this list of species, e.g.:

```
# Generate the list of species that is present in more than one sample
# TODO: What is in column $6? (be explicit)
awk '$6 > 1 {print $6}' \
    < ./midas_output/species/species_prevalence.tsv \
    > ./species.tsv
```

Now, when running `run_snps` or `run_genes` for multiple samples in parallel
(e.g. using a workflow manager or background processes in UNIX),
MIDAS2 will see that the necessary parts of the MIDASDB have already been
downloaded and will not repeat that step.

It is also possible for users to
[construct their own MIDASDB](#advanced-building-your-own-midasdb)
from a custom genome collection (e.g. for metagenome assembled genomes).

# Module: Species Selection

- TODO: Any special details for run_species, merge_species? Necessary flags?
- TODO: How/when to preload the marker genes database
- TODO: Mention that species selection itself happens downstream in the SNV and
  CNV modules.

# Module: Single-Nucleotide Variant Analysis

TODO

# Module: Copy-Number Variant Analysis

TODO

# Advanced: Installation

### Conda

```
conda config --set channel_priority flexible
conda install -c zhaoc1 -c anaconda -c bioconda -c conda-forge -c defaults midas2
```

### Developer Installation

<!--
TODO: I updated this to show where the midas2.yml file came from
as a bonus, this is how I'd suggest developers install it.
Makes the repo editable.
-->
```
git clone https://github.com/czbiohub/MIDAS2.0.git midas2
cd midas2
conda env create -n midas2 -f midas2.yml
cpanm Bio::SearchIO::hmmer --force # Temporary fix for Prokka
# TODO: Temporary? When will it no longer be needed?
pip install .  # Include -e flag for editable installation/development
```

### Docker Image

```
docker pull zhaoc1/midas2:latest
docker run --volume "/home/ubuntu/.aws":"/root/.aws":ro --rm -it midas2:latest
# TODO: Is the above correct? What about binding the DBs and stuff?
```

### Singularity Image

```
# TODO
```

# Advanced: Building Your Own MIDASDB

MIDAS2 users can locally build a new MIDAS DB for a custom collection of
genomes. The target layout of MIDAS DB can refer to [this page](TODO).  This
page is focused specifically on the database construction commands.

To start with, users need to organize the genomes in a specific format and
produce the TOC `genomes.tsv` according to the
[Inputs](https://github.com/czbiohub/MIDAS2.0/wiki/4.-MIDAS-Database#inputs)

We have prepared a toy collections of genomes.

```
git clone https://github.com/czbiohub/MIDAS2.0.git midas2
cd midas2/tests/genomes
```

We will build the new MIDASDB under the directory of `midas2/tests/genomes`.
There are two command-line parameters that we need to pass:

```
--debug: keep the local file after successfully build the database
--force: re-build the database even if already locally exists
```

First, annotate all the genomes.

<!--
TODO: Wrap bash lines
-->
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

# Advanced: MIDAS2 Implementation

## Output Directory Layout

TODO

## MIDASDB Directory Layout

TODO

# Advanced: Contributing / Development Notes

All development is done with git and GitHub

## Installation

Use the [developer installation](#developer-installation) as described;
be sure to run pip with the `-e` or `--editable` flag.

## Run Integration Tests

```
bash tests/midas2_analysis.sh 8
```

## Export Your Conda Environment

<!--
TODO: Moved this from the installation page because it wasn't clear
why users would want to do this.
TODO: Maybe this belongs in the developer notes?
-->

```
conda update --all
conda clean –all
conda env export --no-builds | grep -v "^prefix:" > midas2.yml
```

# End Matter: Glossary

- *Module*: i.e. one of Species, SNPs, Genes (or maybe also MIDAS DB building?)
- *Analysis*: Either the single-sample alignment-based tallying of reads
  belonging to species, SNVs, or genes (depending on which module) OR the
  merging, filtering, and summarization of these results across samples.
- *Workflow*: The overall conceptual order of events from an implicit database
  and pre-requisite shotgun metagenomic samples through species selection and
  SNPs/genes modules, to results that will go into some sort of downstream
  analysis.
- *Reference database*: The upstream UHGG, GTDB, something user supplies, etc.
  that is then pre-processed to make a specific...
- *MIDAS DB*: The pre-processed reference data in a format to be consumed by
  the MIDAS modules, including marker genes, reference genomes, pangenomes,
  etc.
- *Genome collections*: TODO
- *Bowtie2 Index*: Rather than bowtie2 database or some other ambiguous term
- Species-level pangenome refers to the set of non-redundant genes (centroids)
  clustered from all the genomes within one species cluster.

# End Matter: Release Notes

### 0.8

Rename IGGtools to MIDAS2.0.

### 0.7

IGGtools used for simulation and benchmark.

### 0.6

Ability to run subcommands under AWS Batch with status tracking in s3
operations bucket.

### 0.5

Implement `build_pangenomes` subcommand.
[[d0a6f9ede]](https://github.com/czbiohub/iggtools/commit/d0a6f9ede)
