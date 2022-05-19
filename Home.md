# Metagenomic Intra-Species Diversity Analysis, version 2

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
next minor version? MIDAS 2.1? 2.0.2?
-->
MIDAS2 implements the same analyses as the original
(here referred to as [MIDAS1](https://github.com/snayfach/MIDAS)),
but re-engineered to (1) allow for multiple, alternative reference databases
(MIDASDBs), and (2) optimize scaling to collections of thousands of samples.

This documentation includes instructions for [installing software](#installation),
downloading one or more reference databases(#getting-databases),
and running the three standard modules: (species filtering)[#species-filtering-module],
(SNV profiling)[#snv-module], and (CNV profiling)[#cnv-module].
In addition, the philosophy and details of how MIDAS operates are described,
along with instructions for more advanced usage (e.g. (database customization)[#custom-databases]).

# Quickstart

### Basic installation (conda and pip)

These instruction assume that [conda has already been installed](TODO-link-to-conda).

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

### Download test data

<!--
TODO: wget URL
TODO: Download *two* example files so that we can run a full workflow.
-->
```
mkdir -p midas_example
cd midas_example
wget TODO
tar -zxf TODO  # Inflates into a reads/ directory with R1/R2 for each of two samples.
```

### Pre-download marker genes

```
# TODO: midas2 download_db  ... # Download the MIDASDB-UHGG to ./midas_db/uhgg
```

### Identify abundant species

This is designed only to select high-coverage species in each sample.
IMPORTANT: It is not meant to quantify species abundance.

<!--
TODO: It would be nice if there was a file in midasdb_dir that told it the
name of the database. Why does it even need to know the name? Couldn't I
just point it to an arbitrary path with the database?
TODO: I dropped num_cores to 2 from 8. A personal laptop is less likely to have
8 cores...
-->
```
for sample_name in sample1 sample2
do
    midas2 run_species --sample_name ${sample_name} -1 reads/${sample_name}_R1.fastq.gz \
        --midasdb_name uhgg --midasdb_dir ./midas_db/uhgg --num_cores 2 midas_output
done

```

### Run SNVs module

<!--
TODO: Can we drop the extra flags (e.g. chunksize, select_by, etc)?
Are there defaults that can be used for this quickstart?
-->
```
for sample_name in sample1 sample2
do
    midas2 run_snps --sample_name ${sample_name} -1 reads/${sample_name}_R1.fastq.gz \
        --midasdb_name uhgg --midasdb_dir midas_db/uhgg --chunk_size 500000 \
        --select_by median_marker_coverage,unique_fraction_covered --select_threshold=5,0.5 \
        --num_cores 2 midas_output
done
# TODO: midas2 merge_snps
```

The SNV results can be found in TODO,
summarizing the number of sites that have been identified
as polymorphic within and across samples for each species found to be at
sufficient coverage.

### Run CNVs module

<!--
TODO: Can we drop the extra flags (e.g. chunksize, select_by, etc)?
Are there defaults that can be used for this quickstart?
-->
```
for sample_name in sample1 sample2
do
    midas2 run_genes --sample_name ${sample_name} -1 reads/${sample_name}_R1.fastq.gz \
        --midasdb_name uhgg --midasdb_dir midas_db/uhgg \
        --select_by median_marker_coverage,unique_fraction_covered --select_threshold=5,0.5 \
        --num_cores 2 midas_output
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

![Figure 1: MIDAS 2.0 Analysis Modules](https://github.com/czbiohub/MIDAS2.0/blob/master/docs/images/Fig.Modules.png?raw=true "Title")

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

TODO: What the difference between the two high-level DBs?
TODO: DBs are downloaded when needed
TODO: Users _can_ make custom DBs: link to that section
TODO: The DB includes different parts for the different modules.
TODO: DBs can be downloaded and shared across projects

# Overview: Maximizing Performance

TODO: running samples in parallel
TODO: when to use extra threads in the same process?
TODO: pre-downloading databases
TODO: how to choose a chunksize
TODO: mention prebuilt bowtie2 indexes
TODO: etc.

# Overview: The MIDAS Interface

## Common CLI options:

TODO: output dir, sample name, database dir/name, species selection, etc.

## Key input file locations/formats

TODO

## Key output file locations/formats

TODO

# Module: Species Selection

TODO: Any special details for run_species, merge_species? Necessary flags?
TODO: How/when to preload the marker genes database
TODO: Mention that species selection itself happens downstream in the SNV and
CNV modules.

# Module: Single-Nucleotide Variant Analysis

TODO

# Module: Copy-Number Variant Analysis

TODO

# Advanced: Installation

TODO

# Advanced: Database Management

## Downloading Databases

TODO

## Building Your Own MIDASDB

TODO

# Advanced: MIDAS2 Implementation

## Output Directory Layout

TODO

## MIDASDB Directory Layout

TODO

# Advanced: How to Contribute

TODO

# End Matter: Glossary

TODO

# End Matter: Release Notes

TODO
