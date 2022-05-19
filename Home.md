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
next minor version? MIDAS 2.1? MIDAS 2.0 v0.1, 2.0.2? MIDAS2 is also cleaner to my eye and
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
TODO: Revise paths throughout to point to this midas_example root directory
with e.g. midas_example/midasdb as the uhgg download location,
midas_example/reads/sample1_R1.fastq.gz, and midas_example/midas_output/
as the output directory in all example code.
Users would also implicitly already be in midas_example (in fact we'd only
reference it once in the quickstart) and everything else would work
with relative paths from that project root.
-->
```
mkdir -p midas_example
cd midas_example
wget TODO
tar -zxf TODO
# Inflates into a reads/ directory with R1/R2 for each of two samples: sample1 and sample2.
```

3. **Pre-download marker genes.**

```
midas2 database --init --midasdb_name uhgg --midasdb_dir midasdb/uhgg
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
        --midasdb_dir midas_db \
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

![
MIDAS2 Analysis Modules
](static/Fig.Modules.png)

Before running these modules, however, the MIDAS workflow starts by identifying
species at high coverage in each sample (species module).
The results of this analysis are then be used to filter the databases used
in subsequent steps.
The goal is to optimize the Bowtie2 indexes used for each analysis
in order to maximize specificity and sensitivity.
This is important because fundamentally there is a trade-off between:

- "cannibalization" of reads by reference genomes with shared sequence
- incorrect mapping of reads to genomes of interest from different,
  closely related organisms found in the same sample.

# Overview: The MIDASDBs

- TODO: What the difference between the two high-level DBs?
- TODO: DBs are downloaded when needed
- TODO: Users _can_ make custom DBs: link to that section
- TODO: The DB includes different parts for the different modules.
- TODO: DBs can be downloaded and shared across projects

MIDAS2 is a reference-based strain-level genomic variation analysis
pipeline, and it also presuppose a reference database construction step has
already taken place.  "MIDAS reference database (MIDASDB)" refers to a set of
custom files needed for the strain-level metagenomic analysis.

MIDAS1 provided a default bacterial reference databases (see
[Figure 1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5088602/)).
[MIDASDB v1.2](http://lighthouse.ucsf.edu/MIDAS/midas_db_v1.2.tar.gz)
was constructed from a collection of 5952 bacterial species clusters
representing 31,007 high-quality bacterial genome.

In the past few years, the number of sequenced microbial genomes have increased
vastly, in particular with the addition of metagenome-assembled genomes (MAGs)
sequenced from varied habitats. Therefore, it is necessary to update MIDAS
reference database accordingly. On the other hand, processing the large amount
of available genome sequences poses a significant compute challenge. For MIDAS2,
instead of generating the species clusters from scratch, we take advantage
of several published collections of prokaryotic genome databases, and build
MIDAS reference databases individually.


# Overview: Maximizing Performance

- TODO: running samples in parallel
- TODO: when to use extra threads in the same process?
- TODO: pre-downloading databases
- TODO: how to choose a chunksize
- TODO: mention prebuilt bowtie2 indexes
- TODO: etc.

# Overview: Similar Tools

TODO: Link out to e.g. StrainPhlan, GT-Pro, maybe explain briefly what the pros/cons
of other tools are.

# Overview: The MIDAS Interface

## Common CLI options

TODO: output dir, sample name, database dir/name, species selection, etc.

## Output

The three single-sample commands (`run_species`, `run_snps` and `run_genes`), and
share the following
command-line options.

- User-specified root output directory: `midas_outdir=/path/to/root/outdir`.
  This is always passed as a mandatory positional argument to each of the
  MIDAS2 command.

- Unique sample name: `sample_name=/unique/sample/identifier`.

Together, `${midas_outdir}/${sample_name}` constitutes the unique output
directory for single-sample analysis.


## Input

### Single-Sample Analysis

The FASTA/FASTQ file containing single-end or paired-ends sequencing reads:

```
R1=/path/to/forward/reads/R1.fq.gz
R2=/path/to/reverse/reads/R2.fq.gz
```

`${R1}` and/or `${R2}` need to be passed to MIDAS2 analysis commands via
arguments `-1` and `-2` as: `-1 ${R1} -2 ${R2}`

### Across-Samples Analysis

A TSV file, which lists the _sample_name_ and single-sample root output
directory _midas_outdir_, is required for across-samples analyses
in the species, SNV, and CNV modules.
Users
need to pass the local path of the TSV file
(`my_sample_list=/path/to/tsv/file`) to the command-line argument
`sample_list`, e.g. `--sample_list ${my_sample_list}`

A template _sample_list_ is shown here:

```
sample_name       midas_outdir
SRR172902         /home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample
SRR172903         /home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample
```

## MIDAS DB

For all MIDAS2 analysis, users need to choose (1) a valid precomputed MIDAS DB
name (uhgg, gtdb): `my_midasdb_name=uhgg`, and (2) a valid path to local MIDAS
DB: `my_midasdb_dir=/path/to/local/midasdb/uhgg`.
<!--
TODO: I think I've changed my mind and I don't like the shell variable version.
Can we switch to using an "example dir" that we also reference in the
quickstart and also serves as a tutorial of sorts?
-->

MIDAS2 analysis can take two arguments as: `--midasdb_name
${my_midasdb_name} --midasdb_dir ${my_midasdb_dir}`. If the `--midasdb_dir` is
not specified, MIDAS DB will be downloaded to the current directory.
<!--
TODO: This seems dangerous. I expect some users to forget that one flag and
start an enormous download... :-[
TODO: "Can" take two arguments? What if "--midasdb_name" is not given? Doesn't
it *need* those arguments?
-->


### Others Parameters

Users can set the `--num_cores` to the number of physical cores to use: e.g.
`--num_cores 16`.

And all MIDAS2 analysis can print out the full help message and exit by `-h` or
`--help`.


## Input files

TODO: Fastq format. What preprocessing do these reads need?

## Output files

TODO: What are the key tables from CNV and SNV modules that users may want to
put into downstream analyses?

# Module: Database Management

A MIDASDB is comprised of marker genes, representative genomes, and species pangenomes, representative genomes, and
marker genes.
We have already constructed and made available for download MIDASDBs for two
comprehensive public microbial genome collections.

TODO: Link to the websites for each of these.

```
$ midas2 database --list
uhgg 286997 genomes from 4644 species version 1.0
gtdb 258405 genomes from 47893 species version r202
```

<!--
TODO: This is redundant with the "common cli arguments" material.
I kinda like it better here, though. Maybe we drop the common cli arguments
section in favor of introducing flags when they're actually important.
Alternatively, we could mention them in common args, link to here, and then put
all the details about how to use them in this section.
TODO: Right now there's a lot of redundancy across sections.
-->
Most MIDAS commands take as an argument a path to a local mirror of the MIDASDB,
as well as the name of the database.

```
--midasdb_dir: /path/to/the/local/MIDASDB
--midasdb_name: the name prokaryotic genome database (either 'uhgg' or 'gtdb')
```

For the purposes of this documentation we'll generally assume that we're working
with the prebuilt `uhgg` MIDASDB and that the local mirror is in a subdirectory
of the project root: `midasdb`.
However, in practice, users may prefer to have a centralized
MIDASDB so that multiple projects and users can share cached files.

## Preloading species marker genes

We highly recommend that first time users initialize a local MIDASDB.
This downloads the minimal information needed to run the
[species module](#module-species-selection) needed for species selection by the
[SNV](#module-single-nucleotide-variant-analysis) and
[CNV](#module-copy-number-variant-analysis) modules.
<!--
TODO: Add cross-links throughout.
-->

```
midas2 database --init --midasdb_name uhgg --midasdb_dir midasdb
```

## Preloading all species

While it _is_ possible to download an entire MIDASDB using the following
command:

```
midas2 database --download --midasdb_name uhgg --midasdb_dir midasdb -s all
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
species) to a plain text file named `species.list`:

```
102506
122337
```

we can then run the following to preload all of the data needed for these species:

```
midas2 database \
    --download --midasdb_name uhgg --midasdb_dir midasdb \
    --species_list species.list
```

Advanced users may be interested in using the outputs of the `merge_species`
command used in the [species selection module](#module-species-selection)
to generate this list of species, e.g.:

```
# Generate the list of species that is present in more than one sample
# TODO: What is in column $6? (be explicit)
awk '$6 > 1 {print $6}' \
    < midas_output/species/species_prevalence.tsv \
    > species.tsv
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

<!--
TODO: Citations for the below?
-->
Reference-based metagenotyping depends crucially
on the choice of reference sequences and incorrect mapping of
reads is a major problem.
Microbiome data usually contains
hundreds of species in one sample, and an ideal reference database
is both representative and comprehensive in terms of the
species in the sample.
A badly chosen reference may suffer both from ambiguous mapping
of reads to two or more sequences or spurious mapping to incorrect
sequences.
Therefore, a typical MIDAS2 workflow starts
with a species selection step, which
filters the reference database to species believed to be abundant in each
particular sample based on profiling 15 universal single copy marker genes.


```
midas2 run_species \
    --sample_name ${sample_name} \
    -1 reads/${sample_name}_R1.fastq.gz \
    --midasdb_name uhgg \
    --midasdb_dir midas_db \
    --num_cores 2 \
    midas_output
```

<!--
TODO: Full path to species_profile.tsv
TODO: The details of this output, how it's constructed and the meaning of each
column, should go into the outputs reference section under Advanced topics.
-->

The primary output of this command is `TODO/species_profile.tsv`
which describes the coverage of each species' marker genes in the sample.
Species are sorted in decreasing order of `median_marker_coverage`.
Only species with more than two marker genes covered with more than two reads
(a very low bar) are reported.

```
species_id  marker_read_counts  median_marker_coverage  marker_coverage  marker_relative_abundance   unique_fraction_covered
102337      4110                28.48                   28.91            0.30                        1.00
102506      734                 4.98                    4.98             0.05                        0.93
```

* _marker_read_counts_: total mapped read counts
* _median_marker_coverage_: median coverage of the 15 SCGs
* _marker_coverage_: mean coverage of the 15 SCGs
* _marker_relative_abundance_: computed based on _marker_coverage_
* _unique_fraction_covered_: the fraction of uniquely mapped SCGs genes

<!--
TODO: I think this needs to move up to the CLI interface overview.
This is clearly important material, but I don't think it belongs with
`run_species`; you should introduced e.g. `run_snps` first.
-->

In order to use these marker gene profiles to select species for
index building in the `run_snps` and `run_genes` commands,
users pass flags specifying those parameters:
`--select_by` followed by a comma separated list of column names in
`species_profile.tsv` and
`--select_threshold` followed by a comma-separated list of threshold values for
selection.

For most analyses we recommend using the combination
of `median_marker_coverage > 2X` and `unique_fraction_covered > 0.5`:

```
--select_by median_marker_coverage,unique_fraction_covered --select_threshold=2,0.5
```

Some users may wish to genotype low abundant species
and should adjust the parameters accordingly:

```
--select_by median_marker_coverage,unique_fraction_covered --select_threshold=0,0.5
```

Alternatively, users can directly pick a list of species using the
`--species_list` option
It is worth noting that the species in the provided species
list are still subject to the `--select_threshold` restriction.
Users can set `--select_threshold=-1` to escape species selection filters based
on the SCG species profiling:

```
--species_list 102337,102506 --select_threshold=-1
```

Sample-specific rep-genome and/or pan-genome database would be built only for
the species passing the above mentioned filters in the single-sample SNV or CNV
module.

<!--
TODO: We should link here to instructions for how to build the index just once
and use it for all samples.
-->


## Across-Samples Abundant Species

For some analysis, using one comprehensive list of species across samples in
the same study may be desired, in part because it can avoid the time needed to
build a new index for each sample. If taking this approach, users run
`merge_species` subcommands to merge all the single-sample SCGs profiling
results for all the samples listed in the
`my_samples_list=/path/to/list/of/sample/tsv`

### Sample command

```
midas2 merge_species --samples_list ${my_sample_list} ${midas_outdir}
```

### Output files

The primary output of the across-samples species merging analysis is the file `species_prevalence.tsv`:

```
species_id  median_abundance  mean_abundance  median_coverage  mean_coverage  sample_counts
102337      0.186             0.186           16.205           16.205         2
102506      0.035             0.035           2.967            2.967          2
```

* _median_abundance_: median _marker_relative_abundance_ across samples
* _mean_abundance_: mean _marker_relative_abundance_ across samples
* _median_coverage_: median _median_marker_coverge_ across samples
* _mean_coverage_: mean _median_marker_coverge_ across samples
* _sample_counts_: number of samples with _median_marker_coverge_ > 0


<!--
TODO: The below should be removed from this section. Perhaps moved to the
output directory structure advanced section below.
-->
- Each column in the single-sample `species_profile.tsv` are merged across
  samples into a species-by-sample matrix, shown as following:

  - `species_marker_median_coverage.tsv`: species-by-sample _median_marker_coverge_ matrix

```
species_id   SRR172902   SRR172903
102337       3.926       28.484
102506       0.951       4.983
```

  - `species_unique_fraction_covered.tsv`: species-by-sample _unique_fraction_covered_ matrix

```
species_id   SRR172902   SRR172903
102337       1           1
102506       0.92        1
```

  - `species_marker_coverage.tsv`: species-by-sample _marker_coverage_ matrix:

```
species_id   SRR172902   SRR172903
102337       3.926       28.484
102506       0.951       4.983
```

  - `species_marker_read_counts.tsv`: species-by-sample _marker_read_counts_ matrix:

```
species_id   SRR172902   SRR172903
102337       1565        4110
102506       143         734
```

  - `species_relative_abundance.tsv`: species-by-sample _marker_relative_abundance_ matrix:

```
species_id   SRR172902   SRR172903
102337       0.072       0.301
102506       0.019       0.052
```

Having finished the species-selection step, users can now go to the SNV or CNV
modules, depending on their scientific aims.

# Module: Single-Nucleotide Variant Analysis

The SNV module proceeds in two phages: (1) single-sample read pileup (2)
population variants calling across all the samples. The first step can be run
in parallel.  We presuppose users have already completed the
[species selection module](#module-species-selection) and have
`species_profile.tsv` for each sample.
Alternatively, advanced users can pass a
[prebuilt rep-genome database](TODO) ready for the SNV module.
<!--
TODO: Link to how to make this prebuild rep-genome database.
TODO: Consider moving this suggestion to later in the instructions; it's
an advanced topic.
TODO: "rep-genome" database gets very confusing when we're also talking about
MIDASDBs. Can we change that term to "representative genome index"
or "representative genome bowtie2 index"?
That also works with "pangenome bowtie2 index".
-->

Typically, the `run_snps` command proceeds by

1. selecting species based on taxonomic marker gene profiles;
2. building a sample-specific representative genome bowtie2 index;
3. mapping reads with bowtie2 to this index;
4. outputting read mapping results on a per-species basis.

A subsequent call to `merge_snps` combines results from all
samples and filters polymorphic sites based on their
representation across samples.

MIDAS2 purposely holds any filter or species selection upon the
single-sample pileup results until across-samples SNV analysis.

<!--
TODO: What does the below mean? Is it important to be here at the top of
this section?
-->
That being
said, read mapping summary is reported in `snps_summay.tsv`, and
pileup/variants calling results are organized by species. Therefore users can
easily customize species selection on their own.

A typical call to `run_snps` for one sample is:

```
midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
        --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
        --select_by median_marker_coverage,unique_fraction_covered \
        --select_threshold=2,0.5 \
        --num_cores 12 ${midas_outdir}
```

After running all samples in this way, SNV calling is performed across
samples as follows:

```
midas2 merge_snps \
    --samples_list ${my_sample_list} \
    --midasdb_name ${my_midasdb_name} \
    --midasdb_dir ${my_midasdb_dir} \
    --num_cores 32 \
    ${midas_outdir}
```

Where `${my_sample_list}` is a sample manifest file
[described here](https://github.com/czbiohub/MIDAS2.0/wiki/Common-Command-Line-Arguments#across-samples-analysis).

The key outputs of this command are:

TODO

## Advanced SNV Analyses

<!--
TODO: Ideally, the start of each module is a quick summary
of what th module does, followed by
instructions for how to do it the default way.
Next, we describe ONLY the key outputs that they'll want to use in a downstream
analysis.
Then we delve into the details at an intermediate level.
Finally we save truly advanced materials for either a different section entirely "Advanced:"
or we put them at the very bottom of the Module: section.
TODO: Break up this advanced section into subsection. My non-expert understanding
is they should be something like this:
- modifying mapping parameters (quality filters, paired only, fragment length, etc.)
- modified outputs for single-sample analysis
- concepts of how to call SNPs across samples
- passing a prebuilt index (I think this could also go into the maximizing
  performance section or be linked from there)
- Performance improvements for merging samples (this should also go into the
  performance optimization advanced section below.)
-->

Advanced users can adjust post-alignment filters via the following command-line
options (default values indicated):

```
--mapq >= 20: discard read alignments with alignment quality < 20
--mapid >= 0.94: discard read alignments with alignment identity < 0.94
--aln_readq >= 20: discard read alignment with mean quality < 20
--aln_cov >= 0.75: discard read alignment with alignment coverage < 0.75
--aln_baseq >= 30: discard bases with quality < 30
```

- Single-sample pileup for all the species in the restricted species profile
  with paired-ends based post-alignment quality filter.

With `paired_only` and proper `fragment_length`, MIDAS2 will only recruit
properly paired-ends reads for pileup. And all the post-alignment metrics will
be computed based on a read pair, instead of a single read.

```
midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
        --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
        --select_by median_marker_coverage,unique_fraction_covered \
        --select_threshold=2,0.5 \
        --fragment_length 2000 --paired_only \
        --num_cores 12 ${midas_outdir}
```

- Single-sample variant calling for all the species in the restricted
  species profile with paired-ends based post-alignment quality filter.

In recognition of the need for single-sample variant calling, we added an
`advanced` mode to the single-sample SNV analysis in MIDAS2. In the
`advanced` mode, per-species pileup results will also report major allele and
minor allele for all the genomic sites covered by at least two reads, upon
which custom variant calling filter can be applied by the users. Users are
advised to use the setting `ignore_ambiguous` to avoid falsely calling
major/minor alleles for sites with tied read counts.

```
midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
        --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
        --select_by median_marker_coverage,unique_fraction_covered \
        --select_threshold=2,0.5 \
        --fragment_length 2000 --paired_only \
        --advanced --ignore_ambiguous \
        --num_cores 12 ${midas_outdir}
```

- Single-sample pileup for all the species in a prebuilt rep-genome database.

```
midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
        --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
        --prebuilt_bowtie2_indexes ${midas_output}/bt2_indexes/repgenome \
        --prebuilt_bowtie2_species ${midas_output}/bt2_indexes/repgenome.species \
        --select_threshold=-1 --num_cores 12 ${midas_output}
```

### Output files

- `snps_summary.tsv`: summary of read mapping and pileup for all the species in
  the rep-genome database

```
species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered   mean_coverage
102506      5339468        2373275        8045342      468667         224553        0.444              3.390
102337      2749621        2566404        47723458     1479479        1010530       0.933              18.595
```

- _species_id_: six digits species id
- _genome_length_: genome length
- _covered_bases_: number of bases covered by at least two reads
- _total_depth_: total read depth across all _covered_bases_
- _aligned_reads_: total read counts across _covered_bases_ before post-alignment filter
- _mapped_reads_: total read counts across _covered_bases_ after post-alignment filter
- _fraction_covered_: fraction of _covered_bases_ (horizontal genome coverage)
- _mean_coverage_: mean read depth across all _covered_bases_ (vertical genome coverage)


- `{species}.snps.tsv.lz4`: per-species read pileup results

```
ref_id                    ref_pos   ref_allele  depth   count_a  count_c  count_g  count_t
gnl|Prokka|UHGG144544_1   881435    T           11      0        0        0        11
gnl|Prokka|UHGG144544_1   881436    T           13      0        5        0        8
gnl|Prokka|UHGG144544_1   881437    T           12      0        6        0        6
```

- _ref_id_: scaffold/contig id
- _ref_pos_: reference position
- _ref_allele_: reference nucleodie
- _depth_: number of mapped reads
- _count_a_: read counts of A allele
- _count_c_: read counts of C allele
- _count_g_: read counts of G allele
- _count_t_: read counts of T allele

- In the `advanced` mode, per-species pileup results will also include the
  called variants for all the covered genomic sites.

```
ref_id                    ref_pos   ref_allele  depth   count_a  count_c  count_g  count_t  major_allele  minor_allele  major_allele_freq  minor_allele_freq  allele_counts
gnl|Prokka|UHGG144544_1   881435    T           11      0        0        0        11       T             T             1.000              0.000              1
gnl|Prokka|UHGG144544_1   881436    T           13      0        5        0        8        T             C             0.615              0.385              2
gnl|Prokka|UHGG144544_1   881437    T           12      0        6        0        6        C             T             0.500              0.500              2
```

- _major_allele_: the allele with the most read counts
- _minor_allele_: the allele with the 2nd most read counts; same with major_allele if only one allele is observed
- _major_allele_freq_: allele frequency of _major_allele_
- _minor_allele_freq_: allele frequency of _minor_allele_; 0.0 if only one allele is observed
- _allele_counts_: number of alleles observed

## Population SNV Calling Analysis

Having run the single-sample SNV steps for all the samples, users next can
perform the population SNV analysis using the `merge_snps` command.

### Important Concepts

In this section, we will introduce the species and sample filters, the genomic
site filters, the compute of population SNV in MIDAS2, and the chunkified
pileup. Beginner users can skip this section and go straight to [Sample
commands]().

1. **<species, sub-samples-lists> selection**

   Population SNV analysis restricts attention to "sufficiently well" covered
   species in "sufficiently many" samples.

   To be specific, a given <species, sample> pair will only be kept if it has
   more than 40% horizontal genome coverage (`genome_coverage`) and 5X vertical
   genome coverage (`genome_depth`). Furthermore, only "sufficiently prevalent"
   species with "sufficiently many" (`sample_counts`) would be included for the
   population SNV analysis. Therefore, different species may have different
   lists of relevant samples.

2. **<site, relevant sample> selection**

   For each genomic site, a sample is considered to be "relevant" if the
   corresponding site depth falls between the range defined by the input
   arguments `site_depth` and `site_ratio * mean_genome_coverage`; otherwise it
   is ignored for the pooled-SNV compute.

   Therefore, different genomic sites from the same species may have different
   panels of "relevant samples".  And genomic site prevalence can be computed
   as the ratio of the number of relevant samples for the given site over the
   total number of relevant samples for the given species.

3. **Relevant site**

   For each species, a site is considered to be "relevant" if the site
   prevalence meets the range defined by the input arguments `snv_type` and
   `site_prev`. By default, common SNV with more than 90% prevalence are
   reported.

4. **Population SNV Computation**

   There are three main steps to compute and report population SNV in MIDAS2.

   First, for each relevant genomic site, MIDAS2 determines the set of
   alleles present across all relevant samples.  Specifically, for each
   allele (A, C, G, T), `merge_snps` subcommand (1) tallys the sample counts
   (_sc_) of relevant samples containing corresponding allele (`scA:scT`),
   and (2) sums up the read counts (_rc_) of the corresponding allele across
   all the relevant samples (`rc_G:rc_T`).

```
site_id                            rc_A  rc_C  rc_G  rc_T  sc_A  sc_C  sc_G  sc_T
gnl|Prokka|UHGG000587_14|34360|A   26    10    0     0     1     2     0     0
```

Second, population major and minor alleles for a single site can be computed
based on the accumulated read counts or sample counts across all relevant
samples. The population major allele refers to the most abundant/prevalent
allele, and the population minor allele refers to the second most
prevalent/abundant allele.

For example, the population major allele of the site
`gnl|Prokka|UHGG000587_14|34360|A` in the above example is `A` defined by
accumulated read counts and `C` defined by accumulated sample counts.

Third, MIDAS2 collects and reports the sample-by-site matrix of the
corresponding (1) site depth and (2) allele frequency of the above calculated
population minor allele for all the relevant samples. In these two
matrices, MIDAS2 encode `site_depth = 0` and `allele_frequency = -1` with
the special meaning of missing <site, sample> pair.

5. **Chunkified Pileup Implementation**

Both single-sample and across-samples pileup are parallelized on the unit of
chunk of sites, which is indexed by <species_id, chunk_id>. Only when all
chunks from the same species finished processing, chunk-level pileup results
will merged into species-level pileup result.

This implementation makes population SNV analysis across thousands of samples
possible. To compute the population SNV for one chunk, all the pileup results
of corresponding sites across all the samples need to be read into memory. With
the uses of multiple CPUs, multiple chunks can be processed at the same time.
Therefore, for large collections of samples, we recommend higher CPU counts and
smaller chunk size to optimally balance memory and I/O usage, especially for
highly prevalent species. Users can adjust the number of sites per chunk via
`chunk_size` (default value = 1000000). MIDAS2 also has a `robust_chunk`
option, where assigning different chunk sizes to different species based on the
species prevalence.


### Sample commands

- The species, samples, and sites filters, as well as the post-alignment
  filters can be customized with command-line options. For example,

- We can select species with horizontal coverage > 40%, vertical coverage > 3X
  and present in more than 30 relevant samples:

```
--genome_coverage 0.4 --genome_depth 3 --sample_counts 30
```

- We can apply the following site selections: only consider site with minimal
  read depth >= 5, and maximal read depth <= 3 * _mean_coverage_, and the
  minimal allele frequency to call an allele present is 0.05.

```
--site_depth 5 --site_ratio 3 --snp_maf 0.05
```

- We can only report populations SNV meeting the following criteria:
  bi-allelic, common population SNV (present in more than 80% of the
  population) from the protein coding genes based on accumulated sample counts.

```
--snp_type bi --snv_type common --site_prev 0.8 --locus_type CDS --snp_pooled_method prevalence
```

Now we can put all the above-mentioned filters in one `merge_snps` command:

```
midas2 merge_snps --samples_list ${my_sample_list} \
    --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
    --genome_coverage 0.4 --genome_depth 3 --sample_counts 30 \
    --site_depth 5 --site_ratio 3 --snp_maf 0.05 \
    --snp_type bi --snv_type common --site_prev 0.8 --locus_type CDS --snp_pooled_method prevalence \
    --num_cores 32 ${midas_outdir}
```

### Output files

- `snps_summary.tsv`: merged single-sample pileup summary for all the species
  in the restricted species profile. The reported columns
  _covered_bases_:_mean_coverage_ are the same with single-sample pileup
  summary.

```
sample_name  species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered  mean_coverage
SRR172902    100122      2560878        2108551        10782066     248700         207047        0.823             5.113
SRR172903    100122      2560878        2300193        39263110     1180505        820736        0.898             17.069
```

- _sample_name_: unique sample name
- _species_id_: six-digit species id

- `{species_id}.snps_info.tsv.lz4`: per species SNV metadata.

```
site_id                             major_allele  minor_allele  sample_counts  snp_type  rc_A  rc_C  rc_G  rc_T  sc_A  sc_C  sc_G  sc_T  locus_type  gene_id           site_type  amino_acids
gnl|Prokka|UHGG000587_14|34360|A    A             C              2             bi        26    10    0     0     2     2     0     0     CDS         UHGG000587_02083  4D         T,T,T,T
gnl|Prokka|UHGG000587_11|83994|T    G             T              2             bi        0     0     11    45    0     0     2     2     IGR         None              None       None
```

<!--
NOTE: This comment fixes syntax highlighting in vim
-->
- _site_id_: unique site id, composed of f"{ref_id}|{ref_pos}|{ref_allele}"
- _major_allele_: most common/prevalent allele in metagenomes
- _minor_allele_: second most common/prevalent allele in metagenomes
- _sample_counts_: number of relevant samples where metagenomes is found
- _snp_type_: the number of alleles observed at site (mono,bi,tri,quad)
- _rc_A_: accumulated read counts of A allele in metagenomes
- _rc_C_: accumulated read counts of C allele in metagenomes
- _rc_G_: accumulated read counts of G allele in metagenomes
- _rc_T_: accumulated read counts of T allele in metagenomes
- _sc_A_: accumulated sample counts of A allele in metagenomes
- _sc_C_: accumulated sample counts of C allele in metagenomes
- _sc_G_: accumulated sample counts of G allele in metagenomes
- _sc_T_: accumulated sample counts of T allele in metagenomes
- _locus_type_: CDS (site in coding gene), RNA (site in non-coding gene), IGR (site in intergenic region)
- _gene_id_: gene identified if locus type is CDS, or RNA
- _site_type_: indicates degeneracy: 1D, 2D, 3D, 4D
- _amino_acids_: amino acids encoded by 4 possible alleles

- `{species_id}.snps_freq.tsv.lz4`: per species site-by-sample allele frequency
  matrix of population minor allele.

```
site_id                             SRR172902   SRR172903
gnl|Prokka|UHGG000587_11|83994|T    0.692       0.837
gnl|Prokka|UHGG000587_14|34360|A    0.300       0.269
```

- `{species_id}.snps_depth.tsv.lz4`: per species site-by-sample site depth
  matrix. Only accounts for the alleles matching the population major
  and/or minor allele.

```
site_id                             SRR172902   SRR172903
gnl|Prokka|UHGG000587_11|83994|T    13          43
gnl|Prokka|UHGG000587_14|34360|A    10          26
```


# Module: Copy-Number Variant Analysis

Similar to the SNV module, the CNV module also proceeds in two phases: (1)
single-sample pan-genome copy number variants calling (2) merge these results
into a summary across all samples. The first step can be run in parallel. We
presuppose users already follow the [database
customization](https://github.com/czbiohub/MIDAS2.0/wiki/Data-customization)
step, and have `species_profile.tsv` ready for each sample.

## Single-Sample CNV Calling

Species pangenome refers to the set of non-redundant genes (centroids)
clustered from all the genomes within on species cluster. Species in the
restricted species profile are concatenated and used to build the
sample-specific pangenome database, to which reads are aligned using Bowtie2.

Per species per centroid copy numbers are computed in three steps: (1) Per
centroid, read alignment metrics, e.g _mapped_reads_ and _mean_coverage_, are
computed; (2) Per species, median read coverage of all the mapped centroids
corresponding to the 15 universal SCGs are identified; (3) Per centroid, `copy
numbers` are computed and gene presence/absence are further inferred.

### Sample commands

- Single-sample CNV calling for all the species in the restricted species
  profile: `median_marker_coverage > 2` and `unique_fraction_covered` > 0.5.

We presuppose users already [profiling the species
coverage](https://github.com/czbiohub/MIDAS2.0/wiki/Data-customization#species-to-genotype),
and expect `${my_midasdb_dir}/${sample_name}/species/species_profile.tsv`
exists.

```
midas2 run_genes --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
        --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
        --select_by median_marker_coverage,unique_fraction_covered \
        --select_threshold=2,0.5 \
        --num_cores 12 ${midas_outdir}
```

Users can adjust post-alignment quality filter parameters via the command-line
options, and the defaults are:

```
--mapq >= 2: reads aligned to more than one genomic locations equally well are discarded (MAPQ=0,1)
--mapid >= 0.94: discard read alignments with alignment identity < 0.94
--aln_readq >= 20: discard read alignment with mean quality < 20
--aln_cov >= 0.75: discard read alignment with alignment coverage < 0.75
```

### Output files

- `genes_summary.tsv`: read mapping and pan-gene coverage summary for all the
  species in the sample-specific pan-genome database

```
species_id  pangenome_size  covered_genes  fraction_covered  mean_coverage  aligned_reads  mapped_reads   marker_coverage
102337      15578           4468           0.287             16.213         1650361        450353         20.213
102506      731186          4733           0.006             3.803          681335         37272          2.140
```

- _species_id_: six-digit species id
- _pangenome_size_: number of centroids (non-redundant genes) in the species pangenome
- _covered_genes_: number of centroids covered with at least one read
- _fraction_covered_: fraction of _covered_genes_ over _pangenome_size_
- _mean_coverage_: average read depth across _covered_genes_
- _aligned_reads_: total number of aligned reads before post-alignment filter
- _mapped_reads_: total number of aligned reads after post-alignment filter
- _marker_coverage_: average read depth across 15 universal SCGs in the species pangenome

- `{species_id}.genes.tsv.lz4`: per-species gene content profiling summary

```
gene_id              gene_length  aligned_reads  mapped_reads  mean_coverage  fraction_covered  copy_number
UHGG143901_00483     555          14             6             2.961538       0.234234          1.384035
UHGG143901_03589     384          103            57            32.840708      0.294271          15.347667
UHGG143902_04031     207          9              2             1.737500       0.386473          0.811997
```

- _gene_id_: id of centroid in the species pan-genome
- _gene_length_: gene length
- _aligned_reads_: number of aligned reads to _gene_id_ before post-alignment filter
- _mapped_reads_: number of aligned reads to _gene_id_ after post-alignment
  filter
- _mean_coverage_: average read depth of _gene_id_ based on _mapped_reads_
  (total_gene_depth / covered_bases)
- _fraction_covered_: proportion of the _gene_id_ covered by at least one read
  (covered_bases / gene_length)
- _copy_number_: estimated copy number of _gene_id_ based on _mapped_reads_
  (_mean_coverage_ / median_marker_coverage)

## Across-Samples CNV Calling

Take the same [`${my_sample_list}`](#across-samples-analysis):

```
sample_name   midas_outdir
SRR172902     /home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample
SRR172903     /home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample
```

`merge_genes` would expect to locate
`/home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample/SRR172903/genes/genes_summary.tsv`,
generated by `run_genes`.


Having run the single-sample CNV analysis for all the samples listed in the
`${my_sample_list}`, users next can merge the results and product a summary
with `merge_genes` command, and further quantify each pan-gene's
presence/absence by comparing the `copy_number` with the user-defined minimal
gene copy number (`min_copy`).


### Sample commands

- Across-samples CNV analysis using default filters.

```
midas2 merge_genes --samples_list ${my_sample_list} --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} --num_cores 8 ${midas_outdir}
```

By default, merged gene CNV results are reported for genes clustered at 95%
identity. `cluster_pid` and `min_copy` can be customized with the following
command-line options:

```
--genome_depth: filter out species with _mean_coverage_ < 1X.
--min_copy: genes with _copy_number_ >= 0.35 are classified as present.
--cluster_pid: gene CNV results can be reported at various clustering cutoffs {75, 80, 85, 90, 95, 99}.
```

### Output files

- `genes_summary.tsv`: merged single-sample gene content summary. The reported columns _covered_genes_:_marker_coverage_ are the same with single-sample CNV summary.

```
sample_name  species_id  pangenome_size  covered_genes  fraction_covered  mean_coverage  aligned_reads  mapped_reads  marker_coverage
SRR172902    100122      29165           2535           0.087             4.723          263395         53006         1.435
SRR172903    100122      29165           3212           0.110             16.095         1447684        263878        10.713
```
- _sample_name_: unique sample name
- _species_id_: six-digit species id

- `{species_id}.genes_copynum.tsv.lz4`: per species gene-by-sample copy number matrix

```
gene_id            SRR172902     SRR172903
UHGG000587_00401   33.969154     19.891455
UHGG000587_01162   5.703398      2.821237
UHGG000587_00962   2.370930      0.289325
```

- `{species_id}.genes_preabs.tsv.lz4`: per species gene-by-sample presence absence matrix

```
gene_id             SRR172902    SRR172903
UHGG000587_00401    1            1
UHGG000587_01162    1            1
UHGG000587_00962    1            0
```

- `{species_id}.genes_depth.tsv.lz4`: per species gene-by-sample mean coverage matrix

```
gene_id             SRR172902   SRR172903
UHGG000587_00401    48.747945   213.090622
UHGG000587_01162    8.184746    30.222978
UHGG000587_00962    3.402439    3.099448
```

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

# Advanced: Optimizing Index Building

## Population Specific Species Panel

The Bowtie2 rep-genome / pan-genome database can be build upon a list of
customized species across a given panel of samples. Species selection metrics
based on the `species_prevalence.tsv` can be passed to `build_bowtie2eb` via
`--select_by` and `--select_threshold`. For example, to build one rep-genome
database for all the species that is present in more than two samples:

```
midas2 build_bowtie2db \
    --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
    --select_threshold sample_counts --select_by 2 --num_cores 8 \
    --bt2_indexes_name repgenome --bt2_indexes_dir ${midas_outdir}/bt2_indexes
```

  - The generated rep-genome database can be found under the directory `${midas_outdir}/bt2_indexes`
  - The list of customized species can be found at `${midas_outdir}/bt2_indexes/repgenome.species`


If taking this approach, for the single-sample SNV or CNV analysis, users can
pass the pre-built rep-genome to `run_snps` analysis (pan-genome for
`run_genes`), as following:

```
--prebuilt_bowtie2_indexes ${midas_outdir}/bt2_indexes/repgenome \
--prebuilt_bowtie2_species ${midas_outdir}/bt2_indexes/repgenome.species \
--select_threshold=-1
```


# Advanced: Building Your Own MIDASDB

MIDAS2 users can locally build a new MIDASDB for a custom collection of
genomes. The target layout of MIDASDB can refer to [this page](TODO).  This
page is focused specifically on the database construction commands.

To start with, users need to organize the genomes in a specific format and
produce the TOC `genomes.tsv` as described in [the documentation](TODO)

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

### MIDAS2 Results Layout

MIDAS2 writes its outputs to a user-specified root directory, which is always
passed as a mandatory argument to each of the MIDAS2 command. In this
documentation we refer to this directory as:

- `midas_outdir=/path/to/results/root/directory`

Together with the unique sample name:

- `sample_name=/unique/id/of/sample`

`${midas_outdir}/${sample_name}` constitute the unique output directory for
single-sample analysis. All subsequent analysis steps operate within this
directory. For the across-samples SNV or CNV analysis, all subsequent analysis
steps operate within the directory `${midas_outdir}`.


### Single-Sample Results Layout

MIDAS analysis usually starts with species selection which
selects sufficiently abundant species in each sample (subcommand
`run_species`). After completing this step, users can run either of two
strain-level analysis: `run_snps` for single-sample read pileup (SNV module) or
`run_genes` for pan-gene profiling (CNV module).  Here is an example of the
results layout of all single-sample analysis in the local filesystem.

```
Output                                       Producer          Meaning
--------------------------------------------------------------------------------------------------------
{midas_output}/{sample_name}
  |- species
     |- species_profile.tsv                  run_species       Summary of species coverage
     |- markers_profile.tsv                  run_species       Per species marker coverage
  |- snps
     |- snps_summary.tsv                     run_snps          Summary of read mapping to rep-genome
     |- {species}.snps.tsv.lz4               run_snps          Per species pileups
  |- genes
     |- genes_summary.tsv                    run_genes         Summary of read mapping to pan-genome
     |- {species}.genes.tsv.lz4              run_genes         Per species pan-gene coverage

 |- temp
     |- snps
        |- repgenomes.bam                    run_snps          Rep-genome alignment file
        |- {species}/snps_XX.tsv.lz4
     |- genes
        |- pangenome.bam                     run_genes         Pan-genome alignment file
        |- {species}/genes_XX.tsv.lz4
  |- bt2_indexes
     |- snps/repgenomes.*                    run_snps          Sample-specific rep-genome database
     |- genes/pangenomes.*                   run_genes         Sample-specific pan-genome database
```

### Across-Samples Results Layout

For a collection of samples, population SNVs and pan-genome CNVs can be
estimated using subcommands `merge_snps` and `merge_genes`.

```
Output                                             Producer        Meaning
---------------------------------------------------------------------------------------------------------------
{midas_output}
  |- species
    |- species_prevalence.tsv                      merge_species   Per species summary statistics across samples
    |- species/species_read_counts.tsv             merge_species   Species-by-sample read counts matrix
    |- species/species_coverage.tsv                merge_species   Species-by-sample marker coverage matrix
    |- species/species_rel_abundance.tsv           merge_species   Species-by-sample relative abundance matrix
  |- snps
    |- snps_summary.tsv                            merge_snps      Alignment summary statistics per sample
    |- {species}/{species}.snps_info.tsv.lz4       merge_snps      Per species metadata for genomic sites
    |- {species}/{species}.snps_freqs.tsv.lz4      merge_snps      Per species site-by-sample MAF matrix
    |- {species}/{species}.snps_depth.tsv.lz4      merge_snps      Per species site-by-sample read depth matrix
  |-genes
    |- genes_summary.tsv                           merge_genes     Alignment summary statistics per sample
    |- {species}/{species}.genes_presabs.tsv.lz4   merge_genes     Per species gene-by-sample pre-abs matrix
    |- {species}/{species}.genes_copynum.tsv.lz4   merge_genes     Per species gene-by-sample copy number matrix
    |- {species}/{species}.genes_depth.tsv.lz4     merge_genes     Per species gene-by-sample read depth matrix
```

## MIDASDB Directory Layout

To meet the challenge of increased number of available genome sequences, MIDAS2
implemented a new database infrastructure, geared to run on
[AWS Batch](https://aws.amazon.com/batch/) and
[S3](https://aws.amazon.com/s3/), to
achieve [elastic scaling](https://github.com/czbiohub/pairani/wiki) for
building MIDAS reference databases. To be specific, the MIDAS reference
database construction step can be executed in AWS using hundreds of
r5d.24xlarge instances over a period of a couple of days, depositing built
products in S3.  For example, it took ~\$80,000 and a week to build the species
pan-genome for all 47,894 species of GTDB r202.


### Table of Content

The new database infrastructure reads in a table of contents (TOC) file,
containing genome-to-species assignment and a choice of representative genome
for each species cluster.  One TOC file (`genomes.tsv`) per MIDAS reference
database, listing the genome-to-species assignment for all genomes, with per
genome each row. The TOC file has four columns, among which
`genome_is_representative` specify whether the `genome` is the representative
genome for the corresponding `species`. Only one `representative` per
`species`.

```
genome              species      representative        genome_is_representative
GUT_GENOME138501    104351       GUT_GENOME269084      0
GUT_GENOME269084    104351       GUT_GENOME269084      1
```

By default, MIDAS2 inherits the representative genome assignments from
published prokaryotic genome databases. Inspired by the importance of selecting
proper reference genome for accurate template-based SNP calling, this new
infrastructure empowers user the flexibility to dynamically re-assign the
representative genomes, simply by modifying the `genomes.tsv` file accordingly.

**Unified Human Gastrointestinal Genome (UHGG)**: A collection of 286,997 genomes
assembled from metagenomes, isolates and single cells from human stool samples
has been clustered into 4,644 gut-only species in [UHGG 1.0
catalogues](http://ftp.ebi.ac.uk/pub/databases/metagenomics/mgnify_genomes/human-gut/v1.0/).
The collection of all the UHGG genomes were mirrored in
[S3](s3://jason.shi-bucket/IGGdb2.0/clean_set/), which serves as the input to
the database construction. [Six-digit numeric species
ids](s3://jason.shi-bucket/IGGdb2.0/alt_species_ids.tsv) were arbitrarily
assigned. Instead of species name, these `species_id` are used as species
identifier in all the reports generated by MIDAS2.

**Genome Taxonomy Database (GTDB)**:
[GTDB R06-RS202](https://gtdb.ecogenomic.org/stats/r202) contains 45,555
bacterial
and 2,339 archaeal species clusters spanning 258,406 genomes, released on April
27th, 2021. The genome members for each species cluster is specified in the
[sp_clusters_r202.tsv](https://data.ace.uq.edu.au/public/gtdb/data/releases/release202/202.0/auxillary_files/sp_clusters_r202.tsv),
upon which order six-digit numeric species ids are assigned. GTDB only provided
the sequences of the representative genomes, and we downloaded all the genomes
from NCBI genomes repository using
[genome_updater](https://github.com/pirovc/genome_updater).


### Database Target Layout and Construction

MIDAS reference database (MIDASDB) primarily consist of three parts: rep-genome
databases, pan-genome databases, and universal single copy genes (SGC) marker
database. The target layout of any MIDASDB follow the same relative structure,
based on the base directory of the database, both in the S3 bucket and locally.
The following toy example demonstrates the major steps to construct the MIDASDB
and the target layout using a toy collection of genomes with only one species
cluster `species1` with two genomes (`genome1` and `genome2`).


![
MIDAS Reference Database Target Layout and Construction Steps
](static/Fig.DB.Layout.png)


### Inputs

The input collection of genomes need to be organized in the format as
`cleaned_genomes/<species>/<genome>/<genome>.fna`. And the table of content
`genomes.tsv` file needs to be generated accordingly, with randomly
assigned six-digit `species_id`, to replace the species name. The `genome` name
can be kept as it is.

```
genome     species    representative    genome_is_representative
genome1    100001     genome2           0
genome2    100001     genome2           1
```

### Rep-Genome Databases

The genome annotation for all the genomes were done by
[Prokka](https://github.com/tseemann/prokka), and the annotated genes were kept
under the directory of `genes_annotations/<species>/<genome>`. The rep-genome
databases for the SNPs module analysis only included the gene annotations and
sequences for the representative genomes, as specified in the TOC.

```
gene_annotations/100001/genome2/genome2.fna.lz4
gene_annotations/100001/genome2/genome2.ffn.lz4
gene_annotations/100001/genome2/genome2.genes.lz4
```

### SCG Marker Database

Marker genes are defined as universal, single-copy gene families. MIDAS uses
a subset (15) of the [PhyEco gene
families](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0077033).
The pre-computed HMM model of this set of 15 single copy genes (SCGs) are
available at:

```
s3://microbiome-pollardlab/uhgg_v1/marker_gene_models/phyeco/marker_genes.hmm.lz4
s3://microbiome-pollardlab/uhgg_v1/marker_gene_models/phyeco/marker_genes.mapping_cutoffs.lz4
```

For each annotated genome, the homologs of 15 SCGs were identified with
`hmmsearch`, as well as the mapping of gene id to corresponding marker gene id,
under the directory of `marker_genes/phyeco/temp/<species>/<genome>`.

```
marker_genes/phyeco/temp/100001/genome2/genome2.markers.fa
marker_genes/phyeco/temp/100001/genome2/genome2.markers.map
```

For all the representative genomes, the identified marker genes were
concatenated into monolithic `marker_genes.fa`, from which `hs-blastn`
index would be constructed. The indexed `marker_genes.fa` serves as the SCG
marker databases.

```
marker_genes/phyeco/marker_genes.fa
marker_genes/phyeco/marker_genes.fa.sa
marker_genes/phyeco/marker_genes.fa.bwt
marker_genes/phyeco/marker_genes.fa.sequence
```

### Pan-Genome Databases

Species-level pan-genome refers to the set of non-redundant genes that
represent the genetic diversity within one species cluster. In order to
construct the pan-genome database for each species, the first step if to
concatenate the annotated genes from its all genome members into
`pangenomes/100001/genes.ffn`. The second step, which is also the most
time-consuming step, is to cluster the concatenated genes based on 99% percent
identity (PID) using [`vsearch`](https://github.com/torognes/vsearch). Each
cluster was represented by the gene at its center - centroid gene
(`centroids.99.ffn`). The `centroid.99` genes were further on clustered to 95,
90, ..., PID, respectively, and the mapping relationships were listed in
`centroid_info.txt`. The top level `centroids.ffn` file represents the 99
percent identity clusters, and serves as the species pan-genome databases.
Reads are aligned to the pan-genome databases to determine the gene content of
strains in a sample (`midas run_genes`), and reads can optionally aggregated
into gene clusters at any of the lower clustering thresholds across samples
(`midas2 merge_gene`).

```
pangenomes/100001/centroids.ffn
pangenomes/100001/centroid_info.txt
```

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
