# Metagenomic Intra-Species Diversity Analysis 2.0

Metagenomic Intra-Species Diversity Analysis ([MIDAS](https://genome.cshlp.org/content/26/11/1612)) is an integrated pipeline for profiling strain-level genomic variation and gene copy number variations for metagenomic data. Its analysis steps are run against a database of 5,926 bacterial species extracted from 30,000 genomes (MIDAS DB v1.2).

MIDAS 2.0 represent a reimplementation of the same analysis steps as the original [MIDAS tool](https://github.com/snayfach/MIDAS), but able to operate on the more comprehensive MIDAS DB for a larger collection of samples in a fast and scalable manner. Similar to the original MIDAS, MIDAS 2.0 also presuppose a database construction step has already taken place. This page is focused specifically on the analysis steps, and database construction steps can refer [here](https://github.com/czbiohub/MIDAS2.0/wiki/MIDAS-DB). 

***

# MIDAS Results Layout

MIDAS is an all-in-one strain-level metagenomics pipeline. Its scope ranges from building custom MIDAS databases, species profiling, reads alignment, post-alignment filter, to strain-level metagenotyping and pan-gene copy number profiling. There are three modules of MIDAS pipeline: Species, SNPs, Genes. Each module consists of single-sample and across-samples two workflow.

[insert a figure here.]

## Single-sample Results Layout

All single-sample workflow takes as a parameter the path to MIDAS results root directory (`midas_outdir`); and together with `sample_name` constitute the unique output directory {`output_dir`}, i.e.,  `{midas_outdir}/{sample_name}`.  All subsequent analysis steps operate within that directory.

Single-sample analysis start with finding out abundant species in the given sample (`run_species`), followed by `run_snps` for pileup and `run_genes` pan-gene profiling.  Here is an example of layout of the results of all three single-sample modules in the local filesystem.

```
Output                                       Producer             Meaning
------------------------------------------------------------------------------------------------------------
{midas_output}/{sample_name}
  |- species
     |- species_profile.tsv                  midas_run_species    Summary of species coverage
     |- markers_profile.tsv                  midas_run_species    Per species marker coverage
  |- snps
     |- snps_summary.tsv                     midas_run_snps       Summary of reads mapping to rep-genome
     |- {species}.snps.tsv.lz4               midas_run_snps       Per species pileups
  |- genes 
     |- genes_summary.tsv                    midas_run_genes      Summary of reads mapping to pan-geneme
     |- {species}.genes.tsv.lz4              midas_run_genes      Per species pan-gene coverage
 
 |- temp                                                          Temporary Files
     |- snps
        |- repgenomes.bam
        |- {species}/snps_XX.tsv.lz4
     |- genes
        |- pangenome.bam
        |- {species}/genes_XX.tsv.lz4
  |- bt2_indexes                                                  Sample-specific genome database
     |- snps/repgenomes.*
     |- genes/pangenomes.*
```

## Across-samples Results Layout

For a collection of samples, population SNPs and pan-gene copy numbers can be estimated using subcommands `merge_snps` and `merge_genes`. User need to provide the output directory (`{midas_outdir}`) as the input parameter, which is the root of the layout below.

```
Output                                          Producer             Meaning
-----------------------------------------------------------------------------------------------------------
species
  |- species_prevalence.tsv                     midas_merge_species  Per species summary statistics across samples
  |- species/species_read_counts.tsv            midas_merge_species  Species-by-sample read counts matrix
  |- species/species_coverage.tsv               midas_merge_species  Species-by-sample marker coverage matrix
  |- species/species_rel_abundance.tsv          midas_merge_species  Species-by-sample relative abundance matrix
snps
  |- snps_summary.tsv                           midas_merge_snps     Alignment summary statistics per sample
  |- {species}/{species}.snps_info.tsv.lz4      midas_merge_snps     Per species metadata for genomic sites.
  |- {species}/{species}.snps_freqs.tsv.lz4     midas_merge_snps     Per species site-by-sample MAF matrix
  |- {species}/{species}.snps_depth.tsv.lz4     midas_merge_snps     Per species site-by-sample read depth matrix
genes
  |- genes_summary.tsv                          midas_merge_genes    Alignment summary statistics per sample
  |- {species}/{species}.genes_presabs.tsv.lz4  midas_merge_genes    Per species gene-by-sample pre-abs matrix
  |- {species}/{species}.genes_copynum.tsv.lz4  midas_merge_genes    Per species gene-by-sample copy number matrix
  |- {species}/{species}.genes_depth.tsv.lz4    midas_merge_genes    Per species gene-by-sample read depth matrix
```


***

# Single-sample Analysis

## Species Abundance Estimation

For each sample, the analysis begins with a simple species profiling. The goal of the Species flow is to detect abundant species that is present in the sample, which can be used to construct the sample-specific representative genome database (rep-genome) and pangenome database (pan-genome). Raw metagenomic reads were mapped to the 15 universal single copy genes (SCGs). And for each marker gene, uniquely mapped read counts were computed, and ambiguous reads were probabilistically assigned.  The marker coverage is computed as the total alignment length over the gene length. Only species with more than two marker genes with more than two reads are reported. 

### Sample command

  ```
  midas2 run_species --sample_name ${my_sample} -1 /path/to/R1 -2 /path/to/R2 \
         --midasdb_name uhgg --midasdb_dir /path/to/local/midasdb --num_cores 8 ${midas_outdir}
  ```

### Output files

- `species_profile.tsv`: sorted in decreasing order of `median_coverage`. 

   ```
   species_id  marker_read_counts  median_marker_coverage  marker_coverage  marker_relative_abundance  total_covered_marker  unique_covered_marker  ambiguous_covered_marker  total_marker_counts  unique_fraction_covered  total_marker_length
   102470      329                 3.53                    4.16             0.54                       15                    15                     5                         15                   1.00                     9870     
   100039      199                 3.45                    2.45             0.32                       11                    11                     0                         15                   0.73                     10095
   ```

- `markers_profile.tsv`: species-by-marker coverage matrix


   ```
   species_id  marker_id  marker_length  gene_id           total_reads  total_alnbps  coverage  uniq_reads  ambi_reads   uniq_alnbps     ambi_alnbps
   102470      B000032    702            UHGG143484_04271  22           2763          3.94      18          4            2259            504
   102470      B000039    606            UHGG143484_01873  7            856           1.41      3           4            379             477
   ```

Refer to [original MIDAS's estimate species abundance](https://github.com/snayfach/MIDAS/blob/master/docs/species.md) for more details.


## Single Nucleotide Polymorphisms (SNPs) calling

Only species passing the user specific filter based on the above mentioned Species module would be genotyped and construct the sample-specific genome database. Although the default parameters for the species selection is set to detect abundant species in the sample, users can adjust the parameters to suit the purpose of their research. That being said, to genotype low abundant species, we suggested loosen the parameters to `median_marker_coverage` > 0. 

To explore intra-species variations in the sample, raw metagenomic reads were aligned to the sample-specific genome databases. Nucleotide variation for each genomic site was then quantified via pileup and allele count. Reads mapping summary for all the species in the rep-genome database were reported in the per-sample snps summary file. MIDAS doesn't apply any more filters based on the reads mapping at this stage.


### Sample command


- Perform Pileup for all the species in the restricted species profile: `median_marker_coverage` > 2 and `unique_fraction_covered` > 0.5
   
   ```
   midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R2} \
         --midasdb_name uhgg --midasdb_dir /path/to/local/midasdb \
         --num_cores 12 --fragment_length 1000 \
         --select_by median_marker_coverage,unique_fraction_covered \
         --select_threshold=2,0.5 \
          ${midas_outdir}
    ```

- Genotyping with Prebuilt Genome Database

  `--select_threshold=-1`: skip the Species module and pileup for all the species in a prebuilt Bowtie2 genome databases. Use with caution.

   ```
   midas2 run_snps --sample_name ${sample_name} -1 $R1 -2 $R2 \
          --midasdb_name uhgg --midasdb_dir /path/to/local/midasdb \
          --prebuilt_bowtie2_indexes ${bt2_indexes/${bt2_name} \
          --prebuilt_bowtie2_species ${bt2_indexes/${bt2_name}.species \
          --select_threshold=-1 --num_cores 12 ${midas_outdir}
   ``` 


### Output files

- Reads mapping and Pileup summary: `snps_summary.tsv`

   ```
   species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered mean_coverage
   102470      5774847        4379800        23295116     213486         202365        0.758            5.319
   100039      2852528        1764680        25765375     225417         224630        0.619            14.601
   ```

- Per-species Pileup results: `{species}.snps.tsv.lz4`

   ```
   ref_id                    ref_pos   ref_allele  depth  count_a  count_c  count_g  count_t  major_allele  minor_allele  major_allele_freq  minor_allele_freq  allele_counts
   gnl|Prokka|UHGG143484_2   531422    C           5       0        5       0        0        C              C            1.000              1.000              1
   gnl|Prokka|UHGG143484_2   531423    T           6       2        0       0        4        T              A            0.667              0.333              2
   gnl|Prokka|UHGG143484_2   531424    A           6       6        0       0        0        A              A            1.000              1.000              1
   ```

### New Features

- **Chunkified Pileup**

  Single-sample Pileup was parallelized on the unit of chunk of sites, which is indexed by `species_id, chunk_id`. When all chunks from the same species finished processing, chunk-level Pileup results were then merged into species-level Pileup file (`{species}.snps.tsv.lz4`).


- **Re-assign Representative Genome**

  Users can re-select representative genome by modifying the table of content `genomes.tsv` accordingly.

- **Custom MIDAS DB**

  The new infrastructure of MIDAS 2.0 dramatically simplifies the prior knowledge needed to build a custom MIDAS database. The new implementation of MIDAS DB reads in a Table Of Contents (TOC) file containing genome-to-species assignment, and a choice of representative genome for each species.

  First, generate the `{custom_midasdb}/genomes.tsv` in the following format:

     ```
     genome            species representative          genome_is_representative
     GUT_GENOME091053  100001  GUT_GENOME091053        1
     GUT_GENOME091054  100001  GUT_GENOME091053        0
     GUT_GENOME178957  100003  GUT_GENOME178957        1
     GUT_GENOME178958  100003  GUT_GENOME178957        0
     GUT_GENOME178959  100003  GUT_GENOME178957        0
     ```

  Second, collect all the representative genomes as following:

     ```
     ${custom_midasdb}/cleaned_imports/{species_id}/{genome_id}/{genome_id}.fna
     ${custom_midasdb}/cleaned_imports/100001/GUT_GENOME091053/GUT_GENOME091053.fna.lz4
     ${custom_midasdb}/cleaned_imports/100003/GUT_GENOME178957/GUT_GENOME178957.fna.lz4
     ```

  Third, run `run_snps` with `--midasdb_name {custom_midasdb} --midasdb_dir /path/to/midasdb_dir`.


Refer to [MIDAS's call single nucleotide polymorphisms](https://github.com/snayfach/MIDAS/blob/master/docs/cnvs.md) for more details.


## Pangenome Profiling

Similar with the single-sample SNPs module, only abundant species in the restricted species profile would be used to build the sample-specific pangenome database. For each species, the hierarchical compuation for each chunk of genes is: (1) For each gene, compute reads alignment based metrics, e.g. `aligned_reads`, `mapped_reads`, etc; (2) For all the pan-genes, compute the average vertical coverage of the 15 universal SCGs; (3) For each gene, infer the `copy number` normalized to the SGCs coverage, and quantify the pan-gene presence/absence.


### Sample command

   ```
   midas2 run_genes --sample_name ${sample_name} -1 ${R1} -2 ${R2} \
          --midasdb_name uhgg --midasdb_dir /path/to/local/midasdb \
          --select_by median_marker_coverage,unique_fraction_covered \
          --select_threshold=2,0.5 \
          --num_cores 12 ${midas_outdir}
   ```

### Output files

- `genes_summary.tsv`

   ```
   species_id  pangenome_size  covered_genes  fraction_covered  mean_coverage  aligned_reads  mapped_reads   marker_coverage
   102470      87136           6122           0.070             9.519          455045         269656         3.447
   100039      79571           2350           0.030             9.462          206684         113644         0.000
   ```

- `{species_id}.genes.tsv.lz4`: 

   ```
   gene_id              gene_length  aligned_reads  mapped_reads  mean_coverage  fraction_covered  copy_number
   UHGG000947_00080     3138         53             6             2.208723       0.102294          0.640716
   UHGG000947_00081     609          27             25            5.083770       0.940887          1.474722
   ```

Refer to [MIDAS's predict pan-genome gene content](https://github.com/snayfach/MIDAS/blob/master/docs/cnvs.md) for more details.


***

# Across-samples Analysis

MIDAS 2.0 can merge single-sample analysis results across multiple samples, such as performing across-samples core-genome SNPs calling, or merging pan-genome profiling results across samples. 

All three merge subcommands take a `samples_list` as the input argument, which is a TSV file with `sample name` and single-sample output directory `midas_outdir`.  For example, `midas_merge_species` expects to locate `/mnt/cz/hmp-test/SRS011134/species/species_profile.tsv`, generated by `run-species`.

   ```
   sample_name   midas_outdir
   SRS011134     /mnt/cz/hmp-test
   SRS011271     s3://microbiome-chunyu/iggtools-hmp-test
   ```


## Species Abundance Profile

## Sample command

   ```
   midas2 merge_species --samples_list /path/to/sample/lists ${merged_midas_outdir}
   ```

### Output files

- `species_prevalence.tsv`: summary statistics for each species across samples

   ```
   species_id  median_abundance  mean_abundance  median_coverage  mean_coverage  sample_counts
   102470      0.462             0.334           3.533            2.355          2
   100039      0.317             0.320           3.455            3.455          3
   ```

- `species_marker_read_counts.tsv`: species-by-samples marker read counts matrix

   ```
   species_id   sample1   sample2   sample3  
   102470       329       0.000     332
   100039       199       199       196
   ```

- `species_marker_median_coverage.tsv`: species-by-samples median marker coverage matrix

   ```
   species_id   sample1   sample2    sample3  
   102470       3.533      0.000     3.533
   100039       3.455      3.455     3.455
   ```

- `species_marker_coverage.tsv`: species-by-samples mean marker coverage matrix

   ```
   species_id   sample1   sample2    sample3  
   102470       4.161      0.000     4.199
   100039       2.446      2.446     2.408
   ```

- `species_relative_abundance.tsv`: species-by-samples relative abundance matrix

   ```
   species_id   sample1   sample2    sample3  
   102470       0.539     0.000      0.462
   100039       0.317     0.378      0.265
   ```

Refer to [MIDAS's merge species abundance](https://github.com/snayfach/MIDAS/blob/master/docs/merge_species.md) for more details.


## Population SNPs Calling

For each **relevant** genomic site, MIDAS 2.0 determine the set of alleles present for that site across all **relevant** samples. For each allele A, C, G, T we count the samples that are **relevant** for the site and contain that allele, and sum that allele's depths across those samples.  Then we computed the across-samples major alleles for each genomic site, following by collecting the corresponding read depth and allele frequency for each sample.  We use one or the other of these metrics as a proxy for allele frequency, as specified by the `snp_pooled_method ` argument.

- **pipeline details**

  1. <species, sub-sample-lists> selection: the analysis restricts attention to "sufficiently well" covered species in sufficiently many samples. (`genome_depth`, `genome_coverage`, and `sample_counts`). Therefore, despite the provided list of samples by the user, different species may have different lists of relevant samples.

  2. **relevant sample**: determined by per-site sample filters. A sample is considered relevant for a given genomic site when the read depth at that site in that sample falls between the parameters `site_depth` and `site_ratio * genome_depth`. Otherwise, we won't include that sample for the compute of the pooled-SNPs for that site.

  3. **relevant site** across-samples site filters (core presets todo)

More details about the compute can be found at [Cross-Sample SNP Analysis Tools (xsnp)](https://github.com/czbiohub/xsnp/wiki/Data-Schema-And-Computation)


### Sample commands

- Default parameters

   ```
   midas2 merge_snps --samples_list /path/to/sample/lists --num_cores 32 ${merged_midas_outdir}
   ```

### Output files

- `snps_summary.tsv`

   ```
   sample_name species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered  mean_coverage
   sample1     102470      5774847        4379800        23295116     213486         202365        0.758             5.319
   sample1     100039      2852528        1764680        25765375     225497         224630        0.619             14.601
   ```

- `{species_id}/{species_id}.snps_info.tsv.lz4`

   ```
   site_id                          major_allele  minor_allele  sample_counts  snp_type  rc_A  rc_C  rc_G  rc_T  sc_A  sc_C  sc_G  sc_T  locus_type  gene_id           site_type  amino_acids
   gnl|Prokka|UHGG143484_1|810|C    C             T             2              bi        0     10     0     2     0     2    0      2    CDS         UHGG143484_00002  4D         G,G,G,G
   gnl|Prokka|UHGG143484_1|905|T    C             T             2              bi        0     2      0     16    0     2    0      2    CDS         UHGG143484_00002  1D         Y,S,C,F
   ```

- `{species_id}/{species_id}.snps_freq.tsv.lz4`

  ```
  site_id                           sample1   sample3
  gnl|Prokka|UHGG143484_1|810|C	    0.167     0.167
  gnl|Prokka|UHGG143484_1|905|T     0.889     0.889
  ```

- `{species_id}/{species_id}.snps_depth.tsv.lz4`: site-by-sample number of mapped reads, **only** accounts for reads matching either major or minor allele

  ```
  site_id                           sample1  sample3
  gnl|Prokka|UHGG143484_1|810|C	    6        6
  gnl|Prokka|UHGG143484_1|905|T	    9	     9
  ```

Refer to [MIDAS's merge SNPs](https://github.com/snayfach/MIDAS/blob/master/docs/merge_snvs.md) for more details.


## Pan-gene Copy Number


### Sample command

   ```
   midas2 merge_genes --samples_list /path/to/tsv --num_cores 8 ${merged_midas_outdir}
   ```

### Target output files

- `genes_summary.tsv`

   ```
   sample_name  species_id  pangenome_size  covered_genes  fraction_covered  mean_coverage  aligned_reads  mapped_reads  marker_depth
   sample1      100247      87136           6122           0.070             9.519          455045         269656        3.447
   sample3      100247      57785           1889           0.033             4.849          108041         30072         0.000
   ```

- `{species_id}/{species_id}.genes_copynum.tsv.lz4`: copy number

  ```
  gene_id            sample1     sample3
  UHGG000947_00080   0.640716    0.640716
  UHGG000947_00081   1.474722    1.474722
  ```

- `{species_id}/{species_id}.genes_preabs.tsv.lz4`: presence/absence
 
  ```
  gene_id             sample1     sample3
  UHGG000947_00080    1           1
  UHGG000947_00081    1           1
  ```

- `{species_id}/{species_id}.genes_depth.tsv.lz4`: depth

  ```
  gene_id             sample1     sample3
  UHGG000947_00080    2.208723    2.208723
  UHGG000947_00081    5.083770    5.083770
  ```

Refer to [MIDAS's merge gene content](https://github.com/snayfach/MIDAS/blob/master/docs/merge_cnvs.md) for more details.
