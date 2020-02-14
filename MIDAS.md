# Metagenomic Intra-Species Diversity Analysis Subcommands

The original [MIDAS tool](https://github.com/snayfach/MIDAS) is an integrated pipeline for estimating bacterial species abundance and strain-level genomic variation, including pan-gene content and SNPs analysis.   Its analyses steps are run against a database of 5,926 bacterial species extracted from 30,000 genomes.

The MIDAS subcommands in the IGGTOOLS package represent a reimplementation of the same analysis steps as the original [MIDAS tool](https://github.com/snayfach/MIDAS), but able to operate on the more comprehensive UHGG dataset, which consists of 4,644 gut-only species extracted from 286,997 genomes.

Similar to the original MIDAS tool, the IGGTOOLS MIDAS subcommands presuppose a database construction step has already taken place.  This construction step for IGGTOOLS and the UHGG dataset is documented [here](https://github.com/czbiohub/iggtools/wiki).  It was executed in AWS using hundreds of r5d.24xlarge instances over a period of a couple of days, depositing built products in S3.  The commands below implicitly reference the products of that build.  This page is focused specifically on the analysis steps, not the database construction steps.

# Single-sample result layout

For each sample, the analysis begins with a species profiling step.  The identified set of species that are abundant in the sample is then used to perform pan-genome analysis and SNP analysis.  The results of all three steps are laid out in the local filesystem as follows.

```
Output                                             Producer             Meaning
------------------------------------------------------------------------------------------------------------
{sample_name}/species/species_profile.tsv          midas_run_species    List of abundant species in sample
{sample_name}/snps/output/{species_id}.snps.lz4    midas_run_snps       Pileup results for each species_id
{sample_name}/snps/output/summary.txt              midas_run_snps       Summary of the SNPs analysis results
{sample_name}/genes/output/{species_id}.genes.lz4  midas_run_genes      Gene coverage per species_id
{sample_name}/genes/output/summary.txt             midas_run_genes      Pangenome alignment stats
```

Each sample analysis subcommand operates on a single sample.   It takes as a parameter the path to a unique output directory for that sample, which is the root of the layout above.

- `{output_dir}`: output directory unique for the sample, i.e., `{output_root}/{sample_name}`

The first subcommand to run for the sample is `midas_run_species`, and it will create that output directory if it does not exist.  All subsequent analysis steps operate within that directory.


# Pooled samples result layout

Results for multiple samples can be pooled using the corresponding subcommands `midas_merge_species`, `midas_merge_genes`, and `midas_merge_snps`.  The result layout for those looks as follows:

...

All merge subcommands take a list of single-sample output dirs as input.


# midas_run_species: species abundance estimation

## required parameters

-  one or two input fastqs for a single sample

- `{output_dir}`: output directory unique for the sample, i.e., `{output_root}/{sample_name}`;  this will be created if it does not exist

## output files

- **`{output_dir}/{sample_name}/species/species_profile.txt`**: tab-delimited with header, contains all the species with `species_cov` > 0, sorted by decreasing relative abundance.

   ```
   species_id    count_reads     coverage        relative_abundance
   102455        15053           137.635         0.130
   100044        10797           96.509          0.091
   101640        10688           84.354          0.080
   100089        9835            86.637          0.082
   102478        8603            78.545          0.074
   ```

## temp files

- `{output_dir}/{sample_name}/temp/alignment.m8`: alignment files of mapping reads to marker genes using hs-blastn

# midas_run_snps: single nucleotide polymorphism prediction

## required parameters

--  one or two input fastqs for a single sample

-- `{output_dir}`: must match the argument to run_midas_species for the same sample


## optional parameters

-- `species_cov`: select species present in the sample with minimal vertical coverage

-- `species_id`: limit analysys to a single species  **todo**

## inputs

- `{output_dir}/{sample_name}/species/species_profile.txt` --- as produced by run_midas_species

## output files

- `{sample_name}/snps/output/{species_id}.snps.lz4`: per species pileup results

   ```
   ref_id                          ref_pos ref_allele      depth   count_a count_c count_g count_t
   UHGG143505_C0_L5444.9k_H7fb7ad  44696   A               9       9       0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44697   A               9       9       0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44698   G               10      0       0       10      0
   UHGG143505_C0_L5444.9k_H7fb7ad  44699   A               10      10      0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44700   A               10      10      0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44701   C               10      0       10      0       0
   ```

- `{sample_name}/snps/output/summary.txt`: alignment stats for each species

   ```
   species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered mean_coverage
   102478      5444912        4526401        38190009     356763         273537        0.831            8.437
   ```

## temp files

- `{sample_name}/snps/output/temp_sc.{species_cov}/{species_id}/{genome}.fa: the representative genomes of  all the species with `species_cov`> specie_cov.

- `{sample_name}/snps/output/temp_sc.{species_cov}/repgenomes.fa`: the collated representative genome sequences of all the species with `species_cov`> specie_cov.

- `{sample_name}/snps/output/temp_sc.{species_cov}/repgenomes.{1..4, rev.1..2}.bt2`: the Bowtie2 index of the collected representative genomes

- `{sample_name}/snps/output/temp_sc.{species_cov}/repgenomes.{bam, bam.bai}`: the Bowtie2 alignment files of mapping reads to repgenomes.fa


# midas_run_genes: pan genome profiling

## input files

- `{output_dir}/{sample_name}/species/species_profile.txt`

- `species_cov`: select species present in the sample with minimal vertical coverage

- **todo**: add `species_id` option

## output files

- `{genes_output_dir}` := `{output_dir}/{sample_name}/genes`
- `{genes_output_dir}/output_sc.{species_cov}/{species_id}.genes.lz4`: 

   ```
   gene_id              count_reads     coverage        copy_number
   UHGG239769_04714     22              0.571323        0.000000
   UHGG050950_03155     7               0.182088        0.000000
   UHGG221050_00301     7               0.181818        0.000000
   UHGG175788_02287     3               0.078563        0.000000
   ```

- `{genes_output_dir}/output_sc.{species_cov}/summary.txt`: alignment stats for each species

   ```
   species_id  pangenome_size  covered_genes  fraction_covered  mean_coverage marker_coverage aligned_reads   mapped_reads
   102478      704500          145717         0.206837          1.212148      0.000000        1710757         1259775
   ```

## temp files

- `{genes_output_dir}/temp_sc.{species_cov}/{species_id}/centroids.ffn: the centroid genes of all the species with `species_cov`> specie_cov.

- `{genes_output_dir}/temp_sc.{species_cov}/pangenomes.fa`: the collated centroid genes of all the species with `species_cov`> specie_cov.

- `{genes_output_dir}/temp_sc.{species_cov}/pangenomes.{1..4, rev.1..2}.bt2`: the Bowtie2 index of the collected centroid pan genes

- `{genes_output_dir}/temp_sc.{species_cov}/pangenomes.{bam, bam.bai}`: the Bowtie2 alignment files of mapping reads to pan genes

# **merge_interface**

- `{input_files}`: map `sample_name` to its corresponding midas_output_dir

   ```
   sample_name    midas_output_dir
   SRS011134      /mnt/chunyu_6TB/iggtools-hmp-test/SRS011134
   SRS011271      s3://microbiome-chunyu/iggtools-hmp-test/SRS011271
   ```

- {merged_output_dir}: the output directory of merging MIDAS results

# midas_merge_species

## input files

- `sample_depth`: minimum per-sample marker-gene-depth (min_cov) for estimating species prevalence (1.0)

## output files

- `merged_species_outdir`:= `{merged_output_dir}/merged/species`

- `relative_abundance.tsv`: species-by-sample relative abundance matrix

   ```
   species_id   SRS011271   SRS011061   SRS011134  
   102455       0.091       0.130   0.000
   102549       0.000       0.011   0.049
   ```

- `count_reads.tsv`: species-by-samples read counts matrix

   ```
   species_id   SRS011271   SRS011061   SRS011134  
   102455       7352        50348       0.000
   102549       0.000       2           7125
   ```

- coverage.tsv`: species-by-samples genome coverage matrix

   ```
   species_id   SRS011271   SRS011061   SRS011134  
   102455       64.208      443.538     0.000
   102549       0.000       0.011       62.520
   ```

- `species_prevalence.tsv`: summary statistics for each species across samples; species are sorted by descending orders of median_abundance

   ```
   species_id median_abundance mean_abundance median_coverage mean_coverage  prevalence
   102293     0.049            0.134          64.208          169.249        2
   102181     0.034            0.023          38.567          27.530         2
   ```

# midas_merge_snps

The pooled-samples core-genome SNPs calling pipeline can be broken down into the following steps:

- Given list of samples that we are interested in merging the SNPs pileup, we select pairs of <species, sub-sample-list> based on horizontal/vertical pileup coverage, and species prevalence.

- For each species, accumulate A, C, G, T read counts, compute sample counts, and collect read counts for all the samples

- For the second pass, if a genomic site pass the site filters and snp types, then compute the major/minor allele based on accumulated read counts or sample counts (from previous step).

## input files

The input TSV file include two columns: list of samples and the midas_run_snps results path.

## Species and sample filters

Restrict to species that are covered sufficiently well in sufficiently many samples.

A species is covered "sufficiently well" by a sample when both of the following criteria are satisfied

```
- `genome_depth`: [vertical coverage]  per-sample average read depth (5X)

- `genome_coverage`: [horizontal coverage] fraction of reference genome sites covered by at least one read (40%)
```

```
- `sample_counts`: "sufficiently many" samples
```

## Site filters

### Per-sample site filter:

For each site of a given species, we only include <site, sample> pairs that passing the following filters:

```
- `site_depth`: minimum number of reads mapped to genomic site (2)

- `site_ratio`: maximum ratio of site depth over genome depth (5.0)
```

The number of samples passing the per-sample site filters is reported in `count_samples` in the output `snps_info.tsv`.

### Across-sample site filters (core presets):

For the pairs of <site, samples passing per sample site filter>, we calculated the prevalence for each site and only report the core sites (defined by the `site_prevalence` cutoff in the report.

```
- `site_prevalence`: minimum fraction of samples where a genomic site pass the *site filters*

- `site_maf`: minimum pooled minor-allele-frequency of site across samples for calling an allele present (0.1)

- `snp_type`: (for sites > SITE_MAF) specify one or more of the following 
```

- `allele_frequency`: (I don't think this is needed anymore) within sample min_allele_frequency to call a SNP


## output files

- `merged_snps_outdir`:= `{merged_output_dir}/merged/snps`

- 'snps_summary.tsv': 

```
species_id  sample_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered 
 mean_coverage
102293  SRS011271  3612475  2110215  770004185  14515249  9417166  0.584  364.894
102293  SRS011134  3612475  2706546  209982386  3175007  2683395  0.749  77.583
```

For each species, we generate three files:

- `snps_info.tsv`

```
site_id major_allele    minor_allele    sample_counts   snp_type        rc_A    rc_C    rc_G    rc_T    sc_A    sc_C    sc_G    sc_T
UHGG047905_C0_L562.0k_H31cf56|139|C     C       T       2       bi      0       25      0       20      0       1       0       1
UHGG047905_C0_L562.0k_H31cf56|162|C     C       A       2       bi      1       41      0       0       1       2       0       0
```

- `snps_freq.tsv`: site-by-sample minor allele frequency matrix

```
site_id SRS011271       SRS011134
UHGG047905_C0_L562.0k_H31cf56|139|C     1.000   0.000
UHGG047905_C0_L562.0k_H31cf56|162|C     0.000   0.043
```

- `snps_depth.tsv`: site-by-sample number of mapped reads, only accounts for reads matching either major or minor allele

```
site_id SRS011271       SRS011134
UHGG047905_C0_L562.0k_H31cf56|139|C     20      25
UHGG047905_C0_L562.0k_H31cf56|162|C     19      23
```

# midas_merge_genes

## Input files

- `sample_counts`
- `min_copy`
- `cluster_pid`

## Output files

- genes_summary.tsv.lz4

```
sample_id species_id pangenome_size covered_genes fraction_covered mean_coverage marker_coverage
```

For each species (passing the species and sample filter), there are three gene_id-by-sample matrices

- {species_id}/genes_copynum.tsv
- {species_id}/genes_preabs.tsv
- {species_id}/genes_coverage.tsv


- {sa