# Metagenomic Intra-Species Diversity Analysis System 2.0 

MIDAS 2.0 is an integrated pipeline that estimate bacterial species abundance and strain-level genomic variation, including pan-gene content and SNPs analysis, using the UHGG dataset.

# midas_run_species: species abundance estimation

## input parameters

- `sample_name`: unique identifier for each each sample

- `{output_dir}/{sample_name}`: the base output directory of MIDAS 2.0 results

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

## input parameters

- `{output_dir}/{sample_name}/species/species_profile.txt`

- `species_cov`: select species present in the sample with minimal vertical coverage

- **todo**: add `species_id` option

## output files

- `{snps_output_dir}` := `{output_dir}/{sample_name}/snps`

- `{snps_output_dir}/output_sc.{species_cov}/{species_id}.snps.lz4`: per speices pileup results

   ```
   ref_id                          ref_pos ref_allele      depth   count_a count_c count_g count_t
   UHGG143505_C0_L5444.9k_H7fb7ad  44696   A               9       9       0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44697   A               9       9       0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44698   G               10      0       0       10      0
   UHGG143505_C0_L5444.9k_H7fb7ad  44699   A               10      10      0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44700   A               10      10      0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44701   C               10      0       10      0       0
   ```

- `{snps_output_dir}/output_sc.{species_cov}/summary.txt`: alignment stats for each species

   ```
   species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered mean_coverage
   102478      5444912        4526401        38190009     356763         273537        0.831            8.437
   ```

## temp files

- `{snps_output_dir}/temp_sc.{species_cov}/{species_id}/{genome}.fa: the representative genomes of  all the species with `species_cov`> specie_cov.

- `{snps_output_dir}/temp_sc.{species_cov}/repgenomes.fa`: the collated representative genome sequences of all the species with `species_cov`> specie_cov.

- `{snps_output_dir}/temp_sc.{species_cov}/repgenomes.{1..4, rev.1..2}.bt2`: the Bowtie2 index of the collected representative genomes

- `{snps_output_dir}/temp_sc.{species_cov}/repgenomes.{bam, bam.bai}`: the Bowtie2 alignment files of mapping reads to repgenomes.fa


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

The pooled-sample core-genome SNP calling pipeline can be broken down into the following steps:

- first pass: pool nucleotide variants from the input samples

- second pass: (based on the pooled data)
  determine if a genomic site is core: non-zero depth in >= 95% of samples (`site_prevalence`), and further on if the core genomic site is a SNP: major/minor alleles (`allele_frequency`)

## input files

The input TSV file include two columns: lists of samples and the midas_run_snps results path.

## species filters

- `sample_counts`: species with >= MIN_SAMPLES (1)

## sample filters

- `sample_depth`: minimum per-sample average read depth (5X)

- `sample_coverage`: [horizontal_coverage] fraction of reference genome sites covered by at least one read (40%)


## Site filters

### Per-sample site filter:

For each site of a given species, we only include <site, sample> pairs that passing the following filters:

```
- `site_depth`: minimum number of reads mapped to genomic site (2)

- `site_ratio`: maximum ratio of site depth over genome depth (5.0)
```

The number of samples passing the per-sample site filters is reported in `sample_counts` in the output `snps_info.tsv`.

### Across-sample site filters (core presets):

For the pairs of <site, samples passing per sample site filter>, we calculated the prevalence for each site and only report the core sites (defined by the `site_prevalence` cutoff in the report.

- `site_prevalence`: minimum fraction of samples where a genomic site pass the *site filters*

- `site_maf`: minimum average-minor-allele-frequency of site across samples for calling an allele present (0.1)

- `snp_type`: (for sites > SITE_MAF) specify one or more of the following 




- `allele_frequency`: (I don't think this is needed anymore)



## output files

- `merged_snps_outdir`:= `{merged_output_dir}/merged/snps`

- `snps_freq.tsv`: site-by-sample minor allele frequency matrix

- `snps_depth.tsv`: site-by-sample number of mapped reads, only accounts for reads matching either major or minor allele

- `snps_log.tsv`: 


