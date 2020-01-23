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

**todo**: add `species_id` option

## output files

- `{output_dir}/{sample_name}/snps/{output}_sc.{species_cov}/{species_id}.snps.lz4`: per speices pileup results

```
ref_id                          ref_pos ref_allele      depth   count_a count_c count_g count_t
UHGG143505_C0_L5444.9k_H7fb7ad  44696   A               9       9       0       0       0
UHGG143505_C0_L5444.9k_H7fb7ad  44697   A               9       9       0       0       0
UHGG143505_C0_L5444.9k_H7fb7ad  44698   G               10      0       0       10      0
UHGG143505_C0_L5444.9k_H7fb7ad  44699   A               10      10      0       0       0
UHGG143505_C0_L5444.9k_H7fb7ad  44700   A               10      10      0       0       0
UHGG143505_C0_L5444.9k_H7fb7ad  44701   C               10      0       10      0       0
```

- `{output_dir}/{sample_name}/snps/{output}_sc.{species_cov}/summary.txt`: alignment stats for each species

```
species_id      genome_length   covered_bases   total_depth     aligned_reads   mapped_reads    fraction_covered        mean_coverage
102478          5444912         4526401         38190009        356763          273537          0.831         8.437
```


## temp files



## midas_run_genes

## midas_merge_species



## midas_merge_snps

## midas_merge_genes

