# Metagenomic Intra-Species Diversity Analysis Subcommands

The original [MIDAS tool](https://github.com/snayfach/MIDAS) is an integrated pipeline for estimating bacterial species abundance and strain-level genomic variation, including pan-gene content and SNPs analysis.   Its analyses steps are run against a database of 5,926 bacterial species extracted from 30,000 genomes.

The MIDAS subcommands in the IGGTOOLS package represent a reimplementation of the same analysis steps as the original [MIDAS tool](https://github.com/snayfach/MIDAS), but able to operate on the more comprehensive UHGG dataset, which consists of 4,644 gut-only species extracted from 286,997 genomes.

Similar to the original MIDAS tool, the IGGTOOLS MIDAS subcommands presuppose a database construction step has already taken place.  This construction step for IGGTOOLS and the UHGG dataset is documented [here](https://github.com/czbiohub/iggtools/wiki).  It was executed in AWS using hundreds of r5d.24xlarge instances over a period of a couple of days, depositing built products in S3.  The commands below implicitly reference the products of that build.  This page is focused specifically on the analysis steps, not the database construction steps.

# Single-sample result layout

For each sample, the analysis begins with a species profiling step.  The identified set of species that are abundant in the sample is then used to perform pan-genome analysis and representative genome SNP analysis.  The results of all three steps are laid out in the local filesystem as follows.

```
Output                                          Producer            Meaning
------------------------------------------------------------------------------------------------------------
midas_iggdb                                     DB-related files    Mirror s3://miocriombe-igg/2.0/

{sample_name}/species/species_profile.tsv       midas_run_species   List of abundant species in sample

{sample_name}/snps/snps_summary.tsv             midas_run_snps      Summary of the SNPs analysis results
{sample_name}/snps/{species_id}.snps.tsv.lz4    midas_run_snps      Pileup results for each species_id

{sample_name}/genes/genes_summary.tsv           midas_run_genes     Pangenome alignment stats
{sample_name}/genes/{species_id}.genes.tsv.lz4  midas_run_genes     Gene coverage per species_id
```

Each sample analysis subcommand operates on a single sample. It takes as a parameter the path to MIDAS results root directory (`midas_outdir`) and a parameter for the sample name (`sample_name`); together constitute the unique output directory for that sample.

- `{output_dir}`: output directory unique for the sample, i.e., `{midas_outdir}/{sample_name}`

The first subcommand to run for the sample is `midas_run_species`, and it will create that output directory `output_dir` if it does not exist.  All subsequent analysis steps operate within that directory.

# Single-sample analysis

Multiple steps analysis happen for each sample, which usually happen in the following order.

```
{output_dir}
 |- species
 |- snps
 |- genes
 |- temp
 |  |- snps/repgenomes.bam
 |  |- genes/pangenome.bam
 |- bt2_indexes
 |  |- snps/repgenomes.*
 |  |- genes/pangenomes.*
```

## Species abundance estimation

### Example command

  ```
  python -m iggtools midas_run_species \
         --sample_name ${sample_name} -1 /path/to/R1 -2 /path/to/R2 \
         --midas_iggdb /path/to/local/midas/iggdb --num_cores 4 --debug \
         $midas_outdir
  ```

### Target output files

- `species_profile.tsv`: species present in the sample (`coverage` > 0), sorted in decreasing relative abundance. 

   ```
   species_id    read_counts     coverage        rel_abundance
   102455        15053           137.635         0.130
   100044        10797           96.509          0.091
   ```

Refer to [original MIDAS's estimate species abundance](https://github.com/snayfach/MIDAS/blob/master/docs/species.md) for details.


## Single nucleotide polymorphisms calling

To explore within-species variations for the species present in the sample data, metagenomics shotgun reads were aligned to a Bowtie2 indexes of a collection of representative genomes; nucleotide variation for each genomic site was quantify via pileup and alleles count. 

- **Bowtie2 indexes options**

  In addition to the original MIDAS's approach of on-the-fly build the Bowtie2 database for species in the restricted species profile, uses can also provide a prebuilt Bowtie2 indexes, e.g. one Bowtie2 database for all the samples in one study. Given the following three considerations:

    1. Despite the limitation to only abundant species to each sample, build bowtie2 indexes still takes significant amount of CPU time.
    2. For samples from the same study, the microbiome compositions are more or less similar.
    3. Given the multi-mapped reads issues between similar genomes in the Bowtie2 database, per-sample Bowtie2 indexes may cause bias for the pooled-sample core-genome SNP calling.

 Note: we still parse the pileup for abundance species present in each sample, selected by `genome_coverage`.

- **Chunk-of-sites** 

  One chunk-of-sites is indexed by `species_id, chunk_id`, represents `contig_end - contig_start` number of sites for one `contig_id`. For each chunk, pileup counts, mapped reads and aligned reads were computed independently. Chunks were computed in parallel. When all chunks from the same species are processed, the per-chunk pileup results are merged into one pileup file per species (`{species_id}.snps.tsv.lz4`) which avoid the explosion of intermediate files.

- **Q & A**

  1. What if I want to choose a different set representative genomes for given species.
     Ans: you can modify the table-of-content `genomes.tsv` accordingly.

  2. Can I do variant calling for genomes/species not in the UHGG, using the `midas_run_snps`script?
     Ans: Yes. User can build a `mock-midas-iggdb` by following the same file structures. For example, given a set of genomes that you want to perform SNP calling. 

  First, generate the `${mock_midas_iggdb}/genomes.tsv` in the desired format:

  ```
  genome            species representative          genome_is_representative
  GUT_GENOME091053  100001  GUT_GENOME091053        1
  GUT_GENOME212267  100002  GUT_GENOME212267        1
  GUT_GENOME178957  100003  GUT_GENOME178957        1
  ```

  Secondly, set up the representative genomes in the following structure:

  ```
  # cleaned_imports/{species_id}/{genome_id}/{genome_id}.fna
  ${mock_midas_iggdb}/cleaned_imports/100001/GUT_GENOME091053/GUT_GENOME091053.fna.lz4
  ${mock_midas_iggdb}/cleaned_imports/100002/GUT_GENOME212267/GUT_GENOME212267.fna.lz4
  ${mock_midas_iggdb}/cleaned_imports/100003/GUT_GENOME212267/GUT_GENOME212267.fna.lz4
  ```

  Then provide the `${mock_midas_iggdb}` to midas_run_snps by --midas_iggdb.  

  3. When provide midas_run_snps and midas_run_genes with existing Bowtie2 indexes (`--prebuilt_bowtie2_indexes`), MIDAS also needs to know what species were being included during the database build step (`--prebuilt_bowtie2_species`).  
MIDAS will only perform pileup on abundant species that passing the `--marker_depth`. When the provided species of interest (`--species_list`) don't pass the marker_depth filter, or not present in the prebuilt_bowtie2_indexes, then MIDAS won't perform pileup analysis on those species, even if it is provided by `--species_list`.

  
### Example command

- Perform Pileup for all the species with mean vertical genome coverage higher than 10X (need to run midas_run_species before). The following would happen: (1) build per-sample Bowtie2 indexes with all the species (`marker_depth` > 10) (2) align reads to per-sample Bowtie2 indexes (3) pileup and SNPs calling

   ```
   python -m iggtools midas_run_snps --sample_name ${sample_name} \
          -1 $R1 -2 $R2 --marker_depth 10 --num_cores 8 \
          /path/to/midas/output
   ```

- Perform Pileup for all the species with mean vertical genome coverage higher than 10X, and using existing Bowtie2 databases.

  ```
  python -m iggtools midas_run_snps --sample_name ${sample_name} \
         -1 $R1 -2 $R2 \
         --prebuilt_bowtie2_indexes ${bt2_indexes/${bt2_name} \
         --prebuilt_bowtie2_species ${bt2_indexes/${bt2_name}.species \
         --marker_depth 10 --num_cores 8 \
         /path/to/midas/output
   ```

- Special usage of `marker_depth`

  1. `--marker_depth=-1`: skip running species flow and perform pileup for all the species in the Bowtie2 indexes; use with caution, only when you know exactly what you want to do. Can be used with 

  2. `--marker_depth 0`: perform pileup for all the species present in the given sample

- Variant calling for custom mock-midas-iggdb

   ```
   python -m iggtools midas_run_snps --sample_name ${sample_name} \
          -1 $R1 -2 $R2 \
          --midas_iggdb /path/to/mock/midas/iggdb \
          --prebuilt_bowtie2_indexes ${bt2_indexes/${bt2_name} \
          --prebuilt_bowtie2_species ${bt2_indexes/${bt2_name}.species \
          --marker_depth=-1 --num_cores 8 \
          /path/to/midas/output
   ``` 

### Target output files

- `snps_summary.tsv`

   ```
   species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered mean_coverage
   102478      5444912        4526401        38190009     356763         273537        0.831            8.437
   ```

- `{species_id}.snps.tsv.lz4`: per-species pileup results

   ```
   ref_id                          ref_pos ref_allele      depth   count_a count_c count_g count_t
   UHGG143505_C0_L5444.9k_H7fb7ad  44696   A               9       9       0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44697   A               9       9       0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44698   G               10      0       0       10      0
   ```

Refer to [MIDAS's call single nucleotide polymorphisms](https://github.com/snayfach/MIDAS/blob/master/docs/cnvs.md) for more details.


## Pangenome profiling

To quantify the pangenome genes presence/absence for the species(es) of interest in the mNGS data, reads were aligned to the all the centroid_99 genes per species, and gene copy number are estimated.

Similar to above snps flow, genes flow also adopted the prebuilt **Bowtie2 database** and **chunk-of-genes** compute schema. The hierarchical compute for each chunk-of-genes start with: (1) for each gene, aligned reads, mapped reads, read depths and gene length are computed; (2) (after collecting reads alignment information) for all the genes for one species, compute the read depths of the 15 single copy marker genes; (3) for each gene, infer the copy number. 

### Example command

- Profile pangenome for given species, if their mean vertical genome coverage higher than 10X (need to run midas_run_species before). 

   ```
   python -m iggtools midas_run_genes --sample_name ${sample_name} \
          -1 $R1 -2 $R2 --marker_depth 10 --num_cores 8 \
          --species_list 100044,101302,102478 \
          /path/to/midas/output
   ```

### Target output files

- `genes_summary.tsv`

   ```
   species_id  pangenome_size  covered_genes  fraction_covered  mean_coverage marker_coverage aligned_reads   mapped_reads
   102478      704500          145717         0.206837          1.212148      0.000000        1710757         1259775
   ```

- `{species_id}.genes.tsv.lz4`: 

   ```
   gene_id              count_reads     coverage        copy_number
   UHGG239769_04714     22              0.571323        0.000000
   UHGG050950_03155     7               0.182088        0.000000
   ```

Refer to [MIDAS's predict pan-genome gene content](https://github.com/snayfach/MIDAS/blob/master/docs/cnvs.md) for more details.


# Pooled samples analysis

To merge MIDAS single sample results across multiple samples, e.g. generate species abundance matrix across samples, or perform pooled-samples core-genome SNPs calling, or merge pan-genome profiling results across samples, users need to provide a `--samples_list` TSV file, which specifies the `sample_name` and `midas_outdir`:

   ```
   sample_name    midas_outdir
   SRS011134      /mnt/chunyu_6TB/iggtools-hmp-test
   SRS011271      s3://microbiome-chunyu/iggtools-hmp-test
   ```

For example, MIDAS expects to locate `/mnt/chunyu_6TB/iggtools-hmp-test/SRS011134/species/species_profile.tsv`, generated from single-sample species flow. The merged output files are structured as following:

   ```
   ${midas_outdir}
    |- species
    |- snps
    |- genes
    |- temp
    |  |- snps
    |- [bt2_indexes]
   ```

- `{midas_outdir}`: output directory is provided by the user, which is the root of the layout above and below.


# Pooled samples results layout

Results for multiple samples can be pooled using the corresponding subcommands `midas_merge_species`, `midas_merge_genes`, and `midas_merge_snps`.  The result layout for those looks as follows:

```
Output                                      Producer             Meaning
-----------------------------------------------------------------------------------------------------------
species/species_prevalence.tsv              midas_merge_species  Summary statistics per species across samples
species/species_read_counts.tsv             midas_merge_species  Species-by-sample read counts matrix
species/species_coverage.tsv                midas_merge_species  Species-by-sample genome coverage matrix
species/species_rel_abundance.tsv           midas_merge_species  Species-by-sample relative abundance matrix


snps/snps_summary.tsv                        midas_merge_snps     Alignment summary statistics per sample
snps/{sp_id}/{sp_id}.snps_info.tsv.lz4       midas_merge_snps     Per species metadata for genomic sites.
snps/{sp_id}/{sp_id}.snps_freqs.tsv.lz4      midas_merge_snps     Per species site-by-sample MAF matrix
snps/{sp_id}/{sp_id}.snps_depth.tsv.lz4      midas_merge_snps     Per species site-by-sample read depth matrix


genes/genes_summary.tsv                      midas_run_genes      Alignment summary statistics per sample
genes/{sp_id}/{sp_id}.genes_presabs.tsv.lz4  midas_run_genes      Per species gene-by-sample pre-abs matrix
genes/{sp_id}/{sp_id}.genes_copynum.tsv.lz4  midas_run_genes      Per species gene-by-sample copy number matrix
genes/{sp_id}/{sp_id}.genes_depth.tsv.lz4    midas_run_genes      Per species gene-by-sample read depth matrix
```

All merge subcommands take a `samples_list` input argument, which is a TSV file with sample name and single-sample unique output directory.


## Merge species abundance profile

## Example command

User can chose to build `local_bowtie2_indexes` given the merged species profile across samples.

   ```
   python -m iggtools midas_merge_species --samples_list /path/to/sample/lists ${midas_outdir}
   ```

### Target output files

- `species_prevalence.tsv`: summary statistics for each species across samples

   ```
   species_id  median_abundance  mean_abundance  median_coverage  mean_coverage  prevalence
   102293      0.049             0.134           64.208           169.249        2
   102181      0.034             0.023           38.567           27.530         2
   ```

- `species_rel_abundance.tsv`

   ```
   species_id   SRS011271   SRS011061   SRS011134  
   102455       0.091       0.130   0.000
   102549       0.000       0.011   0.049
   ```

- `species_read_counts.tsv`: species-by-samples read counts matrix

   ```
   species_id   SRS011271   SRS011061   SRS011134  
   102455       7352        50348       0.000
   102549       0.000       2           7125
   ```

- `species_coverage.tsv`: species-by-samples genome coverage matrix

   ```
   species_id   SRS011271   SRS011061   SRS011134  
   102455       64.208      443.538     0.000
   102549       0.000       0.011       62.520
   ```

Refer to [MIDAS's merge species abundance](https://github.com/snayfach/MIDAS/blob/master/docs/merge_species.md) for more details.


## Build Bowtie2 indexes with merged species profile

When the microbiome composition for a collection of samples are similar, it is more efficient to build one Bowtie2 indexes with species which are either prevalent (`saample_counts`) or abundant (`mean_coverage`), after running `midas_merge_species`.

Example command to build Bowtie2 database with species present in at least 2 sample:

   ```
   python -m iggtools build_bowtie2_indexes --midas_iggdb /path/to/local/midas/iggdb \
          --species_profile ${merged_midas_outdir}/species/species_prevalence.tsv \
          --select_by sample_counts --select_threshold 2 \
          --num_cores 8 /path/to/bt2_indexes
   ```

## Pooled-samples core-genome SNPs calling

To compute the across-samples pooled-SNPs for each genomic site in the representative genome, `read_counts` and `sample_counts` of A,C,G,T are accumulated across all samples. Then we can compute the across-samples-major-alleles for each site, following by collecting the corresponding read depth and allele frequency for each sample. 

For each relevant site, we determine the set of alleles present for that site across all relevant samples. For each allele A, C, G, T we count the samples that are relevant for the site and contain that allele, and sum that allele's depths across those samples.  We use one or the other of these metrics as a proxy for allele frequency, as specified by the `pool_method` argument.

- **pipeline details**

  1. <species, sub-sample-lists> selection: the analysis restricts attention to "sufficiently well" covered species in sufficiently many samples. (`genome_depth`, `genome_coverage`, and `sample_counts`).

  2. **relevant sample**: per-site sample filters: a sample is considered relevant for a given genomic site when the read depth at that site in that sample falls between the parameters `site_depth` and `site_ratio * genome_depth`. Otherwise, we won't include that sample for the compute of the pooled-SNPs for that site.

  3. **relevant site** across-samples site filters (core presets)

More details about the compute can be found at [Cross-Sample SNP Analysis Tools (xsnp)](https://github.com/czbiohub/xsnp/wiki/Data-Schema-And-Computation)

### Example commands

- Default parameters

   ```
   python -m iggtools midas_merge_snps --samples_list /path/to/tsv --num_cores 8 path/to/merged/midas/outdir
   ```

### Target output files

- `snps_summary.tsv`

   ```
   species_id sample_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered  mean_coverage
   102293     SRS011271  3612475        2110215        770004185    14515249       9417166       0.584             364.894
   102293     SRS011134  3612475        2706546        209982386    3175007        2683395       0.749             77.583
   ```

- `{species_id}/snps_info.tsv`

   ```
   site_id                              major_allele  minor_allele sample_counts snp_type  rc_A    rc_C    rc_G    rc_T    sc_A    sc_C    sc_G    sc_T
   UHGG047905_C0_L562.0k_H31cf56|139|C  C             T            2             bi        0       25      0       20      0       1       0       1
   UHGG047905_C0_L562.0k_H31cf56|162|C  C             A            2             bi        1       41      0       0       1       2       0       0
   ```

- `{species_id}/snps_freq.tsv`

  ```
  site_id                                 SRS011271   SRS011134
  UHGG047905_C0_L562.0k_H31cf56|139|C     1.000       0.000
  UHGG047905_C0_L562.0k_H31cf56|162|C     0.000       0.043
  ```

- `snps_depth.tsv`: site-by-sample number of mapped reads, only accounts for reads matching either major or minor allele

  ```
  site_id                                 SRS011271   SRS011134
  UHGG047905_C0_L562.0k_H31cf56|139|C     20          25
  UHGG047905_C0_L562.0k_H31cf56|162|C     19          23
  ```

Refer to [MIDAS's merge SNPs](https://github.com/snayfach/MIDAS/blob/master/docs/merge_snvs.md) for more details.


## Merge results from pangenome profiling

Considering the number of pan-genome genes is relatively smaller than the genome size, we merge results from pan-genome profiling across samples for all genes per species.

### Example command

   ```
   python -m iggtools midas_merge_genes --samples_list /path/to/tsv --num_cores 8 /path/to/merged/midas/outdir
   ```

### Target output files

- `genes_summary.tsv`

   ```
   species_id  sample_name  pangenome_size  covered_genes  fraction_covered  mean_coverage  aligned_reads mapped_reads marker_depth
   102478      SRS011061    704500          241525         0.343             3.738          8590041       7841889      0.000
   102478      SRS011134    704500          312336         0.443             2.144          6255805       5203378      0.000
   ```

- `{species_id}.genes_copynum.tsv.lz4`

  ```
  gene_id            SRS011061  SRS011134
  UHGG000186_01791   1.099      1.099
  UHGG120544_00344   0.425      0.638
  ```

- `{species_id}.genes_preabs.tsv.lz4`
 
  ```
  gene_id             SRS011061  SRS011134
  UHGG000186_01791    1          1
  UHGG120544_00344    0          1
  ```

- `{species_id}.genes_coverage.tsv.lz4`

  ```
  gene_id            SRS011061  SRS011134
  UHGG000186_01791   8.854      0.000
  UHGG120544_00344   26.730     12.037
  ```

Refer to [MIDAS's merge gene content](https://github.com/snayfach/MIDAS/blob/master/docs/merge_cnvs.md) for more details.
