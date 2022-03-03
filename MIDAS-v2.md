# Metagenomic Intra-Species Diversity Analysis Subcommands

[MIDAS](https://genome.cshlp.org/content/26/11/1612) is an integrated pipeline for profiling strain-level genomic and functional variation for metagenomics data. Its analyses steps are run against a database of 5,926 bacterial species extracted from 30,000 genomes (MIDAS DB v1.2).

The MIDAS subcommands in the IGGTOOLS package represent a reimplementation of the same analysis steps as the original [MIDAS tool](https://github.com/snayfach/MIDAS), but able to operate on the more comprehensive MIDAS DB, in a fast and scalable manner.


Similar to the original MIDAS tool, the IGGTOOLS MIDAS subcommands presuppose a database construction step has already taken place. The construction step for the [UHGG1.0 catalogue](https://www.ebi.ac.uk/metagenomics/genome-catalogues/human-gut-v2-0), which consists of 4,644 gut-only species extracted from 286,997 genomes, is documented [here](https://github.com/czbiohub/iggtools/wiki). It was executed in AWS using hundreds of r5d.24xlarge instances over a period of a couple of days, depositing built products in S3.  The commands below implicitly reference the products of that build.  This page is focused specifically on the analysis steps, not the database construction steps.  MIDAS users can build custom databases for any collections of genomes.


# Single-sample result layout

For each sample, the analysis begins with a species profiling step.  The identified set of species that are abundant in the sample is then used to perform pan-genome analysis and representative genome SNP analysis.  The results of all three steps are laid out in the local filesystem as follows.

```
Output                                          Producer            Meaning
------------------------------------------------------------------------------------------------------------
midas_iggdb                                     DB-related files    Mirror s3://miocriombe-igg/2.0/

{sample_name}/species/species_profile.tsv       midas_run_species   List of species present in sample
{sample_name}/species/markers_profile.tsv       midas_run_species   Marker abundance for each species_id

{sample_name}/snps/snps_summary.tsv             midas_run_snps      Summary of the SNPs analysis results
{sample_name}/snps/{species_id}.snps.tsv.lz4    midas_run_snps      Pileup results for each species_id

{sample_name}/genes/genes_summary.tsv           midas_run_genes     Pangenome alignment stats
{sample_name}/genes/{species_id}.genes.tsv.lz4  midas_run_genes     Gene coverage per species_id
```

Each sample analysis subcommand operates on a single sample. It takes as a parameter the path to MIDAS results root directory (`midas_outdir`) and a parameter for the sample name (`sample_name`); together constitute the unique output directory for that sample.

- `{output_dir}`: output directory unique for the sample, i.e., `{midas_outdir}/{sample_name}`

The first subcommand to run for the sample is `midas_run_species`, to report species present in the sample that we can genotype in the `midas_run_snps` flow.. This command will create the output directory `output_dir` if it does not exist.  All subsequent analysis steps operate within that directory.


# Single-sample analysis

Multiple steps analysis happen for each sample, which usually happen in the following order.

```
{output_dir}
 |- species
 |- snps
    |- snps_summary.tsv 
    |- {species_id}.snps.tsv.lz4
 |- genes
    |- genes_summary.tsv 
    |- {species_id}.genes.tsv.lz4
 |- temp
    |- snps
       |- repgenomes.bam
       |- {species_id}/snps_XX.tsv
    |- genes
       |- pangenome.bam
       |- {species_id}/genes_XX.tsv
 |- bt2_indexes
    |- snps/repgenomes.*
    |- genes/pangenomes.*
```

## Species abundance estimation

### Example command

  ```
  python -m iggtools midas_run_species \
         --sample_name ${my_sample} -1 /path/to/R1 -2 /path/to/R2 \
         --midas_iggdb /path/to/local/midas/iggdb --num_cores 8 --debug \
         ${midas_outdir}
  ```

### Species flow output files

The goal of the species flow is to detect abundant species present in the sample, which can be genotyped later. We mapped the raw reads to 15 single copy Phyeco marker genes using `hs-blastn`. We computed the uniquely mapped read counts for each marker genes, and probabilistically assign the ambiguous reads to marker genes. The coverage of each marker gene is computed as the total alignment length over the marker gene length. Only species with more than marker genes covered by more than two reads are reported. We further call a species is **present** in the sample if the `median_marker_coverage` > 0.

- `species_profile.tsv`: sorted in decreasing order of `median_coverage`. 

   ```
   species_id  read_counts  median_coverage  coverage  relative_abundance total_covered_marker  unique_covered_marker  ambiguous_covered_marker  total_marker_length
   102455      15053        120.05           137.635    0.130             15                    14                     14                        14
   100044      10797        91.20            96.509     0.091             14                    14                     12
   ```

- `markers_profile.tsv`: species-by-marker 


   ```
   species_id  marker_id  marker_length  gene_id          total_reads  total_alnbps  coverage  uniq_reads  ambi_reads   uniq_alnbps     ambi_alnbps
   102455      B000032    1052           195103.peg.1451  1680         2000          11.4      8200        1800         9680            2320
   ```

Refer to [original MIDAS's estimate species abundance](https://github.com/snayfach/MIDAS/blob/master/docs/species.md) for more details.


## Single nucleotide polymorphisms calling

To explore within-species variations for the species present in the sample data, metagenomics shotgun reads were aligned to a Bowtie2 indexes of a collection of representative genomes; nucleotide variation for each genomic site was quantify via pileup and alleles count. 

- **Bowtie2 indexes options**

  In addition to the original MIDAS's approach of on-the-fly build the Bowtie2 database for species in the restricted species profile, uses can also provide a prebuilt Bowtie2 indexes, e.g. one Bowtie2 database for all the samples in one study, or a Bowtie2 database for one particular species across samples. Given the following three considerations:

    1. Despite the limitation to only abundant species to each sample, build bowtie2 indexes still takes significant amount of CPU time.
    2. For samples from the same study, the microbiome compositions are more or less similar.
    3. Given the multi-mapped reads issues between similar genomes in the Bowtie2 database, per-sample Bowtie2 indexes may cause bias for the pooled-sample core-genome SNP calling.

  Despite the species present in varying Bowtie2 database, we only parse the Pileup results for abundant species in the sample, selected by `genome_coverage`.

- **Chunk-of-sites** 

  One chunk-of-sites is indexed by `species_id, chunk_id`, represents `contig_end - contig_start` number of sites for one `contig_id`. For each chunk, pileup counts, mapped reads and aligned reads were computed independently. Chunks were computed in parallel. When all chunks from the same species are processed, the per-chunk pileup results are merged into one pileup file per species (`{species_id}.snps.tsv.lz4`) which avoid the explosion of intermediate files.

- **Q & A**

  1. What if I want to choose a different set representative genomes for given species.
     
     Ans: you can modify the table-of-content `genomes.tsv` accordingly.

  2. Can I do variant calling for genomes/species not in the UHGG, using the `midas_run_snps`script?
     
     Ans: Yes. User can build a `custom_midas_iggdb` by following the same file structures. For example, given a set of genomes that you want to perform SNP calling. 

     First, generate the `${custom_midas_iggdb}/genomes.tsv` in the desired format:

     ```
     genome            species representative          genome_is_representative
     GUT_GENOME091053  100001  GUT_GENOME091053        1
     GUT_GENOME212267  100002  GUT_GENOME212267        1
     GUT_GENOME178957  100003  GUT_GENOME178957        1
     ```

     Secondly, set up the representative genomes in the following structure:

     ```
     ${custom_midas_iggdb}/cleaned_imports/{species_id}/{genome_id}/{genome_id}.fna
     ${custom_midas_iggdb}/cleaned_imports/100001/GUT_GENOME091053/GUT_GENOME091053.fna.lz4
     ${custom_midas_iggdb}/cleaned_imports/100002/GUT_GENOME212267/GUT_GENOME212267.fna.lz4
     ${custom_midas_iggdb}/cleaned_imports/100003/GUT_GENOME212267/GUT_GENOME212267.fna.lz4
     ```

     Then provide the `${custom_midas_iggdb}` to `midas_run_snps` by `--midas_iggdb`.  

  3. When provide `midas_run_snps` and `midas_run_genes` with existing Bowtie2 indexes (`--prebuilt_bowtie2_indexes`), users also need to provide `--prebuilt_bowtie2_species` to specify what species were being included during the database build step.  

  4. MIDAS only perform pileup-based variant calling on abundant species that passing the `--marker_depth`. When the provided species of interest (`--species_list`) don't pass the marker_depth filter, and / or not present in the prebuilt_bowtie2_indexes, then MIDAS won't perform SNPs calling on those species.

  
### Example command

- Perform Pileup for all the species with average vertical genome coverage (todo a more accurate way should be average sgc marker gene vertical coverage) higher than 10X (need to run `midas_run_species beforehand`). 

   ```
   python -m iggtools midas_run_species --sample_name ${sample_name} \
          -1 ${R1} -2 ${R2} --num_cores 4 ${midas_outdir}

   python -m iggtools midas_run_snps --sample_name ${sample_name} \
          -1 ${R1} -2 ${R2} --marker_depth 10 --num_cores 8 ${midas_outdir}
   ```

   The following would happen: (1) build per-sample Bowtie2 indexes with all the species (`marker_depth` > 10) (2) align reads to per-sample Bowtie2 indexes (3) pileup and SNPs calling


- Perform Pileup for all the species with average vertical genome coverage higher than 10X, and using existing Bowtie2 databases.
   
   ```
   python -m iggtools midas_run_snps --sample_name ${sample_name} \
          -1 ${R1} -2 ${R2} --num_cores 8 --marker_depth 10 \
          --prebuilt_bowtie2_indexes ${bt2_indexes/${bt2_name} \
          --prebuilt_bowtie2_species ${bt2_indexes/${bt2_name}.species \
          ${midas_outdir}
    ```

- Special usage of `marker_depth`

  1. `--marker_depth=-1`: skip running species flow and perform pileup for all the species in the Bowtie2 indexes. Use with caution, only when you know exactly what you want to do.  

  2. `--marker_depth 0`: perform pileup for all the species present in the given sample, computed from the `midas_run_species` flow.


- Variant calling for custom-midas-iggdb

   ```
   python -m iggtools midas_run_snps --sample_name ${sample_name} \
          -1 $R1 -2 $R2 \
          --midas_iggdb ${custom_midas_iggdb} \
          --prebuilt_bowtie2_indexes ${bt2_indexes/${bt2_name} \
          --prebuilt_bowtie2_species ${bt2_indexes/${bt2_name}.species \
          --marker_depth=-1 --num_cores 8 ${midas_outdir}
   ``` 

### SNPs flow output files

- Pileup summary: `snps/snps_summary.tsv`

   ```
   species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered mean_coverage
   102478      5444912        4526401        38190009     356763         273537        0.831            8.437
   ```

- Per-species pileup results: `snps/{species_id}.snps.tsv.lz4`

   ```
   ref_id                          ref_pos ref_allele      depth   count_a count_c count_g count_t
   UHGG143505_C0_L5444.9k_H7fb7ad  44696   A               9       9       0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44697   A               9       9       0       0       0
   UHGG143505_C0_L5444.9k_H7fb7ad  44698   G               10      0       0       10      0
   ```

Refer to [MIDAS's call single nucleotide polymorphisms](https://github.com/snayfach/MIDAS/blob/master/docs/cnvs.md) for more details.


## Pangenome profiling

To quantify the pangenome genes presence/absence for the species of interest in the shotgun metagenomics sequencing data, reads were aligned to all the centroid_99 genes per species, and gene copy number are estimated.

Similar to the above snps flow, genes flow also adopted the prebuilt **Bowtie2 database** and **chunk-of-genes** compute schema. The hierarchical compute for each chunk-of-genes start with: (1) For each gene, collect read alignments based metrics, such as `aligned_reads`, `mapped_reads`, `read_depths` and `gene_length` were computed; (2) For all the genes for given species, compute the average read depths of the 15 single copy marker genes; (3) For each gene, infer the copy number compared to the SGC marker gene vertical coverage. 

### Example command

- Profile pangenome for given species, if their average vertical genome coverage higher than 10X (need to run midas_run_species beforehand). 

   ```
   python -m iggtools midas_run_genes --sample_name ${sample_name} \
          -1 ${R1} -2 ${R2} --marker_depth 10 --num_cores 8 \
          --species_list 100044,101302,102478 \
          ${midas_outdir}
   ```

### Genes flow output files

- `genes/genes_summary.tsv`

   ```
   species_id  pangenome_size  covered_genes  fraction_covered  mean_coverage marker_coverage aligned_reads   mapped_reads
   102478      704500          145717         0.206837          1.212148      0.000000        1710757         1259775
   ```

- `genes/{species_id}.genes.tsv.lz4`: 

   ```
   gene_id              count_reads     coverage        copy_number
   UHGG239769_04714     22              0.571323        0.000000
   UHGG050950_03155     7               0.182088        0.100000
   ```

Refer to [MIDAS's predict pan-genome gene content](https://github.com/snayfach/MIDAS/blob/master/docs/cnvs.md) for more details.


# Pooled samples analysis

Iggtools MIDAS can merge single sample analysis results across multiple samples, e.g. generate species abundance matrix across samples, or perform pooled-samples core-genome SNPs calling, or merge pan-genome profiling results across samples. 

All three merge subcommands take a `samples_list` as input argument, which is a TSV file with `sample name` and single-sample unique output directory `midas_outdir`. For example, `midas_merge_species` expects to locate `/mnt/chunyu/iggtools-hmp-test/SRS011134/species/species_profile.tsv`, generated by single-sample species flow.

   ```
   sample_name   midas_outdir
   SRS011134     /mnt/chunyu/iggtools-hmp-test
   SRS011271     s3://microbiome-chunyu/iggtools-hmp-test
   ```

The merged output files are structured as following:

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


## Merge species abundance profile

## Example command

User can chose to build `local_bowtie2_indexes` given the merged species profile across samples.

   ```
   python -m iggtools midas_merge_species --samples_list /path/to/sample/lists ${midas_outdir}
   ```

### Output files

- `species/species_prevalence.tsv`: summary statistics for each species across samples

   ```
   species_id  median_abundance  mean_abundance  median_coverage  mean_coverage  prevalence
   102293      0.049             0.134           64.208           169.249        2
   102181      0.034             0.023           38.567           27.530         2
   ```

- `species/species_rel_abundance.tsv`

   ```
   species_id   SRS011271   SRS011061   SRS011134  
   102455       0.091       0.130   0.000
   102549       0.000       0.011   0.049
   ```

- `species/species_read_counts.tsv`: species-by-samples read counts matrix

   ```
   species_id   SRS011271   SRS011061   SRS011134  
   102455       7352        50348       0.000
   102549       0.000       2           7125
   ```

- `species/species_coverage.tsv`: species-by-samples genome coverage matrix

   ```
   species_id   SRS011271   SRS011061   SRS011134  
   102455       64.208      443.538     0.000
   102549       0.000       0.011       62.520
   ```

Refer to [MIDAS's merge species abundance](https://github.com/snayfach/MIDAS/blob/master/docs/merge_species.md) for more details.


## Build Bowtie2 indexes with merged species profile

When the microbiome composition for a collection of samples are similar, it is more efficient to build one Bowtie2 indexes with species which are either prevalent (`sample_counts`) or abundant (`mean_coverage`), after running `midas_merge_species`, using the `build_bowtie2_indexes` subcommand.

Example command to build Bowtie2 database with species present in at least three samples:

   ```
   python -m iggtools build_bowtie2_indexes --midas_iggdb ${midas_iggdb} \
          --species_profile ${merged_midas_outdir}/species/species_prevalence.tsv \
          --select_by sample_counts --select_threshold 3 \
          --num_cores 8 /path/to/custom_bt2_indexes
   ```

## Pooled-samples core-genome SNPs calling

To compute the across-samples pooled SNPs for each genomic site in the representative genome, `read_counts` and `sample_counts` of A,C,G,T are accumulated across all samples. Then we computed the across-samples major alleles for each genomic site, following by collecting the corresponding read depth and allele frequency for each sample. 

For each **relevant** genomic site, we determine the set of alleles present for that site across all **relevant** samples. For each allele A, C, G, T we count the samples that are **relevant** for the site and contain that allele, and sum that allele's depths across those samples. We use one or the other of these metrics as a proxy for allele frequency, as specified by the `pool_method` argument.

- **pipeline details**

  1. <species, sub-sample-lists> selection: the analysis restricts attention to "sufficiently well" covered species in sufficiently many samples. (`genome_depth`, `genome_coverage`, and `sample_counts`). Therefore, despite the provided list of samples by the user, different species may have different lists of relevant samples.

  2. **relevant sample**: determined by per-site sample filters. A sample is considered relevant for a given genomic site when the read depth at that site in that sample falls between the parameters `site_depth` and `site_ratio * genome_depth`. Otherwise, we won't include that sample for the compute of the pooled-SNPs for that site.

  3. **relevant site** across-samples site filters (core presets todo)

More details about the compute can be found at [Cross-Sample SNP Analysis Tools (xsnp)](https://github.com/czbiohub/xsnp/wiki/Data-Schema-And-Computation)


### Example commands

- Default parameters

   ```
   python -m iggtools midas_merge_snps --samples_list /path/to/tsv --num_cores 8 path/to/merged/midas/outdir
   ```

### Output files

- `snps/snps_summary.tsv`

   ```
   species_id sample_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered  mean_coverage
   102293     SRS011271  3612475        2110215        770004185    14515249       9417166       0.584             364.894
   102293     SRS011134  3612475        2706546        209982386    3175007        2683395       0.749             77.583
   ```

- `snps/{species_id}/snps_info.tsv`

   ```
   site_id                              major_allele  minor_allele sample_counts snp_type  rc_A    rc_C    rc_G    rc_T    sc_A    sc_C    sc_G    sc_T
   UHGG047905_C0_L562.0k_H31cf56|139|C  C             T            2             bi        0       25      0       20      0       1       0       1
   UHGG047905_C0_L562.0k_H31cf56|162|C  C             A            2             bi        1       41      0       0       1       2       0       0
   ```

- `snps/{species_id}/snps_freq.tsv`

  ```
  site_id                                 SRS011271   SRS011134
  UHGG047905_C0_L562.0k_H31cf56|139|C     1.000       0.000
  UHGG047905_C0_L562.0k_H31cf56|162|C     0.000       0.043
  ```

- `snps/snps_depth.tsv`: site-by-sample number of mapped reads, only accounts for reads matching either major or minor allele

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

- `genes/genes_summary.tsv`

   ```
   species_id  sample_name  pangenome_size  covered_genes  fraction_covered  mean_coverage  aligned_reads mapped_reads marker_depth
   102478      SRS011061    704500          241525         0.343             3.738          8590041       7841889      0.000
   102478      SRS011134    704500          312336         0.443             2.144          6255805       5203378      0.000
   ```

- `genes/{species_id}.genes_copynum.tsv.lz4`

  ```
  gene_id            SRS011061  SRS011134
  UHGG000186_01791   1.099      1.099
  UHGG120544_00344   0.425      0.638
  ```

- `genes/{species_id}.genes_preabs.tsv.lz4`
 
  ```
  gene_id             SRS011061  SRS011134
  UHGG000186_01791    1          1
  UHGG120544_00344    0          1
  ```

- `genes/{species_id}.genes_coverage.tsv.lz4`

  ```
  gene_id            SRS011061  SRS011134
  UHGG000186_01791   8.854      0.000
  UHGG120544_00344   26.730     12.037
  ```

Refer to [MIDAS's merge gene content](https://github.com/snayfach/MIDAS/blob/master/docs/merge_cnvs.md) for more details.
