
# SNV Module: Population Single Nucleotide Variants Calling

The SNV module proceeds in two phages: (1) single-sample read pileup (2) population variants calling across all the species. The first steps can be potentially run in parallel.  We presuppose users already run the [database customization](https://github.com/czbiohub/MIDAS2.0/wiki/Data-customization) step, and either have single-sample `species_profile.tsv` or a prebuilt Bowtie2 rep-genome database ready for the SNV module.


## Single-Sample SNV Analysis

In a standard workflow, rep-genome database of the species in the restricted species profile were built for each sample, to which reads were aligned using Bowtie2. Per genomic site read pileup and nucleotide variation for **all** the species in the rep-genome database are reported by MIDAS 2.0.  

MIDAS 2.0 purposely holds any filter or species selection upon the single-sample pileup results until across-samples SNV analysis. That being said, read mapping summary is reported in `snps_summay.tsv`, and pileup/variants calling results are reported by species. Therefore users can easily customize species selection on their own.

### Sample commands

- Single-sample pileup for all the species in the restricted species profile: `median_marker_coverage > 2` and `unique_fraction_covered > 0.5`. 

  We presuppose users already [profiling the species coverage](https://github.com/czbiohub/MIDAS2.0/wiki/Data-customization#species-to-genotype), and expect `${my_midasdb_dir}/${sample_name}/species/species_profile.tsv` exists.

   ```
   midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
         --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
         --select_by median_marker_coverage,unique_fraction_covered \
         --select_threshold=2,0.5 \
         --num_cores 12 ${midas_outdir}         
    ```

  Users can adjust post-alignment quality filter parameters via the command line arguments, and the defaults are:

   ```
   --mapq >= 20: discard read alignments with alignment quality < 20
   --mapid >= 0.94: discard read alignments with alignment identity < 0.94
   --aln_readq >= 20: discard read alignment with mean quality < 20
   --aln_cov >= 0.75: discard read alignment with alignment coverage < 0.75
   --aln_baseq >= 30: discard bases with quality < 30
   ```

- Single-sample pileup for all the species in the restricted species profile with paired-ends based post-alignment quality filter.
  
  Users can recruit only properly paired reads for pileup, by passing the `--paired_only` with proper `--fragment_length`. In this case, the post-alignment metrics will be computed based on a read pair, instead of single read.

   ```
   midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
         --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
         --select_by median_marker_coverage,unique_fraction_covered \
         --select_threshold=2,0.5 \
         --fragment_length 2000 --paired_only \
         --num_cores 12 ${midas_outdir}         
    ```

- Single-sample variant calling for all the species in the restricted species profile with paired-ends based post-alignment quality filter.
  
  In recognition of the need for single-sample variant calling, we added an `advanced` mode to the single-sample SNV analysis in MIDAS 2.0. In the advanced mode, per-species pileup results will also report major allele and minor allele for all the genomic sites covered by at least two reads, upon which custom variant calling filter can be applied by the users. Users are advised to use the setting `--ignore_ambiguous` to avoid falsely calling major/minor alleles for sites with tied read counts.

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
  
  We presuppose users already followed the [database customization](https://github.com/czbiohub/MIDAS2.0/wiki/Data-customization#population-specific-species-panel) step, and have the prebuilt rep-genome database located at `${midas_outdir}/bt2_indexes`. 

   ```
   midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
         --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
          --prebuilt_bowtie2_indexes ${midas_output}/bt2_indexes/repgenome \
          --prebuilt_bowtie2_species ${midas_output}/bt2_indexes/repgenome.species \
          --select_threshold=-1 --num_cores 12 ${midas_output}
   ``` 

### Output files

- `snps_summary.tsv`: the statistics summary of read mapping and pileup for all the species in the rep-genome database

   ```
   species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered   mean_coverage
   102506      5339468        2373275        8045342      468667         224553        0.444              3.390
   102337      2749621        2566404        47723458     1479479        1010530       0.933              18.595
   ```
   - _genome_length_: genome length
   - _covered_bases_: number of bases covered by at least two reads
   - _total_depth_: total read depth across all _covered_bases_
   - _aligned_reads_: total read counts across sites before post-alignment filter
   - _mapped_reads_: total read counts across sites after post-alignment filter
   - _fraction_covered_: fraction of _covered_bases_; horizontal genome coverage
   - _mean_coverage_: mean read depth across all _covered_bases_; vertical genome coverage


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
   
- In the `advanced` mode, per-species pileup file will also include the called variants for **all** the genomic sites. 

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

Having run the single-sample SNV step for all the sample, users next can compute the across-sample SNV analysis using the `merge_snps` command. `merge_snps` requires a TSV file (`${my_sample_list}`) specifying the sample name `sample_name` and root output directory of single-sample SNV results `midas_outdir`. See [this page](https://github.com/czbiohub/MIDAS2.0/wiki/Common-Command-Line-Arguments#across-samples-analysis) for details.


### Important Concepts

In this section, we will introduce the species and sample filters, the genomic site filters, the compute of population SNV in MIDAS 2.0, and the chunkified pileup. Beginner users can skip this section and go straight to [Sample commands](). 

1. **<species, sub-samples-lists> selection**

   Population SNV analysis **restricts attention to "sufficiently well" covered species in "sufficiently many" samples**. 

   To be specific, a given <species, sample> pair will only be kept if it has more than 40% horizontal genome coverage (`--genome_coverage`) and 5X vertical genome coverage (`--genome_depth`). Furthermore, only "sufficiently prevalent" species with "sufficiently many" (`--sample_counts`) would be included for the population SNV analysis. Therefore, different species may have different lists of relevant samples.


2. **<site, relevant sample> selection**

   For each genomic site, a sample is considered to be "relevant" if the corresponding site depth falls between the range defined by the input arguments `site_depth` and `site_ratio * mean_genome_coverage`; otherwise it is ignored for the pooled-SNV compute.  

   Therefore, different genomic sites from the same species may have different panels of "relevant samples".  And genomic site prevalence can be computed as the ratio of the number of relevant samples for the given site over the total number of relevant samples for the given species.

3. **Relevant site** 

   For each species, a site is considered to be "relevant" if the site prevalence meets the range defined by the input arguments `--snv_type` and `--site_prev`. By default, common SNV with more than 90% prevalence are reported.

4. **Population SNV Computation**

   There are three main steps to compute and report population SNV in MIDAS 2.0.

   First, for each **relevant** genomic site, MIDAS 2.0 determines the set of alleles present across **all** relevant samples.  Specifically, for each allele (A, C, G, T), `merge_snps` subcommand (1) tallys the sample counts (_sc_) of **relevant samples** containing corresponding allele (`scA:scT`), and (2) sums up the read counts (_rc_) of the corresponding allele across all the relevant samples (`rc_G:rc_T`).

   ```
   site_id                            rc_A  rc_C  rc_G  rc_T  sc_A  sc_C  sc_G  sc_T
   gnl|Prokka|UHGG000587_14|34360|A   26    10    0     0     1     2     0     0
   ```

   Second, population major and minor alleles for a single site can be computed based on accumulated read counts or sample counts across all relevant samples (specified via `--snp_pooled_method`). The population allele refers to the most abundant/prevalent allele, and the population minor allele refers to the second most prevalent/abundant allele. 

   For example, the population major allele of the site `gnl|Prokka|UHGG000587_14|34360|A` in the above example is `A` defined by accumulated read counts and `C` defined by accumulated sample counts. 

   Third, MIDAS 2.0 collects and reports the sample-by-site matrix of the corresponding (1) site depths and (2) allele frequency of the above calculated **population minor allele** for each sample.

5. **Chunkified Pileup Implementation**

   Both single-sample and across-samples pileup was parallelized on the unit of chunk of sites, which is indexed by <species_id, chunk_id>. Only when all chunks from the same species finished processing, chunk-level pileup results would be merged into species-level pileup file. This implementation makes population SNV analysis across thousands of samples possible. To compute the population SNV for one chunk, all the pileup results of corresponding sites across all the samples need to be read into memory. With the uses of multiple CPUs, multiple chunks can be processed at the same time. Therefore, for large collections of samples, we recommend higher CPU counts and smaller chunk sizes to optimally balance memory and I/O usage, especially for highly prevalent species. Users can adjust the number of sites per chunk via the `--chunk_size`. MIDAS 2.0 also has a `--robust_chunk` option, where adjusting chunk size based on species prevalence.


### Sample commands

- Across-samples SNV calling using default filters.

   ```
   midas2 merge_snps --samples_list ${my_sample_list} --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \ --num_cores 32 ${midas_outdir}
   ```

- Users can customize species, sample, site filters. 

   For example, we can apply the species and sample filters as following:  
   - consider species with horizontal coverage > 40%, vertical coverage > 3X and present in more than 30 samples; only consider genomic site with 
     ```
     --genome_coverage 0.4 --genome_depth 3 --sample_counts 30 
     ```

   For example, we can apply genomic sites filters as following: 
   - include site with minimal read depth >= 5, and maximal read depth <= 3 * _mean_coverage_
     ```
     --site_depth 5 --site_ratio 3 
     ```
   
   For example, we can choose to only report SNV meeting the following criteria:
   - compute all the bi-allelic, common population SNV (present in more than 80% of the population) from the protein coding genes based on accumulated sample counts. The minimal allele frequency to call allele present is 0.05.
     ```
     --snp_type bi --snp_maf 0.05 --locus_type CDS --snp_pooled_method prevalence --site_prev 0.8
     ```
   
   Now we can put all the above-mentioned filters in one `merge_snps` command:
  
   ```
   midas2 merge_snps --samples_list ${my_sample_list} \
        --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
        --genome_coverage 0.4 --genome_depth 3 --sample_counts 30 \
        --site_depth 5 --site_ratio 3 \
        --snp_type bi --snp_maf 0.05 --locus_type CDS --snp_pooled_method prevalence \
        --num_cores 32 ${midas_outdir}
   ```

### Output files

- `snps_summary.tsv`: merged single-sample pileup summary. The reported columns _covered_bases_:_mean_coverage_ are the same with single-sample pileup summary.

   ```
   sample_name  species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered  mean_coverage
   SRR172902    100122      2560878        2108551        10782066     248700         207047        0.823             5.113
   SRR172903    100122      2560878        2300193        39263110     1180505        820736        0.898             17.069
   ```
   - sample_name: unique sample name
   - species_id: six-digit species id

- `{species_id}.snps_info.tsv.lz4`: per species SNV metadata information.  

   ```
   site_id                             major_allele  minor_allele  sample_counts  snp_type  rc_A  rc_C  rc_G  rc_T  sc_A  sc_C  sc_G  sc_T  locus_type  gene_id           site_type  amino_acids
   gnl|Prokka|UHGG000587_14|34360|A    A             C              2             bi        26    10    0     0     2     2     0     0     CDS         UHGG000587_02083  4D         T,T,T,T
   gnl|Prokka|UHGG000587_11|83994|T    G             T              2             bi        0     0     11    45    0     0     2     2     IGR         None              None       None
   ```
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

- `{species_id}.snps_freq.tsv.lz4`: per species site-by-sample allele frequency matrix of **population minor allele**.

  ```
  site_id                             SRR172902   SRR172903
  gnl|Prokka|UHGG000587_11|83994|T    0.692       0.837
  gnl|Prokka|UHGG000587_14|34360|A    0.300       0.269
  ```

- `{species_id}.snps_depth.tsv.lz4`: per species site-by-sample site depth matrix. **Only** accounts for the alleles matching the population major and/or population minor allele.

  ```
  site_id                             SRR172902   SRR172903
  gnl|Prokka|UHGG000587_11|83994|T    13          43
  gnl|Prokka|UHGG000587_14|34360|A    10          26
  ```

