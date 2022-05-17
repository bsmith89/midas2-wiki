
# Population Single Nucleotide Variants (SNVs) Calling

The SNV module proceeds in two phages: (1) single-sample read pileup (2) population variants calling across all the species. The first steps can be potentially in parallel.  We presuppose users already run the [database customization](https://github.com/czbiohub/MIDAS2.0/wiki/Data-customization) step, and either have single-sample `species_profile.tsv` or a prebuilt Bowtie2 rep-genome database ready for the SNV module.


## Single-Sample SNV Analysis

In a standard workflow, rep-genome database of the species in the restricted species profile were built for each sample, to which reads were aligned using Bowtie2. Per genomic site read pileup and nucleotide variation for **all** the species in the rep-genome database are reported by MIDAS 2.0.  

MIDAS 2.0 purposely hold any species selection based on the pileup results until across-samples SNVs analysis. That being said, per species reads mapping summary are reported in `snps_summay.tsv`, and pileup/variants calling results are organized by species. Therefore users can easily filter species accordingly based on the reported `horizontal_coverage` and `vertical_coverage` in the `snps_summary.tsv` on their own.

### Sample command

- Single-sample Pileup for all the species in the restricted species profile: `median_marker_coverage > 2` and `unique_fraction_covered > 0.5`. 

  We presuppose users already [profiling the species coverage](https://github.com/czbiohub/MIDAS2.0/wiki/Data-customization#species-to-genotype), and expect `${my_midasdb_dir}/${sample_name}/species/species_profile.tsv` exists.

   ```
   midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
         --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
         --select_by median_marker_coverage,unique_fraction_covered \
         --select_threshold=2,0.5 \
         --num_cores 12 ${midas_outdir}         
    ```

  The default post-alignment filters are:

   ```
   --mapq >= 20: discard read alignments with alignment quality < 20
   --mapid >= 0.94: discard read alignments with alignment identity < 0.94
   --aln_readq >= 20: discard read alignment with mean quality < 20
   --aln_cov >= 0.75: discard read alignment with alignment coverage < 0.75
   --aln_baseq >= 30: discard bases with quality < 30
   ```

  Users can adjust these post-alignment quality filter parameters via the command line arguments. 


- Single-sample pileup for all the species in the restricted species profile with paired-ends based post-alignment quality filter.
  
  We recommend to use the `--paired_only` with proper `--fragment_length`, to recruiting only properly paired reads for pileup. The post-alignment metrics would be computed based on a reads pair, instead of single read.

   ```
   midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
         --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
         --select_by median_marker_coverage,unique_fraction_covered \
         --select_threshold=2,0.5 \
         --fragment_length 2000 --paired_only \
         --num_cores 12 ${midas_outdir}         
    ```

- Single-sample variant calling for all the species in the restricted species profile with paired-ends based post-alignment quality filter.
  
  In recognition of the need for single-sample SNV calling, we added an `--advanced` mode together with `--ignore_ambiguous`, to the single-sample SNV step in MIDAS 2.0 to report per species major allele and minor allele for all the genomic sites covered by at least two reads, upon which custom variant calling filter should be applied by the users. MIDAS 2.0 ignore ambiguous alleles calling, where genomic sites recruit tied read counts.
`--ignore_ambiguous`   recommend users set the parameter  to ignore ambiguous alleles, 


- Single-sample pileup for all the species in a prebuilt rep-genome database. 
  
  We presuppose users already follow the [database customization](https://github.com/czbiohub/MIDAS2.0/wiki/Data-customization#population-specific-species-panel) step, and have the prebuilt rep-genome database located at `${midas_outdir}/bt2_indexes`. 

   ```
   midas2 run_snps --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
         --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
          --prebuilt_bowtie2_indexes ${midas_output}/bt2_indexes/repgenome \
          --prebuilt_bowtie2_species ${midas_output}/bt2_indexes/repgenome.species \
          --select_threshold=-1 --num_cores 12 ${midas_output}
   ``` 

### Output files

- `snps_summary.tsv`: the statistics of read mapping and pileup summary

   ```
   species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered   mean_coverage
   102506      5339468        2373275        8045342      468667         224553        0.444              3.390
   102337      2749621        2566404        47723458     1479479        1010530       0.933              18.595
   ```
   - genome_length: genome length
   - covered_bases: number of bases covered by at least two reads
   - total_depth: total read-depth across all covered_bases
   - aligned_reads: total read counts across sites
   - mapped_reads: total read counts across sites after quality filter
   - fraction_covered: fraction of sites covered by at least two reads, horizontal genome coverage
   - mean_coverage: 


- Per-species reads pileup results: `{species}.snps.tsv.lz4`

   ```
   ref_id                    ref_pos   ref_allele  depth   count_a  count_c  count_g  count_t
   gnl|Prokka|UHGG144544_1   881435    T           11      0        0        0        11
   gnl|Prokka|UHGG144544_1   881436    T           13      0        5        0        8
   gnl|Prokka|UHGG144544_1   881437    T           12      0        6        0        6
   ```

- In the advanced mode (`--advanced`), the per-species pileup results would also report the called variants for **all** the genomic sites. We recommend setting `--ignore_ambiguous` to avoid falsely calling major/minor allele for site with read counts ties (e.g. ref_position: 881437).

   ```
   ref_id                    ref_pos   ref_allele  depth   count_a  count_c  count_g  count_t  major_allele  minor_allele  major_allele_freq  minor_allele_freq  allele_counts
   gnl|Prokka|UHGG144544_1   881435    T           11      0        0        0        11       T             T             1.000              0.000              1
   gnl|Prokka|UHGG144544_1   881436    T           13      0        5        0        8        T             C             0.615              0.385              2
   gnl|Prokka|UHGG144544_1   881437    T           12      0        6        0        6        C             T             0.500              0.500              2
   ```



# Across-samples Analysis

Population metagenomic SNVs calling and pangenome CNVs estimation can be accomplished by merging single-sample analysis results using the subcommands `merge_snps` and `merge_genes`.  All across-samples subcommands require the input file `samples_list`: a TSV file with `sample_name` and single-sample midas output directory `midas_outdir`. For example, given the following `samples_list`, `merge_snps` would expect to locate `/home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample/SRR172903/snps/snps_summary.tsv`, generated by `run_snps`.

   ```
   sample_name   midas_outdir
   SRR172902     /home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample
   SRR172903     /home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample
   ```


## Population SNVs Calling

1. **<species, sub-samples-lists> selection**

   As mentioned in the single-sample step, MIDAS 2.0 purposely hold any filter upon the single-sample pileup summary until the across-samples step.  The across-samples SNPs analysis restricts attention to "sufficiently well" covered species in "sufficiently many" samples. To be specific, only <species, sample> pair with more than 40% horizontal genome coverage (`fraction_covered`) and 5X vertical genome coverage (`mean_coverage`) would be kept. Furthermore, only "sufficiently prevalent" species with "sufficiently many" (`sample_counts`) would be included for the population SNVs analysis. Therefore, different species may have different lists of relevant samples.


2. **<site, relevant sample> selection**

   For each genomic site, a sample is considered to be "relevant" if the corresponding site read depth falls between the range defined by the input parameters `site_depth` and `site_ratio * genome_depth`; otherwise ignored for the pooled-SNPs compute.  Therefore, different genomic sites from the same species may have different panels of "relevant samples".  And genomic site prevalence can be computed as the ratio of the number of relevant samples for the given site over the total number of relevant samples for the given species.

3. **Relevant site** 

   For each species, a site is considered to be "relevant" if the site prevalence meets the range defined by the input argument `snv_type` and `site_prev`. 


4. **Population SNVs**

   (1) For each **relevant** genomic site, MIDAS 2.0 determines the set of alleles present across **all** relevant samples. Specifically, for each allele (A, C, G, T), `merge_snps` subcommand (1) counts the sample counts (sc) of **relevant samples** containing corresponding allele (`scA:scT`), and (2) sums up the read counts (rc) of the corresponding allele across all the relevant samples (`rc_G:rc_T`).

   ```
   site_id                            rc_A  rc_C  rc_G  rc_T  sc_A  sc_C  sc_G  sc_T
   gnl|Prokka|UHGG000587_14|34360|A   26    10    0     0     2     2     0     0
   gnl|Prokka|UHGG000587_11|83994|T   0     0     11    45    0     0     2     2
   ```

   (2) Across-samples major and minor allele for a single site can be computed by the metric specified via the `snps_pooled_method`: sample counts as in `prevalence` and read counts as in `abundance`. For example, the population major allele of site UHGG143484_1|905|T in the above example is T defined by accumulated read counts and C defined by accumulated sample counts. And the population minor allele refers to the second most prevalent/abundant allele.  

   (3) The sample by site matrix of the corresponding site depths and allele frequency of the above calculated across-samples minor allele for each sample will be collected and reported.

   More details about the compute can be found at [Cross-Sample SNP Analysis Tools (xsnp)](https://github.com/czbiohub/xsnp/wiki/Data-Schema-And-Computation)


### Sample commands

   ```
   midas2 merge_snps --samples_list /path/to/sample/lists --num_cores 32 /path/to/merged/midas/outdir
   ```

### Output files

- Pooled single-sample pileup summary: `snps_summary.tsv`

   ```
   sample_name  species_id  genome_length  covered_bases  total_depth  aligned_reads  mapped_reads  fraction_covered  mean_coverage
   SRR172902    100122      2560878        2108551        10782066     248700         207047        0.823             5.113
   SRR172903    100122      2560878        2300193        39263110     1180505        820736        0.898             17.069
   ```

- Per species SNPs information file: `{species_id}.snps_info.tsv.lz4`.  It contains the computed population major and minor alleles, together with all the metadata.

   ```
   site_id                             major_allele  minor_allele  sample_counts  snp_type  rc_A  rc_C  rc_G  rc_T  sc_A  sc_C  sc_G  sc_T  locus_type  gene_id           site_type  amino_acids
   gnl|Prokka|UHGG000587_14|34360|A    A             C              2             bi        26    10    0     0     2     2     0     0     CDS         UHGG000587_02083  4D         T,T,T,T
   gnl|Prokka|UHGG000587_11|83994|T    G             T              2             bi        0     0     11    45    0     0     2     2     IGR         None              None       None
   ```

- Per species site-by-sample allele frequency matrix of **population minor allele**: `{species_id}.snps_freq.tsv.lz4`.

  ```
  site_id                             SRR172902   SRR172903
  gnl|Prokka|UHGG000587_11|83994|T    0.692       0.837
  gnl|Prokka|UHGG000587_14|34360|A    0.300       0.269
  ```

- Per species site-by-sample site depth matrix: `{species_id}.snps_depth.tsv.lz4`. **Only** accounts for the alleles matching either population major or/and population minor allele.

  ```
  site_id                             SRR172902   SRR172903
  gnl|Prokka|UHGG000587_11|83994|T    13          43
  gnl|Prokka|UHGG000587_14|34360|A    10          26
  ```


