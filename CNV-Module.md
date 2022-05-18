# CNV Module: Population Pangenome Copy Number Variants Calling

Similar to the SNV module, the CNV module also proceeds in two phases: (1) single-sample pan-genome copy number variants calling (2) merge these results into a summary across all samples. The first step can be run in parallel. We presuppose users already follow the [database customization](https://github.com/czbiohub/MIDAS2.0/wiki/Data-customization) step, and have `species_profile.tsv` ready for each sample.

## Single-Sample CNV Calling

Species pangenome refers to the set of non-redundant genes (centroids) clustered from all the genomes within on species cluster. Species in the restricted species profile are concatenated and used to build the sample-specific pangenome database, to which reads are aligned using Bowtie2. 

Per species per centroid **copy numbers** are computed in three steps: (1) Per centroid, read alignment metrics, e.g _mapped_reads_ and _mean_coverage_, are computed; (2) Per species, median read coverage of all the mapped centroids corresponding to the 15 universal SCGs are identified; (3) Per centroid, `copy numbers` are computed and gene presence/absence are further inferred.

### Sample commands

- Single-sample CNV calling for all the species in the restricted species profile: `median_marker_coverage > 2` and `unique_fraction_covered` > 0.5. 

  We presuppose users already [profiling the species coverage](https://github.com/czbiohub/MIDAS2.0/wiki/Data-customization#species-to-genotype), and expect `${my_midasdb_dir}/${sample_name}/species/species_profile.tsv` exists.

   ```
   midas2 run_genes --sample_name ${sample_name} -1 ${R1} -2 ${R1} \
         --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \
         --select_by median_marker_coverage,unique_fraction_covered \
         --select_threshold=2,0.5 \
         --num_cores 12 ${midas_outdir}     
   ```

   Users can adjust post-alignment quality filter parameters via the command-line options, and the defaults are:

   ```
   --mapq >= 2: reads aligned to more than one genomic locations equally well are discarded (MAPQ=0,1)
   --mapid >= 0.94: discard read alignments with alignment identity < 0.94
   --aln_readq >= 20: discard read alignment with mean quality < 20
   --aln_cov >= 0.75: discard read alignment with alignment coverage < 0.75
   ```

### Output files

- `genes_summary.tsv`: read mapping and pan-gene coverage summary for all the species in the sample-specific pan-genome database

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
   - _mapped_reads_: number of aligned reads to _gene_id_ after post-alignment filter
   - _mean_coverage_: average read depth of _gene_id_ based on _mapped_reads_ (total_gene_depth / covered_bases)
   - _fraction_covered_: proportion of the _gene_id_ covered by at least one read (covered_bases / gene_length)
   - _copy_number_: estimated copy number of _gene_id_ based on _mapped_reads_ (_mean_coverage_ / median_marker_coverage)

## Across-Samples CNV Calling

Take the same [`${my_sample_list}`](https://github.com/czbiohub/MIDAS2.0/wiki/Common-Command-Line-Arguments#across-samples-analysis):

   ```
   sample_name   midas_outdir
   SRR172902     /home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample
   SRR172903     /home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample
   ```

`merge_genes` would expect to locate `/home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample/SRR172903/genes/genes_summary.tsv`, generated by `run_genes`.


Having run the single-sample CNV analysis for all the samples listed in the `${my_sample_list}`, users next can merge the results and product a summary with `merge_genes` command, and further quantify each pan-gene's presence/absence by comparing the `copy_number` with the user-defined minimal gene copy number (`min_copy`).


### Sample commands

- Across-samples CNV analysis using default filters.

   ```
   midas2 merge_genes --samples_list ${my_sample_list} --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} --num_cores 8 ${midas_outdir}
   ```

   By default, merged gene CNV results are reported for genes clustered at 95% identity. `cluster_pid` and `min_copy` can be customized with the following command-line options:

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
