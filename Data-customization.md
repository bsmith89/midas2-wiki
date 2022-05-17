
Reference-based metagenotyping pipeline requires users to choose a reference genome(s) as the **template genome database**. Microbiome data usually contains hundreds of species in one sample, and only species with enough read coverage can be used for reliable strain-level analysis. A good reference database should be both representative and comprehensive in terms of the sufficiently abundant species in the sample. Therefore, a typical MIDAS 2.0 workflow starts with a "database customization" step which build **sample-specific reference database** of sufficiently abundant species selected via profiling 15 universal single copy marker genes. 


## Single Sample Abundant Species Detection

Species coverage was estimated via profiling 15 universal single copy genes (SCGs) (`run_species`). The simple 15 universal SCGs serves the purpose of quickly screening the panel of abundant species in the sample, instead of taxonomic profiling.

### Sample command

  ```
  midas2 run_species --sample_name ${sample_name} -1 ${R1} -2 ${R1} --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} --num_cores 8 ${midas_outdir}
  ```

### Output files

- `species_profile.tsv`: the primary output of the abundant species detection analysis. Species are sorted in decreasing order of `median_marker_coverage`. Only species with more than two marker genes covered with more than two reads are reported. 
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

## Species To Genotype

Results from the SCGs profiling are used to identify the panel of species eligible for strain-level genomic variation analysis. In the `run_snps` and `run_genes` analysis, users need to pass the `--select_by` and `--select_threshold` accordingly to select the list of species to genotype. `--select_by` are comma separated columns from above mentioned `species_profile.tsv`, and `--select_threshold` are comma separated corresponding cutoff values.

For sufficiently abundant species selection, we recommend using the combination of `median_marker_coverage > 2X` and `unique_fraction_covered > 0.5`:

```
--select_by median_marker_coverage,unique_fraction_covered --select_threshold=2,0.5 
```

For genotyping low abundant species, users need to adjust the parameters properly:

```
--select_by median_marker_coverage,unique_fraction_covered --select_threshold=0,0.5 
```

An alternative way is to pass a comma separated species of interests to `--species_list`. It is worth noting that the species in the provided species list is still subject to the `--select_threshold` restriction. Users can set `--select_threshold=-1` to escape any filters:

```
--species_list 102337,102506 --select_threshold=-1
```

Sample-specific rep-genome and/or pan-genome database would be built only for the species passing the above mentioned filters in the single-sample SNV or CNV module. 

  
## Across-Samples Abundant Species

For some analysis, using one comprehensive list of species across samples in the same study may be desired, in part because it can avoid the time needed to build a new index for each sample. If taking this approach, users run `merge_species` subcommands to merge all the single-sample SCGs profiling results for all the samples listed in the `my_samples_list=/path/to/list/of/sample/tsv`

### Sample command

   ```
   midas2 merge_species --samples_list ${my_sample_list} ${midas_outdir}
   ```

### Output files

- `species_prevalence.tsv`: the primary output of the across-samples species merging analysis. 

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

- Each column in the single-sample `species_profile.tsv` are merged across samples into a species-by-sample matrix, shown as following:
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

## Population Specific Species Panel

The Bowtie2 rep-genome / pan-genome database can be build upon a list of customized species across a given panel of samples. Species selection metrics based on the `species_prevalence.tsv` can be passed to `build_bowtie2eb` via `--select_by` and `--select_threshold`. For example, to build one rep-genome database for all the species that is present in more than two samples:

```
build_bowtie2db \
    --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} \ 
    --select_threshold sample_counts --select_by 2 --num_cores 8 \
    --bt2_indexes_name repgenome --bt2_indexes_dir ${midas_outdir}/bt2_indexes
```

  - The generated rep-genome database can be found under the directory `${midas_outdir}/bt2_indexes`
  - The list of customized species can be found at `${midas_outdir}/bt2_indexes/repgenome.species`


If taking this approach, for the single-sample SNV or CNV analysis, users can pass the pre-built rep-genome to `run_snps` analysis (pan-genome for `run_genes`), as following:

```
--prebuilt_bowtie2_indexes ${merge_midas_outdir}/bt2_indexes/repgenome \
--prebuilt_bowtie2_species ${merge_midas_outdir}/bt2_indexes/repgenome.species \
--select_threshold=-1
```

Having finished the database customization step, users can now go to [SNV](https://github.com/czbiohub/MIDAS2.0/wiki/SNV-Module) or [CNV](https://github.com/czbiohub/MIDAS2.0/wiki/CNV-Module) modules, depending on the scientific aims. 

