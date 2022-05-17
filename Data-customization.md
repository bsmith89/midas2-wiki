
Reference-based metagenotyping pipeline requires users to choose a reference genome(s) as the **template genome database**. Microbiome data usually contains hundreds of species in one sample, and only species with enough read coverage can be used for reliable strain-level analysis. A good reference database should be both representative and comprehensive in terms of the sufficiently abundant species in the sample. Therefore, a typical MIDAS 2.0 workflow starts with a "database customization" step which build **sample-specific reference database** of sufficiently abundant species selected via profiling 15 universal single copy marker genes. 


## Abundant Species Detection

Species coverage was estimated via profiling 15 universal single copy genes (SCGs) (`run_species` subcommand). The simple 15 universal SCGs serves the purpose of quickly screening the panel of abundant species in the sample, instead of taxonomic profiling.

### Sample command

  ```
  midas2 run_species --sample_name ${sample_name} -1 ${R1} -2 ${R1} --midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir} --num_cores 8 ${midas_outdir}
  ```

### Output files

- The primary output of the abundant species detection analysis is: `species_profile.tsv`. 

Only species with more than two marker genes covered with more than two reads are reported in the `species_profile.tsv`. And the species are sorted in decreasing order of `median_marker_coverage`. 

   ```
   species_id  marker_read_counts  median_marker_coverage  marker_coverage  marker_relative_abundance   unique_fraction_covered
   102337      4110                28.48                   28.91            0.30                        1.00
   102506      734                 4.98                    4.98             0.05                        0.93
   ```


## Species To Genotype

The panel of species eligible for strain-level genomic variation analysis is determined by the SCG profiling results. Sample-specific rep-genome and pan-genome database would be built only for the species passing the filter.  We recommend using the combination of `median_marker_coverage` and `unique_fraction_covered` as the metrics to determine the list of abundant species. However, if the goal of the experiment is to genotype low abundant species, then users need to set the parameters properly, e.g. `median_marker_coverage > 0`.  The two subcommands `run_snps` and `run_genes` share the same command line arguments for filtering species:

```
  --select_by SELECT_BY
                        Comma separated columns from species_profile to filter
                        species.
  --select_threshold CHAR
                        Comma separated cutoff correspond to select_by to filter
                        species (> XX) (2, )
```
  
## Species Merge and Fit in species_list

An alternative way to select the species for downstream analysis is to directly provide the list of species. The species in the provided species list is still subject to the `select_threshold` restriction. Users can set `--select_threshold=-1` to escape any filters.

```
  --species_list CHAR   Comma separated list of species ids
```

