
## Output

The two single-sample strain-level analysis (`run_snps` and `run_genes`), and the database customization analysis (`run_species`) share the following command-line options.

- User-specified root output directory: `midas_outdir=/path/to/root/outdir`. This is always passed as a mandatory positional argument to each of the MIDAS 2.0 command.

- Unique sample name: `sample_name=/unique/sample/identifier`.

Together, `${midas_outdir}/${sample_name}` constitutes the unique output directory for single-sample analysis. 


## Input

### Single-Sample Analysis

The FASTA/FASTQ file containing single-end or paired-ends sequencing reads:

- `R1=/path/to/forward/reads/R1.fq.gz`

- `R2=/path/to/reverse/reads/R2.fq.gz`

`${R1}` and/or `${R2}` need to be passed to MIDAS 2.0 analysis commands via arguments `-1` and `-2` as: ```-1 ${R1} -2 ${R2}```

### Across-Samples Analysis

A TSV file, which lists the _sample_name_ and single-sample root output directory _midas_outdir_, is required for all across-samples analysis. Users need to pass the local path of the TSV file (`my_sample_list=/path/to/tsv/file`) to the command-line argument `sample_list`, e.g. `--sample_list ${my_sample_list}`

A template _sample_list_ is shown here:

  ```
   sample_name       midas_outdir
   SRR172902         /home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample
   SRR172903         /home/ubuntu/hmp_mock/midas2_output_uhgg/single_sample
   ```

## MIDAS DB

For all MIDAS 2.0 analysis, users need to choose (1) a valid precomputed MIDAS DB name (uhgg, gtdb): `my_midasdb_name=uhgg`, and (2) a valid path to local MIDAS DB: `my_midasdb_dir=/path/to/local/midasdb`.

MIDAS 2.0 analysis can takes in two arguments as: `--midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir}`. If the `--midasdb_dir` is not specified, MIDAS DB will be downloaded to the current directory.


### Others Parameters

Users can set the `--num_cores` to the number of physical cores to use: `--num_cores 16`.

And all MIDAS 2.0 analysis can print out the full help message and exit by `-h` or `--help`.

