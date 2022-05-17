
## Single-Sample Analysis

### Input

The FASTA/FASTQ file containing single-end or paired-ends sequencing reads:

- `R1=/path/to/forward/reads/R1.fq.gz`

- `R2=/path/to/reverse/reads/R2.fq.gz`

`${R1}` and `${R2}` needs to be provided to MIDAS 2.0 subcommands via arguments `-1` and `-2` as: `-1 ${R1} -2 ${R2}` 

### Output

The two single-sample strain-level analysis subcommands (`run_snps` and `run_genes`), and the database customization subcommand (`run_species`) share the following command line arguments.

- User-specified root output directory: `midas_outdir=/path/to/root/outdir`. This is always passed as a mandatory positional argument to each of the MIDAS 2.0 command.

- Unique sample name: `sample_name=/unique/sample/identifier`.

Together, `${midas_outdir}/${sample_name}` constitute the unique output directory for single-sample analysis. 


## MIDAS DB

For all MIDAS 2.0 subcommands, users need to choose (1) a valid precomputed MIDAS DB name (uhgg, gtdb):

- `my_midasdb_name=UHGG`

and (2) a valid path to local MIDAS DB:

- `my_midasdb_dir=/path/to/local/midasdb`.

MIDAS 2.0 subcommands can takes in two arguments as: `--midasdb_name ${my_midasdb_name} --midasdb_dir ${my_midasdb_dir}`. If the `--midasdb_dir` is not specified, MIDAS DB would be downloaded to the current directory.


### Others

Users can set the `--num_cores` to the number of physical cores to use: `--num_cores 16`.

And all MIDAS subcommands can print out the full help message and exit by `-h` or `--help`.

