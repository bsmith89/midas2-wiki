
# Single-sample Analysis

Reference-based variants calling and copy number estimation requires choosing a reference genome(s) as the template. Microbiome data usually contains hundreds of species in one sample, and only species with enough reads coverage (horizontal and vertical) can be used for reliable strain-level analysis. A good reference database should be both representative and comprehensive in terms of the abundant species in the given sample. Therefore, single-sample analysis of MIDAS 2.0 start with building a sample-specific reference database of species selected via profiling 15 universal single copy marker genes (Species module).


## Common Command Line Arguments

All three single-sample analysis subcommands of MIDAS 2.0 share the following command line arguments.

```
positional arguments:
  midas_outdir          Path to base directory to store results. 

optional arguments:
  --sample_name SAMPLE_NAME   
                        Unique sample identifier. The output results will be stored at 
                        midas_outdir/sample_name.  
  -1 R1                 FASTA/FASTQ file containing 1st mate if using paired-end 
                        reads. Otherwise FASTA/FASTQ containing unpaired reads.
  -2 R2                 FASTA/FASTQ file containing 2nd mate if using paired-end reads.
  
  --midasdb_name {uhgg,gtdb,testdb}
                        MIDAS Database name. If not specificed, default option is uhgg.
  --midasdb_dir         Path to Local MIDASDB. If not specified, MIDASDB would be 
                        on-demand downloaded to current directory (.).

  --num_cores INT       Number of physical cores to use (8)
  -h, --help            Show help message and exit.
  -f, --force           Force rebuild of pre-existing outputs.
  -g, --debug           Debug mode: skip cleanup on error, extra prints.
```
