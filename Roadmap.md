# Overall plan of attack

 A)  [Boris, Chunyu]   Build a small DB, based on a few important species, using new AWS batch flow.

 B)  [Boris]   Make MIDAS work on that small DB.

 C)  [Chunyu]  Validate results.

 D)  [Boris]   Extend DB to all species for GTPro validation.

 E)  [Chunyu]  Ensure results still good.

 F)  [Boris]   Build entire DB

#  Task A:  Build a small DB, based on a few important species, using new AWS batch flow.

1.  [DONE] [Chunyu]  Document and test MIDAS DB build steps (prokka, vsearch, hmmsearch).

2.  [DONE] [Chunyu, Boris]  Blueprint for DB layout in S3.

3.  [DONE] [Boris]  Build docker image with all tools.

4.  [DONE] [Boris]  Validate local NVME init/sharing from Batch containers.

5.  [IN PROGRESS]  [Boris]  Wrap container/instance init code in the python iggtools package.

6.  [IN PROGRESS]  [Boris]  iggtools subcommands for prokka, vsearch, hmmsearch with appropriate use of S3 and NVME to stage results

