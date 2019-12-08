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

5.  [DONE]  [Boris]  Wrap container/instance init code in the python iggtools package.

6.  [IN PROGRESS]  [Boris]  iggtools subcommands for prokka, vsearch, hmmsearch with appropriate use of S3 and NVME to stage results

#  Technical debt 

This section is a light-weight form issue tracking.

## Support AWS instance types other than r5.12xlarge

This should be easy with the recently released update to aegea.  Involves removing the magic numbers 838 and 1715518 from aws_batch_init.

## Handle clusters that do not include their own centroids.

Let cX and cY be 99% clusters with centroids X and Y, respectively.   Normally X is an element of cX and does not belong to any other 99% clusters.  In some rare degenerate cases, X is also a member of cY.  Subsequent  coarser reclustering at 95, 90, ... ANI for the elements of cX would then produce incorrect results.  We need to modify the reclustering assignments to handle this case correctly.

There is a hypothesis this case occurs primarily when contig IDs clash between genomes, so ongoing work to prevent those clashes could address this problem. 