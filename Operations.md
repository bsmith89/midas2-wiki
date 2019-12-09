# Prerequisites

Python 3.7 or above (as installed in Dockerfile).

AWS credentials/roles/permissions configured appropriately for aws cli to be able to access S3, ECR, Batch, etc.


# System-wide install

```
pip3 install 'git+git://github.com/czbiohub/iggtools' --upgrade

iggtools --version
```

# Lint and run locally for development

```
cd /path/to/iggtools

pylint iggtools

python3 -m iggtools --version
```

# Build and push docker image for batch

Essential commands to update the docker container for Batch jobs.

```
make changes to Dockerfile as needed
docker build .
docker tag <locally_built_image> 423543210473.dkr.ecr.us-west-2.amazonaws.com/iggtools:latest
bash -c "`aws ecr get-login --region us-west-2 --no-include-email`"
docker push 423543210473.dkr.ecr.us-west-2.amazonaws.com/iggtools:latest
```

For more background, see the [PairANI instructions for managing container images](https://github.com/czbiohub/pairani/wiki#managing-container-images).


# Submit subcommands to run in AWS Batch

```
iggtools aws_batch_submit --batch_command "iggtools build_pangenome -s 1:100"
```
will build 1% of pangenomes in AWS, the ones that have IDs equal to 1 modulo 100.

To build all genomes using 100 parallel jobs:
```
for i in {0...99}; do
    iggtools aws_batch_submit --batch_command "iggtools build_pangenome -s ${i}:100"
done
```
This could be expensive, so make sure it works on 1% first.

Operational status updates will appear under 
```
s3://microbiome-igg/2.0/operations/<utc_date>/<unix_time>__<event_type>__<job_id>.json
```
