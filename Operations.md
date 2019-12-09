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
iggtools aws_batch_submit --batch_command "iggtools build_pangenome -s 13:100"
```
will build 1% of pangenomes in AWS, those whose species IDs are congruent to 13 modulo 100.

To build all pangenomes via 100 concurrent jobs:
```
for i in {0...99}; do
    iggtools aws_batch_submit --batch_command "iggtools build_pangenome -s ${i}:100"
done
```
This could be expensive, so make sure it works on 1% first.

To fully realize 100-way parallelism, the AWS Batch compute environment max_vcpus and EC2 instance limit for r5.12x instances must be raised to 4800 and 100, respectively.  Otherwise, 100 jobs will be queued, but a smaller number (as permitted by those limits) would execute in parallel, with the rest waiting behind.

# Tracking operational events

Submitting a job creates a unique record under
```
s3://microbiome-igg/2.0/operations/<utc_date>/<unix_time>__<event_type>__<job_id>.json
```
where event_type is `aws_batch_submit`.  The record points to the job's various logs and results.  In the above example, where 100 jobs are submitted, there would be 100 distinct records.