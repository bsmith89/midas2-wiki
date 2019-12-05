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
