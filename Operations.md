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

Follow the instructions for [managing container images in pairani](https://github.com/czbiohub/pairani/wiki#managing-container-images)

The repository is called `iggtools`.

```
docker build .
docker tag <locally_built_image> 423543210473.dkr.ecr.us-west-2.amazonaws.com/iggtools:latest
aws ecr get-login --region us-west-2 --no-include-email | xargs -Icmd bash -c "cmd"
docker push 423543210473.dkr.ecr.us-west-2.amazonaws.com/iggtools:latest
```