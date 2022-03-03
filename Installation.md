
# Install with Conda


Add proper channel order in the `.condarc` file, when downloading [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) for the first time.
```
conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge
```

Create Conda environment.

```
conda env create --name=iggtools --quiet --file iggtools.yml
cpanm Bio::SearchIO::hmmer --force # Temporary fix for Prokka
```

Install Iggtools.

```
git clone https://github.com/czbiohub/iggtools && cd iggtools && python3 setup.py build && python3 setup.py install

# Unit test
bash tests/run_midas.sh 8
```

Export environment.

```
# After make some custom changes 
conda update --all 
conda clean –all

conda env export --no-builds | grep -v "^prefix:" > iggtools.update.yml
```

# Install with Docker

Run.

```
docker pull zhaoc1/midas:latest
sudo docker run --volume "/home/ubuntu/.aws":"/root/.aws":ro --rm -it midas:latest
```

Build.

```
docker image build -t midas:v0.7 -t midas:latest -f Dockerfile .

docker image tag midas:latest zhaoc1/midas:v0.7
docker image tag midas:latest zhaoc1/midas:latest
```

Push.

```
docker image push zhaoc1/midas:v0.7
docker image push zhaoc1/midas:latest
```



