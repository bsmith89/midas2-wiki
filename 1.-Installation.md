
# Install with Conda


Add proper channel order in the `.condarc` file, when downloading [Miniconda3](https://docs.conda.io/en/latest/miniconda.html) for the first time.
```
conda config --add channels anaconda
conda config --add channels bioconda
conda config --add channels conda-forge
```

Create Conda environment.

```
conda env create --name=midas2.0 --quiet --file midas2.0.yml
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
docker pull zhaoc1/midas2:latest
sudo docker run --volume "/home/ubuntu/.aws":"/root/.aws":ro --rm -it midas2:latest
```

Build.

```
docker image build -t midas2:v0.8 -t midas2:latest -f Dockerfile .

docker image tag midas2:latest zhaoc1/midas2:v0.8
docker image tag midas2:latest zhaoc1/midas2:latest
```

Push.

```
docker image push zhaoc1/midas2:v0.8
docker image push zhaoc1/midas2:latest
```



