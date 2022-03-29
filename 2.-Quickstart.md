
After installing MIDAS 2.0, users can test the main functionality of MIDAS 2.0 with the provided unit test data (`tests`). Users doesn't need to pre-download the databases, and it will be downloaded in an on-demand manner.

```
cd MIDAS2.0
bash tests/run_midas.sh 8
```


We also provided the scripts to build custom databases for any given collection of geomes. 

```
cd MIDAS2.0
bash tests/build_midasdb.sh 8
```