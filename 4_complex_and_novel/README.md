To run this benchmark, after activating the master conda environment as per `../README.md` run the Snakemake
```
PYTHONPATH=.. snakemake -c 32 --configfile ../1_novel_strains/config-8threads.yaml --use-conda
```
