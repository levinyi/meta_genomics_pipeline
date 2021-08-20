this is the manual of meta genomics data analysis pipeline

step by step:

cornavirus fasta need pre-build reference index before using this pipeline.

or using nt database for blastn

create a new prokka_env
```
conda activate prokka_env
```

### 
# note that prokka must using under prokka_env:
conda create -n prokka_env -c conda-forge -c bioconda prokka
conda activate prokka_env
###
# some databases need pre-build reference index before using it. for example:
1) cornavirus fasta
2) kraken2 database
3) prokka database
4) meghit
5) humann3 database
6) metaphlan3 database
7) blastn to kraken2db
example:
```
cd example
# configuring data.list and config.txt then run:
sh work.sh

```
