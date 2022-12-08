## Run scRNA-Seq pipeline

Firstly, we need to setup environments:

```shell
source env.sh
```

How to run workflow:

```shell
run_toil.sh workflow/rhapsody_wta_1.8/main.cwl data/3DPI-MOCK-2.yml
```

Restart last workflow without losing any progress:

```shell
run_toil.sh -r workflow/rhapsody_wta_1.8/main.cwl data/3DPI-MOCK-2.yml
```

Start the same workflow from scratch:

```shell
run_toil.sh -c workflow/rhapsody_wta_1.8/main.cwl data/3DPI-MOCK-2.yml
```

Run toil leader within a named screen session:

```shell
run_toil.sh -s workflow/rhapsody_wta_1.8/main.cwl data/3DPI-MOCK-2.yml
```

## Run LCM-seq pipeline

The dependencies of LCM-Seq pipeline will be managed by Snakemake, so there's no need to create a conda environment by youself.

```shell
source env.sh
run_lcmseq.sh all
```

## Develop Environment

### Setup Conda

```shell
# Install miniconda
wget https://repo.anaconda.com/miniconda/Miniconda3-py39_4.12.0-Linux-x86_64.sh
bash Miniconda3-py39_4.12.0-Linux-x86_64.sh

# Add channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict

# Create environments
conda create -n renv python=3.8
conda activate renv
conda install -n renv -c conda-forge r-base=4.2
```

The following dependencies are required for building R packages from sources, you can also install them using system package manager if you have root privilege.

```shell
conda install -n renv -c conda-forge pkg-config gmp libgit2 fftw pandoc
```

### Setup renv

Python package `velocyto.py` imports `numpy` and `Cython` in the file `setup.py`, so we need to install two of them manually prior to restoring dependencies.

```shell
source ./renv/python/virtualenvs/renv-python-3.8/bin/activate
pip install numpy==1.21.6 Cython==0.29.30
```

Restore R packages:

Why `R.utils` not recorded in renv.lock?

```R
renv::restore()
```
