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
conda create -n rdev python=3.8
conda activate rdev
conda install -n rdev -c conda-forge r-base=4.2
```

The following dependencies are required for building R packages from sources, you can also install them using system package manager if you have root privilege.

```shell
conda install -n rdev -c conda-forge pkg-config gmp libgit2 fftw pandoc hdf5
```

### Setup renv

Python package `velocyto.py` imports `numpy` and `Cython` in the file `setup.py`, so we need to install two of them manually prior to restoring dependencies.

```shell
source ./renv/python/virtualenvs/renv-python-3.8/bin/activate
pip install numpy==1.21.6 Cython==0.29.30
```

Restore R packages:

```R
renv::restore()
```

**Note**: If you have problems compiling packages such as `Cairo` from sources, you can install `r-cairo` from conda with the same version to the record in `renv.lock`. Some of the following packages may require hydration: `XML (r-xml)`, `clock (r-clock)` and `Cairo (r-cairo)`.

```shell
conda install -n rdev -c conda-forge r-cairo=1.6_0
```

Copy the cache from conda using `hydrate()`:

```R
renv::hydrate(packages = "Cairo", sources="~/miniconda3/envs/rdev/lib/R/library")
```

Then continue to `renv::restore()`.
