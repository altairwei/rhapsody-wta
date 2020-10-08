# Unified Single-Cell Pipeline by BD Genomics

The Unified Single-Cell Pipeline (UP) was created by the  BD Genomics Bioinformatics Group and is also known 
as "MIST": Molecular Index-based Sequencing Analysis Tools.

## Installation


### Docker-based

The pipeline is distributed via docker. To dockerize, run:

```bash
make docker
```

The resulting CWL in `dist/internal_Targeted/MIST.cwl` or `dist/internal_WTA/MIST.cwl` can be run via `cwltool`.

To deploy on Sevenbridges Genomics, run:

```bash
make deploy
```

### Creating a development environment

MIST relies on `pip` and `miniconda` for automated dependency resolution. To install locally, create a conda environment,
allow it to resolve the build dependencies for your machine and, finally, install mist itself. Note, this only works on
Linux.

```bash
conda create -n mist --file conda-env.txt
conda activate mist
python setup.py install
```

By default, this installs the cpython version of mist. To force a cython installation, run the following instead:

```bash
export CYTHON_ENABLED_FOR_MIST_PIPELINE=1 \
&& python setup.py install
```

Currently, it is impossible to install the MIST dependencies on macOS. This is because the ClusterAnalysis node relies 
on the incompatible `tsne` package. Discussions are planned to remove the ClusterAnalysis node, which relies on it.

The other non-conda packages, which may be useful for testing, can be installed on macOS thusly:

```bash
pip install $(grep -v tsne < dependency_resolution/requirements.txt)
pip install -r dependency_resolution/requirements-test.txt
```

## Updating dependencies
The files in the `dependency_resolution` folder define the specification for the build system. Miniconda builds the
non-python extensions deterministically by installing pre-compiled binaries listed in `conda-spec-file.txt`. These 
binaries solve the specification defined in `conda-env.txt`.'

To update the dependencies managed by conda, specify the desired change in `conda-env.txt` and run 
`make update_conda_lockfile`.

There are additional python-specific dependencies. Some of those simply lack a conda build recipe (e.g. `tsne`); these
are therefore defined in `requirements.txt` and are to be installed by `pip`. The dependencies for the test suite are also 
defined in `requirements-testing.txt`. These are not installed in the final docker image but, rather, are installed to 
a temporary stage. Jenkins rescues and uses this temporary stage to run the test suite.