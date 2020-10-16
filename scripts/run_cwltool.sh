#!/usr/bin/env bash
#PBS -l mem=192gb,nodes=1:ppn=24,walltime=48:00:00
#PBS -q fat
#PBS -N cwltool-wta
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

WORKFLOW=$1
WORKFLOW_INPUT=$2

cd $PBS_O_WORKDIR
source activate rhapsody
RESULTS_FOLDER=results/$(date +%Y%m%d_%H%M%S)
mkdir -p $RESULTS_FOLDER
mkdir -p tmp/docker_tmp

export CWL_SINGULARITY_CACHE="$HOME/src/rhapsody-wta/dockerImages"
export TMPDIR=$(pwd)/tmp/docker_tmp

cwltool \
  --singularity --parallel --outdir $RESULTS_FOLDER \
  $WORKFLOW $WORKFLOW_INPUT