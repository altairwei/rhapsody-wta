#!/usr/bin/env bash
#PBS -l mem=192gb,nodes=1:ppn=24,walltime=48:00:00
#PBS -q fat
#PBS -N altairwei-rhapsody-wta
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

WORKFLOW=$1
WORKFLOW_INPUT=$2

cd $PBS_O_WORKDIR
source activate wta
cwltool \
  --udocker --parallel --debug --outdir results \
  $WORKFLOW $WORKFLOW_INPUT