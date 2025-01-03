#!/usr/bin/env bash
#PBS -l mem=20gb,nodes=1:ppn=1,walltime=24:00:00
#PBS -q batch
#PBS -N convert_to_mtx
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR

source activate renv
source env.sh

INPUT_DIR=$1

SeuratAnalysis.R -c "${INPUT_DIR}"