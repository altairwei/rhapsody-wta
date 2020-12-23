#!/usr/bin/env bash
#PBS -l mem=20gb,nodes=1:ppn=1,walltime=24:00:00
#PBS -q batch
#PBS -N seurat_integrate
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR

source activate renv
source env.sh

OUTPUT_DIR="$1"
shift
INPUT_DIRS="$@"

SeuratAnalysis.R -m --plot -O "${OUTPUT_DIR}" --produce-cache --integrate ${INPUT_DIRS}