#!/usr/bin/env bash
#PBS -l mem=128gb,nodes=1:ppn=6,walltime=2400:00:00
#PBS -q fat
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

SeuratAnalysis.R --mem-size 128 -m -d -c --process 6 --integrate --group-rule "\\dDPI-(MOCK|PNR2|TR4)" -O "${OUTPUT_DIR}" ${INPUT_DIRS}