#!/usr/bin/env bash
#PBS -l mem=64gb,nodes=1:ppn=6,walltime=24:00:00
#PBS -q batch
#PBS -N seurat_analysis
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR

source activate renv
source env.sh

SeuratAnalysis.R --mem-size 64 --compress -d --process 6 -c "$@"