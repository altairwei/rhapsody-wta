#!/usr/bin/env bash
#PBS -l mem=6gb,nodes=1:ppn=12,walltime=48:00:00
#PBS -q batch
#PBS -N bigwig
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR

source activate rhapsody

BAM_FILE=$1
BIGWIG_FILE="${BAM_FILE}.bigwig"

bamCoverage -p 12 -b "${BAM_FILE}" -o "${BIGWIG_FILE}"