#!/usr/bin/env bash
#PBS -l mem=6gb,nodes=1:ppn=6,walltime=48:00:00
#PBS -q fat
#PBS -N samtools-index
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR

source activate wta

BAM_FILE=$1

samtools \
  index \
  -@ 6 \
  -c \
  -- "$BAM_FILE"
