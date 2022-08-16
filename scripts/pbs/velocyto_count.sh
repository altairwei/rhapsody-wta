#!/usr/bin/env bash
#PBS -l mem=24gb,nodes=1:ppn=6,walltime=24:00:00
#PBS -q batch
#PBS -N velocyto_count
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR

# RunVelocyte.py depends on samtools
source activate wta
source renv/python/virtualenvs/renv-python-3.8/bin/activate
source env.sh

INPUT_DIR=$1
GTF_FILE=$2

RunVelocyto.py -R --samtools-threads 6 --samtools-memory 4096 ${INPUT_DIR} ${GTF_FILE}