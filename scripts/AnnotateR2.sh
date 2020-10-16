#!/usr/bin/env bash
#PBS -l mem=48gb,nodes=1:ppn=6,walltime=48:00:00
#PBS -q batch
#PBS -N AnnotateR2
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

# Annotation files with GTF format
GTF=$1
# STAR generated genome index tarball
INDEX=$2
# R2 of PE reads
R2_READS=$3

cd $PBS_O_WORKDIR
source activate rhapsody

export CORES_ALLOCATED_PER_CWL_PROCESS=6

mist_annotate_R2.py --index "$GTF,$INDEX" --R2 "$R2_READS"