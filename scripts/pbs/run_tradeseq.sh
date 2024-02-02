#!/usr/bin/env bash
#PBS -N tradeSeq_fitGAM
#PBS -o logs/tradeSeq
#PBS -l mem=96gb,nodes=1:ppn=1,walltime=120:00:00
#PBS -q fat
#PBS -j oe

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

## setup Conda env
cd $PBS_O_WORKDIR
source activate rdev

logfile=$(mktemp -p logs/tradeSeq memory.XXXXXX)
nohup ./scripts/JobMemRecord.sh $PBS_JOBID $logfile &

## run R
Rscript $1
