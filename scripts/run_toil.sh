#!/usr/bin/env bash
#PBS -l mem=192gb,nodes=1:ppn=24,walltime=48:00:00
#PBS -q fat
#PBS -N toil-wta
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

export CWL_SINGULARITY_CACHE="$HOME/src/rhapsody-wta/dockerImages"

WORKFLOW=$1
WORKFLOW_INPUT=$2

cd $PBS_O_WORKDIR
source activate toil
ID_STRING=$(date +%Y%m%d_%H%M%S)
RESULTS_FOLDER=results/$ID_STRING

mkdir -p $RESULTS_FOLDER
mkdir -p tmp/$ID_STRING
mkdir -p tmp/$ID_STRING/tmpwork

export TMPDIR=tmp/$ID_STRING/tmpwork

toil-cwl-runner \
  --jobStore file:tmp/$ID_STRING/rhapsody-wta-job-store \
  --singularity \
  --maxCores 24 \
  --defaultMemory 16G \
  --defaultCores 8 \
  --defaultDisk 120G \
  --outdir $RESULTS_FOLDER \
  --workDir tmp/$ID_STRING/tmpwork \
  --writeLogs logs \
  --logFile logs/cwltoil_$ID_STRING.log \
  --noStdOutErr \
  --disableCaching \
  --logLevel INFO \
  --retryCount 2 \
  --maxLogFileSize 20000000000 \
  --stats \
  $WORKFLOW $WORKFLOW_INPUT