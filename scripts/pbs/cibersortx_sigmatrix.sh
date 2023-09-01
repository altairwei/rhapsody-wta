#!/usr/bin/env bash
#PBS -l mem=64gb,nodes=1:ppn=1,walltime=24:00:00
#PBS -q batch
#PBS -N cibersortx_sigmatrix
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR
source secrets

INPUT_DIR="$1"

singularity --version

singularity exec --no-home -c \
  -B $INPUT_DIR:/src/data \
  -B $INPUT_DIR:/src/outdir \
  dockerImages/cibersortx_fractions:latest.sif \
  /src/CIBERSORTxFractions \
  --username $USERNAME \
  --token $TOKEN \
  --verbose TRUE \
  --single_cell TRUE \
  --fraction 0 \
  --refsample sc_mock_counts.tsv

mv $INPUT_DIR/CIBERSORTx_sc_mock_counts_inferred_phenoclasses.CIBERSORTx_sc_mock_counts_inferred_refsample.bm.K999.txt \
   $INPUT_DIR/sc_mock_counts_sigmatrix.tsv
