#!/usr/bin/env bash
#PBS -l mem=96gb,nodes=1:ppn=1,walltime=24:00:00
#PBS -q fat
#PBS -N cibersortx_fractions
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR
source secrets

singularity --version

singularity exec --no-home -c \
  -B $PWD/results/Deconv:/src/data \
  -B $PWD/results/Deconv:/src/outdir \
  dockerImages/cibersortx_fractions:latest.sif \
  /src/CIBERSORTxFractions \
  --username $USERNAME \
  --token $TOKEN \
  --refsample sc_mock_counts.tsv \
  --mixture lcm_tpm_matrix.tsv \
  --sigmatrix sc_mock_counts_sigmatrix.tsv \
  --perm 100 \
  --rmbatchSmode TRUE \
  --verbose TRUE
