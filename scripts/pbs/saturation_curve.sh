#!/usr/bin/env bash
#PBS -l mem=2gb,nodes=1:ppn=1,walltime=24:00:00
#PBS -q batch
#PBS -N saturation_curve
#PBS -o logs
#PBS -e logs

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR

source activate renv
source env.sh

INPUT_DIR=$1
BAM_FILE="${INPUT_DIR}/*_final.BAM"
QC_DIR="${INPUT_DIR}/QualityControl"

mkdir -p "${QC_DIR}"

saturation.py ${BAM_FILE} | tee "${QC_DIR}/QC_Saturation_Curve.csv" | PlotSaturation.R -o "${QC_DIR}/QC_Saturation_Curve.png"