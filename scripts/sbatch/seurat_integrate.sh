#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=360G
#SBATCH -p fat
#SBATCH -J seurat_integrate
#SBATCH -o logs/seurat_integrate-%j.out
#SBATCH -e logs/seurat_integrate-%j.err

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd ${SLURM_SUBMIT_DIR}

source activate renv
source env.sh

OUTPUT_DIR="$1"
shift
INPUT_DIRS="$@"

CMD="SeuratAnalysis.R \
  --reduction rpca \
  --anchors 20 \
  --use-matrix \
  --produce-cache \
  --integrate \
  --group-rule \\dDPI-(MOCK|PNR2|TR4) \
  --mem-size 360 \
  --process 64 \
  --data-list results/datalist.csv \
  -O ${OUTPUT_DIR} ${INPUT_DIRS}"

echo "$CMD" >&2

${CMD}