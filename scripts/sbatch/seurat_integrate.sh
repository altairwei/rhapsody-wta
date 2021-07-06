#!/usr/bin/env bash
#SBATCH -N 1
#SBATCH -n 64
#SBATCH --mem=128G
#SBATCH -p fat
#SBATCH -J seurat_integrate
#SBATCH -o logs/seurat_integrate-%j.out
#SBATCH -e logs/seurat_integrate-%j.err

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $SLURM_SUBMIT_DIR

source activate renv
source env.sh

OUTPUT_DIR="$1"
shift
INPUT_DIRS="$@"

REF="1,3,5,7,9,11,13,15,17,20"

CMD="SeuratAnalysis.R \
  -r rpca --anchors 5 \
  -m -d -c \
  --integrate \
  --group-rule \\dDPI-(MOCK|PNR2|TR4) \
  --mem-size 128 --process 64 \
  -O ${OUTPUT_DIR} ${INPUT_DIRS}"

echo "$CMD" >&2

${CMD}