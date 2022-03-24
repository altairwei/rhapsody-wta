#!/usr/bin/env bash

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

source activate lcmseq

RUN_ID=$(date +%Y%m%d_%H%M%S)
LOG_DIR=logs/$(date +%Y%m%d)/${RUN_ID}

mkdir -p ${LOG_DIR}

snakemake \
    --use-conda \
    --config Log_Dir=${LOG_DIR} \
    --cluster "workflow/lcm_seq/profile/pbs-torque/pbs-submit.py -e ${LOG_DIR} -o ${LOG_DIR}" \
    --cluster-status "workflow/lcm_seq/profile/pbs-torque/pbs-status.py" \
    --jobscript "workflow/lcm_seq/profile/pbs-torque/pbs-jobscript.sh" \
    --jobs 5000 \
    --latency-wait 10 \
    --verbose \
    --notemp \
    --snakefile workflow/lcm_seq/Snakemake \
    "$@"