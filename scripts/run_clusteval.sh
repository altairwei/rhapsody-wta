#!/usr/bin/env bash

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

RUN_ID=$(date +%Y%m%d_%H%M%S)
LOG_DIR=logs/$(date +%Y%m%d)/${RUN_ID}

mkdir -p ${LOG_DIR}

snakemake \
    --use-conda \
    --config Log_Dir=${LOG_DIR} conda_env=rdev \
    --cluster "workflow/profile/pbs-torque/pbs-submit.py -e ${LOG_DIR} -o ${LOG_DIR}" \
    --cluster-status "workflow/profile/pbs-torque/pbs-status.py" \
    --cluster-cancel "qdel" \
    --jobscript "workflow/profile/pbs-torque/pbs-jobscript.sh" \
    --jobname "{name}.{jobid}" \
    --jobs 100 \
    --retries 10 \
    --latency-wait 10 \
    --notemp \
    --snakefile workflow/clusteval/Snakemake \
    "$@"
