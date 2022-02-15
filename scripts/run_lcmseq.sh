#!/usr/bin/env bash

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

snakemake \
    --use-conda \
    --profile workflow/lcm_seq/profile/pbs-torque \
    --snakefile workflow/lcm_seq/Snakemake \
    "$@"