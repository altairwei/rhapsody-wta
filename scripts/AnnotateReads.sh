#!/usr/bin/env bash
#PBS -l mem=32gb,nodes=1:ppn=1,walltime=48:00:00
#PBS -q batch
#PBS -N AnnotateReads
#PBS -o logs/AnnotateReads_${PBS_JOBID}_stdout.log
#PBS -e logs/AnnotateReads_${PBS_JOBID}_stderr.log

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

cd $PBS_O_WORKDIR
source activate rhapsody

function CollectFiles {
    local folder=$1
    local suffix=$2

    local result
    while read line; do
        result="${result+${result},}${line}"
    done < <(find ${folder} -name "*${suffix}")

    printf "${result}"
}

INPUT_DIR=$1

FILTERING_STATS="$(CollectFiles ${INPUT_DIR} _read_quality.csv.gz)"
ANNOT_R1="$(CollectFiles ${INPUT_DIR} _Annotation_R1.csv.gz)"
ANNOT_R2="$(CollectFiles ${INPUT_DIR} _Triticum_aestivum_WTA_Annotation_R2.csv.gz)"
R2_QUALITY_METRICS="$(CollectFiles ${INPUT_DIR} _picard_quality_metrics.csv.gz)"
REFERENCE_PANEL_NAMES="$(find ${INPUT_DIR} -name reference_panel_names.json -print -quit)"

printf "Found %s filtering_stats\n" "$(echo ${FILTERING_STATS} | tr "," "\n" | wc -l)" 1>&2
printf "Found %s annot_r1\n" "$(echo ${ANNOT_R1} | tr "," "\n" | wc -l)" 1>&2
printf "Found %s annot_r2\n" "$(echo ${ANNOT_R2} | tr "," "\n" | wc -l)" 1>&2
printf "Found %s r2_quality_metrics\n" "$(echo ${R2_QUALITY_METRICS} | tr "," "\n" | wc -l)" 1>&2
printf "Found %s reference_panel_names\n" "$(echo ${REFERENCE_PANEL_NAMES} | tr "," "\n" | wc -l)" 1>&2

source env.sh

AnnotateReads.py \
  --filtering-stats "${FILTERING_STATS}" \
  --annotR1 "${ANNOT_R1}" \
  --annotR2 "${ANNOT_R2}" \
  --r2-quality-metrics "${R2_QUALITY_METRICS}" \
  --label-version 2 \
  --index data/annotations/Triticum_aestivum.IWGSC.48.gtf,data/reference_sequences/Triticum_aestivum.IWGSC.48.tar.gz \
  --reference-panel-names "${REFERENCE_PANEL_NAMES}"