#!/usr/bin/env bash

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

print_help() {
  echo "Script usage: $(basename $0) [-r] [-s] [-c] [-w] [-h] workflow sample" >&2
}

while getopts 'rscw:h' OPTION; do
  case "$OPTION" in
    r)
      echo "Restarting the workflow..."
      RESTART_OPT="--restart"
      ;;
    s)
      echo "Run toil leader within a named screen session..."
      USE_SCREEN=true
      ;;
    c)
      echo "Clean job store before running workflow.."
      CLEAN_JOBSTORE=true
      ;;
    w)
      echo "Override CWL workflow requirements.."
      OVERRIDES_OPT="--overrides ${OPTARG}"
      ;;
    h)
      print_help
      exit
      ;;
    ?)
      print_help
      exit 1
      ;;
  esac
done
shift "$(($OPTIND -1))"

# Positional arguments
WORKFLOW=$1
WORKFLOW_INPUT=$2

# Infer the sample name from the input
SAMPLE_NAME=$(basename "${WORKFLOW_INPUT}")
SAMPLE_NAME=${SAMPLE_NAME%.*}
RUN_ID=$(date +%Y%m%d_%H%M%S)

# Create folders for workflow
RESULTS_FOLDER=results/$SAMPLE_NAME
LOG_FOLDER=logs/$(date +%Y%m%d)
export TMPDIR=tmp/$SAMPLE_NAME/tmpwork
mkdir -p $RESULTS_FOLDER
mkdir -p $LOG_FOLDER
mkdir -p $TMPDIR

# Arguments for toil
JOBSTORE=tmp/$SAMPLE_NAME/jobstore
LOG_FILE=$LOG_FOLDER/cwltoil_${SAMPLE_NAME}_${RUN_ID}.log

# Environment variables for toil
export TOIL_TORQUE_ARGS="-q batch"
export TOIL_TORQUE_REQS="walltime=72:00:00"
export CWL_SINGULARITY_CACHE="$PWD/dockerImages"

TOIL_CMD="toil-cwl-runner ${RESTART_OPT} ${OVERRIDES_OPT} \
  --maxLocalJobs 100 \
  --batchSystem torque \
  --singularity \
  --defaultDisk 140G \
  --jobStore file:${JOBSTORE} \
  --outdir $RESULTS_FOLDER \
  --workDir $TMPDIR \
  --writeLogs logs \
  --logFile $LOG_FILE \
  --disableCaching \
  --logLevel INFO \
  --retryCount 3 \
  --maxLogFileSize 20000000000 \
  --stats \
  $WORKFLOW $WORKFLOW_INPUT"

echo "Options:"
echo "- USE_SCREEN: ${USE_SCREEN-false}"
echo "- WORKFLOW: ${WORKFLOW}"
echo "- WORKFLOW_INPUT: ${WORKFLOW_INPUT}"
echo "- SAMPLE_NAME: ${SAMPLE_NAME}"
echo "- RUN_ID: ${RUN_ID}"
echo "- RESTART_OPT: ${RESTART_OPT-None}"
echo "- OVERRIDES_OPT: ${OVERRIDES_OPT-None}"
echo "- RESULTS_FOLDER: ${RESULTS_FOLDER}"
echo "- JOBSTORE: ${JOBSTORE}"
echo "- LOG_FILE: ${LOG_FILE}"
echo "- CMD: ${TOIL_CMD}"

source activate rhapsody

if [ "$CLEAN_JOBSTORE" = true ]; then
  toil clean ${JOBSTORE}
fi

if [ "$USE_SCREEN" = true ]; then
  screen -S ${SAMPLE_NAME} -d -m ${TOIL_CMD}
  echo "You can enter the screen session by command: screen -r ${SAMPLE_NAME}"
else
  ${TOIL_CMD}
fi