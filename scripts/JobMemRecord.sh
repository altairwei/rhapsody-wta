#!/bin/bash

# Check arguments
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <job_id> <log_file>"
    exit 1
fi

JOB_ID=$1
LOG_FILE=$2

# Define qmem function
function qmem() {
  local jobid=${1}
  local mem=$(qstat -f -1 ${jobid} | sed -rn '5s/.*= ([[:digit:]]+)kb/\1/p')
  numfmt --from iec --to iec "${mem}K"
}

# Check job status
function check_job_status() {
  local status=$(qstat -f ${JOB_ID} | grep "job_state" | awk '{print $3}')
  echo $status
}

# Monitor memory usage every 10 minutes
while true; do
    job_status=$(check_job_status)

    # If job is not running or completed, exit the loop
    if [ "$job_status" != "R" ]; then
        if [ "$job_status" == "C" ]; then
            echo "Job ${JOB_ID} has completed. Stopping memory usage monitoring." >> ${LOG_FILE}
        else
            echo "Job ${JOB_ID} is no longer in a running state. Stopping memory usage monitoring." >> ${LOG_FILE}
        fi
        break
    fi

    # Get memory usage
    MEMORY_USAGE=$(qmem ${JOB_ID})
    # Get current time
    CURRENT_TIME=$(date '+%Y-%m-%d %H:%M:%S')

    # Log time and memory usage to file
    echo "${CURRENT_TIME}: ${MEMORY_USAGE}" >> ${LOG_FILE}

    # Wait for 10 minutes
    sleep 600
done
