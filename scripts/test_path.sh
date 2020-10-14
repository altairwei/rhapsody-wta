#!/usr/bin/env bash
#PBS -q low
#PBS -N test-path
#PBS -o /public/home/weiwanqian/src/rhapsody-wta/toil_workflow_57b6a086-6d74-4578-8941-0ce2401e1c48_job_8_batch_torque_${PBS_JOBID}_std_output.log
#PBS -e /public/home/weiwanqian/src/rhapsody-wta/toil_workflow_57b6a086-6d74-4578-8941-0ce2401e1c48_job_8_batch_torque_${PBS_JOBID}_std_error.log

set -e -o pipefail
shopt -s failglob
export LC_ALL=C

echo "~~~~~~~~~~~~~~~~~~ Hello World ~~~~~~~~~~~~~~~~~~~"
echo "/tmp/toil_workflow_57b6a086-6d74-4578-8941-0ce2401e1c48_job_8_batch_torque_${PBS_JOBID}_std_output.log"
echo "/tmp/toil_workflow_57b6a086-6d74-4578-8941-0ce2401e1c48_job_8_batch_torque_${PBS_JOBID}_std_error.log"
echo "$PBS_JOBID"