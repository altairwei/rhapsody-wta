#!/usr/bin/env python

import sys
import argparse
from argparse import ArgumentError
from typing import Dict

from toil.jobStores.fileJobStore import FileJobStore
from toil.jobGraph import JobGraph

def bytes_fmt(num, suffix='B'):
    for unit in ['','K','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, 'Y', suffix)


def print_job(job: JobGraph):
    state = ", ".join([
        "Cores: " + str(job.cores),
        "Memory: " + bytes_fmt(job.memory),
        "Disk: " + bytes_fmt(job.disk),
        "Retry: " + str(job.remainingRetryCount)
    ])
    print("ID: " + job.jobStoreID)
    print("Name: " + job.jobName)
    print(state)
    print()


def parseLocator(locator):
    if locator[0] in '/.' or ':' not in locator:
        return 'file', locator
    else:
        try:
            name, rest = locator.split(':', 1)
        except ValueError:
            raise RuntimeError('Invalid job store locator syntax.')
        else:
            return name, rest


def update_job(job: JobGraph, attr: Dict) -> JobGraph:
    job._cores = attr.get("cores", job._cores)
    job._memory = attr.get("memory", job._memory)
    job._disk = attr.get("disk", job._disk)
    return job


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update jobs within toil job store.")
    parser.add_argument("-j", "--jobStore", dest="job_store", type=str, metavar='STRING', required=True)
    parser.add_argument("-n", "--jobName", dest="job_name", type=str, metavar='STRING')
    parser.add_argument("-i", "--jobID", dest="job_id", type=str, metavar='STRING')

    provgroup = parser.add_argument_group("Options for updating jobs' attibutes.")
    provgroup.add_argument('-c', '--cores', dest='cores', type=int, default=None, metavar='FLOAT')
    provgroup.add_argument('-m', '--memory', dest='memory', type=int, default=None, metavar='INT')
    provgroup.add_argument('-d', '--disk', dest='disk', type=int, default=None, metavar='INT')

    options = parser.parse_args()

    name, rest = parseLocator(options.job_store)
    if name != "file":
        raise ArgumentError("--jobStore", "Support for local file job store only.")
    jobstore = FileJobStore(rest)

    attr = {}
    if options.cores:
        attr["cores"] = options.cores
    if options.memory:
        attr["memory"] = options.memory
    if options.disk:
        attr["disk"] = options.disk

    if options.job_id:
        job = jobstore.load(options.job_id)
        job = update_job(job, attr)
        jobstore.update(job)
        print("Job updated:", job.jobStoreID)

    if options.job_name:
        for job in jobstore.jobs():
            if job.jobName == options.job_name:
                jobstore.update(update_job(job, attr))
                print("Job updated:", job.jobStoreID)

    # The default behavior is to show the job states
    if not options.job_name and not options.job_id:
        for job in jobstore.jobs():
            print_job(job)
