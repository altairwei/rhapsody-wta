#!/usr/bin/env python

import os
import sys
import json
import argparse
from argparse import ArgumentError, Namespace
from typing import Dict

from toil.jobStores.fileJobStore import FileJobStore
from toil.jobGraph import JobGraph
from toil.job import Job
from toil.common import Toil, Config
from toil.utils.toilDebugFile import recursiveGlob


def bytes_fmt(num, suffix="B"):
    for unit in ["", "K", "M", "G", "T", "P", "E", "Z"]:
        if abs(num) < 1024.0:
            return "%3.1f%s%s" % (num, unit, suffix)
        num /= 1024.0
    return "%.1f%s%s" % (num, "Y", suffix)


def print_job(job: JobGraph):
    print("ID:", job.jobStoreID)
    print("Name:", job.jobName)
    print("Children:", len(job.stack))
    print("Command:", job.command)
    print(
        ", ".join(
            [
                "Cores: " + str(job.cores),
                "Memory: " + bytes_fmt(job.memory),
                "Disk: " + bytes_fmt(job.disk),
                "Retry: " + str(job.remainingRetryCount),
            ]
        )
    )
    print()


def update_job(job: JobGraph, attr: Dict) -> JobGraph:
    job._cores = attr.get("cores", job._cores)
    job._memory = attr.get("memory", job._memory)
    job._disk = attr.get("disk", job._disk)
    return job


def update_config(conf: Config, new_obj: Dict) -> None:
    for name in new_obj:
        setattr(conf, name, new_obj[name])


def show_action(options: Namespace) -> None:
    jobstore = Toil.resumeJobStore(options.job_store)
    if options.print_config:
        json.dump(vars(jobstore.config), sys.stdout, indent=2)
        print()
    elif options.job_id:
        job_graph = jobstore.load(options.job_id)
        job = Job._loadJob(job_graph.command, jobstore)
    elif options.listfiles:
        files = recursiveGlob(jobstore.jobFilesDir, "*")
        files = list(filter(lambda f: not f.endswith('.new'), files))
        try:
            for file_item in files:
                print(file_item)
        except BrokenPipeError:
            pass
    else:
        for job in jobstore.jobs():
            print_job(job)


def update_action(options: Namespace) -> None:
    jobstore = Toil.resumeJobStore(options.job_store)

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
        print("Job updated:", job.jobStoreID, file=sys.stderr)

    if options.job_name:
        for job in jobstore.jobs():
            if job.jobName == options.job_name:
                jobstore.update(update_job(job, attr))
                print("Job updated:", job.jobStoreID, file=sys.stderr)

    if options.config_file:
        with open(options.config_file, "r", encoding="UTF-8") as f:
            new_conf = json.load(f)
            update_config(jobstore.config, new_conf)
            jobstore.writeConfig()
            print("Jobstore config updated.", file=sys.stderr)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Update jobs within toil job store.")
    subparsers = parser.add_subparsers()

    show_cmd = subparsers.add_parser("show", help="Show jobs or config.")
    show_cmd.add_argument("job_store", type=str, metavar="JOB_STORE")
    show_cmd.add_argument("-i", "--jobID", dest="job_id", type=str, metavar="STRING")
    show_cmd.add_argument(
        "--print-config",
        dest="print_config",
        action="store_true",
        default=False,
        help="Print the config of workflow.",
    )
    show_cmd.add_argument(
        "-l",
        "--list-files",
        dest="listfiles",
        action="store_true",
        default=False,
        help="List all files within jobstore.",
    )
    show_cmd.set_defaults(func=show_action)

    update_cmd = subparsers.add_parser("update")
    update_cmd.add_argument("job_store", type=str, metavar="JOB_STORE")
    update_cmd.add_argument(
        "-n", "--jobName", dest="job_name", type=str, metavar="STRING"
    )
    update_cmd.add_argument("-i", "--jobID", dest="job_id", type=str, metavar="STRING")
    update_cmd.add_argument(
        "--update-config",
        dest="config_file",
        type=str,
        metavar="CONFIG_FILE",
        default=None,
        help="Update workflow config using a json file",
    )
    update_cmd.add_argument(
        "-c",
        "--cores",
        dest="cores",
        type=int,
        default=None,
        metavar="FLOAT",
        help="Update cores requirement",
    )
    update_cmd.add_argument(
        "-m",
        "--memory",
        dest="memory",
        type=int,
        default=None,
        metavar="INT",
        help="Update memory requirement",
    )
    update_cmd.add_argument(
        "-d",
        "--disk",
        dest="disk",
        type=int,
        default=None,
        metavar="INT",
        help="Update disk requirement",
    )

    options = parser.parse_args()
    options.func(options)
