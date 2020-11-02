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
from toil.jobStores.abstractJobStore import AbstractJobStore

try:
    import cPickle as pickle
except ImportError:
    import pickle


class PrintEncoder(json.JSONEncoder):
    def default(self, o):
        return repr(o)


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


def serialise_job(jobsotre: AbstractJobStore, job: Job, graph: JobGraph) -> None:
    with jobsotre.batch():
        with jobsotre.writeFileStream(graph.jobStoreID, cleanup=True) as (
            fileHandle,
            fileStoreID,
        ):
            pickle.dump(job, fileHandle, pickle.HIGHEST_PROTOCOL)
        userScript = job.getUserScript().globalize()
        graph.command = " ".join(("_toil", fileStoreID) + userScript.toCommand())
        jobsotre.update(graph)


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
        if options.attribute:
            attr_chain = options.attribute.split(".")
            p = job
            for attr in attr_chain:
                try:
                    p = getattr(p, attr)
                except AttributeError:
                    p = p[attr]
            print(p)
        else:
            print(vars(job))
    elif options.listfiles:
        files = recursiveGlob(jobstore.jobFilesDir, "*")
        files = list(filter(lambda f: not f.endswith(".new"), files))
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

    graph_attr = {}
    if options.cores:
        graph_attr["cores"] = options.cores
    if options.memory:
        graph_attr["memory"] = options.memory
    if options.disk:
        graph_attr["disk"] = options.disk

    job_attr = {}

    if options.job_id:
        job_graph = jobstore.load(options.job_id)
        if job_graph.command:
            job = Job._loadJob(job_graph.command, jobstore)
            for key_chain, val in options.attribute_pairs:
                attr_chain = key_chain.split(".")
                last_attr = attr_chain.pop()
                p = job
                for attr in attr_chain:
                    try:
                        p = getattr(p, attr)
                    except AttributeError:
                        p = p[attr]
                if val.startswith("json:"):
                    val = json.loads(val[len("json:") :])
                try:
                    old_val = getattr(p, last_attr)
                    setattr(p, last_attr, val)
                    new_val = getattr(p, last_attr)
                except AttributeError:
                    old_val = p[last_attr]
                    p[last_attr] = val
                    new_val = p[last_attr]
                print(
                    "Update job.{0} from {1} to {2}".format(
                        key_chain, old_val, new_val
                    ),
                    file=sys.stderr,
                )
            serialise_job(jobstore, job, job_graph)
        job_graph = update_job(job_graph, graph_attr)
        jobstore.update(job_graph)
        print("Job updated:", job_graph.jobStoreID, file=sys.stderr)

    if options.job_name:
        for job_graph in jobstore.jobs():
            if job_graph.jobName == options.job_name:
                jobstore.update(update_job(job_graph, graph_attr))
                print("Job updated:", job_graph.jobStoreID, file=sys.stderr)

    if options.config_file:
        with open(options.config_file, "r", encoding="UTF-8") as f:
            new_conf = json.load(f)
            update_config(jobstore.config, new_conf)
            jobstore.writeConfig()
            print("Jobstore config updated.", file=sys.stderr)


def export_action(options: Namespace) -> None:
    jobstore = Toil.resumeJobStore(options.job_store)
    if options.filename and options.output:
        localDirPath = os.path.dirname(options.output)
        if not localDirPath:
            options.output = os.path.join(os.getcwd(), options.output)
        jobstore.readFile(options.filename, options.output)


if __name__ == "__main__":

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument("job_store", type=str, metavar="JOB_STORE")

    parser = argparse.ArgumentParser(description="Update jobs within toil job store.")
    subparsers = parser.add_subparsers()

    show_cmd = subparsers.add_parser(
        "show", help="Show jobs or config.", parents=[parent_parser]
    )
    show_cmd.set_defaults(func=show_action)
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
    show_cmd.add_argument(
        "-t",
        "--attribute",
        dest="attribute",
        type=str,
        default=None,
        metavar="KEY",
        help="Show job with given attribute.",
    )

    update_cmd = subparsers.add_parser(
        "update", help="Update the jobstore.", parents=[parent_parser]
    )
    update_cmd.set_defaults(func=update_action)
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
    update_cmd.add_argument(
        "--cwlfile",
        dest="cwlfile",
        type=str,
        default=None,
        help="Update CWL job with a CWL file.",
    )
    update_cmd.add_argument(
        "-t",
        "--attribute-pair",
        dest="attribute_pairs",
        type=str,
        action="append",
        nargs=2,
        default=None,
        metavar=("KEY", "VAL"),
        help="Update job with given attributes.",
    )

    export_cmd = subparsers.add_parser(
        "export", help="Export file from jobstore.", parents=[parent_parser]
    )
    export_cmd.set_defaults(func=export_action)
    export_cmd.add_argument(
        "-f",
        "--filename",
        dest="filename",
        type=str,
        default=None,
        help="File name to export",
    )
    export_cmd.add_argument(
        "-o",
        "--output",
        dest="output",
        type=str,
        default=None,
        help="Output file name.",
    )

    options = parser.parse_args()
    options.func(options)
