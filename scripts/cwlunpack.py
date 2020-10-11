#!/usr/bin/env python3

import pathlib
import argparse
import json
import re
import sys
import os
import urllib.parse
from typing import (
    Optional,
    MutableMapping,
    MutableSequence
)

import cwltool
import cwltool.main
import cwltool.pack
import cwltool.resolver
import cwltool.load_tool
import cwltool.context
import cwltool.workflow

import ruamel.yaml
from ruamel.yaml.comments import CommentedMap, CommentedSeq


def unpack_cwl(cwlfile, outname, loadingContext):
    uri, tool_file_uri = cwltool.load_tool.resolve_tool_uri(
        cwlfile,
        resolver=loadingContext.resolver,
        fetcher_constructor=loadingContext.fetcher_constructor,
    )
    loadingContext, workflowobj, uri = cwltool.load_tool.fetch_document(uri, loadingContext)
    loadingContext, uri = cwltool.load_tool.resolve_and_validate_document(
        loadingContext,
        workflowobj,
        uri,
        preprocess_only=True
    )
    # Responsible for parsing the entire document
    # For those `packed` CWL documents, we will get MutableSequence
    # In fact, the content of `$graph` is MutableSequence.
    # Metadata is everything except for `$graph` for `packed` CWL 
    processobj, metadata = loadingContext.loader.resolve_ref(uri)

    yaml = ruamel.yaml.YAML()
    idx = {} # This is used to store parsed documents.
    if isinstance(processobj, MutableMapping):
        # This is a normal handwritten CWL document.
        idx[processobj["id"]] = CommentedMap(processobj.items())
    elif isinstance(processobj, MutableSequence):
        # This is a `packed` CWL document.
        _, frag = urllib.parse.urldefrag(uri)
        for po in processobj:
            if not frag:
                if po["id"].endswith("#main"):
                    # Find the main workflow
                    uri = po["id"]
            idx[po["id"]] = CommentedMap(po.items())
        idx[metadata["id"]] = CommentedMap(metadata.items())
    
    # Find all the ids
    def loadref(base: Optional[str], lr_uri: str):
        lr_loadingContext = loadingContext.copy()
        lr_loadingContext.metadata = {}
        lr_loadingContext, lr_workflowobj, lr_uri = cwltool.load_tool.fetch_document(
            lr_uri, lr_loadingContext
        )
        lr_loadingContext, lr_uri = cwltool.load_tool.resolve_and_validate_document(
            lr_loadingContext, lr_workflowobj, lr_uri
        )
        if lr_loadingContext.loader is None:
            raise Exception("loader should not be None")
        return lr_loadingContext.loader.resolve_ref(lr_uri, base_url=base)[0]

    ids = set()  # type: Set[str]
    cwltool.pack.find_ids(processobj, ids)

    runs = {uri}
    cwltool.pack.find_run(processobj, loadref, runs)

    for f in runs:
        find_ids(loadingContext.loader.resolve_ref(f)[0], ids)
    
    mainpath, _ = urllib.parse.urldefrag(uri)

    # Rewrite id back to their original names.
    rewrite = {}
    def rewrite_id(id_name: str, mainuri: str) -> None:
        if id_name == mainuri:
            rewrite[id_name] = "#main"
        elif id_name.startswith(mainuri) and id_name[len(mainuri)] in ("#", "/"):
            if id_name[len(mainuri) :].startswith("#main/"):
                rewrite[id_name] = "#" + uniquename(id_name[len(mainuri) + 1 :], names)
            else:
                rewrite[id_name] = "#" + uniquename("main/" + id_name[len(mainuri) + 1 :], names)
        else:
            path, frag = urllib.parse.urldefrag(id_name)
            if path == mainpath:
                rewrite[id_name] = "#" + uniquename(frag, names)
            else:
                if path not in rewrite:
                    rewrite[path] = "#" + uniquename(shortname(path), names)

    for id_name in sorted(ids):
        rewrite_id(id_name, uri)

    # Update workflow


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Unpack cwl.")
    parser.add_argument("cwlfile")
    parser.add_argument("outname")
    parser.add_argument("-C", "--outdir", type=str, default=os.getcwd(),
        help="Output folder for the unpacked CWL files.")
    options = parser.parse_args()

    loading_context = cwltool.context.LoadingContext(vars(options))
    loading_context.construct_tool_object = cwltool.workflow.default_make_tool
    loading_context.resolver = cwltool.resolver.tool_resolver
    unpack_cwl(options.cwlfile, options.outname, loading_context)
    
