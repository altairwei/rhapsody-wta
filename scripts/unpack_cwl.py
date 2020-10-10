#!/usr/bin/env python3

import pathlib
import argparse
import json
import re
import sys

import cwltool
import cwltool.main
import cwltool.resolver
import cwltool.load_tool
import cwltool.context
import cwltool.workflow


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
    processobj, metadata = loadingContext.loader.resolve_ref(uri)
    #with open(outname, "w", encoding="UTF-8") as f:
    #    f.write(cwltool.main.print_pack(loadingContext, uri))
    print(processobj)
    for po in processobj:
        print(po["id"])


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Unpack cwl.")
    parser.add_argument("cwlfile")
    parser.add_argument("outname")
    options = parser.parse_args()

    loading_context = cwltool.context.LoadingContext(vars(options))
    loading_context.construct_tool_object = cwltool.workflow.default_make_tool
    loading_context.resolver = cwltool.resolver.tool_resolver
    unpack_cwl(options.cwlfile, options.outname, loading_context)
    