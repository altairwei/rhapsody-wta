#!/usr/bin/env python3

import pathlib
import argparse
import json
import re
import sys
import os
import urllib.parse
from typing import (
    Union,
    Dict,
    List,
    Any,
    Set,
    Optional,
    MutableMapping,
    MutableSequence,
    Callable
)

import cwltool
import cwltool.main
import cwltool.pack
import cwltool.resolver
import cwltool.load_tool
import cwltool.context
import cwltool.workflow

from cwltool.utils import CWLObjectType, CWLOutputType

import ruamel.yaml
from ruamel.yaml import YAML
from ruamel.yaml.comments import CommentedMap, CommentedSeq



def traverse_cwl_objects(
    d: Union[CWLObjectType, CWLOutputType, MutableSequence[CWLObjectType], None],
    callback: Callable[[MutableMapping], None],
) -> None:
    if isinstance(d, MutableSequence):
        for s in d:
            traverse_cwl_objects(s, callback)
    elif isinstance(d, MutableMapping):
        callback(d)
        for s2 in d.values():
            traverse_cwl_objects(s2, callback)


def convert_id(id_name: str, mainuri: str) -> str:
    mainpath, _ = urllib.parse.urldefrag(mainuri)
    if id_name == mainuri:
        return "main"
    elif id_name.startswith(mainuri) and id_name[len(mainuri)] in ("#", "/"):
        # 处理 main 文档内部的关系
        path, frag = urllib.parse.urldefrag(id_name)
        parts = frag.split("/")
        if parts[0] == "main":
            if len(parts) > 2:
                return "/".join(parts[2:])
            else:
                return "/".join(parts[1:])
        else:
            return frag
    else:
        path, frag = urllib.parse.urldefrag(id_name)
        # 处理打包成的外部 cwl 文件
        if path == mainpath:
            if frag.endswith(".cwl"):
                # 主文档 id
                return frag[:-4]
            else:
                parts = frag.split("/")
                if parts[0].endswith(".cwl"):
                    if len(parts) > 2:
                        return "/".join(parts[2:])
                    else:
                        return "/".join(parts[1:])
                else:
                    return frag


def convert_run(run_name: str, mainuri: str):
    _, frag = urllib.parse.urldefrag(run_name)
    if frag == "main":
        return "main.cwl"
    else:
        return frag


def replace_ids(
    d: Union[CWLObjectType, CWLOutputType, MutableSequence[CWLObjectType], None],
    id_map: Dict[str, str],
) -> None:
    def callback(doc: MutableMapping):
        for i in ("id", "name"):
            if i in doc and isinstance(doc[i], str):
                new_id = id_map[doc[i]]
                doc[i] = new_id
    traverse_cwl_objects(d, callback)


def replace_runs(
    d: Union[CWLObjectType, CWLOutputType, MutableSequence[CWLObjectType], None],
    id_map: Dict[str, str],
) -> None:
    def callback(doc: MutableMapping):
        if "run" in doc and isinstance(doc["run"], str):
            new_id = id_map[doc["run"]]
            doc["run"] = new_id
    traverse_cwl_objects(d, callback)


def find_refs(
    d: Union[CWLObjectType, CWLOutputType, MutableSequence[CWLObjectType], None],
    refs: Set[str],
    stem: str
) -> None:
    def callback(doc: MutableMapping):
        for key, val in doc.items():
            if key not in ("id", "name", "run") and \
                isinstance(val, str) and val.startswith(stem):
                refs.add(val)
    traverse_cwl_objects(d, callback)


def replace_refs(d: Any, rewrite: Dict[str, str], stem: str, newstem: str) -> None:
    if isinstance(d, MutableSequence):
        for s, v in enumerate(d):
            if isinstance(v, str):
                if v in rewrite:
                    d[s] = rewrite[v]
                elif v.startswith(stem):
                    d[s] = newstem + v[len(stem) :]
                    rewrite[v] = d[s]
            else:
                replace_refs(v, rewrite, stem, newstem)
    elif isinstance(d, MutableMapping):
        for key, val in d.items():
            if isinstance(val, str):
                if val in rewrite:
                    d[key] = rewrite[val]
                elif val.startswith(stem):
                    id_ = val[len(stem) :]
                    # prevent appending newstems if tool is already packed
                    if id_.startswith(newstem.strip("#")):
                        d[key] = "#" + id_
                    else:
                        d[key] = newstem + id_
                    rewrite[val] = d[key]
            elif key == "out" and isinstance(val, MutableSequence):
                new_out = []
                for out_ref in val:
                    if out_ref.startswith(stem):
                        _, frag = urllib.parse.urldefrag(out_ref)
                        new_out.append(frag.split("/")[-1])
                if new_out:
                    d[key] = new_out
            else:
                replace_refs(val, rewrite, stem, newstem)


def yaml_dump(obj: CWLObjectType, filename: str):
    yaml=YAML()
    yaml.default_flow_style = False
    with open(filename, "w", encoding="UTF-8") as f:
        formated_obj = json.loads(json.dumps(obj, indent=4))
        yaml.dump(formated_obj, f)


def json_dump(obj: CWLObjectType, filename: str):
        with open(filename, "w", encoding="UTF-8") as fh:
            json.dump(obj, fh, indent=4)


def unpack_cwl(cwlfile, loadingContext) -> Dict[str, CWLObjectType]:
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

    mainpath, _ = urllib.parse.urldefrag(uri)

    ids = set()  # type: Set[str]
    cwltool.pack.find_ids(processobj, ids)

    runs = {uri}
    cwltool.pack.find_run(processobj, loadref, runs)

    for f in runs:
        cwltool.pack.find_ids(loadingContext.loader.resolve_ref(f)[0], ids)

    refs = set()
    find_refs(processobj, refs, mainpath)

    id_map = {}
    for id_name in sorted(ids):
        id_map[id_name] = convert_id(id_name, uri)

    run_map = {}  
    for run_name in sorted(runs):
        run_map[run_name] = convert_run(run_name, uri)

    # Save each CWL object into separated file
    results = {}
    for obj in processobj:
        old_id = obj["id"]
        new_id = id_map[old_id]
        replace_ids(obj, id_map)
        replace_runs(obj, run_map)

        for old_ref in list(id_map.keys()):
            new_ref = id_map[old_ref]
            # Skip CWL obj main id
            if old_ref != old_id:
                replace_refs(obj, id_map, old_ref, new_ref)

        if "http://commonwl.org/cwltool#original_cwlVersion" in obj:
            del obj["http://commonwl.org/cwltool#original_cwlVersion"]

        obj["cwlVersion"] = metadata["cwlVersion"]

        if new_id == "main":
            obj.update(metadata)
            del obj["id"]

        results[new_id] = obj
    
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Unpack cwl.")
    parser.add_argument("cwlfile")
    parser.add_argument("-f", "--output-format", choices=["json", "yaml"],
        type=str, default="json", help="Decide the output cwl file format.")
    parser.add_argument("-C", "--outdir", type=str, default=os.getcwd(),
        help="Output folder for the unpacked CWL files.")
    options = parser.parse_args()

    loading_context = cwltool.context.LoadingContext(vars(options))
    loading_context.construct_tool_object = cwltool.workflow.default_make_tool
    loading_context.resolver = cwltool.resolver.tool_resolver
    cwl_obj_map = unpack_cwl(options.cwlfile, loading_context)

    os.makedirs(options.outdir, exist_ok=True)
    for obj_id in cwl_obj_map:
        cwl_obj_filename = os.path.join(options.outdir, obj_id + ".cwl")
        cwl_obj = cwl_obj_map[obj_id]
        if options.output_format == "yaml":
            yaml_dump(cwl_obj, cwl_obj_filename)
        elif options.output_format == "json":
            json_dump(cwl_obj, cwl_obj_filename)
        else:
            raise Exception("Unknown output format.")

    
