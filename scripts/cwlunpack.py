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
    Any,
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

from cwltool.utils import CWLObjectType, CWLOutputType

import ruamel.yaml
from ruamel.yaml.comments import CommentedMap, CommentedSeq


def replace_ids(
    d: Union[CWLObjectType, CWLOutputType, MutableSequence[CWLObjectType], None],
    id_map: Dict[str, str],
) -> None:
    if isinstance(d, MutableSequence):
        for s in d:
            replace_ids(s, id_map)
    elif isinstance(d, MutableMapping):
        for i in ("id", "name"):
            if i in d and isinstance(d[i], str):
                new_id = id_map[d[i]]
                d[i] = new_id
        for s2 in d.values():
            replace_ids(s2, id_map)


def replace_runs(
    d: Union[CWLObjectType, CWLOutputType, MutableSequence[CWLObjectType], None],
    id_map: Dict[str, str],
) -> None:
    if isinstance(d, MutableSequence):
        for s in d:
            replace_runs(s, id_map)
    elif isinstance(d, MutableMapping):
        if "run" in d and isinstance(d["run"], str):
            new_id = id_map[d["run"]]
            d["run"] = new_id
        for s2 in d.values():
            replace_runs(s2, id_map)


def unpack_cwl(cwlfile, loadingContext):
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
        cwltool.pack.find_ids(loadingContext.loader.resolve_ref(f)[0], ids)
    
    mainpath, _ = urllib.parse.urldefrag(uri)

    #print("main uri:", uri)
    #for i in ids:
    #    print(i)
    # Rewrite id back to their original names.
    rewrite = {}
    names = set()  # type: Set[str]
    def rewrite_id(id_name: str, mainuri: str) -> None:
        if id_name == mainuri:
            print("Condition 1", id_name)
            # 重命名主文档，最外层的那个文档
            rewrite[id_name] = "#main"
            print("==========>", rewrite[id_name])
        elif id_name.startswith(mainuri) and id_name[len(mainuri)] in ("#", "/"):
            # 在主 cwl 文档内部，处理内部关系
            if id_name[len(mainuri) :].startswith("#main/"):
                # 假设用户自定义了一个 “main” id，那么就要把它转换成一个名字
                print("Condition 2-1", id_name)
                rewrite[id_name] = "#" + cwltool.process.uniquename(id_name[len(mainuri) + 1 :], names)
                print("============>", rewrite[id_name])
            else:
                # 没有 pack 的文档大部分都不会以 main 作为 id
                print("Condition 2-2", id_name)
                # 给原来的步骤冠上 #main 前缀。
                rewrite[id_name] = "#" + cwltool.process.uniquename("main/" + id_name[len(mainuri) + 1 :], names)
                print("============>", rewrite[id_name])
        else:
            # 外部的资源
            path, frag = urllib.parse.urldefrag(id_name)
            if path == mainpath:
                # 同一个 cwl 文档内，不同的 cwl 对象，即 id 不同。可能用于 packed 文档，那些 $graph 中的对象。
                # 他们处在同一个 cwl 文档内，但不属于同一个 cwl 对象。
                print("Condition 3-1", id_name)
                rewrite[id_name] = "#" + cwltool.process.uniquename(frag, names)
                print("============>", rewrite[id_name])
            else:
                # 引用外部 cwl 文件
                print("Condition 3-2", id_name)
                if path not in rewrite:
                    # 为什么这里不更改外部 cwl 的内部 id 呢？应该是通过 replace_refs 函数完成的。
                    rewrite[path] = "#" + cwltool.process.uniquename(cwltool.process.shortname(path), names)
                    print("============>", rewrite[path])

    rewrite_id_map = {}
    def split_id(id_name: str, mainuri: str, main_id: str = "main") -> None:
        if id_name == mainuri:
            # TODO: 改写成用户指定的 outname
            rewrite_id_map[id_name] = main_id
        elif id_name.startswith(mainuri) and id_name[len(mainuri)] in ("#", "/"):
            # 处理 main 文档内部的关系
            path, frag = urllib.parse.urldefrag(id_name)
            parts = frag.split("/")
            if parts[0] == "main":
                rewrite_id_map[id_name] = "/".join(parts[1:])
            else:
                rewrite_id_map[id_name] = frag
        else:
            path, frag = urllib.parse.urldefrag(id_name)
            # 处理打包成的外部 cwl 文件
            if path == mainpath:
                if frag.endswith(".cwl"):
                    # 主文档 id
                    rewrite_id_map[id_name] = frag[:-4]
                else:
                    parts = frag.split("/")
                    if parts[0].endswith(".cwl"):
                        rewrite_id_map[id_name] = "/".join(parts[1:])
                    else:
                        rewrite_id_map[id_name] = frag


    for id_name in sorted(ids):
        split_id(id_name, uri)

    rewrite_run_map = {}
    def rewrite_run(run_name: str, mainuri: str):
        path, frag = urllib.parse.urldefrag(run_name)
        if frag == "main":
            pass
        else:
            rewrite_run_map[run_name] = frag
    
    for run_name in sorted(runs):
        rewrite_run(run_name, uri)

    # Save each CWL object into separated file
    for obj in processobj:
        old_id = obj["id"]
        new_id = rewrite_id_map[old_id]
        replace_ids(obj, rewrite_id_map)
        replace_runs(obj, rewrite_run_map)

        for old_ref in list(rewrite_id_map.keys()):
            new_ref = rewrite_id_map[old_ref]
            # Skip CWL obj main id
            if old_ref != old_id:
                cwltool.pack.replace_refs(obj, rewrite_id_map, old_ref, new_ref)

        if new_id == "main":
            obj.update(metadata)
        with open(new_id + ".cwl", "w", encoding="UTF-8") as f:
            json.dump(obj, f, indent=4)


    #json.dump(rewrite_id_map, sys.stdout, indent=2)

    # Update workflow


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Unpack cwl.")
    parser.add_argument("cwlfile")
    parser.add_argument("-C", "--outdir", type=str, default=os.getcwd(),
        help="Output folder for the unpacked CWL files.")
    options = parser.parse_args()

    loading_context = cwltool.context.LoadingContext(vars(options))
    loading_context.construct_tool_object = cwltool.workflow.default_make_tool
    loading_context.resolver = cwltool.resolver.tool_resolver
    unpack_cwl(options.cwlfile, loading_context)
    
