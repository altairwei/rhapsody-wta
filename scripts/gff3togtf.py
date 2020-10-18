#!/usr/bin/env python

import argparse
import sys
import json
from typing import Dict, List

import HTSeq
from HTSeq import GenomicFeature


def attr_to_string(attrs: Dict):
    attr_list = []
    for id_attr in ("gene_id", "transcript_id"):
        if id_attr in attrs:
            attr_list.append("{} \"{}\"".format(id_attr, attrs[id_attr]))
            del attrs[id_attr]
    attr_list.extend(
        "{} \"{}\"".format(str(key), str(val)) for (key, val) in attrs.items()
    )
    return "; ".join(attr_list)


def get_transcript_parent(gff3_file):
    print("Start to get transcript parent.", file=sys.stderr)
    transcript_parent = {}
    gff3 = HTSeq.GFF_Reader(options.gff3_file)
    i = 0
    for feature in gff3:
        if "Parent" in feature.attr:
            parent_type, parent_id = feature.attr["Parent"].split(":")
            if parent_type == "gene" and "transcript_id" in feature.attr:
                transcript_parent[feature.attr["transcript_id"]] = parent_id
        i += 1
        if i % 100000 == 0:
            print("%d GFF lines processed." % i, file=sys.stderr)
    return transcript_parent


def get_gtf_line(
    feature: GenomicFeature,
    transcript_parent: Dict[str, str],
    id_prefix: List[str],
    type_mapping: Dict[str, str]
):
    attr_dict: Dict[str, str] = {}
    attr_dict.update(feature.attr)

    for prefix in id_prefix:
        if prefix + "_id" in attr_dict:
            if not attr_dict[prefix + "_id"].startswith(prefix):
                attr_dict[prefix + "_id"] = prefix + ":" + attr_dict[prefix + "_id"]

    # Fill other necessary attributes
    if "Parent" in attr_dict:
        parent_type, parent_id = attr_dict["Parent"].split(":")
        if parent_type in id_prefix:
            full_id = attr_dict["Parent"]
        else:
            full_id = parent_id
        if parent_type == "transcript":
            attr_dict["transcript_id"] = full_id
        elif parent_type == "gene":
            attr_dict["gene_id"] = full_id
            if "transcript_id" in attr_dict:
                transcript_parent[attr_dict["transcript_id"]] = full_id

    if "gene_id" not in attr_dict:
        if feature.type != "chromosome":
            attr_dict["gene_id"] = transcript_parent[attr_dict["transcript_id"]]

    if feature.type in type_mapping:
        feature.type = type_mapping[feature.type]

    if feature.type == "exon":
        # Check for gene_id and transcript_id exists
        if "gene_id" not in attr_dict or "transcript_id" not in attr_dict:
            raise Exception("Exon must contain both 'gene_id' and 'transcript_id'")

    return "\t".join([
        feature.iv.chrom,
        feature.source,
        feature.type,
        str(feature.iv.start + 1), # See https://htseq.readthedocs.io/en/master/genomic.html#HTSeq.GenomicInterval
        str(feature.iv.end), # See https://htseq.readthedocs.io/en/master/genomic.html#HTSeq.GenomicInterval
        feature.score,
        feature.iv.strand,
        str(feature.frame),
        attr_to_string(attr_dict)
    ]) + "\n"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Converts Ensembl's favored GFF3 to GTF.")
    parser.add_argument("gff3_file", help="GFF3 file obtained from Ensembl.")
    parser.add_argument("gtf_file", help="Converted file in GTF format.")
    parser.add_argument("-p", "--reserve-id-prefix",
        dest="id_prefix", action="append", default=[],
        help="The id prefix to be reserved. Example `-p transcript -p gene`"
             " will produce 'transcript:G9200.1' and 'gene:G9200' as ids.")
    parser.add_argument("-t", "--transform-feature-type",
        dest="type_mapping", action="append", default=[],
        help="Rewrite the type of a feature to another one. Example `-t mRNA:transcript`"
             " will change the type of mRNA feature to transcript type.")
    options = parser.parse_args()

    type_mapping = {}
    for type_aes in options.type_mapping:
        old_type, new_type = type_aes.split(":")
        type_mapping[old_type] = new_type

    with open(options.gtf_file, "w", encoding="UTF-8") as fh:
        gff3 = HTSeq.GFF_Reader(options.gff3_file)
        i = 0
        transcript_parent = {}
        feature_type = {}
        for feature in gff3:
            feature_type[feature.type] = feature_type.get(feature.type, 0) + 1
            line = get_gtf_line(
                feature, transcript_parent, options.id_prefix, type_mapping)
            fh.write(line)
            i += 1
            if i % 100000 == 0:
                print("%d GFF lines processed." % i, file=sys.stderr)
        json.dump(feature_type, sys.stderr, indent=2)
        print("", file=sys.stderr)