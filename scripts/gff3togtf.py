#!/usr/bin/env python

import argparse
import sys
from typing import Dict

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


def get_gtf_line(feature: GenomicFeature, transcript_parent: Dict):
    attr_dict = {}
    attr_dict.update(feature.attr)

    # Fill other necessary attributes
    if "Parent" in attr_dict:
        parent_type, parent_id = attr_dict["Parent"].split(":")
        if parent_type == "transcript":
            attr_dict["transcript_id"] = parent_id
        elif parent_type == "gene":
            attr_dict["gene_id"] = parent_id
            if "transcript_id" in attr_dict:
                transcript_parent[attr_dict["transcript_id"]] = parent_id

    if "gene_id" not in attr_dict:
        try:
            if feature.type != "chromosome":
                attr_dict["gene_id"] = transcript_parent[attr_dict["transcript_id"]]
        except KeyError:
            print(transcript_parent)
            raise

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
    parser = argparse.ArgumentParser("Convert GFF3 to GTF.")
    parser.add_argument("gff3_file")
    parser.add_argument("gtf_file")
    options = parser.parse_args()

    transcript_parent = get_transcript_parent(options.gff3_file)
    with open(options.gtf_file, "w", encoding="UTF-8") as fh:
        print("Start to write GTF.", file=sys.stderr)
        gff3 = HTSeq.GFF_Reader(options.gff3_file)
        i = 0
        for feature in gff3:
            fh.write(get_gtf_line(feature, transcript_parent))
            i += 1
            if i % 100000 == 0:
                print("%d GFF lines processed." % i, file=sys.stderr)