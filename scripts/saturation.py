#!/usr/bin/env python

import sys
import argparse
import time
import random
import csv

import pysam


def get_alignment_count(bamfile):
    stats_output = pysam.idxstats(bamfile).strip()
    total_align = 0
    for line in stats_output.split("\n"):
        segs = line.split("\t")
        total_align += (int(segs[2]) + int(segs[3]))
    return total_align


def get_raw_reads_count(bamfile):
    raw_reads = 0
    for align in samfile.fetch(until_eof=True):
        if not align.is_duplicate and not align.is_secondary:
            raw_reads += 1
    return raw_reads


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("bamfile", type=str, metavar="BAM_FILE")
    parser.add_argument("-v", "--verbose", action="store_true")
    options = parser.parse_args()

    samfile = pysam.AlignmentFile(options.bamfile, "rb")

    random.seed(round(time.time()))

    n_depth = 100
    buckets = [{"read_count": 0, "genes": set()} for i in range(n_depth)]
    i = 0
    for alignment in samfile.fetch(until_eof=True):
        # Select a bucket for each reads.
        bucket_idx = random.randrange(n_depth)
        buckets[bucket_idx]["read_count"] += 1

        if alignment.has_tag("XF") and alignment.has_tag("CN"):
            gene = alignment.get_tag("XF")
            cell = alignment.get_tag("CN")
            if cell == "T":
                buckets[bucket_idx]["genes"].add(gene)

        if options.verbose:
            i += 1
            if i % 1000000 == 0:
                sys.stderr.write("Processed %i alignments\r" % i)

    output_writer = csv.writer(sys.stdout)
    output_writer.writerow(["depths", "detected_genes"])

    depth = 0
    detected_genes = set()
    for buck in buckets:
        depth += buck["read_count"]
        detected_genes = detected_genes.union(buck["genes"])
        #output_writer.writerow([buck["read_count"], len(buck["genes"])])
        output_writer.writerow([depth, len(detected_genes)])
