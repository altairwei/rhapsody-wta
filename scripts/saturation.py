#!/usr/bin/env python

import sys
import argparse
import time
import random
import csv

import pysam


class ProgressBar(object):

    def __init__(self, total, prefix="", suffix="", ncol=60, file=sys.stderr):
        self.count = total
        self.prefix = prefix
        self.suffix = suffix
        self.ncol = ncol
        self.file = file
        self.done = 0
        self.last_time = time.time()

    def update(self, amount):
        self.done += amount
        x = int(self.ncol*self.done/self.count)
        self.file.write("%s[%s%s] %i/%i %s\r" % (
            self.prefix, "#"*x, " "*(self.ncol-x), self.done, self.count,
            self.suffix))
        self.file.flush()

    def __enter__(self):
        self.update(0)
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.file.write("%s[%s%s] %i/%i %s\r\n" % (
            self.prefix, "#"*self.ncol, " "*0, self.done, self.count,
            self.suffix))
        self.file.flush()


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
    options = parser.parse_args()

    samfile = pysam.AlignmentFile(options.bamfile, "rb")
    if not samfile.check_index():
        print("Error: can not find BAM index.")
        sys.exit(1)

    random.seed(round(time.time()))

    total = get_alignment_count(options.bamfile)
    n_depth = round(total/10000000)
    buckets = [{"num": 0, "genes": set()} for i in range(n_depth)]
    with ProgressBar(total, "Processing: ") as bar:
        for alignment in samfile.fetch(until_eof=True):
            try:
                gene = alignment.get_tag("XF")
            except KeyError:
                continue
            if not gene.startswith("__"):
                bucket_idx = random.randrange(n_depth)
                buckets[bucket_idx]["num"] += 1
                buckets[bucket_idx]["genes"].add(gene)
            bar.update(1)

    output_writer = csv.writer(sys.stdout)
    output_writer.writerow(["Depths", "Detected_Genes"])

    depth = 0
    detected_genes = set()
    for buck in buckets:
        depth += buck["num"]
        detected_genes = detected_genes.union(buck["genes"])
        output_writer.writerow([depth, len(detected_genes)])
