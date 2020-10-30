#!/usr/bin/env python

import argparse
import sys
import json

import HTSeq


def processGFF(gff_filename, feature_type):
    gff = HTSeq.GFF_Reader(gff_filename)
    features = HTSeq.GenomicArrayOfSets("auto", True)
    counts = {}
    i = 0

    try:
        for f in gff:
            # print(f.type)
            if f.type == feature_type:
                try:
                    feature_id = f.attr["gene_name"]
                except KeyError:
                    try:
                        feature_id = f.attr["gene_id"]
                    except KeyError:
                        raise ValueError(
                            "Feature %s does not contain a gene_name or gene_id attribute"
                            % (f.name)
                        )

                if f.iv.strand == ".":
                    raise ValueError(
                        "Feature %s at %s does not have strand information but you are "
                        % (f.name, f.iv)
                    )

                features[f.iv] += feature_id
                counts[feature_id] = 0

            i += 1
            if i % 100000 == 0:
                print("%d GFF lines processed." % i, file=sys.stderr)
    except:
        print(
            "Error occured when processing GFF file (%s):\n"
            % gff.get_line_number_string(),
            file=sys.stderr,
        )
        raise
        print("%d GFF lines processed.\n" % i, file=sys.stderr)

    if len(counts) == 0:
        print(
            "Warning: No features of type '%s' found.\n" % feature_type, file=sys.stderr
        )

    return features, counts


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Process GFF")
    parser.add_argument("gff_filename")
    parser.add_argument("feature_type")

    options = parser.parse_args()

    features, counts = processGFF(options.gff_filename, options.feature_type)
    print(
        "Total %s of feature type '%s'" % (len(counts), options.feature_type),
        file=sys.stderr,
    )
    json.dump(counts, sys.stdout, indent=2)
