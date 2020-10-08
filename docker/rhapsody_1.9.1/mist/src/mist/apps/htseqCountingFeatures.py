#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
09/07/2016: Modified from htseq-counting to annotate 'no feature' reads
and to include default htseq settings to include headers in SAM files
Internal pipeline for Resolve Analysis WTA v0.1.1 as of 11/25/2015 (on VM cellinux01)

'''

import os
from collections import defaultdict
import HTSeq
import logging


class UnknownChrom( Exception ):
    pass

def invert_strand( iv ):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2

def getFeatures(iv_seq, features):
    fs = set()
    for iv in iv_seq:
        if iv.chrom not in features.chrom_vectors:
            return "UnknownChrom"
        for iv2, fs2 in features[ iv ].steps():
            fs = fs.union( fs2 )
    return fs

def getAlignment(r, stranded):
    if stranded != "reverse":
        iv_seq = ( co.ref_iv for co in r.cigar if co.type == "M" and co.size > 0 )
    else:
        iv_seq = ( invert_strand( co.ref_iv ) for co in r.cigar if co.type == "M" and co.size > 0 )
    return iv_seq


def count_reads_in_features(bamout, bam_filename, gff_filename, total_counts_file):
    groups = defaultdict().fromkeys(['empty', 'ambiguous', 'notaligned', 'lowqual', 'nonunique'], 0)
    # subgroup "no feature" reads under default htseqcounts
    groups['nofeature'] = {k:0 for k in ['__introns_','__antisense_exons_','__antisense_introns_','__intergenic']}

    def getAssignment(r):
        try:
            for subgr in ['htseq', '__introns_', '__antisense_exons_', '__antisense_introns_']:
                if subgr == 'htseq':
                    iv_seq = getAlignment(r, 'yes')
                    fs = getFeatures(iv_seq, features_exon)
                elif subgr == '__introns_':
                    iv_seq = getAlignment(r, 'yes')
                    fs = getFeatures(iv_seq, features_gene)
                elif subgr == '__antisense_exons_':
                    iv_seq = getAlignment(r, 'reverse')
                    fs = getFeatures(iv_seq, features_exon)
                elif subgr == '__antisense_introns_':
                    iv_seq = getAlignment(r, 'reverse')
                    fs = getFeatures(iv_seq, features_gene)

                if fs == "UnknownChrom": raise UnknownChrom
                if fs is not None and len(fs) >= 1: break

            if (fs is None or len(fs) == 0):
                writeout = "__intergenic"
                groups['nofeature']['__intergenic'] += 1
            elif len(fs) > 1:
                if subgr == 'htseq':
                    writeout = "__ambiguous[" + '+'.join(fs) + "]"
                    groups['ambiguous'] += 1
                else:
                    # ambigously mapped subgrouped reads
                    writeout = "__ambiguous" + subgr[2:] + "[" + '+'.join(fs) + "]"
                    groups['nofeature'][subgr] += 1
            else:
                counts_gene[list(fs)[0]] += 1
                if subgr == 'htseq':
                    writeout = list(fs)[0]
                    counts_exon[list(fs)[0]] += 1
                else:
                    writeout = subgr + list(fs)[0]
                    groups['nofeature'][subgr] += 1
        except UnknownChrom:
            writeout = "__no_feature"
            groups['empty'] += 1

        return writeout

    def processGFF(gff_filename, feature_type):
        gff = HTSeq.GFF_Reader(gff_filename)
        features = HTSeq.GenomicArrayOfSets("auto", True)
        counts = {}
        i = 0

        try:
            for f in gff:
                if f.type == feature_type:
                    try:
                        feature_id = f.attr['gene_name']
                    except KeyError:
                        try:
                            feature_id = f.attr['gene_id']
                        except KeyError:
                            raise ValueError("Feature %s does not contain a gene_name or gene_id attribute" % (f.name))

                    if f.iv.strand == ".":
                        raise ValueError("Feature %s at %s does not have strand information but you are " % (f.name, f.iv))

                    features[f.iv] += feature_id
                    counts[feature_id] = 0

                i += 1
                if i % 100000 == 0:
                    logging.info("%d GFF lines processed.\n" % i)
        except:
            logging.info("Error occured when processing GFF file (%s):\n" % gff.get_line_number_string() )
            raise
        logging.info("%d GFF lines processed.\n" % i)

        if len( counts ) == 0:
            logging.info( "Warning: No features of type '%s' found.\n" % feature_type )

        return features, counts

    logging.info( "Process GFF for feature type exon.\n")
    logging.info(os.path.isfile(gff_filename))
    features_exon, counts_exon = processGFF(gff_filename, 'exon')
    features_gene, counts_gene = processGFF(gff_filename, 'gene')

    try:
        read_seq = HTSeq.BAM_Reader(bam_filename)
    except:
        logging.info( "Error occured when reading beginning of BAM file.\n" )
        raise

    bam_writer = HTSeq.BAM_Writer.from_BAM_Reader(bamout, read_seq)

    i = 0
    for r in read_seq:
        if i > 0 and i % 100000 == 0:
            logging.info( "%d SAM alignment records processed.\n" % i)

        i += 1
        if i == 1:
            if r.paired_end: # if pe mode
                raise ValueError("paired alignment.\n")

        if not r.aligned:
            groups['notaligned'] += 1
            r.optional_fields.append(('XF', '__not_aligned'))
            bam_writer.write(r)
            continue
        writeout = getAssignment(r)
        r.optional_fields.append(('XF', writeout))
        bam_writer.write(r)


    totaloutfile = open(total_counts_file, "w")
    for fn in sorted( counts_exon.keys() ):
        totaloutfile.write("%s\t%d" % ( fn, counts_exon[fn] ) + "\n")
    totaloutfile.close()

    logging.info("__no_feature\t%d" % groups['empty'])
    logging.info("__ambiguous\t%d" % groups['ambiguous'])
    logging.info("__not_aligned\t%d" % groups['notaligned'])
    logging.info("__alignment_not_unique\t%d" % groups['nonunique'])
    logging.info("__nofeature_introns\t%d" % groups['nofeature']['__introns_'])
    logging.info("__nofeature_antisense_exons\t%d" % groups['nofeature']['__antisense_exons_'])
    logging.info("__nofeature_antisense_introns\t%d" % groups['nofeature']['__antisense_introns_'])
    logging.info("__nofeature_intergenic\t%d" % groups['nofeature']['__intergenic'])

    return bamout, total_counts_file