#!/usr/bin/env python

import sys
import os
import glob
import click
import logging
import pysam
from typing import *
from velocyto.commands._run import _run


logging.basicConfig(stream=sys.stdout, format='%(asctime)s - %(levelname)s - %(message)s', level=logging.DEBUG)


def rhapsody_wta_correct(bamfilepath: str, corrected_output: str):
    parentpath, bamfilename = os.path.split(bamfilepath)
    logging.info(f"Correcting {bamfilename}")
    infile = pysam.AlignmentFile(bamfilepath, mode="rb")
    if corrected_output is None:
        bam_out_path = os.path.join(parentpath, f"correct_{bamfilename}")
    else:
        bam_out_path = corrected_output
    outfile = pysam.AlignmentFile(bam_out_path, mode="wb", template=infile)

    i = 0
    for read in infile:
        if read.has_tag("CN"):
            is_putative_cell = read.get_tag("CN")
            if is_putative_cell == "T":
                RSEC_adjusted_UMI = read.get_tag("MA")
                read.set_tag("UB", RSEC_adjusted_UMI, value_type="Z")
                outfile.write(read)
        i += 1
        if i % 100000000 == 0:
            print("%d sligned segment processed." % i, file=sys.stderr)

    infile.close()
    outfile.close()

    logging.info("Done")


@click.command(short_help="Runs the velocity analysis for a Rhapsody WTA sample")
@click.argument("samplefolder",
                type=click.Path(exists=True,
                                file_okay=False,
                                dir_okay=True,
                                readable=True,
                                writable=True,
                                resolve_path=True))
@click.argument("gtffile",
                type=click.Path(exists=True,
                                file_okay=True,
                                dir_okay=False,
                                readable=True,
                                resolve_path=True))
@click.option("--metadatatable", "-s",
              help="Table containing metadata of the various samples (csv fortmated rows are samples and cols are entries)",
              default=None,
              type=click.Path(resolve_path=True,
                              file_okay=True,
                              dir_okay=False,
                              readable=True))
@click.option("--mask", "-m",
              help=".gtf file containing intervals to mask",
              default=None,
              type=click.Path(resolve_path=True,
                              file_okay=True,
                              dir_okay=False,
                              readable=True))
@click.option("--logic", "-l",
              help="The logic to use for the filtering (default: Default)",
              default="Default")
@click.option("--multimap", "-M",
              help="""Consider not unique mappings (not reccomended)""",
              default=False,
              is_flag=True)
@click.option("--samtools-threads", "-@",
              help="The number of threads to use to sort the bam by cellID file using samtools",
              default=16)
@click.option("--samtools-memory",
              help="The number of MB used for every thread by samtools to sort the bam file",
              default=2048)
@click.option("--dtype", "-t",
              help="The dtype of the loom file layers - if more than 6000 molecules/reads per gene per cell are expected set uint32 to avoid truncation (default run_10x: uint16)",
              default="uint16")
@click.option("--dump", "-d",
              help="For debugging purposes only: it will dump a molecular mapping report to hdf5. --dump N, saves a cell every N cells. If p is prepended a more complete (but huge) pickle report is printed (default: 0)",
              default="0")
@click.option('--reuse-corrected-bam', '-R',
              help="Reuse corrected bamfile if exists.",
              default=False,
              is_flag=True)
@click.option('--verbose', '-v',
              help="Set the vebosity level: -v (only warinings) -vv (warinings and info) -vvv (warinings, info and debug)",
              count=True, default=1)
def run(samplefolder: str, gtffile: str,
           metadatatable: str, mask: str, logic: str, multimap: bool,
           samtools_threads: int, samtools_memory: int, dtype: str, dump: str,
           reuse_corrected_bam : bool, verbose: str) -> None:
    """Runs the velocity analysis for a Rhapsody WTA Sample

    SAMPLEFOLDER specifies the Rhapsody WTA pipline sample folder

    GTFFILE genome annotation file
    """

    # Check that the Rhapsody WTA analysis was run successfully
    metricmatches = glob.glob(os.path.join(samplefolder, "Logs", "*-mist_metrics.log"))
    if len(metricmatches) == 0:
        logging.error("The outputs of Rhapsody WTA are not ready")
        sys.exit(1)
    metricfilename = metricmatches[0]
    if "Completed Metrics" not in open(metricfilename).read():
        logging.error("The outputs of Rhapsody WTA are not ready")
        sys.exit(1)

    outputfolder = os.path.join(samplefolder, "Velocyto")
    if not os.path.exists(outputfolder):
        os.mkdir(outputfolder)

    # Check and correct bam file
    bammatches = glob.glob(os.path.join(samplefolder, "*_final.BAM"))
    if len(bammatches) == 0:
        logging.error("The bam file not exist!")
        sys.exit(1)
    bamfilename = bammatches[0]

    corrected = os.path.join(outputfolder, f"correct_{os.path.basename(bamfilename)}")
    if os.path.exists(corrected) and reuse_corrected_bam:
        logging.info("Using existing corrected bamfile %s" % corrected)
    else:
        rhapsody_wta_correct(bamfilename, corrected)
    bamfilename = corrected

    # Check valid cell names
    bcmatches = glob.glob(os.path.join(samplefolder, "*_Expression_Matrix.mtx.colnames"))
    if len(bcmatches) == 0:
        logging.error("Please convert expression csv file into MatrixMarket format first.")
        sys.exit(1)
    bcfile = bcmatches[0]

    sampleid = os.path.basename(samplefolder.rstrip("/").rstrip("\\"))
    assert not os.path.exists(os.path.join(outputfolder, f"{sampleid}.loom")), "The output already exist. Aborted!"
    additional_ca = {}

    return _run(bamfile=(bamfilename, ), gtffile=gtffile, bcfile=bcfile, outputfolder=outputfolder,
                sampleid=sampleid, metadatatable=metadatatable, repmask=mask, onefilepercell=False,
                logic=logic, without_umi=False, umi_extension="no", multimap=multimap, test=False, samtools_threads=samtools_threads,
                samtools_memory=samtools_memory, dump=dump, loom_numeric_dtype=dtype, verbose=verbose, additional_ca=additional_ca)

if __name__ == "__main__":
    run()
