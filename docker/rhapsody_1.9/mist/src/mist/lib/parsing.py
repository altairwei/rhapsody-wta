# from typing import List, Union  # disabled until update to python 2.7.11 or python 3.x
import argparse
import logging
import pprint


seq_stats = '--seq-stats'
annot_mol_file = '--annot-mol-file'


def add_parsing(parser, required_arguments=None, optional_arguments=None):
    """MetadataParser helps avoid duplicated code for parsing metadata"""
    if required_arguments is None:
        required_arguments = []
    if optional_arguments is None:
        optional_arguments = []
    arguments = required_arguments + optional_arguments

    if seq_stats in arguments:
        # TODO: instead of passing along the sequencing metrics file, which is a file
        #       that contains the metadata incidentally, pass along the metadata itself
        parser.add_argument(seq_stats,
                            required=seq_stats in required_arguments,
                            help='SeqMetrics generated in AnnotateReads')
    if annot_mol_file in arguments:
        parser.add_argument(annot_mol_file,
                            required=annot_mol_file in required_arguments,
                            help='Metrics-files archive generated by GDT.')
    return parser


def log_node_arguments(args):
    # type: (Union[argparse.Namespace, dict]) -> None
    cli_logger = logging.getLogger('cli')
    try:
        args_dict = args.__dict__
    except AttributeError:
        args_dict = args
    cli_logger.debug('Running with options:\n{}'.format(pprint.pformat(args_dict)))
