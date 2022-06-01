# Author: Roll 820

import logging
import os
import time
import argparse
import numpy as np
from Bio import SeqIO
import pandas as pd

from utils import amend_des, load_seq
from domain import run_hmmscan
from chain import get_block_seq
from align import align_block, concatenate_align_block, refine_align_edge
from evaluate import run_fastSP


def run_DAMPSA(seqs, fp_DAMPSA_aln_out, fp_std_aln_out='', fp_dom_out='',
         fp_ref_aln='', fp_evl_out='', focus_clan=None, cache_dom=None,
         standard_app='clustalo', domain_app='clustalo', linker_app='clustalo',
         do_edge_refine=True, check_linker_len=False, n_thread=7):
    """Main function of DAMPSA.

    Args:
        seqs (dict): dictionary storing Bio.SeqRecords.
        fp_std_aln_out (str): path to store alignment output from standard algorithm.
        fp_DAMPSA_aln_out (str): path to store DAMPSA alignment output, in fasta format.
        fp_dom_out (str): path to store domain prediction output.
        fp_ref_aln (str): path to reference alignment file; if empty, no evaluation.
        fp_evl_out (str): path to store evaluation (FastSP) output; if empty, no evaluation.
        focus_clan (list): list of Pfam Clan (domain superfamily) ID. If supplied, will only consider these clans.
        cache_dom (str): file path to cached domain annotation table (TSV format).
        standard_app (str): ``clustalo`` or ``mafft`` for aligning entire sequence as a standard reference.
        domain_app (str): ``clustalo`` or ``mafft`` for aligning domain blocks.
        linker_app (str): ``clustalo`` or `mafft` for aligning linker blocks.
        do_edge_refine (bool): whether to refine alignment block edge.
        check_linker_len (bool): whether check if the linker is too long (likely domain mis-detection).
        n_thread (int): thread number for computation, default 7.
    """
    logger = logging.getLogger('DAMPSA.main')

    # Hmmscan (slow)
    if cache_dom is None:
        domain_table = run_hmmscan(seqs, n_thread=n_thread)
    else:
        domain_table = pd.read_table(cache_dom, sep='\t')

    # filter clans
    if focus_clan is not None:
        domain_table = domain_table[domain_table.clan.isin(focus_clan)]

    # extract sequences
    chain, domain, linker, raw_seq, faulty_seq, domain_table = get_block_seq(seqs, domain_table, get_mask=False,
                                                                             check_linker_len=check_linker_len)
    logger.info('Domain chain: {}.'.format(' )=( '.join(chain)))
    if len(list(np.unique(domain_table['seq_id']))) == 1:
        raise NotImplementedError('Invalid input: only one sequence has detected domain.')

    if fp_dom_out:
        domain_table.to_csv(fp_dom_out, sep='\t', index=False)

    # align with standard method
    if fp_std_aln_out:
        logger.info('Aligning with the standard clustalo...')
        raw_seq_sorted = {}
        keys = sorted(raw_seq.keys())
        for k in keys:
            raw_seq_sorted[k] = raw_seq[k]
        raw_align = align_block(raw_seq_sorted, mode='single', app=standard_app)
        _ = SeqIO.write(raw_align, fp_std_aln_out, 'fasta')

    # align domain and linker
    logger.info('Aligning blocks...')
    domain_align = {}
    for k in domain.keys():
        if domain[k]:
            domain_align[k] = align_block(domain[k], mode='single', app=domain_app)
    linker_align = {}
    for k in linker.keys():
        if linker[k]:
            linker_align[k] = align_block(linker[k], mode='single', app=linker_app)
    all_id = list(np.setdiff1d(np.array(list(raw_seq.keys())),
                               np.array(list(faulty_seq.keys()))))
    final_align = concatenate_align_block(chain, domain_align, linker_align, all_id)

    # refine alignment
    if do_edge_refine:
        logger.info('Refining aligned block edges...')
        try:
            final_align = refine_align_edge(final_align,
                                            domain_table[~ domain_table['seq_id'].isin(list(faulty_seq.keys()))])
        except:
            logger.info('Failed to refine aligned block edges, '
                        'maybe no protein contains all domains detected in the group.')

    # align faulty sequences separately
    if faulty_seq:
        final_align = align_block(faulty_seq, final_align, mode='profile', app='clustalo')
        final_align.sort()

    final_align = amend_des(raw_seqs=seqs, alignment=final_align)
    _ = SeqIO.write(final_align, fp_DAMPSA_aln_out, "fasta")

    # evaluate the output (DAMPSA vs. clustalo)
    if fp_ref_aln != '' and fp_std_aln_out != '' and fp_evl_out != '':
        logger.info("Evaluating alignment...")
        evl = run_fastSP(ref_file=fp_ref_aln, test_files=[fp_std_aln_out, fp_DAMPSA_aln_out], test_file_dir=None)
        evl.to_csv(fp_evl_out, sep='\t', index=False)
    return


def get_parser():
    """Setup commandline parser.
    """
    parser = argparse.ArgumentParser(description="DAMPSA input arguments.")
    parser.add_argument(
        "--input",
        type=str,
        default=None,
        help="Path to the input .fasta file.",
    )
    parser.add_argument(
        "--output",
        type=str,
        default=None,
        help="Path to the alignment .fasta file output.",
    )
    parser.add_argument(
        "--domain-out",
        type=str,
        default=None,
        help="Path to output domain annotation results.",
    )
    parser.add_argument(
        "--refine-edge",
        action='store_true',
        help="Refine alignments at the edge between domain and linker segments.",
    )
    parser.add_argument(
        "--no-check-linker",
        action='store_true',
        help="Not to check if the linker is too long - likely contains mis-detected domains.",
    )
    parser.add_argument(
        "--focus-clan",
        type=str,
        default=None,
        help="Only consider the specified Clan IDs (domain superfamily) - comma separated.",
    )
    parser.add_argument(
        "--cache-dom",
        type=str,
        default=None,
        help="Skip hmmscan, use supplied filepath to cached domain table (TSV-like).",
    )
    parser.add_argument(
        "--domain-app",
        type=str,
        default='clustalo',
        help="Aligner for domain segments (clustalo or mafft), default clustalo.",
    )
    parser.add_argument(
        "--linker-app",
        type=str,
        default='clustalo',
        help="Aligner for linker segments (clustalo or mafft), default clustalo.",
    )
    parser.add_argument(
        "--log",
        action='store_true',
        help="Store log file in the same folder as the alignment",
    )
    parser.add_argument(
        "--n-thread",
        type=int,
        default=7,
        help="Thread number for running hmmscan (domain annotation), default 7.",
    )
    return parser


def setup_logger(file_out=''):
    """Setup logger.

    Args:
        file_out (str): File path to write (append) log; if empty, will not generate log file.
    """
    logger = logging.getLogger('DAMPSA')
    logger.setLevel(logging.INFO)
    formator = logging.Formatter(fmt="%(asctime)s %(filename)s %(levelname)s %(message)s",
                                 datefmt="%Y/%m/%d %X")

    str_hdl = logging.StreamHandler()
    str_hdl.setFormatter(formator)
    logger.addHandler(str_hdl)

    if file_out:
        file_hdl = logging.FileHandler(file_out, encoding="UTF-8")
        file_hdl.setFormatter(formator)
        logger.addHandler(file_hdl)
    return


if __name__ == '__main__':
    # setup
    start_time = time.time()
    args = get_parser().parse_args()
    if args.log:
        setup_logger(file_out=os.path.join(os.path.dirname(args.input), 'log.txt'))
    else:
        setup_logger(file_out='')
    logger = logging.getLogger('DAMPSA')
    logger.info(vars(args))
    logger.info('Start processing {}.'.format(os.path.basename(args.input)))

    seqs = load_seq(args.input)
    if args.domain_out is None:
        args.domain_out = ''
    if args.focus_clan is not None:
        args.focus_clan = list(map(lambda x:x.strip(), args.focus_clan.split(',')))

    # run application
    run_DAMPSA(seqs, fp_DAMPSA_aln_out=args.output, fp_dom_out=args.domain_out,
               domain_app=args.domain_app, linker_app=args.linker_app,
               do_edge_refine=args.refine_edge, check_linker_len=not args.no_check_linker,
               n_thread=args.n_thread, focus_clan=args.focus_clan, cache_dom=args.cache_dom,
               # below for evaluation only, not implemented in standard run
               fp_ref_aln='', fp_evl_out='', standard_app='', fp_std_aln_out='')
    logger.info('Finished processing {} in {:.2f}s.'.format(os.path.basename(args.input),
                                                            time.time() - start_time))
    logger.info('='*50)
