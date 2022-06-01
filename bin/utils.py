# Author: Roll 820

import pandas as pd
import numpy as np
import re
import os
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment


def get_overlap_pair(segms):
    """Get every pair of overlapping segment in a list.

    Args:
        segms (list): a list of segments :code:`[start, end]`.

    Returns:
        A list of paired segments that overlap.
    """

    ovl = []
    for i in range(len(segms)):
        for j in range(i+1, len(segms)):
            a = segms[i]
            b = segms[j]
            if not (a[1] < b[0] or b[1] < a[0]):
                ovl.append((a,b))
    return ovl


def get_within_pair(segms):
    """Get every segment that is contained by another segment.

    Args:
        segms (list): a list of segments :code:`[start, end]`.

    Returns:
        A list of paired segments that have a "within" relationship.
    """

    ovl = []
    for i in range(len(segms)):
        for j in range(i+1, len(segms)):
            a = segms[i]
            b = segms[j]
            if (a[1] < b[1] and a[0] > b[0]):
                ovl.append(a)
            elif (b[1] < a[1] and b[0] > a[0]):
                ovl.append(b)
    return ovl


def load_seq(fp, allow_suffix=None):
    """Read all sequences together (may take much memory)
    
    Args:
        fp (str): filepath to sequences or the folder containing them.
        allow_suffix (list): allow other suffix than :file:`.fasta`.
    Returns:
        A dictionary with keys: sequence ID and values: :class:`Bio.SeqRecord` object.
    """
    if allow_suffix is None:
        allow_suffix = []

    if os.path.isdir(fp):
        ipt_ls = os.listdir(fp)
        record = {}
        for i in ipt_ls:
            for sfx in allow_suffix + ['fasta']:
                if re.search(sfx, i):
                    rd = SeqIO.parse(os.path.join(fp, i), 'fasta')
                    rd = SeqIO.to_dict(rd)
                    record.update(rd)
                    break
    else:
        record = SeqIO.parse(fp, 'fasta')
        record = SeqIO.to_dict(record)
    
    return record


def load_Uniprot_series(series_fp, tab_fp=None, collapse_org=True, max_PE=3, output=None, **kwargs):
    """Load a series of protein sequences from the Uniprot with pre-processing.

    Args:
        series_fp (str): filepath to protein series.
            - either :file:`.fasta` file with multiple proteins, e.g. BLAST output, or one directory with multiple
             :file:`.fasta` files).
        collapse_org (bool): collapse multiple sequences belonging to the same organism or not.
            - If `True`, then collapse based on the following criteria:

                - If :code:`tab_fp` supplied, pick "reviewed" sequence of the organism;
                otherwise, pick the sequence with the lowest PE.

                - If still exists multiple sequences, pick the longest one.

        tab_fp (str): use "reviewed" information in Uniprot tab summary if supplied.
        max_PE (int): maximum protein evidence level, default 3.
        **kwargs: other parameters for loading protein series.
    
    Returns:
        UniProt protein sequences.
    """
    logger = logging.getLogger('DAMPSA.loadUniProt')
    kwargs['fp'] = series_fp
    seqs = load_seq(**kwargs)
    to_rm = []
    organism_id = {}
    for key in seqs.keys():
        seq = seqs[key]
        mth = re.search('OS=(.+) OX=.+ PE=(.+) SV',seq.description)
        if int(mth.group(2))>max_PE:
            to_rm.append(key)
        else:
            if mth.group(1) not in organism_id.keys():
                organism_id[mth.group(1)] = [key]
            else:
                organism_id[mth.group(1)].append(key)
    logger.info('Filtering {}/{} out by evidence level PE={}.'.format(len(to_rm), len(seqs.keys()), max_PE))
    
    for i in to_rm:
        del seqs[i]

    if tab_fp is not None:
        tab = pd.read_csv(tab_fp, sep='\t')
    else:
        tab = None

    if collapse_org:
        to_rm = []
        for key in organism_id.keys():
            orgnm = organism_id[key]
            if len(orgnm) > 1:
                orgnm_short = list(map(lambda x:x.split('|')[1], orgnm))
                flag = False if tab is None else True
                if tab is not None:
                    sub_tab = tab.loc[orgnm_short]
                    rvd = sub_tab[sub_tab['Status']=='reviewed'].Entry
                    if rvd.shape[0]!=0:
                        # pick reviewed
                        sel = rvd.index[0]
                    else:
                        flag = False
                if not flag:
                    # pick lowest PE
                    seq_PE = list(map(lambda x:
                        int(re.search('PE=(.+) SV',seqs[x].description).group(1)), orgnm))
                    min_idx = np.where(seq_PE == np.min(seq_PE))[0]
                    sel = list(np.array(orgnm)[min_idx])
                    if len(sel)>1:
                        seq_len = list(map(lambda x:len(seqs[x]), sel))
                        # pick longest, the first one
                        sel = sel[np.where(seq_len == np.max(seq_len))[0][0]]
                    else:
                        sel = sel[0]
                for o in orgnm:
                    if o != sel:
                        to_rm.append(o)
        logging.info('Filtering {}/{} out by collapsing the same organism.'.format(len(to_rm), len(seqs.keys())))
    
        for i in to_rm:
            del seqs[i]
    
    logging.info('Remaining {} sequences.'.format(len(seqs)))

    # change sequence ID to UniProt Entry
    new_seqs = {}
    for key in seqs.keys():
        sq = seqs[key]
        sq.id = key.split('|')[1]
        new_seqs[key.split('|')[1]] = sq

    if output:
        SeqIO.write(new_seqs.values(), output, 'fasta')
    
    return new_seqs


def amend_des(raw_seqs, alignment):
    """Flow description in raw sequences into alignment.
    
    Args:
        raw_seqs (dict): sequence dictionary.
        alignment (Bio.Align.MultipleSeqAlignment): alignment amend with metadata.
    """
    des = {}
    aln_seqs = []
    for k in raw_seqs.keys():
        des[k] = raw_seqs[k].description
    for k in range(len(raw_seqs)):
        sq = alignment[k,:]
        sq.description = des[sq.id]
        # sq.id = sq.id.replace('_', '')
        aln_seqs.append(sq)
    return MultipleSeqAlignment(aln_seqs)


def convert_aln(alignment):
    """Replace :code:`.` in :file:`.fasta` alignment with :code:`-`. Also uppercase each char.
    """
    new_ref = []
    for rf in alignment:
        rf.seq = Seq(rf.seq.replace('.','-').upper())
        new_ref.append(rf)
    return MultipleSeqAlignment(new_ref)
