# Author: Roll 820

from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Application import ApplicationError
from Bio.Align import MultipleSeqAlignment, AlignInfo
from Bio import AlignIO, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import uuid
import copy
import subprocess
import numpy as np
import re
import pandas as pd

from external.alignment.alignment.align import multi_sequence_alignment


def align_block(seqs, profile=None, mode='single', app='clustalo'):
    """Align a block of sequences with selected aligner.

    Args:
        seqs (dict or list): dictionary or list of :class:`Bio.SeqRecord`.
        profile (dict): :class:`Bio.SegRecord` used to generate a profile, use in ``profile`` mode.
        mode (str): either ``single`` or ``profile``.
        app (str): either ``clustalo`` or ``mafft``. Only applied to ``single`` mode.

    Returns:
        Aligned sequence block.

    Note:
        - In the ``profile`` mode, you can only use `clustalo`.
        - In the ``single`` mode, you can choose to use `clustalo` or `mafft l-ens-i`.
    """
    align = None
    if mode == 'profile' and profile is None:
        raise ValueError('With profile mode, must supply sequences to generate the profile.')
    
    if isinstance(seqs, dict):
        seqs = list(seqs.values())

    if len(seqs) == 1 and mode == 'single':
        # only one sequence, no need to align
        return MultipleSeqAlignment(seqs)

    filename = uuid.uuid4().hex
    in_file = '.' + filename + '.fasta'
    SeqIO.write(seqs, in_file, 'fasta')
    
    if mode == 'profile':
        if len(seqs) == 1:
            flag = True
        else:
            flag = False
        profile_file = '.'+filename+'_pfl.fasta'
        SeqIO.write(profile, profile_file, 'fasta')
        
    out_file = '.' + filename + '_aln.fasta'
    skip = False
    try:
        if mode == 'single':
            if app == 'mafft':
                cmd = 'linsi --quiet {} > {}'.format(in_file, out_file) # MAFFT aligner
                _ = subprocess.run(cmd, shell=True)
            else:
                clustalomega_cline = ClustalOmegaCommandline(infile=in_file,
                                                   outfile=out_file, verbose=False, auto=True)
                clustalomega_cline()
        elif mode == 'profile':
            if not flag:
                cmd = 'clustalo -i {} --p1 {} --is-profile -o {}'.format(in_file, 
                                                        profile_file, out_file)
            else:
                cmd = 'clustalo --p1 {} --p2 {} --is-profile -o {}'.format(in_file, 
                                                        profile_file, out_file)
            _ = subprocess.run(cmd, shell=True)
    except ApplicationError as e:
        if re.search('forcing Viterbi', str(e)):
            to_align = list(map(lambda x:list(str(x.seq)), seqs))
            aligned = multi_sequence_alignment(to_align, gap_penalty=1)
            aligned = np.array(pd.DataFrame(aligned.alignment))
            aln_dic = {}
            for i in range(aligned.shape[0]):
                aln_dic[''.join(aligned[i,:]).replace('_','')] = ''.join(aligned[i,:]).replace('_','-')
            for i in range(len(to_align)):
                seqs[i].seq = Seq(aln_dic[str(seqs[i].seq)])
            align = MultipleSeqAlignment(seqs)
            skip = True
    finally:
        _ = subprocess.run('rm ' + in_file, shell=True)
        if mode == 'profile':
            _ = subprocess.run('rm ' + profile_file, shell=True)

    if not skip:
        align = AlignIO.read(out_file, "fasta")
        _ = subprocess.run('rm ' + out_file, shell=True)

    return align


def concatenate_align_block(chain, domain_align, linker_align, all_id, escape='null'):
    """Connect alignment blocks

    Args:
        chain (list): order of alignment blocks.
        domain_align (dict): alignment of each domain block.
        linker_align (dict): alignment of each linker (and N/C overhang) block.
        all_id (list): list of all sequence ids.
        escape (str): empty alignment block that cannot be chained, default :code:'null'.

    Returns:
        Concatenated alignment block.
    """

    def _check_aln(align_set, all_id, escape='null'):
        """Check missing sequences in the alignment block

        Args:
            align_set (dict): a group of alignment blocks.
            all_id (list): all valid sequences involved in this round of alignment.
            escape (str): the name of alignment blocks no need to check.

        Returns:
            A group of alignment blocks with missing sequences padded.
        """
        align_sets = copy.deepcopy(align_set)
        for k in align_sets.keys():
            aln_id = []
            aln = align_sets[k]
            if k != escape:
                for i in aln:
                    aln_id.append(i.id)
                to_add = list(np.setdiff1d(all_id, aln_id))
                if len(to_add) > 0:
                    len_aln = aln.get_alignment_length()
                    for i in to_add:
                        aln.extend([SeqRecord('-' * len_aln, id=i, name='_'.join([k, i]))])
                aln.sort()
            align_sets[k] = aln
        return align_sets

    domain_align = _check_aln(domain_align, all_id)
    linker_align = _check_aln(linker_align, all_id)
    final_align = []
    for i in chain:
        if i == 'null':
            continue
        if i[0] == '_':
            if i in linker_align.keys():
                final_align.append(linker_align[i])
        else:
            final_align.append(domain_align[i])
    
    if len(final_align) == 1:
        return final_align[0]
    else:
        for i in range(1,len(final_align)):
            final_align[0] += copy.deepcopy(final_align[i])

    return final_align[0]


def refine_align_edge(alignment, domain_table, search_range=10):
    """Refine the edge between linker and domains

    Args:
        alignment (Bio.Align.MultipleSeqAlignment): alignment matrix to refine.
        domain_table (pandas.DataFrame): table storing domain information.
        search_range (int): only backward search domain end for 10 residues

    Returns:
         Edge-refined alignments.
    """
    # alignment = AlignIO.read('../data/src/DAMPSA_aln.fasta', 'fasta')
    # domain_table = pd.read_table('../data/src/src_filtered_dom.txt', sep='\t)
    
    mtx = np.array(alignment)
    # get the consensus sequence
    summary_align = AlignInfo.SummaryInfo(alignment)
    consensus = str(summary_align.dumb_consensus())
    
    # identify sequences that have all domain/linkers
    ndom = len(np.unique(domain_table.clan))
    mid = []
    for i in np.unique(domain_table.seq_id):
        if domain_table[domain_table.seq_id == i].shape[0] == ndom:
            mid.append(i)
    if len(mid) == 0:
        raise ValueError('Cannot find a master sequence.')
    dt = domain_table[domain_table.seq_id.isin(mid)]
    ids = {}
    count = 0
    for aln in alignment:
        ids[aln.id] = count
        count += 1

    # identify the coordinate in alignment matrix
    dom_end = {}
    for dm in np.unique(dt.clan):
        dt_sub = dt[dt.clan == dm]
        j = 0
        master_seq = alignment[ids[dt_sub.seq_id.iloc[0]]]
        ed = dt_sub.env_end.iloc[0]
        all_ids = []
        for i in dt_sub['seq_id']:
            all_ids.append(ids[i])
        while True:
            for i in range(len(master_seq)):
                if master_seq[i] != '-':
                    j += 1
                if j == ed:
                    dom_end[dm] = i
                    break
            # calculate coordinate in all seqs
            len_cur = np.apply_along_axis(lambda 
                                          x:len(''.join(x).replace('-','')),
                                          1, mtx[all_ids,:i])
            len_cur = np.transpose(np.array([dt_sub.env_end, len_cur+1]))
            beyond_seq = list(np.where(len_cur[:,0] > len_cur[:,1])[0])
            if beyond_seq:
                master_seq = alignment[ids[dt_sub.seq_id.iloc[beyond_seq[0]]]]
                ed = dt_sub.env_end.iloc[beyond_seq[0]]
                j = 0
            else:
                break
    
    # for each domain end site, search if anything in the linker can fill gaps.
    for dm in dom_end.keys():
        ed = dom_end[dm]
        for i in range(mtx.shape[0]):
            if mtx[i,ed] == '-':
                # if there is any gap, try to fill it
                inv_seq = ''.join(mtx[i,ed-search_range+1:ed+1])[::-1]
                t = 0
                while inv_seq[t] == '-':
                    t += 1
                    if t == search_range:
                        break
                # get first non gap residue
                p = 1
                flg = False
                while mtx[i, ed + p] == '-':
                    p += 1
                    if ed + p == len(consensus):
                        flg = True
                        break
                if flg:
                    continue
                # modify if match consensus
                o = q = 1
                flg = False
                while q <= t:
                    if consensus[ed - t + q] == mtx[i, ed + p + o - 1]:
                        flg = True
                        mtx[i, ed - t + q] = mtx[i, ed + p + o - 1]
                        mtx[i, ed + p + o - 1] = '-'
                        o += 1
                    elif flg:
                        break
                    q += 1
                    if ed + p + o - 1 == len(consensus):
                        break
                    
    new_aln = []
    for i in range(mtx.shape[0]):
        new_aln.append(SeqRecord(Seq(''.join(mtx[i,:])), id=alignment[i].id,
                                 name=alignment[i].name, 
                                 description=alignment[i].description))
    new_aln = MultipleSeqAlignment(new_aln)
    
    # AlignIO.write(new_aln, 'refine.fasta', 'fasta')
    
    return new_aln
