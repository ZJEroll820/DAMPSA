# Author: Roll 820

import logging
import numpy as np
import copy
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd
import re

from align import align_block
from external.alignment.alignment.align import multi_sequence_alignment


def get_block_seq(seqs, domain_table, dcol='clan',
                  start_nm='env_start', end_nm='env_end',
                  get_mask=False, check_linker_len=False):
    """Get sequence for specific domain.

    Args:
        seqs (dict): dictionary of :class:`Bio.SeqRecord`.
        domain_table (pandas.DataFrame): must have the following columns:
            :code:`seq_id`: sequence ID;
            :code:`start_nm`: start of domain;
            :code:`end_nm`: end of domain;
            :code:`dcol`: column refer to the domain.
        dcol (str): config the name of the column referring to the domain ID.
        start_nm (str): config the name of the column referring to the domain start.
        end_nm (str): config the name of the column referring to the domain end.
        get_mask (bool): whether mask domain with * and output sequences.
        check_linker_len (bool): whether to check overlong linker.

    Returns:
        Sequences corresponding to each domain block (dict).
    """
    logger = logging.getLogger('DAMPSA.domain')

    domain_table = copy.deepcopy(domain_table)
    seqs = copy.deepcopy(seqs)
    miss_count = 0
    all_table = list(domain_table['seq_id'])
    coll = {}
    coll_link = {}

    chain, domain_table, gap_dom = get_chain_by_align(domain_table, dcol=dcol, start_nm=start_nm, gap_penalty_init=0.4)
    if len(gap_dom) > 0:
        logger.warning('{} sequences seem to have gapped domains. Align as a whole later.'.format(len(gap_dom)))
        logger.debug(gap_dom)
    # chain = get_chain(domain_table, dcol=dcol, start_nm=start_nm, end_nm=end_nm) # TODO
    for i in np.unique(domain_table[dcol]):
        coll[i] = {}
    for i in np.setdiff1d(chain, np.unique(domain_table[dcol])):
        coll_link[i] = {}

    # Sort by chain domains
    dom = chain[1::2]
    domain_table[dcol] = domain_table[dcol].astype('category')
    domain_table[dcol] = domain_table[dcol].cat.set_categories(dom)
    domain_table = domain_table.sort_values(['seq_id', dcol])

    # Extract domain and linker sequences
    # !!! Assume no gapped domain chain (e.g. A-B-C and A-C)
    # a lookup for linker between domains
    lnk_up = {}
    for i in range(2,len(chain)-1,2):
        lnk_up['='.join([chain[i-1], chain[i+1]])] = chain[i]
    to_rm = []
    ovl_rm = []
    for i in seqs.keys():
        if i not in all_table:
            miss_count += 1
            to_rm.append(i)
        elif i in gap_dom:
            continue
        else:
            sub = domain_table[domain_table['seq_id'] == i].copy()
            sub.index = range(0, sub.shape[0])
            for j in range(sub.shape[0]):
                st = sub.loc[j, start_nm] - 1
                ed = sub.loc[j, end_nm]
                if j == 0:
                    sub_seq = seqs[i].seq[:st]
                    if len(sub_seq) > 0:
                        coll_link['_N'][i] = SeqRecord(sub_seq, id=i, name='_'.join(['_N', i]))
                if j == sub.shape[0] - 1:
                    sub_seq = seqs[i].seq[ed:]
                    if len(sub_seq) > 0:
                        coll_link['_C'][i] = SeqRecord(sub_seq, id=i, name='_'.join(['_C', i]))
                else:
                    next_st = sub.loc[j + 1, start_nm] - 1
                    if next_st < ed:
                        logger.warning('Overlapping domain detected, align {} as a whole later.'.format(i))
                        ovl_rm.append(i)
                        break
                    link_nm = lnk_up['='.join([sub.loc[j, dcol], sub.loc[j+1, dcol]])]
                    if ed != next_st:
                        sub_seq = seqs[i].seq[ed:next_st]
                        coll_link[link_nm][i] = SeqRecord(sub_seq, id=i, name='_'.join([link_nm, i]))

                sub_seq = seqs[i].seq[st:ed]

                if get_mask:
                    # mask domain with gap,
                    # change in place,
                    # assume that domain does not overlap!!!
                    sub_seq_ls = list(seqs[i].seq)
                    sub_seq_ls[st:ed] = ['*' for _ in range(ed - st)]
                    sub_seq_ls = ''.join(sub_seq_ls)
                    seqs[i].seq = Seq(sub_seq_ls)

                coll[sub.loc[j, dcol]][i] = SeqRecord(sub_seq, id=i, name='_'.join([sub.loc[j, dcol], i]))

    if len(to_rm) > 0:
        logger.warning('{} sequences missing in the domain annotation table. Align as a whole later.'.format(miss_count))
        logger.warning(to_rm)

    # remove information of domain-overlapping seq
    if ovl_rm:
        for i in coll.keys():
            for j in ovl_rm:
                if j in coll[i].keys():
                    del coll[i][j]
        for i in coll_link.keys():
            for j in ovl_rm:
                if j in coll_link[i].keys():
                    del coll_link[i][j]

    faulty_seq = []
    if check_linker_len:
        # check if there is any overlong linker
        faulty_seq, coll, coll_link, chain = refine_linker(coll, coll_link, chain, domain_table)
    if faulty_seq:
        logger.warning('{} sequences have unreliable domain annotation. Align as a whole later.'.format(len(faulty_seq)))
        logger.warning(faulty_seq)

    fss = {}
    for fs in faulty_seq + to_rm + gap_dom + ovl_rm:
        fss[fs] = seqs[fs]
    if len(seqs) < 2:
        raise ValueError('No or only one sequence remained after quality control.')

    return chain, coll, coll_link, seqs, fss, domain_table


def add_linker_to_chain(chain):
    """Add linker annotation like :code:`_link_1` (in order from one) to the domain chain.
    """
    chain_link = []
    ct = 1
    if len(chain) > 1:
        for i in range(len(chain)):
            chain_link.append(chain[i])
            chain_link.append('_link_' + str(ct))
            ct += 1
        chain_link = ['_N'] + chain_link[:-1] + ['_C']  # two overhang
    else:
        chain_link = ['_N'] + chain + ['_C']
    return chain_link


def get_chain_by_align(domain_table, dcol='clan', start_nm='env_start', gap_penalty_init=0.4):
    """Get Clan chain from domain annotation using progressive alignment.

    Args:
        gap_penalty_init (float): Initial gap penalty (0 to 1) to try progressive alignment on the chain.
            If the alignment cannot separate each domain or linker in a unique column, will reduce gap penalty by half.
            If failed for five rounds, will raise an error.

    See :code:`get_block_seq` for other parameters.
    """
    if gap_penalty_init <= 0 or gap_penalty_init > 1:
        raise ValueError('Parameter gap_penalty_init should be between 0 and 1.')

    all_dom = list(np.unique(domain_table[dcol]))
    if len(all_dom) == 1:
        chain_link = add_linker_to_chain(all_dom)
        return chain_link, domain_table, []

    to_align = []
    lk = [] # lookup for unique pattern
    for i in np.unique(domain_table['seq_id']):
        sub = domain_table[domain_table.seq_id == i].sort_values(start_nm)
        dlist = list(sub[dcol])
        if '='.join(dlist) not in lk:
            lk.append('='.join(dlist))
            to_align.append(list(sub[dcol]))

    if len(to_align) == 1:
        chain_link = add_linker_to_chain(to_align[0])
        return chain_link, domain_table, []

    # check if one column only contains one domain, if yes, divide gap penalty by half
    MAX_ITER = 5
    itr = 1
    flg = False
    while itr <= MAX_ITER:
        aligned = multi_sequence_alignment(to_align, gap_penalty=gap_penalty_init)
        aligned = np.array(pd.DataFrame(aligned.alignment))
        flg2 = False
        for i in range(aligned.shape[1]):
            blk = list(np.setdiff1d(np.unique(aligned[:, i]), '_'))
            if len(blk) > 1:
                gap_penalty_init /= 2
                flg2 = True
                break
        if not flg2:
            flg = True
            break
        itr += 1
    if not flg:
        raise NotImplementedError('Two different domains in one column, failed to find a good gap penalty')

    chain = []
    dup_check = []
    for i in range(aligned.shape[1]):
        blk = list(np.setdiff1d(np.unique(aligned[:,i]), '_'))
        if blk[0] in chain and blk[0] not in dup_check:
            # duplicated domains, check later
            dup_check.append(blk[0])
        chain.append(blk[0])

    dup = {}
    aln = copy.deepcopy(aligned)
    if dup_check:
        # if has duplicated domains, label as 1,2,3...; identify all pattern that contains it
        for d in dup_check:
            count = 1
            for i in range(aligned.shape[1]):
                nm = d + '_' + str(count)
                flg = False
                for j in range(aligned.shape[0]):
                    if aligned[j,i] == d:
                        flg = True
                        aligned[j,i] = nm
                if flg:
                    chain[chain.index(d)] = nm
                    count += 1

        # duplicated domain pattern substitution
        for i in range(aligned.shape[0]):
            a = '='.join(aln[i, aln[i,:] != '_'])
            b = '='.join(aligned[i, aligned[i,:] != '_'])
            if a != b:
                dup[a] = b
    
    # generate a lookup table, exclude "gapped" pattern (domain A-C rather than A-B-C)
    faulty_pattern = []
    for i in range(aligned.shape[0]):
        a = '='.join(aln[i, aln[i,:] != '_'])
        full_chain = ''.join(aln[i,:]).replace('_', ' ').lstrip().rstrip()
        if full_chain.replace(' ', '') != full_chain:
            faulty_pattern.append(a)

    # update domain table
    faulty_seq = []
    dom_table = pd.DataFrame()
    for i in np.unique(domain_table['seq_id']):
        sub = domain_table[domain_table.seq_id == i].sort_values(start_nm)
        intersect = list(set(sub[dcol]) & set(dup_check))
        pt = '='.join(list(sub[dcol]))
        if pt in faulty_pattern:
            faulty_seq.append(i)
        if intersect:
            sub[dcol] = dup[pt].split('=')
        dom_table = pd.concat([dom_table, sub], ignore_index=True)
    dom_table = pd.DataFrame(dom_table)

    # check if chain still valid after removing gapped seq
    if faulty_seq:
        sub_dom_table = dom_table[~ dom_table['seq_id'].isin(faulty_seq)]
        extra = list(np.setdiff1d(np.array(chain), np.array(sub_dom_table[dcol])))
        new_chain = []
        for ch in chain:
            if ch not in extra:
                new_chain.append(ch)
        chain = new_chain
    # also check domain table
    dom_table = dom_table[dom_table['clan'].isin(chain)]

    # add chain linker
    chain_link = add_linker_to_chain(chain)
    return chain_link, dom_table, faulty_seq


def get_chain(domain_table, dcol='clan', start_nm='env_start', end_nm='env_end'):
    """(deprecated) Get Clan chain from domain annotation.

    Note:
        Domain sorted according to mean center across all sequences.
        Use IQR to filter outlier (if sequences containing domain <=3, don't do).
    """
    df = []
    for i in np.unique(domain_table[dcol]):
        sub = domain_table[domain_table[dcol] == i]
        center = np.mean(sub[[start_nm, end_nm]], axis=1)
        if len(center) > 2:
            Q1, Q3 = np.percentile(center, 25), np.percentile(center, 75)
            IQR = Q3 - Q1
            upper, lower = center >= (Q3 + 1.5 * IQR), center <= (Q1 - 1.5 * IQR)
            mn = np.mean(center[np.logical_not(np.logical_or(upper, lower))])
        else:
            mn = np.mean(center)
        df.append({'group': i, 'mean_pos': mn})
    df = pd.DataFrame(df)
    df = df.sort_values('mean_pos')
    chain = list(df['group'])
    chain_link = []
    ct = 1
    if len(chain) > 1:
        for i in range(len(chain)):
            chain_link.append(chain[i])
            chain_link.append('_link_' + str(ct))
            ct += 1
        chain_link = ['_N'] + chain_link[:-1] + ['_C']  # two overhang
    else:
        chain_link = ['_N'] + chain + ['_C']
    return chain_link


def check_linker_in_domain(dom, test_seqs, direction='start'):
    """Check if the linker has potential sequence in the domain.
    
    Args:
        dom (list): domain set(s) to test.
        seqs (dict): set of linker segment sequences to test.
        direction (str): either start or end of the test sequence may overlap
            with the domain sequences.
    """
    MIN_ALN = 10  # minimum aligned count define linker-domain overlap
    if direction not in ['start', 'end']:
        raise ValueError('Direction must be either start or end.')
    
    profile = align_block(dom, mode='single')
    idx = profile[:,0].index(re.search('(\w)', profile[:,0]).group(1))
    master_seq = profile[idx,:]
    
    # get the profile block in the profile2seq alignment
    aln = align_block(test_seqs, profile, mode='profile')
    if direction == 'end':
        aln = aln[::-1]
    for i in range(len(aln)):
        if aln[i,:].id == master_seq.id:
            idx = i
    aln_start = str(aln[idx,:].seq).index(str(master_seq.seq)[0])
    
    test_key = list(test_seqs.keys())
    test_idx = {}
    for i in range(len(aln)):
        for j in range(len(test_key)):
            if aln[i,:].id == test_key[j]:
                test_idx[test_key[j]] = i
                break
    
    # check if the linker really overlaps with the domain
    out = []
    for k in test_key:
        test_in_dom = str(aln[test_idx[k],aln_start:aln_start+profile.
                              get_alignment_length()].seq)
        if len(test_in_dom) - test_in_dom.count('-') >= MIN_ALN:
            b4_start = aln[test_idx[k],:aln_start]
            idx = len(b4_start) - str(b4_start.seq).count('-')
            if direction == 'end':
                idx = len(str(aln[test_idx[k],:].seq).replace('-','')) - idx
            out.append({'id':k, 'idx':idx})
    return out


def refine_linker(domain_set, linker_set, chain, domain_table):
    """Check if the linker is reasonable. Identify undetected domains.
    Then re-partition the domain and linker relationship.
    
    Args:
        domain_set (dict): set of domain segment sequences.
        linker_set (dict): set of linker segment sequences.
        chain (list): domain chain.
        domain_table (pandas.DataFrame): domain annotation table.
    """

    faulty_seq = [] # align later
    for k in linker_set.keys():
        ls = linker_set[k]
        lk_len = list(map(lambda x:len(ls[x]), ls.keys()))

        # If the linker is too long, then the domain annotation may not be reliable.
        if len(lk_len) > 2:
            Q1, Q3 = np.percentile(lk_len, 25), np.percentile(lk_len, 75)
            IQR = Q3 - Q1
            upper = lk_len >= (Q3 + 1.5 * IQR)
            outlier = list(np.array(list(ls.keys()))[upper])

            # If in the chain, the adjacent block is detected for these outliers, then the seq is still reliable.
            outlier_val = []
            for o in outlier:
                cidx = chain.index(k)
                check_block = []
                if cidx != 0:
                    check_block.append(chain[cidx-1])
                if cidx != (len(chain)-1):
                    check_block.append(chain[cidx+1])

                detected_dom = list(domain_table[domain_table['seq_id'] == o].clan)
                if len([j for j in check_block if j in detected_dom]) != len(check_block):
                    outlier_val.append(o)

            outlier = outlier_val.copy()
            faulty_seq.extend(outlier)

            # TODO: alternatively, refine domain and linker?
            """
            if len(outlier) > 0:
                idx = chain.index(k)
                seq_to_test = {}
                for ol in outlier:
                    seq_to_test[ol] = ls[ol]
                if idx != len(chain)-1:
                    # overlap with the next domain
                    fwd = check_linker_in_domain(domain_set[chain[idx+1]], 
                                                 seq_to_test, direction='start')
                if idx != 0:
                    # overlap with the last domain
                    back = check_linker_in_domain(domain_set[chain[idx-1]], 
                                                  seq_to_test, direction='end')
            """

    # remove unreliable sequences from the alignment blocks
    for fs in faulty_seq:
        for dom in domain_set.keys():
            if fs in domain_set[dom].keys():
                del domain_set[dom][fs]
        for lnk in linker_set.keys():
            if fs in linker_set[lnk].keys():
                del linker_set[lnk][fs]

    # check if there is any chain or linker missing after removing faulty_seq
    for i in range(len(chain)):
        if chain[i][0] == '_':
            if len(linker_set[chain[i]]) == 0:
                del linker_set[chain[i]]
                chain[i] = 'null'
        else:
            if len(domain_set[chain[i]]) == 0:
                del domain_set[chain[i]]
                chain[i] = 'null'

    return faulty_seq, domain_set, linker_set, chain
