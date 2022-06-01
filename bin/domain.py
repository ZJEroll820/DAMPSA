# Author: Roll 820

import uuid
import logging
import os
import subprocess
import requests
import re
import numpy as np
import pandas as pd
from Bio import SeqIO

from utils import get_within_pair


def run_hmmscan(seqs, n_thread=4):
    """Run Hmmscan for input sequences.

    Args:
        seqs (dict): dictionary of :class:`Bio.SeqRecord` protein sequences.
        n_thread (int): number of thread for running :command:`hmmscan`.
    """
    logger = logging.getLogger('DAMPSA.domain')
    path_to_msa = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

    if not isinstance(seqs, dict):  # single record
        seqs = {seqs.id: seqs}
    filename = uuid.uuid4().hex
    fnm = '.' + filename + '.fasta'
    fnm_out = '.' + filename + '.hmmscan'
    _ = SeqIO.write(seqs.values(), os.path.join(path_to_msa, fnm), 'fasta')
    cmd = 'hmmscan --domtblout {} --cut_ga --cpu {} --noali {} {}'.format(
    #cmd = 'hmmscan --domtblout {} --domE 0.01 --cpu {} --noali {} {}'.format(
        os.path.join(path_to_msa, fnm_out),
        n_thread,
        os.path.join(path_to_msa, 'data/Pfam_scan_db/Pfam-A.hmm'),
        os.path.join(path_to_msa, fnm))

    logger.info('Running hmmscan...')
    try:
        _ = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE)
    finally:
        _ = subprocess.run('rm ' + os.path.join(path_to_msa, fnm), shell=True)
    try:
        tb = read_hmm_scan(os.path.join(path_to_msa, fnm_out))
    finally:
        _ = subprocess.run('rm ' + os.path.join(path_to_msa, fnm_out), shell=True)
    if tb.shape[0] == 0:
        raise ValueError('No domain annotation from hmmscan!')

    logger.info('Annotating domain...')
    # Refine hmmscan output (annotate clan, collapse same clan)
    tb = ann_clan(tb, os.path.join(path_to_msa, 'data/Pfam_scan_db/Pfam-A.hmm.dat'))
    tb = collapse_clan(tb)
    return tb


def read_hmm_scan(fp_hmm_scan):
    """Read :command:`hmmscan` output file in :option:`domtbout` format.

    Args:
        fp_hmm_scan (str): filepath to the :command:`hmmscan` output.

    Returns:
        Domain annotation table.
    """
    # tb = pd.read_table(fp_hmm_scan, sep='\s+', comment='#', header=None)
    tb = pd.read_table(fp_hmm_scan, comment='#', header=None)
    tb_ls = []
    for i in tb.index:
        tb_ls.append(re.split('\s+', tb.loc[i][0])[:22])
    tb = pd.DataFrame(tb_ls)
    thead = ['target name', 'accession', 'tlen', 'query name', 'accession2', 'qlen',
             'E-value_full', 'score_full', 'bias_full', '#', 'of', 'c-Evalue', 'i-Evalue', 'score_dom', 'bias_dom',
             'from_hmm', 'to_hmm', 'from_ali', 'to_ali', 'from_env', 'to_env', 'acc', 'description of target']
    thead = thead[:22]
    # tb = tb.iloc[:,:22]
    tb.columns = thead
    return tb


def ann_clan(hmm_scan_tb, fp_pfam_dat):
    """Annotate pfam domain to the Clan information

    Args:
        hmm_scan_tb (pandas.DataFrame): hmmscan output table, returned from :code:`read_hmm_scan`.
        fp_pfam_dat (str): file path to the :file:`Pfam-A.hmm.dat`.

    Returns:
        Domain table annotated with Clan information.
    """
    all_record = []
    with open(fp_pfam_dat, encoding='utf-8') as f:
        dom = {}
        while True:
            line = f.readline()
            if line:
                if line == '//\n':
                    if 'CL' not in dom.keys():
                        dom['CL'] = None
                    all_record.append(dom)
                    dom = {}
                elif re.search('=', line):
                    info = re.split('\s\s+', line.strip().split('=')[1])
                    dom[info[0].split(' ')[1]] = info[1]
            else:
                break
    all_record = pd.DataFrame(all_record)
    merged = pd.merge(left=hmm_scan_tb, right=all_record[['AC', 'DE', 'CL']], how='left', left_on='accession',
                      right_on='AC')

    # for those without clan info, use PF domain ID directly
    for i in merged.index:
        cl = merged.loc[i].CL
        if cl != cl or cl is None:
            merged.loc[i].CL = merged.loc[i].accession

    return merged


def collapse_clan(hmm_scan_tb_clan, rm_within=True):
    """Collapse domains that belong to the same Clan family.

    Args:
        hmm_scan_tb_clan (pandas.DataFrame): :command:`hmmscan` output table annotated with clan info.
        rm_within (bool): whether remove regions within other regions, like sub-domain.
    Returns:
        Domain table with domains grouped by Clan.
    """

    out = []
    for seq in np.unique(hmm_scan_tb_clan['query name']):
        sub = hmm_scan_tb_clan[hmm_scan_tb_clan['query name'] == seq]
        clan_count = np.unique(sub.CL, return_counts=True)
        dup_clan = clan_count[0][clan_count[1] != 1]
        uni_clan = clan_count[0][clan_count[1] == 1]
        uni_clan = sub[sub['CL'].isin(uni_clan)]
        for i in range(uni_clan.shape[0]):
            out.append({'seq_id': seq,
                        'clan': uni_clan.iloc[i].CL,
                        # Use `env` range (.from_env) - envelope or the full alignment
                        # Use `aln` range (.from_ali) - high confidence alignment
                        'env_start': uni_clan.iloc[i].from_env,
                        'env_end': uni_clan.iloc[i].to_env,
                        'PF': uni_clan.iloc[i].accession})
        for dc in dup_clan:
            dc_sub = sub[sub.CL == dc][['from_env', 'to_env', 'accession']]
            new_from = np.min(np.array(dc_sub['from_env'], dtype=int))
            new_end = np.max(np.array(dc_sub['to_env'], dtype=int))
            out.append({'seq_id': seq, 'clan': dc,
                        'env_start': new_from, 'env_end': new_end,
                        'PF': '_'.join(dc_sub.accession)})
    out = pd.DataFrame(out)
    out['env_start'] = out['env_start'].astype('int')
    out['env_end'] = out['env_end'].astype('int')

    if rm_within:
        to_rm = []
        for i in np.unique(out['seq_id']):
            ssub = out[out['seq_id'] == i][['env_start', 'env_end']]
            ssub = np.transpose(ssub).to_dict()
            check = []
            for s in ssub.keys():
                check.append(tuple(ssub[s].values()))
            wp = get_within_pair(check)
            for w in wp:
                to_rm.append(list(ssub.keys())[check.index(w)])
        out = out.drop(labels=to_rm)
    return out


def get_pfam_from_Uniprot(uniprot_entry):
    """(deprecated) Get Pfam information from UniProt database (cross-reference).

    Args:
        uniprot_entry (str): UniProt Entry (ID) of the protein.

    Returns:
        A list of Pfam domain IDs.
    """
    payload = {'query': 'id:'+uniprot_entry, 'columns': 'database(Pfam)', 'format':'tab'}
    response = requests.get('https://www.uniprot.org/uniprot', params=payload)
    if response.status_code == 200:
        pfams = response.text.split('\n')[1].split(';')
        pfams = pfams[:-1]
        return pfams
    else:
        raise ConnectionError('Failed to connect to UniProt or invalid UniProt entry!')


def read_pfam_scan(fp_pfam_scan):
    """(deprecated) Read pfamscan.pl output.
    """
    df = pd.read_table(fp_pfam_scan, comment='@', skiprows=27, header=None, sep='\s+')
    df = df.iloc[:, 0:15]
    df.columns = ['seq_id', 'aln_start', 'aln_end', 'env_start', 'env_end', 'hmm_acc', 'hmm_name', 'type', 'hmm_start',\
                  'hmm_end', 'hmm_length', 'bit_score',
                  'E_value', 'significance', 'clan']
    df = df[df['type'] == 'Domain']
    df = df.astype({'seq_id': str})
    return df
