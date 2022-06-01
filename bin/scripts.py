# Author: Roll 820

import os
os.chdir('/Users/jefft/Desktop/BioinformaticsProject/MSA/bin')
import re
from Bio import AlignIO
from Bio.Seq import Seq
import numpy as np
from tqdm import trange
from main import setup_logger
import logging
import pandas as pd

from utils import load_seq, load_Uniprot_series, convert_aln
from main import run_DAMPSA
os.chdir('../')

#%% Run BAliBASE

# load BAliBASE sequence
# single
# fps = ['data/bb_test/test_tfa/BB11018.tfa']
# escape = ['BB12001'] # not run for these prefix

# batch
home = 'data/bb_test/test_tfa'
fps = os.listdir(home)
fps = list(map(lambda x:os.path.join(home, x), fps))
escape = []
# escape = list(map(lambda x:x.split('_')[0], os.listdir('data/bb_test/test_domain_out')))

for i in trange(len(fps)):
    fp = fps[i]
    if re.search('DS_Store', fp):
        continue
    prefix = os.path.basename(fp).split('.')[0]
    if prefix in escape:
        continue
    print('=' * 20)
    print('Processing ' + prefix + '.')

    basedir = os.path.dirname(os.path.dirname(fp))
    fp_std_out = os.path.join(basedir, 'test_aln_out', prefix + '_clustalo_aln.fasta')
    fp_aln_out = os.path.join(basedir, 'test_aln_out', prefix + '_DAMPSA_aln.fasta')
    fp_ref_aln = os.path.join(basedir, 'test_msf', prefix + '.msf')
    fp_dom_out = os.path.join(basedir, 'test_domain_out', prefix + '_dom.txt')
    fp_evl_out = os.path.join(basedir, 'evl_score_out', prefix + '_score.txt')
    seqs = load_seq(fp)
    
    error = []
    try:
        run_DAMPSA(seqs, fp_aln_out, fp_std_out, fp_dom_out, fp_ref_aln, fp_evl_out, 
                   do_edge_refine=False, check_linker_len=True,
                   domain_app='clustalo', linker_app='clustalo', standard_app='clustalo')
    except Exception as e:
        print(e)
        error.append(fp)
        

#%% Run UniProt

# load human src seqeuence series
series_fp = 'data/src/src_euky.fasta'
tab_fp = 'data/src/src_euky.tab'  # optional, use if want to consider "review" info when collapsing organisms.
seqs = load_Uniprot_series(series_fp=series_fp, tab_fp=tab_fp)

# run main pipeline
# fp_std_out = 'data/src/src_filtered_clustalo.fasta'
fp_std_out = ''
fp_aln_out = 'data/src/src_filtered_DAMPSA_aln.fasta'
fp_dom_out = 'data/src/src_filtered_dom.txt'
fp_ref_aln = ''
fp_evl_out = ''

run_DAMPSA(seqs, fp_aln_out, fp_std_out, fp_dom_out, fp_ref_aln, fp_evl_out, 
           do_edge_refine=True, check_linker_len=True,
           domain_app='clustalo', linker_app='clustalo', standard_app='clustalo')


#%% Run HomFam
base = 'data/homfam-20110613-25'
out_base = 'data/homfam-out'

pfx = []
for i in os.listdir(base):
    if i != '.DS_Store':
        pfx.append('_'.join(i.split('_')[:-1]))
pfx = list(np.unique(pfx))


for t in trange(len(pfx)):
    prefix = pfx[t]
    print('Processing {}.'.format(prefix))
    
    seqs = AlignIO.read(os.path.join(base, prefix+'_ref.vie'), 'fasta')
    _ = AlignIO.write(seqs, os.path.join(out_base, 'ref', 
                                         prefix+'_ref.txt'),'fasta')
    raw_seqs = {}
    for i in seqs:
        i.seq = Seq(str(i.seq).replace('-',''))
        raw_seqs[i.id] = i
    
    fp_ref_aln = os.path.join(base, prefix+'_ref.vie')
    fp_dom_out = os.path.join(out_base, 'dom', prefix+'_dom.txt')
    fp_evl_out = os.path.join(out_base, 'evl', prefix+'_score.txt')
    fp_aln_out = os.path.join(out_base, 'align', prefix+'_DAMPSA.txt')
    fp_std_out = os.path.join(out_base, 'align', prefix+'_clustalo.txt')
    
    try:
        run_DAMPSA(raw_seqs, fp_aln_out, fp_std_out, fp_dom_out, 
                   fp_ref_aln, fp_evl_out, do_edge_refine=True,
                   domain_app='clustalo', linker_app='clustalo', 
                   standard_app='clustalo')
    except:
        print(pfx[t] + ' cannot be processed.')


#%% Run Prefab
from evaluate import evaluate_Prefab
import traceback

base = 'data/prefab4/in'
out_base = 'data/prefab4-out'

setup_logger(file_out='log.txt')
logger = logging.getLogger('DAMPSA')
pfx = os.listdir(base)
pfx = list(np.unique(pfx))
pfx_ftd = pfx # pick random 100? TODO


for t in trange(len(pfx_ftd)):
    prefix = pfx_ftd[t]
    if prefix not in pfx_ftd: # check
        continue
    print('Processing {}.'.format(prefix))
    
    seqs = load_seq(os.path.join(base, prefix))
    ref = AlignIO.read(os.path.join(os.path.dirname(base), 'ref', prefix), 'fasta')
    new_ref = convert_aln(ref)
    _ = AlignIO.write(new_ref, os.path.join(out_base, 'ref', 
                                            prefix+'_ref.txt'),'fasta')
    
    # Prefab needs a distinct evaluation pipeline
    fp_ref_aln = os.path.join(out_base, 'ref', prefix+'_ref.txt')
    fp_evl_out = os.path.join(out_base, 'evl', prefix+'_score.txt')
    fp_dom_out = os.path.join(out_base, 'dom', prefix+'_dom.txt')
    fp_aln_out = os.path.join(out_base, 'align', prefix+'_DAMPSA.txt')
    fp_std_out = os.path.join(out_base, 'align', prefix+'_clustalo.txt')
    
    '''
    run_DAMPSA(seqs, fp_aln_out, fp_std_out, fp_dom_out, 
               fp_ref_aln='', fp_evl_out='',
               do_edge_refine=False, check_linker_len=True,
               domain_app='clustalo', linker_app='clustalo', 
               standard_app='clustalo')

    tp = list(map(lambda x:prefix+x, ['_clustalo', '_DAMPSA']))
    evl = evaluate_Prefab(ref_file=fp_ref_aln, test_prefix=tp,
                          test_files=[fp_std_out, fp_aln_out])
    evl.to_csv(fp_evl_out, index=False, sep='\t')
    '''
    
    
    try:
        run_DAMPSA(seqs, fp_aln_out, fp_std_out, fp_dom_out, 
                   fp_ref_aln='', fp_evl_out='',
                   do_edge_refine=False, check_linker_len=True,
                   domain_app='clustalo', linker_app='clustalo', 
                   standard_app='clustalo')

        tp = list(map(lambda x:prefix+x, ['_clustalo', '_DAMPSA']))
        evl = evaluate_Prefab(ref_file=fp_ref_aln, test_prefix=tp,
                              test_files=[fp_std_out, fp_aln_out])
        evl.to_csv(fp_evl_out, index=False, sep='\t')
    except:
        logging.error(str(traceback.format_exc()))
        logging.info(pfx_ftd[t] + ' cannot be processed.')



