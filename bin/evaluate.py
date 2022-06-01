# Author: Roll 820

import subprocess
from Bio.Seq import Seq
from Bio import AlignIO, SeqIO
import uuid
import os
import re
import pandas as pd
import numpy as np

from utils import convert_aln


def fasta2msf(in_file, out_file):
    """(deprecated) Convert MSA in fasta format to msf format, require EMBOSS :command:`aligncopy`.
    """
    _ = subprocess.run('aligncopy ' + in_file + ' ' + out_file, shell=True)
    return


def msf2fasta(in_file, out_file):
    """Convert MSA in msf format to fasta format.

    Args:
        in_file (str): Input file path, must end with :file:`.msf`.
        out_file (str): Output file path, must end with :file:`.fasta`.
    """
    aln = AlignIO.read(in_file, 'msf')
    _ = AlignIO.write(aln, out_file, 'fasta')
    return


def run_fastSP(ref_file, test_files, test_file_dir=None, fastSP_jar_path=None):
    """Run fastSP to evaluate MSA performance.
    Args:
        ref_file (str): path to reference alignment file.
            If reference is .msf file, will convert to fasta temp file.
        test_files (list): list of alignment fasta files to test on reference alignment.
        test_file_dir (str): home directory of test files. If None, directly use `test_files`.
        fastSP_jar_path (str): specify path to :file:`FastSP.jar`.
    """
    ref_end = os.path.basename(ref_file).split('.')[-1]
    ref_dir = os.path.dirname(ref_file)
    if ref_end == 'msf':
        filename = uuid.uuid4().hex
        out_file = '.' + filename + '_ref.fasta'
        ref_file_old = ref_file
        ref_file = os.path.join(ref_dir, out_file)
        msf2fasta(in_file=ref_file_old, out_file=ref_file)

    if test_file_dir is not None:
        test_files = list(map(lambda x:os.path.join(test_file_dir, x), test_files))

    try:
        if fastSP_jar_path is not None:
            PATH_TO_FASTSP_JAR = fastSP_jar_path
        else:
            PATH_TO_FASTSP_JAR = 'bin/external/FastSP/FastSP.jar'

        df_coll = pd.DataFrame()
        for tf in test_files:
            tf_name = os.path.basename(tf).split('.')[0]
            cmd = 'java -jar {} -r {} -e {}'.format(PATH_TO_FASTSP_JAR, ref_file, tf)
            out = subprocess.getoutput(cmd)
            out = out.split('\n')
            idx = 0
            for i in range(len(out)):
                if re.search('SP-Score', out[i]):
                    idx = i
                    break
            out = out[idx:-1]
            df = []
            for i in range(len(out)):
                ro = re.search('^(\D+) ([\d\\.]+)$', out[i])
                if ro is not None:
                    df.append({'metric': ro.group(1), 'value': float(ro.group(2))})
                else:
                    ro = re.search('^(\D+) NaN', out[i])
                    df.append({'metric': ro.group(1), 'value': 'NA'})
            df = pd.DataFrame(df)
            df['query'] = tf_name
            df_coll = pd.concat([df_coll, df], ignore_index=True)
    finally:
        if ref_end == 'msf':
            _ = subprocess.run('rm ' + ref_file, shell=True)

    return df_coll


def evaluate_Prefab(ref_file, test_files, test_prefix):
    """Evaluate performance of Prefab.
    
    According to original publication (Edgar 2004, NAR), the accuracy is
    evaluated on the original pair, after align them together with blast hits.
    """
    
    ref = AlignIO.read(ref_file, 'fasta')
    ref = convert_aln(ref)
    ref_nm = []
    for r in ref:
        ref_nm.append(r.id)
    
    new_test_files = [] # temp files with reference only.
    count = 0
    for tf in test_files:
        test = AlignIO.read(tf, 'fasta')
        new_test = []
        new_test_str = []
        for i in test:
            if i.id in ref_nm:
                new_test.append(i)
                new_test_str.append(list(str(i.seq)))
        if len(new_test) != len(ref_nm):
            raise ValueError('Reference not found in test set.')

        # collapse gaps in all test sequences
        new_test_str = np.array(new_test_str)
        collapse_test = [True for _ in range(new_test_str.shape[1])]
        for i in range(new_test_str.shape[1]):
            allchar = list(np.unique(new_test_str[:,i]))
            if len(allchar) == 1 and allchar[0] == '-':
                collapse_test[i] = False
        new_test_str = new_test_str[:,collapse_test]
        for i in range(new_test_str.shape[0]):
            sq = new_test[i]
            sq.seq = Seq(''.join(new_test_str[i,:]))
            new_test[i] = sq

        filename = uuid.uuid4().hex
        out_file = test_prefix[count] + '.' + filename + '_ref.fasta'
        SeqIO.write(new_test, out_file, 'fasta')
        new_test_files.append(out_file)
        count += 1
    
    try:
        df = run_fastSP(ref_file, new_test_files)
    except Exception as e:
        print(e)
        return None
    finally:
        for f in new_test_files:
            _ = subprocess.run('rm ' + f, shell=True)

    return df
