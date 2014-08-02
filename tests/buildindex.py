from __future__ import print_function

import tempfile
import os
import subprocess
import shutil

from functools import partial 
import numpy as np

from util import parseFasta
from xmfa import parseXMFA


DIR = os.path.dirname(os.path.realpath(__file__))
PROGMAUVE_DIR = os.path.join(os.path.dirname(DIR), 'bin', 
                             'progressiveMauveStatic')


class TempDirs(object):

    def __init__(self, num_dirs):
        self.num_dirs = num_dirs

    def __enter__(self):
        self.temp_dirs = [tempfile.mkdtemp() for _ in range(self.num_dirs)]
        return self.temp_dirs

    def __exit__(self, *exc):
        for tmp_dir in self.temp_dirs:
            try:
                shutil.rmtree(tmp_dir)
            except OSError as exc:
                if exc.errno != errno.ENOENT:
                    raise 
        if exc[0] is not None:
            raise exc[1]


def runMauve(fa1, fa2, flags={}):
    with TempDirs(2) as (temp_dir1, temp_dir2):
        temp_dirs_args = ['--scratch-path-1', temp_dir1,
                          '--scratch-path-2', temp_dir2]
        args = [PROGMAUVE_DIR, fa1, fa2] + temp_dirs_args + \
               [str(el) for pair in flags.items() for el in pair]
        proc = subprocess.Popen(args, stdout=subprocess.PIPE, 
                                stdin=subprocess.PIPE, cwd=temp_dir1)
        proc.wait()


def getGenomeLength(genome_fa):
    data = parseFasta(genome_fa)
    return len(data[0][1])


def _processBackboneFile(backbone_fp):
    with open(backbone_fp) as fd:
        fd.readline()  # Throw away the header
        raw_backbone_data = [line.strip().split('\t') for line in fd]
    # Convert values to integers and convert to 0-based indexing
    return map(partial(map, lambda x: int(x) - 1), raw_backbone_data)

def lookupSubAlignment(seq_num, sub_alignment_group):
    sa = None
    try:
        sa = filter(lambda a: a.seq_num == seq_num, sub_alignment_group.alignments)[0]
    except IndexError:
        pass
    return sa

def buildIndex(genome_fa, ref_genome_fa):
    # Generate absolute filepaths
    genome_fp = os.path.abspath(genome_fa)
    ref_genome_fp = os.path.abspath(ref_genome_fa)
    with TempDirs(1) as (temp_dir,):
        output_fp = os.path.join(temp_dir, 'mauveout.xmfa')
        runMauve(genome_fp, ref_genome_fp, flags={'--output': output_fp})
        sub_alignment_groups = parseXMFA(output_fp)
    genome_length = getGenomeLength(genome_fp)
    # LUT to map indices from the engineered genome to the reference genome
    idx_lut_arr = np.empty(genome_length, dtype=int)
    # Initialize all values to -1 (un-mapped indices will be -1)
    idx_lut_arr[:] = -1
    for sub_alignment_group in sub_alignment_groups:
        genome_sa = lookupSubAlignment(1, sub_alignment_group)
        ref_genome_sa = lookupSubAlignment(2, sub_alignment_group)
        if genome_sa is not None and ref_genome_sa is not None:
            # Note that we have to convert to 0-based indexing
            genome_idx = genome_sa.start_idx - 1
            ref_genome_idx = ref_genome_sa.start_idx - 1 
            print(genome_idx, ref_genome_idx)
            for idx, base in enumerate(genome_sa.seq):
                ref_genome_base = ref_genome_sa.seq[idx]
                if base != '-':
                    genome_idx += 1
                if ref_genome_base != '-':
                    ref_genome_idx += 1
                if base == ref_genome_base:
                    try:
                        idx_lut_arr[genome_idx] = ref_genome_idx
                    except IndexError:
                        # TODO: Figure out why we're overrunning the array...
                        print('Array overrun: ', genome_idx, ref_genome_idx)

    return idx_lut_arr

    
index_lut = buildIndex('test_1swap.fa', 'baseseq.fa')

genome1 = parseFasta('test_1swap.fa')[0][1]
genome2 = parseFasta('baseseq.fa')[0][1]

def _idxLookup(idx):
    lut_idx = index_lut[idx]
    return genome2[lut_idx] if lut_idx != -1 else '-'

diff_func = lambda (g1, g2): ' ' if g1 == g2 else '*'

for i in range(len(genome1) / 70):
    genome_1_line = ''.join([genome1[i * 70 + y] for y in range(70)])
    genome_2_line = ''.join([_idxLookup([i * 70 + y]) for y in range(70)])
    diff_line = ''.join(map(diff_func, zip(genome_1_line, genome_2_line)))
    print(i * 70)
    print(genome_1_line)
    print(genome_2_line)
    print(diff_line)
    print()
