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
            for idx, base in enumerate(genome_sa.seq):
                ref_genome_base = ref_genome_sa.seq[idx]
                if base != '-':
                    genome_idx += 1
                if ref_genome_base != '-':
                    ref_genome_idx += 1
                if base == ref_genome_base and genome_idx < genome_length:
                    idx_lut_arr[genome_idx] = ref_genome_idx

    return idx_lut_arr


def cleanIndex(idx_arr, genome, ref_genome, cleaning_radius=4):
    ''' Remove scars in the index mapping that result from mismatches
    and / or border regions surrounding inserts / duplications.

    '''
    for idx, mapped_idx in enumerate(idx_arr):
        if mapped_idx == -1: 
            lower_idx = idx - 1
            while lower_idx > 0 and ((idx - lower_idx) < cleaning_radius):
                if idx_arr[lower_idx] == -1:
                    lower_idx -= 1
                    continue
                for upper_idx in range(idx, idx + cleaning_radius + 1):
                    if upper_idx > idx_arr.shape[0]:
                        break
                    if idx_arr[upper_idx] == -1:
                        continue
                    idx_range = upper_idx - lower_idx
                    if (idx_arr[upper_idx] - idx_arr[lower_idx] == 
                            (idx_range)):
                        idx_arr[idx] = idx_arr[lower_idx] + (idx - lower_idx)
                        break
                    # Handle "edge" cases
                    lower_bound_edge_idx = idx_arr[lower_idx] + (idx-lower_idx)
                    if (ref_genome[lower_bound_edge_idx] == genome[idx]):
                        idx_arr[idx] = lower_bound_edge_idx
                        break
                    upper_bound_edge_idx = idx_arr[upper_idx] + (upper_idx-idx)
                    if (ref_genome[upper_bound_edge_idx] == genome[idx]):
                        idx_arr[idx] = upper_bound_edge_idx
                        break
                lower_idx -= 1

            # l_area = idx_arr[idx-4:idx+5]  # 9 el arr centered on idx
            # if l_area.shape[0] == 9:
            #     # Large area smoothing
            #     if l_area[8] - l_area[0] == 8:
            #         for i in range(9):
            #             l_area[i] = l_area[0] + i
            #     # Local area smoothing
            #     elif l_area[5] - l_area[3] == 2:
            #         l_area[4] = l_area[3] + 1


index_lut = buildIndex('test_mut_dup.fa', 'baseseq.fa')

genome1 = parseFasta('test_mut_dup.fa')[0][1]
genome2 = parseFasta('baseseq.fa')[0][1]

cleanIndex(index_lut, genome1, genome2)

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
