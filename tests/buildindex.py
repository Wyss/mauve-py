from __future__ import print_function

import tempfile
import os
import subprocess
import shutil

from functools import partial 
import numpy as np

from util import parseFasta


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


def buildIndex(genome_fa, ref_genome_fa):
    # Generate absolute filepaths
    genome_fp = os.path.abspath(genome_fa)
    ref_genome_fp = os.path.abspath(ref_genome_fa)
    with TempDirs(1) as (temp_dir,):
        output_fn = os.path.join(temp_dir, 'mauveout.data')
        runMauve(genome_fp, ref_genome_fp, flags={'--output': output_fn})
        backbone_data = _processBackboneFile(output_fn + '.backbone')
    genome_length = getGenomeLength(genome_fp)
    # LUT to map indices from the engineered genome to the reference genome
    idx_lut_arr = np.empty(genome_length, dtype=int)
    # Initialize all values to -1 (un-mapped indices will be -1)
    idx_lut_arr[:] = -1
    for (g1_start_idx, g1_end_idx, g2_start_idx, g2_end_idx) in backbone_data:
        mapping = zip(range(g1_start_idx, g1_end_idx), 
                      range(g2_start_idx, g2_end_idx))
        for idx1, idx2 in mapping:
            idx_lut_arr[idx1] = idx2
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

print(index_lut[98])