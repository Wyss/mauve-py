from __future__ import print_function

import tempfile
import os
import subprocess
import shutil

from functools import partial
import numpy as np
from itertools import starmap

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
        sa = list(filter(lambda a: a.seq_num == seq_num, sub_alignment_group.alignments))[0]
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
    idx_lut_arr = np.empty(genome_length, dtype=np.int32)
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
    # Remove .sslist files, if created
    for fp in [genome_fp, ref_genome_fp]:
        try:
            os.remove(os.path.splitext(fp)[0] + '.sslist')
        except OSError:
            pass

    return idx_lut_arr


def findGaps(idx_arr):
    ''' Find and print gaps in the idx_arr
    '''
    gap_arr = np.zeros(idx_arr.shape[0], dtype=np.int32)
    try:
        for idx, mapped_idx in enumerate(idx_arr[1:], start=1):
            lower_idx = idx - 1
            lower_mapped_idx = idx_arr[idx-1]
            if mapped_idx == -1 and lower_mapped_idx != -1:
                success = False
                for upper_idx, upper_mapped_idx in enumerate(idx_arr[idx+1:],
                                                             start=idx+1):
                    if (upper_idx-idx > 300):
                        break
                    if ((upper_mapped_idx - lower_mapped_idx) ==
                       (upper_idx - lower_idx)):
                        gap_size = upper_idx - lower_idx
                        idx_arr[lower_idx:upper_idx+1] = np.arange(lower_mapped_idx, upper_mapped_idx+1, dtype=int)
                        gap_arr[idx] = gap_size
                        success = True
                        break
                if not success:
                    gap_arr[idx] = -1
    except KeyboardInterrupt:
        print(idx)
        raise
    return gap_arr

def findEdges(idx_arr):
    ''' Find and print edges in the idx_arr
    return a 1 on the 5 prime most index of a segment
    uses slicing
    '''
    edge_arr = np.zeros(idx_arr.shape[0], dtype=np.uint8)
    edge_arr_view = edge_arr[1:]    # slice not a copy
    delta_map_idxs = ind_arr[1:] - ind_arr[:-1]
    edge_arr_view[delta_map_idxs != 1] = 1
    edge_arr[0] = 1     # 0 index is always an edge
    return edge_arr
#end def

def findMismatches(idx_arr, genome, ref_genome):
    mismatch_arr = np.zeros(idx_arr.shape[0], dtype=np.uint8)
    for idx, base in enumerate(genome):
        if base != ref_genome(idx_arr[idx]):
            mismatch_arr[idx] = 1
    return mismatch_arr
# end def


# def cleanIndex(idx_arr, genome, ref_genome, cleaning_radius=4):
#     ''' Remove scars in the index mapping that result from mismatches
#     and / or border regions surrounding inserts / duplications.

#     '''
#     if idx_arr[0] == -1:
#         for upper_idx in range(1, cleaning_radius + 1):
#             if idx_arr[upper_idx] != -1:
#                 upper_bound_edge_idx = (idx_arr[upper_idx] -
#                                         upper_idx)
#             if (ref_genome[upper_bound_edge_idx] == genome[0]):
#                 idx_arr[0] = upper_bound_edge_idx

#     for idx, mapped_idx in enumerate(idx_arr[1:], start=1):

#         if mapped_idx == -1:

#             lower_idx = idx - 1
#             if idx_arr[lower_idx] != -1:
#                 lower_bound_edge_idx = idx_arr[lower_idx] + 1
#                 if (ref_genome[lower_bound_edge_idx] == genome[idx]):
#                     idx_arr[idx] = lower_bound_edge_idx
#                     continue

#                 for upper_idx in range(idx+1, idx + cleaning_radius + 1):
#                     if upper_idx > idx_arr.shape[0]:
#                         break
#                     if idx_arr[upper_idx] == -1:
#                         continue
#                     idx_range = upper_idx - lower_idx
#                     if (idx_arr[upper_idx] - idx_arr[lower_idx] ==
#                             (idx_range)):
#                         idx_arr[idx] = idx_arr[lower_idx] + 1
#                         break

#             upper_bound_edge_idx = idx_arr[idx+1] - 1
#             if (idx_arr[idx+1] != -1 and
#                     ref_genome[upper_bound_edge_idx] == genome[idx]):
#                 idx_arr[idx] = upper_bound_edge_idx
#                 continue


def test():
    import sys

    OUTPUT_TO_FN = False

    # GENOME = 'test_mut_dup.fa'
    # REF_GENOME = 'baseseq.fa'
    # OUTPUT_FN = 'test_mut_dup.out'

    GENOME = os.path.join(DIR, '2014_02_18_gen9_54_failed_seg_fixed_with_gc_fixes.fa')
    REF_GENOME = os.path.join(DIR, 'mds42_full.fa')
    OUTPUT_FN = os.path.join(DIR, '54_failed_seg_test-gapfill.out')

    if OUTPUT_TO_FN:
        print_fd = open(OUTPUT_FN, 'w')
    else:
        print_fd = sys.stdout

    # index_lut = buildIndex(GENOME, REF_GENOME)

    # np.save('index_lut.temp', index_lut)

    index_lut = np.load('index_lut.temp.npy')

    genome1 = parseFasta(GENOME)[0][1]
    genome2 = parseFasta(REF_GENOME)[0][1]

    # cleanIndex(index_lut, genome1, genome2, cleaning_radius=11)

    def _idxLookup(idx):
        lut_idx = index_lut[idx]
        return genome2[lut_idx] if lut_idx != -1 else '-'

    def _findEdge(idx):
        character = ' '
        try:
            if index_lut[idx] == -1:
                pass
            elif index_lut[idx] != index_lut[idx+1] - 1:
                character = '|'
            elif idx > 0 and index_lut[idx] != index_lut[idx-1] + 1:
                character = '|'
        except IndexError:
            pass
        finally:
            return character

    findGaps(index_lut)

    diff_func = lambda g1, g2: ' ' if g1 == g2 else '*'

    # for i in range(len(genome1) // 70):
    #     genome_1_line = ''.join([genome1[i * 70 + y] for y in range(70)])
    #     genome_2_line = ''.join([_idxLookup([i * 70 + y]) for y in range(70)])
    #     diff_line = ''.join(starmap(diff_func, zip(genome_1_line, genome_2_line)))
    #     edge_line = ''.join(map(_findEdge, [i * 70 + y for y in range(70)]))
    #     print(i * 70,file=print_fd)
    #     print(genome_1_line,file=print_fd)
    #     print(genome_2_line,file=print_fd)
    #     print(diff_line,file=print_fd)
    #     print(edge_line,file=print_fd)
    #     print(file=print_fd)

    if OUTPUT_TO_FN:
        print_fd.close()

    for idx in range(353780, 353851):
        print(idx, index_lut[idx])


if __name__ == '__main__':
    # import cProfile as profile
    # profile.run('test()')
    test()
