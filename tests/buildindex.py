from __future__ import print_function

import errno
import tempfile
import os
import subprocess
import shutil

import numpy as np
from itertools import starmap

from util import parseFasta
from xmfa import parseXMFA


DIR = os.path.dirname(os.path.realpath(__file__))

MAUVE_DIR = (os.environ.get('MAUVE_DIR') or 
            os.path.join(os.path.dirname(DIR), 'bin'))

PROG_MAUVE_STATIC_FP = os.path.join(MAUVE_DIR, 'progressiveMauveStatic')

if not (os.path.isfile(PROG_MAUVE_STATIC_FP) and 
        os.access(PROG_MAUVE_STATIC_FP, os.X_OK)):
    raise IOError('Could not find progressiveMauveStatic in {}. Environmental '
                  'variable `MAUVE_DIR` is defined as: {}.'.format(
                   MAUVE_DIR, os.environ.get('MAUVE_DIR')))


# Non-mask array data type
# 32-bit supports genomes up to 4 billion bases in length
IDX_ARRAY_DTYPE = np.int32
# Mask array data type (should only have values of 0 + 1 so unit8 is good)
MASK_ARRAY_DTYPE = np.uint8


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


def runMauve(fasta_files, flags={}):
    ''' Run progressiveMauveStatic with the fasta files specifies in the list
    `fasta_files`, and command line flags `flags`.

        ex. runMauve(['k12_modded.fa', 'k12_ref.fa'], 
                     flags={'output'=k12_diff.xmfa'})
    '''
    abs_paths = [os.path.abspath(fp) for fp in fasta_files]

    with TempDirs(2) as (temp_dir1, temp_dir2):
        flags.update({'--scratch-path-1': temp_dir1,
                      '--scratch-path-2': temp_dir2})
        args = ([os.path.join(MAUVE_DIR, 'progressiveMauveStatic')] + 
                abs_paths + [str(el) for pair in flags.items() for el in pair])
        proc = subprocess.Popen(args, stdout=subprocess.PIPE,
                                stdin=subprocess.PIPE, cwd=temp_dir1)
        proc.wait()

    # Remove .sslist files, if created
    for fp in abs_paths:
        try:
            os.remove(os.path.splitext(fp)[0] + '.sslist')
        except OSError:
            pass


def getGenomeLength(genome_fa):
    ''' Return the genome length from a fasta file containing a single
    header / sequence
    '''
    data = parseFasta(genome_fa)
    return len(data[0][1])


def lookupSubAlignment(seq_num, sub_alignment_group):
    sa = None
    try:
        sa = list(filter(lambda a: a.seq_num == 
                  seq_num, sub_alignment_group.alignments))[0]
    except IndexError:
        pass
    return sa


def fillGaps(idx_lut, max_gap_width=300):
    ''' Fill gaps in the alignment that are up to `max_gap_width` in size. A 
    gap is only filled if flanking indices can be found such that:

        upper_idx - lower_idx == upper_mapped_idx - lower_mapped_idx

    '''
    for idx, mapped_idx in enumerate(idx_lut[1:], start=1):
        if mapped_idx == -1:
            lower_mapped_idx = idx_lut[idx-1]
            if lower_mapped_idx != -1:
                lower_idx = idx - 1
                for upper_idx, upper_mapped_idx in enumerate(idx_lut[idx+1:],
                                                             start=idx+1):
                    if (upper_idx-idx > max_gap_width):
                        break
                    if ((upper_mapped_idx - lower_mapped_idx) ==
                            (upper_idx - lower_idx)):
                        idx_lut[lower_idx:upper_idx+1] = np.arange(
                            lower_mapped_idx, upper_mapped_idx+1, 
                            dtype=IDX_ARRAY_DTYPE)
                        break


def smoothEdges(idx_lut, smoothing_radius=20):
    ''' Corrects issues that result in messy / incorrect alignments in areas
    of high mismatch density. Mauve sometimes chokes a little bit when there 
    are multiple mismatches in close proximity to one another, so we can 
    apply similar principles to those used in fillGaps to "smooth" out these
    areas.

                            TATGTGGATGTCAGGTCAAAGCG
                            TATGTGGAGGTGAGCTCAAAGCG
    Mismatches:                     *  *  *        
    Index mapping edges:            |||||||   

    '''
    prev_mapped_idx = idx_lut[0]
    for idx, mapped_idx in enumerate(idx_lut[1:], start=1):
        if mapped_idx != prev_mapped_idx + 1:
            lower_mapped_idx = idx_lut[idx-1]
            if lower_mapped_idx != -1:
                lower_idx = idx - 1
                for upper_idx, upper_mapped_idx in enumerate(idx_lut[idx+1:],
                                                             start=idx+1):
                    if (upper_idx-idx > smoothing_radius):
                        break
                    if ((upper_mapped_idx - lower_mapped_idx) ==
                            (upper_idx - lower_idx)):
                        idx_lut[lower_idx:upper_idx+1] = np.arange(
                            lower_mapped_idx, upper_mapped_idx+1, 
                            dtype=IDX_ARRAY_DTYPE)
                        break


def buildIndex(genome_fa, ref_genome_fa, fill_gaps=True, max_gap_width=300,
               smooth_edges=True, smoothing_width=20):
    ''' Builds and returns an index mapping from the genome in `genome_fa` to
    the genome in `ref_genome_fa`. If `fillGaps` is True, gaps up to 
    `max_gap_width` will be filled if flanking indices can be found with 
    commensurate map spacing (see `fillGaps` docstring). If `smooth_edges` is
    True, areas of high mismatch frequency with gap-penalty-related problems
    will be smoothed.
    '''
    # Generate absolute filepaths
    genome_fp = os.path.abspath(genome_fa)
    ref_genome_fp = os.path.abspath(ref_genome_fa)
    # Run mauve and parse output (all output files are cleaned up after run / 
    # parsing)
    with TempDirs(1) as (temp_dir,):
        output_fp = os.path.join(temp_dir, 'mauveout.xmfa')
        runMauve(genome_fp, ref_genome_fp, flags={'--output': output_fp})
        sub_alignment_groups = parseXMFA(output_fp)
    genome_length = getGenomeLength(genome_fp)

    # LUT to map indices from the engineered genome to the reference genome
    idx_lut = np.empty(genome_length, dtype=IDX_ARRAY_DTYPE)
    # Initialize all values to -1 (un-mapped indices will be -1)
    idx_lut[:] = -1

    # Parse sub-alignments into index mapping
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
                    idx_lut[genome_idx] = ref_genome_idx
    if fill_gaps:
        fillGaps(idx_lut, max_gap_width=max_gap_width)
    if smooth_edges:
        smoothEdges(idx_lut, smoothing_width=smoothing_width)
    return idx_lut


def findGaps(idx_lut):
    ''' Find and print gaps in the idx_lut
    '''
    gap_arr = np.zeros(idx_lut.shape[0], dtype=np.int32)
    try:
        for idx, mapped_idx in enumerate(idx_lut[1:], start=1):
            lower_idx = idx - 1
            lower_mapped_idx = idx_lut[idx-1]
            if mapped_idx == -1 and lower_mapped_idx != -1:
                success = False
                for upper_idx, upper_mapped_idx in enumerate(idx_lut[idx+1:],
                                                             start=idx+1):
                    if (upper_idx-idx > 300):
                        break
                    if ((upper_mapped_idx - lower_mapped_idx) ==
                       (upper_idx - lower_idx)):
                        gap_size = upper_idx - lower_idx
                        idx_lut[lower_idx:upper_idx+1] = np.arange(lower_mapped_idx, upper_mapped_idx+1, dtype=int)
                        gap_arr[idx] = gap_size
                        success = True
                        break
                if not success:
                    gap_arr[idx] = -1
    except KeyboardInterrupt:
        print(idx)
        raise
    return gap_arr


def findEdges(idx_lut):
    ''' Find and print edges in the idx_lut
    return a 1 on the 5 prime most index of a segment
    uses slicing
    '''
    edge_arr = np.zeros(idx_lut.shape[0], dtype=MASK_ARRAY_DTYPE)
    edge_arr_view = edge_arr[1:]    # slice not a copy
    delta_map_idxs = idx_lut[1:] - idx_lut[:-1]
    edge_arr_view[delta_map_idxs != 1] = 1
    edge_arr[0] = 1     # 0 index is always an edge
    return edge_arr


def findMismatches(idx_lut, genome, ref_genome):
    mismatch_arr = np.zeros(idx_lut.shape[0], dtype=MASK_ARRAY_DTYPE)
    for idx, base in enumerate(genome):
        if idx_lut[idx] != -1:
            if base != ref_genome[idx_lut[idx]]:
                mismatch_arr[idx] = 1
    return mismatch_arr


def findDuplicateMappings(idx_lut):
    ''' Find duplicate mappings in the idx_lut. If everything is working
    properly there should be 1:1, unambiguous mappings.
    '''
    arr_cpy = np.copy(idx_lut)
    arr_cpy = np.clip(arr_cpy, 0, 4294967295)
    mapped_idx_counts = np.bincount(arr_cpy)
    # -1 as 0 will show up here as we are clipping unmapped basses to 0
    num_duplicate_mappings = sum(mapped_idx_counts > 1) - 1
    return num_duplicate_mappings


def test():
    import sys

    OUTPUT_TO_FN = True

    # GENOME = 'test_mut_dup.fa'
    # REF_GENOME = 'baseseq.fa'
    # OUTPUT_FN = 'test_mut_dup.out'

    GENOME = os.path.join(DIR, '2014_02_18_gen9_54_failed_seg_fixed_with_gc_fixes.fa')
    REF_GENOME = os.path.join(DIR, 'mds42_full.fa')
    OUTPUT_FN = os.path.join(DIR, '54_failed_seg_test-testcleanup.out')

    if OUTPUT_TO_FN:
        print_fd = open(OUTPUT_FN, 'w')
    else:
        print_fd = sys.stdout

    # Cache the numpy array
    try:
        idx_lut = np.load('index_lut.npy')
    except IOError:
        idx_lut = buildIndex(GENOME, REF_GENOME, fill_gaps=False, 
                             smooth_edges=False)
        np.save('index_lut.temp', idx_lut)

    genome1 = parseFasta(GENOME)[0][1]
    genome2 = parseFasta(REF_GENOME)[0][1]

    # Test cleaning
    print('pre-cleaning duplicates:', findDuplicateMappings(idx_lut))
    fillGaps(idx_lut)
    smoothEdges(idx_lut)
    print('post-cleaning duplicates:', findDuplicateMappings(idx_lut))


    def _idxLookup(idx):
        lut_idx = idx_lut[idx]
        return genome2[lut_idx] if lut_idx != -1 else '-'

    def _findEdge(idx):
        character = ' '
        try:
            if idx_lut[idx] == -1:
                pass
            elif idx_lut[idx] != idx_lut[idx+1] - 1:
                character = '|'
            elif idx > 0 and idx_lut[idx] != idx_lut[idx-1] + 1:
                character = '|'
        except IndexError:
            pass
        finally:
            return character

    diff_func = lambda g1, g2: ' ' if g1 == g2 else '*'

    for i in range(len(genome1) // 70):
        genome_1_line = ''.join([genome1[i * 70 + y] for y in range(70)])
        genome_2_line = ''.join([_idxLookup([i * 70 + y]) for y in range(70)])
        diff_line = ''.join(starmap(diff_func, zip(genome_1_line, genome_2_line)))
        edge_line = ''.join(map(_findEdge, [i * 70 + y for y in range(70)]))
        print(i * 70,file=print_fd)
        print(genome_1_line,file=print_fd)
        print(genome_2_line,file=print_fd)
        print(diff_line,file=print_fd)
        print(edge_line,file=print_fd)
        print(file=print_fd)

    if OUTPUT_TO_FN:
        print_fd.close()


if __name__ == '__main__':
    # import cProfile as profile
    # profile.run('test()')
    test()
