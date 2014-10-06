from __future__ import print_function

import errno
import tempfile
import os
import sys
import subprocess
import shutil

import bitarray
import numpy as np
from itertools import starmap

MAUVE_DIR = os.environ.get('MAUVE_DIR')
PROG_MAUVE_STATIC_FP = os.path.join(MAUVE_DIR, 'progressiveMauveStatic')

if not (os.path.isfile(PROG_MAUVE_STATIC_FP) and
        os.access(PROG_MAUVE_STATIC_FP, os.X_OK)):
    raise IOError('Cache could not find progressiveMauveStatic in {}. Environmental '
                  'variable `MAUVE_DIR` is defined as: {}.'.format(
                   MAUVE_DIR, os.environ.get('MAUVE_DIR')))

from .util import parseFasta
from .xmfa import parseXMFA
import mauve.indexutils as indexutils

def getSeqFromFile(seq_fp):
    _, ext = os.path.splitext(seq_fp)
    if ext in ('.fa', '.fasta'):
        # Just returns the first sequence in the file for now --
        # this may be undesirable as a general behavior
        return fasta.parseFasta(seq_fp)[0][1]
    else:
        raise IOError("unknown format")

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


def getGenomeLength(genome_fp):
    ''' Return the genome length from a fasta file containing a single
    header / sequence
    '''
    return len(getSeqFromFile(genome_fp))


def buildIndex(genome_fp, ref_genome_fp, genome_seq=None, ref_genome_seq=None, 
               fill_gaps=True, max_gap_width=300, smooth_edges=True, 
               smoothing_radius=20):
    ''' Builds and returns an index mapping from the genome in `genome_fp` to
    the genome in `ref_genome_fp`. If `fillGaps` is True, gaps up to
    `max_gap_width` will be filled if flanking indices can be found with
    commensurate map spacing (see `fillGaps` docstring). If `smooth_edges` is
    True, areas of high mismatch frequency with gap-penalty-related problems
    will be smoothed.
    '''
    # Generate absolute filepaths
    genome_fp = os.path.abspath(genome_fp)
    ref_genome_fp = os.path.abspath(ref_genome_fp)
    # Run mauve and parse output (all output files are cleaned up after run /
    # parsing)
    with TempDirs(1) as (temp_dir,):
        output_fp = os.path.join(temp_dir, 'mauveout.xmfa')
        runMauve([genome_fp, ref_genome_fp], flags={'--output': output_fp})
        sub_alignment_groups = parseXMFA(output_fp)
    genome_length = getGenomeLength(genome_fp)

    # LUT to map indices from the engineered genome to the reference genome
    idx_lut = np.empty(genome_length, dtype=IDX_ARRAY_DTYPE)
    # Initialize all values to -1 (un-mapped indices will be -1)
    idx_lut[:] = -1
    # Parse sub-alignments into index mapping
    for sub_alignment_group in sub_alignment_groups:
        genome_sa = indexutils.lookupSubAlignment(1, sub_alignment_group)
        ref_genome_sa = indexutils.lookupSubAlignment(2, sub_alignment_group)
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
        genome_seq = genome_seq or getSeqFromFile(genome_fp)
        ref_genome_seq = ref_genome_seq or getSeqFromFile(ref_genome_seq)
        indexutils.fixZeroIdx(idx_lut, genome_seq, ref_genome_seq)
        indexutils.fillGaps(idx_lut, max_gap_width=max_gap_width)
    if smooth_edges:
        indexutils.smoothEdges(idx_lut, smoothing_radius=smoothing_radius)
    return idx_lut


def test():
    import sys

    # OUTPUT_TO_FN = True

    # GENOME = 'test_mut_dup.fa'
    # REF_GENOME = 'baseseq.fa'
    # OUTPUT_FN = 'test_mut_dup.out'

    GENOME = os.path.join(LOCAL_DIR, 'mds42_recoded.fa')
    REF_GENOME = os.path.join(LOCAL_DIR, 'mds42_full.fa')
    # OUTPUT_FN = os.path.join(LOCAL_DIR, '54_failed_seg_test-testcleanup.out')

    # if OUTPUT_TO_FN:
    #     print_fd = open(OUTPUT_FN, 'w')
    # else:
    #     print_fd = sys.stdout

    # Cache the numpy array
    try:
        idx_lut = np.load(os.path.join(LOCAL_DIR, 'test_index_lut.npy'))
    except IOError:
        idx_lut = buildIndex(GENOME, REF_GENOME, fill_gaps=False,
                             smooth_edges=False)
        fillGaps(idx_lut)
        smoothEdges(idx_lut)
        np.save(os.path.join(LOCAL_DIR, 'test_index_lut'), idx_lut)

    genome = parseFasta(GENOME)[0][1]
    ref_genome = parseFasta(REF_GENOME)[0][1]

    return idx_lut, genome, ref_genome

if __name__ == '__main__':
    # import cProfile as profile
    # profile.run('test()')
    test()
