import bitarray
import numpy as np

cimport cython
cimport numpy as np


# Non-mask array data type
# 32-bit supports genomes up to 4 billion bases in length
IDX_ARRAY_DTYPE = np.int32
# Mask array data type (should only have values of 0 + 1 so unit8 is good)
MASK_ARRAY_DTYPE = np.uint8

ctypedef np.int32_t IDX_ARRAY_DTYPE_t


def lookupSubAlignment(seq_num, sub_alignment_group):
    sa = None
    try:
        sa = list(filter(lambda a: a.seq_num ==
                  seq_num, sub_alignment_group.alignments))[0]
    except IndexError:
        pass
    return sa


def fixZeroIdx(np.ndarray[IDX_ARRAY_DTYPE_t, ndim=1] idx_lut, genome, ref_genome):
    """ only map the zero index on an exact match
    it's problematic
    """
    cdef int mapped_idx, next_idx, ref_idx
    mapped_idx = idx_lut[0]
    if mapped_idx == -1:
        next_idx = 1
        while (idx_lut[next_idx] == -1) or (next_idx < 10):
            next_idx +=1
        if idx_lut[next_idx] == -1:
            return
        ref_idx = idx_lut[next_idx] - next_idx
        if ref_idx > -1:
            if genome[0] == ref_genome[ref_idx]:
                idx_lut[0] = ref_idx


@cython.wraparound(False)
@cython.boundscheck(False)
def fillGaps(np.ndarray[IDX_ARRAY_DTYPE_t, ndim=1] idx_lut, int max_gap_width=300):
    ''' Fill gaps in the alignment that are up to `max_gap_width` in size. A
    gap is only filled if flanking indices can be found such that:

        upper_idx - lower_idx == upper_mapped_idx - lower_mapped_idx
    '''
    cdef int idx, upper_idx, lower_idx, mapped_idx, upper_mapped_idx, lower_mapped_idx
    # for idx, mapped_idx in enumerate(idx_lut[1:], start=1):
    arange = np.arange
    for idx in range(1, len(idx_lut)):
        mapped_idx = idx_lut[idx]
        if mapped_idx == -1:
            lower_mapped_idx = idx_lut[idx-1]
            if lower_mapped_idx != -1:
                lower_idx = idx - 1
                # for upper_idx, upper_mapped_idx in enumerate(idx_lut[idx+1:],
                #                                              start=idx+1):
                for upper_idx in range(idx+1, len(idx_lut)):
                    upper_mapped_idx = idx_lut[upper_idx]
                    if (upper_idx-idx > max_gap_width):
                        break
                    if ((upper_mapped_idx - lower_mapped_idx) ==
                            (upper_idx - lower_idx)):
                        idx_lut[lower_idx:upper_idx+1] = arange(
                            lower_mapped_idx, upper_mapped_idx+1,
                            dtype=IDX_ARRAY_DTYPE)
                        break


@cython.wraparound(False)
@cython.boundscheck(False)
def smoothEdges(np.ndarray[IDX_ARRAY_DTYPE_t, ndim=1] idx_lut, int smoothing_radius=20):
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
    # for idx, mapped_idx in enumerate(idx_lut[1:], start=1):
    cdef int idx, upper_idx, mapped_idx, lower_mapped_idx, lower_idx, upper_mapped_idx
    for idx in range(1, len(idx_lut)):
        mapped_idx = idx_lut[idx]
        if mapped_idx != idx_lut[idx-1] + 1:
            lower_mapped_idx = idx_lut[idx-1]
            if lower_mapped_idx != -1:
                lower_idx = idx - 1
                for upper_idx in range(idx+1, len(idx_lut)):
                    upper_mapped_idx = idx_lut[upper_idx]
                    if (upper_idx-idx > smoothing_radius):
                        break
                    if ((upper_mapped_idx - lower_mapped_idx) ==
                            (upper_idx - lower_idx)):
                        idx_lut[lower_idx:upper_idx+1] = np.arange(
                            lower_mapped_idx, upper_mapped_idx+1,
                            dtype=IDX_ARRAY_DTYPE)
                        break


def findGaps(np.ndarray[IDX_ARRAY_DTYPE_t, ndim=1] idx_lut):
    ''' Find and print gaps in the idx_lut
    '''
    cdef int idx, upper_idx, lower_mapped_idx, upper_mapped_idx, gap_size, mapped_idx
    gap_arr = np.zeros(idx_lut.shape[0], dtype=np.int32)
    try:
        for idx in range(1, len(idx_lut)):
            mapped_idx  = idx_lut[idx]
            lower_idx = idx - 1
            lower_mapped_idx = idx_lut[idx-1]
            if mapped_idx == -1 and lower_mapped_idx != -1:
                success = False
                for upper_idx in range(idx+1, len(idx_lut)):
                    upper_mapped_idx = idx_lut[upper_idx]
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


def findEdges(np.ndarray[IDX_ARRAY_DTYPE_t, ndim=1] idx_lut):
    ''' Find and print edges in the idx_lut
    return a 1 on the 5 prime most index of a segment
    uses slicing
    '''
    edge_arr = np.zeros(idx_lut.shape[0], dtype=MASK_ARRAY_DTYPE)
    # edge_arr = bitarray.bitarray(idx_lut.shape[0])
    # edge_arr.setall(0)
    edge_arr_view = edge_arr[1:]    # slice not a copy
    delta_map_idxs = idx_lut[1:] - idx_lut[:-1]
    edge_arr_view[delta_map_idxs != 1] = 1
    edge_arr[0] = 1     # 0 index is always an edge
    edge_bitarray = bitarray.bitarray()
    edge_bitarray.pack(edge_arr.tostring())
    return edge_bitarray


def findMismatches(np.ndarray[IDX_ARRAY_DTYPE_t, ndim=1] idx_lut, genome, ref_genome):
    cdef int idx, mapped_idx
    cdef int len_genome = len(genome)
    cdef char base

    cdef char* genome_string = genome
    cdef char* ref_genome_string = ref_genome

    mismatch_arr = bitarray.bitarray(idx_lut.shape[0])
    mismatch_arr.setall(0)
    for idx in range(len_genome):
        base = genome_string[idx]
        mapped_idx = idx_lut[idx]
        if mapped_idx != -1:
            if base != ref_genome_string[mapped_idx]:
                mismatch_arr[idx] = 1
    return mismatch_arr


def findDuplicateMappings(np.ndarray[IDX_ARRAY_DTYPE_t, ndim=1] idx_lut):
    ''' Find duplicate mappings in the idx_lut. If everything is working
    properly there should be 1:1, unambiguous mappings.
    '''
    cdef np.ndarray[IDX_ARRAY_DTYPE_t, ndim=1] arr_cpy = np.clip(idx_lut, 0, (2**31)-1) # clip to 0, 2**31-1
    cdef np.ndarray[np.int_t, ndim=1] mapped_idx_counts = np.bincount(arr_cpy)
    # -1 as 0 will show up here as we are clipping unmapped basses to 0
    cdef int num_duplicate_mappings = np.sum(mapped_idx_counts > 1) - 1
    return num_duplicate_mappings
