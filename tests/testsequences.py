'''
Generate some test genomes for Mauve alignments. All of the genomes will be
derived from the 'mds42_full.fa' genome.

This script's output is non-deterministic.

'''

import random

from util import parseFasta

def moveSection(genome_str, start_idx, end_idx, insertion_point_idx):
    section = genome_str[start_idx:end_idx]
    indexing_offset = 0
    if insertion_point_idx > start_idx:
        indexing_offset = -len(section)
    new_genome_str = genome_str[:start_idx] + genome_str[end_idx:]
    new_genome_str = (new_genome_str[:insertion_point_idx + indexing_offset] +
                      section +
                      new_genome_str[insertion_point_idx + indexing_offset:])
    return new_genome_str


def duplicateSection(genome_str, start_idx, end_idx,
                     insertion_point_idx=None, add_point_mutations=False,
                     mutation_prob=0.01):
    section = genome_str[start_idx:end_idx]
    insertion_point_idx = insertion_point_idx or end_idx
    if add_point_mutations:
        addPointMutations(section, mutation_prob=0.01)
    new_genome_str = (genome_str[:insertion_point_idx] +
                      section +
                      genome_str[insertion_point_idx:])
    return new_genome_str

_RAND_BASE_LUT = {
    'A': 'TGC',
    'T': 'AGC',
    'G': 'ATC',
    'C': 'ATG'
}

def addPointMutations(genome_str, start_idx=None, end_idx=None,
                      mutation_prob=0.01):
    start_idx = start_idx or 0
    end_idx = end_idx or len(genome_str)
    genome_list = list(genome_str)
    for base_idx in range(start_idx, end_idx):
        if random.randint(0,1) < mutation_prob:
            genome_list[base_idx] = random.choice(
                _RAND_BASE_LUT[genome_list[base_idx]])
    return ''.join(genome_list)


if __name__ == '__main__':

    # Read in base sequence (length = 10kb)
    BASE_SEQUENCE = parseFasta('baseseq.fa')[0][1]
    FASTA_BASE_STRING = '>baseseq_{}\n'

    # Sequence with a single rearragement
    swapped_seq = moveSection(BASE_SEQUENCE, 99, 199, 499)
    with open('test_1swap.fa', 'w') as fd:
        fd.write(FASTA_BASE_STRING.format('1swap (bases 100-200 moved to 500)'))
        fd.write(swapped_seq)
        fd.write('\n')

    # Sequence with a single duplication
    dup_seq = duplicateSection(BASE_SEQUENCE, 99, 199)
    with open('test_1dup.fa', 'w') as fd:
        fd.write(FASTA_BASE_STRING.format('1dup (bases 100-200 duplicated)'))
        fd.write(dup_seq)
        fd.write('\n')

    # Sequence with a bunch of mutations
    mut_seq = addPointMutations(BASE_SEQUENCE, mutation_prob=0.1)
    with open('test_muts.fa', 'w') as fd:
        fd.write(FASTA_BASE_STRING.format('muts (10% probability of mutations)'))
        fd.write(mut_seq)
        fd.write('\n')

