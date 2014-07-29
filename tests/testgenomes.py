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
    print len(new_genome_str)
    return new_genome_str


def duplicateSection(genome_str, start_idx, end_idx, 
                     insertion_point_idx=None, add_point_mutations=False,
                     mutation_prob=0.01):
    section = genome_str[start_idx:end_idx]
    insertion_point_idx = insertion_point_idx or end_idx
    if add_point_mutations:
        addPointMutations(section, mutation_prob=0.01)
    new_genome_str = genome_str[:start_idx] + genome_str[end_idx:]
    new_genome_str = (new_genome_str[:insertion_point_idx] + 
                      section + 
                      new_genome_str[insertion_point_idx:])  
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

    # Read in master genome (length = 3,976,195)
    MASTER_GENOME = parseFasta('mds42_full.fa')[0][1]
    FASTA_BASE_STRING = '>mds42_full_{}\n'

    # Genome with two major section rearrangements
    swapped_genome1 = moveSection(MASTER_GENOME, 15000, 20000, 180000)
    swapped_genome1 = moveSection(swapped_genome1, 2500000, 2550000, 1500)
    with open('mds42_full_2swaps.fa', 'w') as fd:
        fd.write(FASTA_BASE_STRING.format('2swaps'))
        fd.write(swapped_genome1)
        fd.write('\n')

    # Genome with two duplications, one with mutations and one without
    dupl_genome1 = duplicateSection(MASTER_GENOME, 15000, 20000)
    dupl_genome1 = duplicateSection(dupl_genome1, 850000, 900000,
                                      add_point_mutations=True,
                                      mutation_prob=0.3)
    with open('mds42_full_2duplications.fa', 'w') as fd:
        fd.write(FASTA_BASE_STRING.format('2duplications'))
        fd.write(dupl_genome1)
        fd.write('\n')    

    # Genome with 3 swaps and 3 duplications
    swapdupl_genome = moveSection(MASTER_GENOME, 15000, 20000, 180000)
    swapdupl_genome = moveSection(swapdupl_genome, 450000, 454000, 500000)
    swapdupl_genome = moveSection(swapdupl_genome, 100, 4200, 1200000)
    swapdupl_genome = duplicateSection(swapdupl_genome, 2500000, 8000,
                                      add_point_mutations=True,
                                      mutation_prob=0.3)
    swapdupl_genome = duplicateSection(swapdupl_genome, 2500000, 8000,
                                      add_point_mutations=True,
                                      mutation_prob=0.3)
    swapdupl_genome = duplicateSection(swapdupl_genome, 2900000, 2902000)    

    with open('mds42_full_3swaps3duplications.fa', 'w') as fd:
        fd.write(FASTA_BASE_STRING.format('3swaps3duplications'))
        fd.write(swapdupl_genome)
        fd.write('\n')    

    # Same as above but with added mutations throughout
    mut_swapdupl_genome = addPointMutations(swapdupl_genome, 
                                            mutation_prob=0.08)
    with open('mds42_full_3swaps3duplications_muts.fa', 'w') as fd:
        fd.write(FASTA_BASE_STRING.format('3swaps3duplications_muts'))
        fd.write(mut_swapdupl_genome)
        fd.write('\n')     
