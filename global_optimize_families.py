#!/usr/bin/env python
import re
import sys
from collections import defaultdict
from itertools import combinations

def check_locus(m, f):
    '''
    checks the genotypes of a cross m x f to see if there is a 
    homozygote and a heterozygote present -> most informative if we are 
    genotyping F1 offspring to check for mendelian inheritance
    '''
    if ((m == 'AA' and f == 'AB') or (m == 'BB' and f == 'AB') or
        (m == 'AB' and f == 'AA') or (m == 'AB' and f == 'BB')):
        return(True)
    else:
        return(False)

def read_genotypes(file_name, pairs):
    '''
    reads in a genotyping file from the affy genotyping suite, for each
    line (locus) checking each pair to see which pairs have informative 
    genotypes.  returns a list of the locus names and a dictionary of 
    tupple -> set()
    (m, f) -> loci
    '''
    informative = defaultdict(lambda: set())
    loci = []
    with open(file_name) as fi:
        for l in fi:
            if l.startswith('#'):
                #comments
                pass
            elif l.startswith('Probe Set ID'):
                #header line
                header = l.rstrip().split('\t')
                sample_names = ([re.sub('.AxiomGT1.chp Call Codes', '', x) 
                                 for x in header[1:-3]])
            else:
                #locus line
                tmp = l.rstrip().split('\t')
                locus_name = tmp[0]
                loci.append(locus_name)
                locus = dict(zip(sample_names, tmp[1:-3]))
                #decide if locus is informative in each pair 
                for pair in pairs:
                    if check_locus(locus[pair[0]], locus[pair[1]]):
                        informative[pair].add(locus_name)
    return informative, loci

def union_len(pair_combo, informative):
    '''
    returns the length of the union of the informative loci of 
    all the crosses passed in
    '''

    uni = set()
    for p in pair_combo:
        uni.update(informative[p])
    return len(uni)


def global_optimize_crosses(max_k, pairs, informative):
    '''
    global optimization, run through each possible combination of each number of crosses

    '''
    best_sets = []
    for k in range(2,max_k+1):
        best_num_loci = 0
        for pair_combo in combinations(pairs, k):
            current_num_loci = union_len(pair_combo, informative)
            if best_num_loci < current_num_loci:
                best_combo = pair_combo
                best_num_loci = current_num_loci
        best_sets.append(best_combo)

    return best_sets


if __name__ == '__main__':
    #business time.

    #read in the possible mating pairs to test as tuples
    pairs = []
    with open(sys.argv[2]) as fi:
        for l in fi:
            pairs.append(tuple(l.split()))

    #get all locus names, loci informative for each cross
    informative, loci = read_genotypes(sys.argv[1], pairs)

    #get the loci captured at each k, the number of added loci,
    #and the pairs decending order of how many more loci they add
    best_combos = global_optimize_crosses(6, pairs, informative)

    for combo in best_combos:
        print(combo, ":", union_len(combo, informative))
