#!/usr/bin/env python
import re
import sys
from collections import defaultdict
import matplotlib.pyplot as plt
from itertools import combinations
from matplotlib_venn import venn3_unweighted 

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

def union_len(things):
    uni = set()
    for x in things:
        uni.unioon(x)
    return len(uni)


def global_optimize_crosses(max_k, pairs, informative):
    '''
    global optimization, run through each possible combination of each number of crosses

    '''
    best_sets = []
    for k in range(2,max_k+1):
        best_num_loci = 0
        for st in  combinations(pairs, k):
            current_num_loci = union_len(st)
            if best_num_loci < current_num_loci:
                best_set = st
                best_num_loci = current_num_loci
        best_sets.append(best_set)

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
    best_sets = global_optimize_crosses(6, pairs, informative)

    #labels for pairs
    best_pairs_labels=["{0} x {1}".format(*x) for x in best_sets]
    '''
    #plot results
    fig = plt.figure(figsize = (8,8))
    ax1 = plt.subplot2grid((3,3), (0,0), colspan = 3, rowspan = 2)
    ax1.set_title('Shared Loci in Top Three Pairs', fontsize = 26)
    #venn diagram on top
    venn3_unweighted([informative[x] for x in best_pairs[:3]], 
                     set_labels = best_pairs_labels[:3], ax=ax1)
    #bar chart on bottom
    ax2 = plt.subplot2grid((3,3), (2,0), colspan=3)
    ax2.bar(range(len(added_loci[:8])), added_loci[:8])
    ax2.set_ylabel('Number of Loci')
    ax2.set_title('Added Informative Loci Per Pair', fontsize = 26)
    ax2.set_xticklabels(best_pairs_labels[:8], rotation = 17 )
    #add a bit of extra space
    fig.subplots_adjust(hspace=.44)
    fig.savefig('summmary.pdf')
    plt.show()
    '''
