#!/usr/bin/env python
import re
import sys
from collections import defaultdict
import matplotlib.pyplot as plt
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

def _opt_crosses(max_k, current_k, pairs, informative, opt_loci, 
                 added_loci, opt_pairs):
    '''
    recursive method to find optimal sets of pairs such that each added pair
    adds the most new loci to the overall set of informative ones
    '''
    current_k += 1
    best_diff = set()
    next_pair = ('','')

    if current_k == max_k or len(pairs) == 0:
        #finished!
        return opt_loci, added_loci, opt_pairs

    elif current_k == 1:
        #to start just choose the set with the most informative loci
        for p in pairs:
            if len(best_diff) < len(informative[p]):
                best_diff = informative[p]
                next_pair = p
        opt_loci.append(informative[next_pair])
        added_loci.append(len(informative[next_pair]))
    else:
        #iteratively find the next best pair to add, based on its added
        #contribution of loci
        best_diff = set()
        for p in pairs:
            tmp_diff = informative[p].difference(opt_loci[-1])
            if len(best_diff) < len(tmp_diff):
                best_diff = tmp_diff
                next_pair = p
        opt_loci.append(opt_loci[-1].union(informative[next_pair]))
        added_loci.append(len(best_diff))

    pairs.remove(next_pair)
    opt_pairs.append(next_pair)
    #recursion!
    return _opt_crosses(max_k, current_k, pairs, informative, opt_loci, 
                        added_loci, opt_pairs)

def optimize_crosses(max_k, pairs, informative):
    '''
    wrapper for _opt_crosses to hide the recursive variable passing
    '''
    pairs_copy = set(pairs)
    return _opt_crosses(max_k, 0, pairs_copy, informative, [], [], [])

def main():
    #business time.

    #read in the possible mating pairs to test as tuples
    pairs = []
    with open(sys.argv[2]) as fi:
        for l in fi:
            pairs.append(tuple(l.split()))

    #get all locus names, loci informative for each cross
    informative, loci = read_genotypes(sys.argv[1], pairs)

    shared_loci, added_loci, best_pairs = optimize_crosses(6, pairs, 
                                                           informative)

    best_pairs_labels=["{0} x {1}".format(*x) for x in best_pairs]

    fig = plt.figure(figsize = (8,8))
    ax1 = plt.subplot2grid((3,3), (0,0), colspan = 3, rowspan = 2)
    ax1.set_title('Shared Loci in Top Three Pairs', fontsize = 26)

    venn3_unweighted(shared_loci[:3], set_labels = best_pairs_labels[:3], 
                     ax=ax1)
    ax2 = plt.subplot2grid((3,3), (2,0), colspan=3)
    ax2.bar(range(len(added_loci)), added_loci)
    ax2.set_ylabel('Number of Loci')
    ax2.set_title('Added Informative Loci Per Pair', fontsize = 26)
    ax2.set_xticklabels(best_pairs_labels, rotation = 17 )

    fig.subplots_adjust(hspace=.44)
    fig.savefig('summmary.pdf')
    plt.show()

if __name__ == '__main__':
    main()
