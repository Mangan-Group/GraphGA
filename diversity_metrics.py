import numpy as np
from itertools import chain
from GA import compare_circuit, is_equal

#genotype (topology)
def geno_diversity(X):
    dupes_dict = compare(X)
    return len(dupes_dict)

def compare(X):
    dupes_dict = dict()
    for i, g in enumerate(X):
        if i not in list(chain.from_iterable(dupes_dict.values())):
            dupes_dict.update({i: [j for j in range(i+1, len(X)) if is_equal(g, X[j])]})
    return dupes_dict

#phenotype (objective)
def pheno_diversity(F, delta=0.5):
    n, bins = grouping(F, delta)
    return len(set(n))

def grouping(F, delta=0.5):
    bins = np.arange(np.min(F), np.max(F)+delta, delta)
    n = np.full_like(F, 0)
    for i in range(len(F)):
        n[i] = np.argwhere(F[i] >= bins).flatten()[-1]
    return n, bins

def first_seen(progression):

    looking = True
    gen_num = 0

    # for each value in reversed list
    # of minimum objective function for
    # each generation
    for gen in reversed(progression):
        # if the value does not equal the
        # final value in reversed list
        # return length of progression
        # (number of generations) - gen_num
        # (number of generations to end)
        if progression[-1] != gen:
            return len(progression) - gen_num
        gen_num += 1

    return 0

