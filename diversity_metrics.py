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

