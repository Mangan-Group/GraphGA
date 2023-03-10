import numpy as np
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from rankcrowding import RankAndCrowding
from GA import (
    crossover,
    mutate 
)
from diversity_metrics import (
    geno_diversity,
    pheno_diversity
)

def single_obj_GA(
        problem: object,
        n_gen: int,
        population: np.ndarray,
        num_circuits: int, 
        obj: np.ndarray,
        metrics: bool =False
):

    obj_min = np.zeros(n_gen + 1)

    # create list to store circuits with min obj function for 
    # initial population and all generations
    all_obj = []
    all_obj.append(obj)
    ind_min = np.argmin(obj)
    obj_min[0] = obj[ind_min]
    circuit_min = []
    circuit_min.append(population[ind_min])

    geno = None
    pheno = None
    if metrics:
        geno = np.zeros(n_gen+1)
        geno[0] = geno_diversity(population)

        pheno = np.zeros_like(geno)
        pheno[0] = pheno_diversity(obj)

    for gen in range(n_gen):
        children = crossover(population, obj)
        mutate(problem, children, 1.)
        obj_children = np.asarray([problem.func(g[0]) for g in children])
        all_obj.append(obj_children)
        obj = np.append(obj, obj_children)
        population = np.vstack((population, children))
        S = np.lexsort([obj])
        obj = obj[S[:num_circuits]]
        population = population[S[:num_circuits], :]

        ind_min = np.argmin(obj)

        obj_min[gen + 1] = obj[ind_min]
        circuit_min.append(population[ind_min])

        if metrics:
             geno[gen+1] = geno_diversity(population)
             pheno[gen+1] = pheno_diversity(obj)

    all_obj = np.asarray(all_obj).reshape(num_circuits*(1 + n_gen), 1)

    return all_obj, obj_min, circuit_min, geno, pheno


def multi_obj_GA(
        problem: object, 
        n_gen: int,
        population: np.ndarray,
        num_circuits: int,
        obj: np.ndarray
):
        
    all_obj = []
    all_obj.append(obj)
    nds = RankAndCrowding()
    for gen in range(n_gen):
        _, rank_dict = nds.do(obj, num_circuits, return_rank=True)
        children = crossover(population, obj, rank_dict)
        mutate(problem, children, 1.)
        obj_children = np.asarray([problem.func(g[0]) for g in children])
        all_obj.append(obj_children)

        obj = np.vstack((obj, obj_children))
        population = np.vstack((population, children))

        S = nds.do(obj, num_circuits)
        obj = obj[S]
        population = population[S, :]

    fronts = NonDominatedSorting().do(obj)
    all_obj = np.asarray(all_obj).reshape(num_circuits*(1 + n_gen), 2)

    return fronts, all_obj