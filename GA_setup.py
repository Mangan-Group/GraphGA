import numpy as np 
from copy import deepcopy
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

# set up GA for single objective
def single_obj_GA(
        problem: object,
        n_gen: int,
        population: np.ndarray,
        num_circuits: int, 
        obj: np.ndarray,
        probability_crossover: float,
        probability_mutation: float,
        mutate_dose: bool,
        metrics: bool =False
):
    # create list to store min obj function, 
    # circuits with min obj function, and all
    # obj functions for initial population 
    # and all generations 
    obj_min = np.zeros(n_gen + 1)


    # all_obj = []
    # all_obj.append(obj)
    # index of list that contains min obj function
    ind_min = np.argmin(obj)
    obj_min[0] = obj[ind_min]
    circuit_min = []
    circuit_min.append(population[ind_min])
    # save all circuits from intial population
    # and top num_circuits from each generation
    top_circuits_each_gen = []
    top_circuits_each_gen.append(population)
    top_obj_each_gen = []
    top_obj_each_gen.append(obj)

    geno = None
    pheno = None
    # if storing metrics, create lists for those
    # metrics and store initial population value
    if metrics:
        geno = np.zeros(n_gen+1)
        geno[0] = geno_diversity(population)

        pheno = np.zeros_like(geno)
        pheno[0] = pheno_diversity(obj)

    for gen in range(n_gen):
        # perform crossover to generate new
        # population (children) from parent 
        # circuits if randomly generated float  
        # is less than probability_crossover
        if np.random.uniform() < probability_crossover:
            children = crossover(population, obj)
        else:
            children = deepcopy(population)

        # perform mutation on children if 
        # randomly generated float is less 
        # than probability_mutation (used 
        # in mutate function)
        mutate(problem, children, 
               probability_mutation, 
               dose=mutate_dose
        )
        ### add all children here###

        # simulate topology and calculate obj
        # function for each circuit in children
        # and append to obj array
        obj_children = np.asarray([problem.func(g[0]) for g in children])
        # all_obj.append(obj_children)
        obj = np.append(obj, obj_children)
        ### add all obj here ###

        # add children to population array
        population = np.vstack((population, children))
        # return array of indices that would sort
        # obj
        S = np.lexsort([obj])
        # select top num_circuits obj from obj array 
        # and top num_circuits from population
        # (initial population + children of each gen)
        obj = obj[S[:num_circuits]]
        population = population[S[:num_circuits], :]
        top_circuits_each_gen.append(population)
        top_obj_each_gen.append(obj)

        # return index of minimum obj function
        # from obj
        ind_min = np.argmin(obj)

        # add min obj to obj_min and circuit with
        # min obj to circuit_min
        obj_min[gen + 1] = obj[ind_min]
        circuit_min.append(population[ind_min])

        # calculate metrics for population
        if metrics:
             geno[gen+1] = geno_diversity(population)
             pheno[gen+1] = pheno_diversity(obj)

    # reshape all_obj to be 1 column array
    # all_obj = np.asarray(all_obj).reshape(num_circuits*(1 + n_gen), 1)
    top_circuits_each_gen = np.asarray(
        top_circuits_each_gen
        ).reshape(num_circuits*(1 + n_gen), 1)
    top_obj_each_gen = np.asarray(
        top_obj_each_gen
        ).reshape(num_circuits*(1 + n_gen), 1)

    return  (obj_min, circuit_min, top_circuits_each_gen,
             top_obj_each_gen ,geno, pheno) #all_obj


def multi_obj_GA(
        problem: object, 
        n_gen: int,
        population: np.ndarray,
        num_circuits: int,
        obj: np.ndarray,
        probability_crossover: float,
        probability_mutation: float,
        mutate_dose: bool,
):
    # create list to store all obj 
    # functions for initial population and 
    # all generations 
    all_obj = []
    all_obj.append(obj)
    # save all circuits from intial population
    # and top num_circuits from each generation
    top_circuits_each_gen = []
    top_circuits_each_gen.append(population)
    top_obj_each_gen = []
    top_obj_each_gen.append(obj)

    # create class instance of non-dominated
    # sorting class (to sort multi-objective
    # and determine pareto front)
    nds = RankAndCrowding()
    for gen in range(n_gen):
        # sort objectives using non-dominated
        # sorting algorithm and return ranks 
        # for each circuit index in population
        _, rank_dict = nds.do(obj, num_circuits, return_rank=True)

        # perform crossover to generate new
        # population (children) from parent 
        # circuits if randomly generated float  
        # is less than probability_crossover
        if np.random.uniform() < probability_crossover:
            children = crossover(population, obj, rank_dict)
        else:
            children = deepcopy(population)

        # perform mutation on children if 
        # randomly generated float is less 
        # than probability_mutation (used 
        # in mutate function)
        mutate(problem, children, 
                probability_mutation, 
                dose=mutate_dose
        )

        # simulate topology and calculate obj
        # function for each circuit in children
        # and append to all_obj list and obj
        # array
        obj_children = np.asarray([problem.func(g[0]) for g in children])
        all_obj.append(obj_children)

        obj = np.vstack((obj, obj_children))
        # add children to population array
        population = np.vstack((population, children))

        # sort objectives using non-dominated
        # sorting algorithm and return indices
        # of num_circuits highest rank
        S = nds.do(obj, num_circuits)

        # select top num_circuits obj from obj array 
        # and top num_circuits from population
        # (initial population + children of each gen)
        obj = obj[S]
        population = population[S, :]
        top_circuits_each_gen.append(population)
        top_obj_each_gen.append(obj)


    fronts = NonDominatedSorting().do(obj)
    # reshape all_obj and top_circuits_each_gen
    # to be 1 column arrays
    all_obj = np.asarray(
        all_obj).reshape(num_circuits*(1 + n_gen), 2)
    top_circuits_each_gen = np.asarray(
        top_circuits_each_gen
        ).reshape(num_circuits*(1 + n_gen), 1)
    top_obj_each_gen = np.asarray(
        top_obj_each_gen
        ).reshape(num_circuits*(1 + n_gen), 1)

    return (fronts, obj, all_obj, top_circuits_each_gen,
            top_obj_each_gen)