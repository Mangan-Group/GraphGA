import numpy as np 
from copy import deepcopy
from pymoo.indicators.hv import HV
from rankcrowding import RankAndCrowding
from GA import (
    sampling,
    crossover,
    mutate 
)
from diversity_metrics import (
    first_seen
)

def single_obj_GA(
        problem: object,
        seed: int
):
    """Sets up the GA for single objective
    optimization."""

    #set seed and create initial population
    np.random.seed(seed)
    population = sampling(
        problem.promo_node,
        problem.num_dict,
        problem.min_dose,
        problem.max_dose,
        problem.dose_interval,
        inhibitor=problem.inhibitor
    )
    num_circuits = len(population)
    # calculate objective for each circuit 
    # in initial population
    obj = np.asarray(
        [problem.func(g[0]) for g in population])

    # create list to store min obj function, 
    # circuits with min obj function for each
    # generation
    obj_min = np.zeros(problem.n_gen + 1)

    # index of obj list that contains min 
    # obj function
    ind_min = np.argmin(obj)
    obj_min[0] = obj[ind_min]
    circuit_min = []
    circuit_min.append(population[ind_min])

    for gen in range(problem.n_gen):
        # perform crossover to generate new
        # population (children) from parent 
        # circuits if randomly generated float  
        # is less than probability_crossover
        if np.random.uniform() < problem.prob_crossover:
            children = crossover(population, obj)
        else:
            children = deepcopy(population)

        # perform mutation on children if 
        # randomly generated float is less 
        # than probability_mutation (used 
        # in mutate function)
        mutate(
            problem, children, 
            problem.prob_mutation, 
            dose=problem.mutate_dose
        )

        # simulate topology and calculate obj
        # function for each circuit in children
        # and append to obj array
        obj_children = np.asarray(
            [problem.func(g[0]) for g in children])
        obj = np.append(obj, obj_children)

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

        # return index of minimum obj function
        # from obj
        ind_min = np.argmin(obj)

        # add min obj to obj_min and circuit with
        # min obj to circuit_min
        obj_min[gen + 1] = obj[ind_min]
        circuit_min.append(population[ind_min])
    
    # find in which gen the min obj first appeared
    gen_converged = first_seen(circuit_min)
    return [obj_min[-1], gen_converged]

def multi_obj_GA(
        problem: object, 
        seed: int
):
    """Sets up the GA for single objective
    optimization."""

    #set seed and create initial population
    np.random.seed(seed)
    population = sampling(
        problem.promo_node,
        problem.num_dict,
        problem.min_dose,
        problem.max_dose,
        problem.dose_interval,
        inhibitor=problem.inhibitor
    )
    num_circuits = len(population)
    # calculate objective for each circuit 
    # in initial population
    obj = np.asarray(
        [problem.func(g[0]) for g in population])
    
    # define reference point and class
    # instance of hypervolume calculator
    if "t_pulse" in '\t'.join(problem.obj_labels):
        if len(problem.obj_labels) == 3:
            ref_point = np.array([problem.max_time, 0, 0])
        else:
            ref_point = np.array([problem.max_time, 0])
    else:
        ref_point = np.array([0, 0])
    hv = HV(ref_point=ref_point)

    # store the progression of hypervolumes
    hypervolumes = []
    
    # create class instance of non-dominated
    # sorting class (to sort multi-objective
    # and determine pareto front)
    nds = RankAndCrowding()
    for gen in range(problem.n_gen):
        # sort objectives using non-dominated
        # sorting algorithm and return ranks 
        # for each circuit index in population
        _, rank_dict = nds.do(obj, num_circuits, return_rank=True)

        # perform crossover to generate new
        # population (children) from parent 
        # circuits if randomly generated float  
        # is less than probability_crossover
        if np.random.uniform() < problem.prob_crossover:
            children = crossover(population, obj, rank_dict)
        else:
            children = deepcopy(population)

        # perform mutation on children if 
        # randomly generated float is less 
        # than probability_mutation (used 
        # in mutate function)
        mutate(problem, children, 
                problem.prob_mutation, 
                dose=problem.mutate_dose
        )

        # simulate topology and calculate obj
        # function for each circuit in children
        # and append to all_obj list and obj
        # array
        obj_children = np.asarray(
            [problem.func(g[0]) for g in children])

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

        # append hypervolume to list
        hypervolumes.append(hv(obj))

    # find in which gen the max/final hv first appeared
    gen_converged = first_seen(hypervolumes)
    return [-hypervolumes[-1], gen_converged]