import numpy as np
import pickle
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator
from rankcrowding import RankAndCrowding
from diversity_metrics import (
    geno_diversity,
    pheno_diversity
)
from GA import (
    sampling,
    crossover,
    mutate 
)
# from define_circuit import Topo

# To-Do:
# add metrics for single obj opt
# separate file with function for single_obj_GA and multi_obj_GA
# save files: add results path attribute to testcases


#amp: 10 1 part, 20 2 part

def set_testcase(
    case: object,
    settings: tuple,
):
    [promo_node, 
    # max_part, #no need for this anymore?
    min_dose, #pack into [min, max, interval]
    max_dose, 
    dose_interval, 
    inhibitor,
    DsRed_inhibitor,
    num_dict, # defines max part
    n_gen, 
    pop,
    num_processes] = settings

    testcase = case(
        promo_node,
        # max_part,
        min_dose, 
        max_dose, 
        dose_interval, 
        inhibitor,
        DsRed_inhibitor,
        num_dict, 
        n_gen, 
        pop,
        num_processes, 
    )

    return testcase

def run(
        testcase: object,
        metrics: bool =False
    ):
    '''Run the genetic algorithm for a given test case class instance'''

    population = sampling(
        testcase.promo_node,
        testcase.num_dict,
        testcase.min_dose,
        testcase.max_dose,
        testcase.dose_interval,
        testcase.inhibitor
    )
    num_circuits = len(population)
    # raise an exception if num_circuits is odd
    if num_circuits % 2 != 0:
        raise Exception("Population size must be an even number")
    # set # generations according to testcase attribute; store as array
    n_gen = testcase.n_gen
    generations = np.arange(n_gen + 1)

    # calculate objective for each circuit in initial population
    obj = np.asarray([testcase.func(g[0]) for g in population])
    
    if isinstance (obj[0], np.ndarray):
        print('objective is an array- using multi-objective optimization')
        all_obj = []
        all_obj.append(obj)
        nds = RankAndCrowding()
        for gen in range(n_gen):
            _, rank_dict = nds.do(obj, num_circuits, return_rank=True)
            children = crossover(population, obj, rank_dict)
            mutate(testcase, children, 1.)
            obj_children = np.asarray([testcase.func(g[0]) for g in children])
            all_obj.append(obj_children)

            obj = np.vstack((obj, obj_children))
            population = np.vstack((population, children))

            S = nds.do(obj, num_circuits)
            obj = obj[S]
            population = population[S, :]

        fronts = NonDominatedSorting().do(obj)
        all_obj = np.asarray(all_obj).reshape(num_circuits*(1 + n_gen), 2)

        return fronts, all_obj

    else:
        print('objective is not an array- using single objective optimization')
        # create array to store min obj function value for initial
        # population and all generations
        obj_min = np.zeros(n_gen + 1)

        # create list to store circuits with min obj function for 
        # initial population and all generations
        circuit_min = []
        ind_min = np.argmin(obj)
        obj_min[0] = obj[ind_min]
        circuit_min.append(population[ind_min])

        for gen in range(n_gen):
            children = crossover(population, obj)
            mutate(testcase, children, 1.)
            obj_children = np.asarray([testcase.func(g[0]) for g in children])
            obj = np.append(obj, obj_children)
            population = np.vstack((population, children))
            S = np.lexsort([obj])
            obj = obj[S[:num_circuits]]
            population = population[S[:num_circuits], :]

            ind_min = np.argmin(obj)

            obj_min[gen + 1] = obj[ind_min]
            circuit_min.append(population[ind_min])

        return obj, obj_min, circuit_min, generations

def run_combinitorial(
        testcase: object,
        topo_path: str
):
    with open(topo_path, "rb") as fid:
        topologies = pickle.load(fid)
    
    objectives = [testcase.func(g) for g in topologies]
    
    return objectives

testcase_amp = set_testcase(
    Amplifier,
    ['P1',
    # 1,
    150,
    150,
    5,
    False,
    False,
    {1: 10, 2: 30},
    30,
    False,
    None,
])

testcase_sc = set_testcase(
    SignalConditioner,
    ['P1',
    # 1,
    150,
    150,
    5,
    False,
    False,
    {1: 10, 2: 30},
    30,
    False,
    None,
])

run(testcase_sc)

# obj_min, topology_min, gens = run(testcase_amp)

# edge_list = [('P1', 'Z1'), ('Z1', 'Z6'), ('Z6', 'Rep')]
# dose_list = {'Z1': 75, 'Z6': 10, 'Rep': 9}
# promo_node = 'P1'

# obj_val = testcase_amp.func(Topo(edge_list, dose_list, promo_node))
# print(isinstance(obj_val, list))

# print(obj_min[-1], topology_min[-1][0].edge_list)