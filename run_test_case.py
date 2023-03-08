import numpy as np
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator
from GA import (
    sampling,
    crossover,
    mutate 
)

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
        num_dict, 
        n_gen, 
        pop,
        num_processes, 
    )

    return testcase

def run(testcase: object):
    population = sampling(
        testcase.promo_node,
        testcase.num_dict,
        testcase.min_dose,
        testcase.max_dose,
        testcase.dose_interval
    )

    num_circuits = len(population)

    # raise an exception if num_circuits is odd
    if num_circuits % 2 != 0:
        raise Exception("Population size must be an even number")

    # set # generations according to testcase attribute
    n_gen = testcase.n_gen

    # calculate objective for each circuit in initial population
    obj = np.asarray([testcase.func(g[0]) for g in population])

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
        # if np.random.uniform() < prob:
        children = crossover(population, obj)
        # print(children)
        # print(np.shape(children))

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

    generations = np.arange(n_gen + 1)

    return obj_min, circuit_min , generations

def run_combinitorial(
        testcase: object,
        topo_path: str
):
    # topologies = 
    
    return

testcase_amp = set_testcase(
    Amplifier,
    ['P1',
    # 1,
    150,
    150,
    5,
    False,
    {1: 10, 2: 30},
    30,
    False,
    None,
])

obj_min, topology_min, gens = run(testcase_amp)

print(obj_min[-1], topology_min[-1][0].edge_list)