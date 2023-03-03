import numpy as np
from amplifier_problem import Amplifier
from GA import (
    sampling,
    crossover,
    mutate 
)

def set_testcase(
    case: object,
    settings: tuple,
):
    [promo_node, 
    max_part, #no need for this anymore?
    min_dose, #pack into [min, max, interval]
    max_dose, 
    dose_interval, 
    inhibitor,
    num_dict, # defines max part
    n_gen, 
    pop,
    num_processes, 
    combinatorial] = settings

    testcase = case(
        promo_node,
        max_part,
        min_dose, 
        max_dose, 
        dose_interval, 
        inhibitor,
        num_dict, 
        n_gen, 
        pop,
        num_processes, 
        combinatorial
    )

    return testcase

testcase = set_testcase(
    Amplifier,
    ['P1',
    1,
    150,
    150,
    5,
    False,
    {1: 3, 2: 0},
    3,
    False,
    None,
    False
])

def run(testcase: object):
    population = sampling(
        testcase.promo_node,
        testcase.num_dict,
        testcase.min_dose,
        testcase.max_dose,
        testcase.dose_interval
    )

    num_circuits = len(population)
    n_gen = testcase.n_gen
    obj = np.asarray([testcase.func(g[0]) for g in population])
    obj_min = np.zeros(n_gen + 1)
    circuit_min = []
    ind_min = np.argmin(obj)
    obj_min[0] = obj[ind_min]
    circuit_min.append(population[ind_min])

    for gen in range(n_gen):
        # if np.random.uniform() < prob:
        children = crossover(population, obj)
        # else:
        # children = deepcopy(population)
        mutate(testcase, children, 1.)
        obj_children = np.asarray([-testcase.simulate(g[0]) / testcase.ref[g[0].promo_node]['on'] for g in children])
        obj = np.append(obj, obj_children)
        population = np.vstack((population, children))
        S = np.lexsort([obj])
        obj = obj[S[:num_circuits]]
        population = population[S[:num_circuits], :]

        ind_min = np.argmin(obj)

        obj_min[gen + 1] = obj[ind_min]
        circuit_min.append(population[ind_min])

    generations = np.arange(n_gen + 1)

    # return