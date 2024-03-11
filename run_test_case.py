import timeit
import numpy as np
import pickle
import pandas as pd
from multiprocessing import Pool
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from saving import make_main_directory
from GA import sampling
from GA_setup import (
    single_obj_GA,
    multi_obj_GA
)
from plot_search_results import plot_obj_distribution
from load_Z_mat_samples import Z_mat_list


def run(
        testcase: object,
        settings: dict,
    ):
    '''Run the genetic algorithm for a given test case'''

    folder_path = make_main_directory(settings)

    problem = testcase(
        promo_node=settings["promo_node"],
        dose_specs=settings["dose_specs"],
        max_part=settings["max_part"],
        inhibitor=settings["inhibitor"],
        DsRed_inhibitor=settings["DsRed_inhibitor"],
        num_dict=settings["num_dict"],
        n_gen=settings["n_gen"],
        probability_crossover=settings["probability_crossover"],
        probability_mutation=settings["probability_mutation"],
        mutate_dose=settings["mutate_dose"],
        pop=settings["pop"],
        num_processes=settings["num_processes"],
        obj_labels=settings["obj_labels"],
        max_time=settings["max_time"]
    )
    
    population = sampling(
        problem.promo_node,
        problem.num_dict,
        problem.min_dose,
        problem.max_dose,
        problem.dose_interval,
        inhibitor=problem.inhibitor
    )

    # save initial population
    file_name = "initial_population.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(population, fid)

    num_circuits = len(population)

    # raise an exception if num_circuits is odd
    if num_circuits % 2 != 0:
        raise Exception("Population size must be an even number")
    
    # calculate objective for each circuit 
    # in initial population
    if settings["pop"]:
        topologies = [g[0] for g in population]
        with Pool(settings["num_processes"]) as pool:
            obj_list = pool.imap(problem.func, topologies)

            pool.close()
            pool.join()
        obj_list = list(obj_list)
        obj = np.asarray(obj_list)
    else:
        obj = np.asarray(
            [problem.func(g[0]) for g in population])
    # print(obj)
    # run multi-objective GA if multiple objectives
    # (determined by type of first first objective 
    # in array- will be array if more than one)
    if isinstance (obj[0], np.ndarray):
        print('objective is an array- using '+
              'multi-objective optimization')
        multi_obj_GA(
            folder_path,
            problem,
            population,
            num_circuits,
            obj,
            settings["plot"]
        )
        # print("final objectives: ", obj)
        
    else:
        print('objective is not an array- using '+
              'single objective optimization')
        single_obj_GA(
            folder_path,
            problem,
            population,
            num_circuits,
            obj,
        )


# make this a .json file (config_Amplifier, config_SignalConditioner, config_PulseGenerator)
settings = {
    "test_case": "SignalConditioner",
    "promo_node": "P1",
    "dose_specs": [5, 75, 5],
    "max_part": 2,
    "inhibitor": True,
    "DsRed_inhibitor": False,
    "num_dict": {1: 46, 2: 122},
    "n_gen": 50,
    "probability_crossover": 0.32, #0.32, increased to 0.5, then 0.75
    "probability_mutation": 0.57, #0.57, increased to 0.75, then 1.0
    "mutate_dose": True,
    "pop": True,
    "num_processes": 8,
    "obj_labels": ["ON_rel", "FI_rel"],
    "max_time": 42,
    "plot": False,
    "seed": 0,
    "repository_path": "/Users/kdreyer/Documents/Github/GraphGA/",
    "folder_name": "Signal_Cond_pop_inhibitor_ZF1_ZF2_new_dose_terms"
}

if __name__ == "__main__":
    # np.random.seed(0)
    # run_sc_pop_GA = '''run(SignalConditioner, settings)'''
    # n = 5
    # result = timeit.timeit(stmt=run_sc_pop_GA, globals=globals(), number=n)
    # print(f"Execution time is {result / n} seconds")

    if settings["test_case"] == "Amplifier":
        test_case = Amplifier
    elif settings["test_case"] == "SignalConditioner":
        test_case = SignalConditioner
    elif settings["test_case"] == "PulseGenerator":
        test_case = PulseGenerator
    else:
        raise Exception("Error: test case not defined")

    for seed in range(0, 1):
        np.random.seed(seed)
        settings["seed"] = seed
        settings["folder_name"] = settings["folder_name"] + "_seed_" + str(seed)

        run(test_case, settings)
        print("seed "+str(seed)+" complete")

# run_combinitorial_pop_samples(SignalConditioner, settings,
#                               Z_mat_list)