import timeit
import numpy as np
import pickle
import networkx
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from GA import sampling
from GA_setup import (
    single_obj_GA,
    multi_obj_GA
)
from load_Z_mat_samples import Z_mat_list
# from define_circuit import Topo

def run(
        testcase: object,
        settings: dict,
        metrics: bool =False
    ):
    '''Run the genetic algorithm for a given test case'''
    # look at GAMES code- create directory and cd into directory
    # based on settings 
    problem = testcase(
        settings["promo_node"],
        settings["dose_specs"],
        settings["max_part"],
        settings["inhibitor"],
        settings["DsRed_inhibitor"],
        settings["num_dict"],
        settings["n_gen"],
        settings["pop"],
        settings["num_processes"],
    )
    
    population = sampling(
        problem.promo_node,
        problem.num_dict,
        problem.min_dose,
        problem.max_dose,
        problem.dose_interval,
        problem.inhibitor
    )
    #add code to save population
    num_circuits = len(population)
    # raise an exception if num_circuits is odd
    if num_circuits % 2 != 0:
        raise Exception("Population size must be an even number")
    
    # set # generations according to testcase attribute; store as array
    n_gen = problem.n_gen
    generations = np.arange(n_gen + 1)

    # calculate objective for each circuit in initial population
    obj = np.asarray([problem.func(g[0]) for g in population])
    
    if isinstance (obj[0], np.ndarray):
        print('objective is an array- using multi-objective optimization')
        fronts, all_obj = multi_obj_GA(
            problem,
            n_gen,
            population,
            num_circuits,
            obj,
            settings["probability_crossover"],
            settings["probability_mutation"],
        )
        results = {
            "settings": settings,
            "fronts": fronts,
            "all_obj": all_obj
        }

        with open(
            settings["results_path"]+
            "GA_results/"+
            settings["file_name"],
            "wb"
        ) as fid:
            pickle.dump(results, fid)

    else:
        print('objective is not an array- using single objective optimization')
        all_obj, obj_min, circuit_min, geno, pheno = single_obj_GA(
            problem,
            n_gen,
            population,
            num_circuits,
            obj,
            settings["probability_crossover"],
            settings["probability_mutation"],
            metrics
        )
        results = {
            "settings": settings,
            "all_obj": all_obj,
            "obj_min": obj_min,
            "circuit_min": circuit_min,
            "genotype": geno,
            "phenotype": pheno
        }

        if settings["pop"]:
            with open(
                settings["results_path"]+
                "GA_results/"+
                settings["file_name"]+
                "_pop"+
                ".pkl",
                "wb"
            ) as fid:
                pickle.dump(results, fid)

        else:
            with open(
                settings["results_path"]+
                "GA_results/"+
                settings["file_name"]+
                ".pkl",
                "wb"
            ) as fid:
                pickle.dump(results, fid)

    print(results["obj_min"])
    print(results["circuit_min"][-1][0].edge_list)

    return results

def run_combinitorial(
        testcase: object,
        settings: dict,
        topo_path: str
):
    
    problem = testcase(
        settings["promo_node"],
        settings["dose_specs"],
        settings["max_part"],
        settings["inhibitor"],
        settings["DsRed_inhibitor"],
        settings["num_dict"],
        settings["n_gen"],
        settings["pop"],
        settings["num_processes"],
    )
    
    with open(topo_path, "rb") as fid:
        topologies = pickle.load(fid)
    
    obj = np.asarray([problem.func(g) for g in topologies])

    if isinstance (obj[0], np.ndarray):
        print('objective is an array- using non-dominated sorting')

    else:
        print('objective is not an array- sorting by single objective')
        sorted_idx = np.argsort(obj)
        obj_sorted = obj[sorted_idx]
        topo_sorted = np.asarray(topologies)[sorted_idx]

        results = {
            "objectives": obj_sorted,
            "topologies": topo_sorted
        }

        if settings["pop"]:
            with open(
                settings["results_path"]+
                "Combinatorial_results/"+
                settings["file_name"]+
                "_pop"+
                ".pkl",
                "wb"
            ) as fid:
                pickle.dump(results, fid)

        else:
            with open(
                settings["results_path"]+
                "Combinatorial_results/"+
                settings["file_name"]+
                ".pkl",
                "wb"
            ) as fid:
                pickle.dump(results, fid)

        print(results["objectives"])
        print(results["topologies"][0].edge_list)
        return results

def run_combinitorial_pop_samples(
        testcase: object,
        settings: dict,
        path: str,
        topo_fname: str,
        sc_obj_fname: str,
        Z_mat_list: list,
        obj_threshold: float
):
    
    problem = testcase(
        settings["promo_node"],
        settings["dose_specs"],
        settings["max_part"],
        settings["inhibitor"],
        settings["DsRed_inhibitor"],
        settings["num_dict"],
        settings["n_gen"],
        settings["pop"],
    )
    
    with open(path+topo_fname, "rb") as fid:
        topologies = pickle.load(fid)

    with open(path+sc_obj_fname, "rb") as fid:
        sc_obj = pickle.load(fid)
    

    idx_list = np.argwhere(sc_obj >= obj_threshold)
    topologies_list = topologies[idx_list]
    Z_mat_sampling = {}
    # Z_mat_sampling["Z_mat_list"] = Z_mat_list
    # Z_mat_sampling["topology_list"] = topologies_list

    for i, topology in enumerate(topologies_list):
        topology_dict = {}
        obj_list = []
        for Z_mat in Z_mat_list:
            problem.Z = Z_mat
            obj = problem.func(topology[0])
            obj_list.append(obj)
        
        topology_dict["objectives_range"] = max(obj_list) - min(obj_list)
        topology_dict["objectives_mean"] = np.mean(obj_list)
        topology_dict["objectives [min, max]"] = [min(obj_list), max(obj_list)]
        topology_dict["objectives_list"] = obj_list
        topology_dict["edge_list"] = topology[0].edge_list
        Z_mat_sampling["topologies_list["+str(i)+"]"] = topology_dict

        with open(
                settings["results_path"]+
                "Combinatorial_results/"+
                settings["file_name"],
                "wb"
            ) as fid:
                pickle.dump(Z_mat_sampling, fid)

    return Z_mat_sampling


# add results folder name 
# make this a .json file (config_Amplifier, config_SignalConditioner, config_PulseGenerator)
settings = {
    "promo_node":"P1",
    "dose_specs": [75, 75, 5],
    "max_part": 2,
    "inhibitor": False,
    "DsRed_inhibitor": False,
    "num_dict": {1: 10, 2: 20},
    "n_gen": 1,
    "pop": True,
    "probability_crossover": 1.0,
    "probability_mutation": 1.0,
    "num_processes": 1,
    "results_path": "/Users/kdreyer/Desktop/Github/GraphGA/Results/",
    "file_name": "230403_Amplifier_pop_sampling3.pkl"
}


# topo_path_1 = "Amplifier/Amplifier_combo_1.pkl"
# topo_path_2 = "Amplifier/Amplifier_combo_2.pkl"

# run_amp_pop_GA = '''run_combinitorial(Amplifier, settings, topo_path_2)'''

# n = 1
# result = timeit.timeit(stmt=run_amp_pop_GA, globals=globals(), number=n)
# print(f"Execution time is {result / n} seconds")







# topo_path_1 = "Amplifier/Amplifier_combo_1.pkl"
# topo_path_2 = "Amplifier/Amplifier_combo_2.pkl"
# results = run_combinitorial(Amplifier, settings, topo_path_2)
# print(results["objectives"][0])
# print(results["topologies"][0].edge_list)

# results = run(
#     Amplifier,
#     settings,
#     metrics=True
# )
# print(results["obj_min"])
# print(results["circuit_min"][0].edge_list)


# results = run(Amplifier, settings, metrics=True)
# print(results["obj_min"][-1])
# print(results["circuit_min"][-1][0].edge_list)
# print(results["circuit_min"][-1][0].dose)

# results["circuit_min"][-1][0].plot_graph()

# print(settings)


Z_mat_sampling = run_combinitorial_pop_samples(
    Amplifier, settings, "/Users/kdreyer/Desktop/Github/GraphGA/Amplifier/",
    "Amplifier_topos_all.pkl", "Amplifier_ON_rel_all.pkl",
    Z_mat_list, 64.86873478791355)

print(Z_mat_sampling)
