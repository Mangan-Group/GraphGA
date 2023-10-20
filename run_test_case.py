import timeit
import numpy as np
import pickle
import networkx as nx
import matplotlib.pyplot as plt
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
# from load_Z_mat_samples import Z_mat_list
# from define_circuit import Topo

seed = 1
np.random.seed(seed)

def run(
        testcase: object,
        settings: dict,
        metrics: bool =False
    ):
    '''Run the genetic algorithm for a given test case'''

    folder_path = make_main_directory(settings)

    problem = testcase(
        settings["promo_node"],
        settings["dose_specs"],
        settings["max_part"],
        settings["inhibitor"],
        settings["DsRed_inhibitor"],
        settings["num_dict"],
        settings["n_gen"],
        settings["pop"],
        num_processes=settings["num_processes"],
    )
    
    population = sampling(
        problem.promo_node,
        problem.num_dict,
        problem.min_dose,
        problem.max_dose,
        problem.dose_interval,
        inhibitor=problem.inhibitor
    )
    file_name = "initial_population.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(population, fid)

    num_circuits = len(population)
    # raise an exception if num_circuits is odd
    if num_circuits % 2 != 0:
        raise Exception("Population size must be an even number")
    
    # set # generations according to testcase attribute
    n_gen = problem.n_gen

    # calculate objective for each circuit in initial population
    obj = np.asarray([problem.func(g[0]) for g in population])
    
    if isinstance (obj[0], np.ndarray):
        print('objective is an array- using'+
              'multi-objective optimization')
        (_, obj, _, top_circuits_each_gen,
         top_obs_each_gen) = multi_obj_GA(
            problem,
            n_gen,
            population,
            num_circuits,
            obj,
            settings["probability_crossover"],
            settings["probability_mutation"],
            settings["mutate_dose"]
        )
        print("final objectives: ", obj)

        file_name = "final_objectives.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(obj, fid)

        file_name = "top_circuits_each_gen.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(top_circuits_each_gen, fid)

        file_name = "top_objs_each_gen.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(top_objs_each_gen, fid)
    
        
    else:
        print('objective is not an array- using'+
              'single objective optimization')
        (obj_min, circuit_min, top_circuits_each_gen,
        top_objs_each_gen, _, _) = single_obj_GA(
            problem,
            n_gen,
            population,
            num_circuits,
            obj,
            settings["probability_crossover"],
            settings["probability_mutation"],
            settings["mutate_dose"],
            metrics
        )

        print("minimum objectives:", obj_min)
        # unique objectives in top num_circuit topologies
        ### could add all circuits from all gens (initial
        ### pop + children from each gen and get unique)
        unique_obj, unique_indices = np.unique(top_objs_each_gen,
                                               return_index=True)
        print(len(top_circuits_each_gen[unique_indices]))

        graph_file_name = "circuit_with_min_obj"
        plot_graph(circuit_min[-1][0], 
                   folder_path + "/" + graph_file_name)

        file_name = "minimum_obj_all_gens.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(obj_min, fid)

        file_name = "min_obj_circuit_all_gens.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(circuit_min, fid)

        file_name = "top_num_circuit_obj_each_gen.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(top_objs_each_gen, fid)

        file_name = "top_num_circuit_circuits_each_gen.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(top_circuits_each_gen, fid)


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

def plot_graph(topology, file_name):
    plt.figure()
    plt.tight_layout()
    nx.draw_networkx(topology.graph, arrows=True, arrowsize=15, node_size=600, node_shape='s')
    plt.savefig(file_name+".svg")


# add results folder name 
# make this a .json file (config_Amplifier, config_SignalConditioner, config_PulseGenerator)
settings = {
    "promo_node":"P1",
    "dose_specs": [75, 75, 5],
    "max_part": 2,
    "inhibitor": True,
    "DsRed_inhibitor": False,
    "num_dict": {1: 26, 2: 26},
    "n_gen": 1,#40,
    "pop": False,
    "probability_crossover": 1.0,
    "probability_mutation": 1.0,
    "mutate_dose": False,
    "num_processes": 1,
    "repository_path": "/Users/kdreyer/Documents/Github/GraphGA/",
    "folder_name": "Amplifier_single_cell_test",
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


# Z_mat_sampling = run_combinitorial_pop_samples(
#     Amplifier, settings, "/Users/kdreyer/Desktop/Github/GraphGA/Amplifier/",
#     "Amplifier_topos_all.pkl", "Amplifier_ON_rel_all.pkl",
#     Z_mat_list, 64.86873478791355)

# print(Z_mat_sampling)

run(Amplifier, settings)
