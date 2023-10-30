import timeit
import numpy as np
import pickle
import pandas as pd
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

# seed = settings["seed"]
# np.random.seed(seed)

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
        settings["probability_crossover"],
        settings["probability_mutation"],
        settings["mutate_dose"],
        settings["pop"],
        settings["CI"],
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
    obj = np.asarray(
        [problem.func(g[0]) for g in population])
    
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
            settings["get_unique"],
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
            settings["get_unique"],
            metrics
        )

# outdated and will need updating to run
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

# need to update for new CI method
def run_combinitorial_pop_samples(
        testcase: object,
        settings: dict,
        Z_mat_list: list,
):
    topo_path = settings["single_cell_results_path"]
    folder_path = make_main_directory(settings)
    
    problem = testcase(
        settings["promo_node"],
        settings["dose_specs"],
        settings["max_part"],
        settings["inhibitor"],
        settings["DsRed_inhibitor"],
        settings["num_dict"],
        settings["n_gen"],
        settings["probability_crossover"],
        settings["probability_mutation"],
        settings["mutate_dose"],
        settings["pop"],
    )
    
    with open(topo_path, "rb") as fid:
        topologies = pickle.load(fid)

    Z_mat_sampling = {}

    if testcase == Amplifier:
        for i, topology in enumerate(topologies):
            topology_dict = {}
            obj_list = []
            for Z_mat in Z_mat_list:
                problem.Z = Z_mat
                obj = problem.func(topology[0])
                obj_list.append(obj)
            
            topology_dict["objectives_range"] = max(obj_list) - min(obj_list)
            topology_dict["objectives_mean"] = np.mean(obj_list)
            topology_dict["objectives [min, max]"] = [[min(obj_list), max(obj_list)]]
            topology_dict["objectives_list"] = [obj_list]
            topology_dict["edge_list"] = [topology[0].edge_list]
            Z_mat_sampling["topologies_list["+str(i)+"]"] = topology_dict
            print("topology " + str(i) + " done")

        Z_mat_sampling_df = pd.DataFrame(Z_mat_sampling)
        Z_mat_sampling_df = Z_mat_sampling_df.transpose()
        print(Z_mat_sampling_df)

        file_name = "Z_mat_sampling_all_GA_circuits.pkl"
        Z_mat_sampling_df.to_pickle(folder_path + "/" + file_name)

        ON_rel_range = Z_mat_sampling_df["objectives_range"].tolist()
        ON_rel_range_mean = round(np.mean(ON_rel_range), 4)
        fig_text = "mean = " + str(ON_rel_range_mean)
        fig_name = "ON_rel_range_distribution.svg"
        fig_path = folder_path + "/" + fig_name
        plot_obj_distribution(fig_path, ON_rel_range, 
                        "ON_rel range", fig_text)

    else:
        for i, topology in enumerate(topologies[:1000]):
            topology_dict = {}
            obj_list = []
            for Z_mat in Z_mat_list:
                problem.Z = Z_mat
                obj = problem.func(topology[0])
                obj_list.append(obj)
            obj_list = np.asarray(obj_list)
            topology_dict["objectives[0]_range"] = max(obj_list[:, 0]) - min(obj_list[:, 0])
            topology_dict["objectives[0]_mean"] = np.mean(obj_list[:, 0])
            topology_dict["objectives[0] [min, max]"] = [[min(obj_list[:, 0]), max(obj_list[:, 0])]]
            topology_dict["objectives[0]_list"] = obj_list[:, 0]
            topology_dict["objectives[1]_range"] = max(obj_list[:, 1]) - min(obj_list[:, 1])
            topology_dict["objectives[1]_mean"] = np.mean(obj_list[:, 1])
            topology_dict["objectives[1] [min, max]"] = [[min(obj_list[:, 1]), max(obj_list[:, 1])]]
            topology_dict["objectives[1]_list"] = obj_list[:, 1]
            topology_dict["edge_list"] = [topology[0].edge_list]
            Z_mat_sampling["topologies_list["+str(i)+"]"] = topology_dict
            print("topology " + str(i) + " done")

        Z_mat_sampling_df = pd.DataFrame(Z_mat_sampling)
        Z_mat_sampling_df = Z_mat_sampling_df.transpose()
        print(Z_mat_sampling_df)

        file_name = "Z_mat_sampling_all_GA_circuits.pkl"
        Z_mat_sampling_df.to_pickle(folder_path + "/" + file_name)

        ON_rel_range = Z_mat_sampling_df["objectives[0]_range"].tolist()
        ON_rel_range_mean = round(np.mean(ON_rel_range), 4)
        fig_text1 = "mean = " + str(ON_rel_range_mean)
        fig_name1 = "ON_rel_range_distribution.svg"
        fig_path1 = folder_path + "/" + fig_name1
        plot_obj_distribution(fig_path1, ON_rel_range, 
                        "ON_rel range", fig_text1)
        
        FI_rel_range = Z_mat_sampling_df["objectives[1]_range"].tolist()
        FI_rel_range_mean = round(np.mean(FI_rel_range), 4)
        fig_text2 = "mean = " + str(FI_rel_range_mean)
        fig_name2 = "FI_rel_range_distribution.svg"
        fig_path2 = folder_path + "/" + fig_name2
        plot_obj_distribution(fig_path2, FI_rel_range, 
                        "FI_rel range", fig_text2)



# make this a .json file (config_Amplifier, config_SignalConditioner, config_PulseGenerator)
settings = {
    "promo_node":"P1",
    "dose_specs": [75, 75, 5],
    "max_part": 2,
    "inhibitor": True,
    "DsRed_inhibitor": False,
    "num_dict": {1: 26, 2: 26},
    "n_gen": 2,
    "probability_crossover": 0.55,
    "probability_mutation": 1.0,
    "mutate_dose": False,
    "pop": False,
    "CI": None,
    "num_processes": 1,
    "get_unique": False,
    "plot": False,
    "seed": 0,
    "repository_path": "/Users/kdreyer/Documents/Github/GraphGA/",
    "folder_name": "Amplifier_single_cell_seed_0"
}


# run_amp_pop_GA = '''run_combinitorial(Amplifier, settings, topo_path_2)'''
# n = 1
# result = timeit.timeit(stmt=run_amp_pop_GA, globals=globals(), number=n)
# print(f"Execution time is {result / n} seconds")


for seed in range(0, 1):
    np.random.seed(seed)
    settings["seed"] = seed
    settings["folder_name"] = "Amplifier_single_cell_test_seed_" + str(seed)

    run(Amplifier, settings)
    print("seed "+str(seed)+" complete")

# run_combinitorial_pop_samples(SignalConditioner, settings,
#                               Z_mat_list)