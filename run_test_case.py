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
# from define_circuit import Topo

# To-Do:
# save files: add results path key:val to settings dict

def run(
        testcase: object,
        settings: dict,
        metrics: bool =False
    ):
    '''Run the genetic algorithm for a given test case class instance'''

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
            obj
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

        with open(
            settings["results_path"]+
            "GA_results/"+
            settings["file_name"],
            "wb"
        ) as fid:
            pickle.dump(results, fid)

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

        with open(
            settings["results_path"]+
            "Combinatorial_results/"+
            settings["file_name"],
            "wb"
        ) as fid:
            pickle.dump(results, fid)

        return results


settings = {
    "promo_node":"P1",
    "dose_specs": [75, 75, 5],
    "max_part": 2,
    "inhibitor": False,
    "DsRed_inhibitor": False,
    "num_dict": {1: 10, 2: 20},
    "n_gen": 30,
    "pop": True,
    "num_processes": 1,
    "results_path": "/Users/kdreyer/Desktop/Github/GraphGA/Results/",
    "file_name": "230309_Amplifier_2part.pkl"
}
# topo_path = "Amplifier/Amplifier_combo_2.pkl"

# results = run_combinitorial(Amplifier, settings, topo_path)
results = run(Amplifier, settings, metrics=True)
print(results["obj_min"][-1])
print(results["circuit_min"][-1][0].edge_list)
print(results["circuit_min"][-1][0].dose)

results["circuit_min"][-1][0].plot_graph()

# print(settings)

# if __name__ == '__main__':
#     run_amp_pop_GA = '''run(Amplifier, settings)'''

#     n = 5
#     result = timeit.timeit(stmt=run_amp_pop_GA, globals=globals(), number=n)
#     print(f"Execution time is {result / n} seconds")

