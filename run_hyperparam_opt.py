import numpy as np
import json
import pickle
import pandas as pd
from statistics import mean
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from parmoo.extras.libe import libE_MOOP
from parmoo import MOOP
from parmoo.searches import LatinHypercube
from parmoo.surrogates import GaussRBF
from parmoo.acquisitions import RandomConstraint
from parmoo.optimizers import LocalGPS
from libensemble.tools import parse_args
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator
from GA import sampling
from GA_setup_hyperparameter_opt import (
    single_obj_GA,
    multi_obj_GA
)

_, is_manager, _, _ = parse_args()
seed = 20

settings = {
"test_case": "SignalConditioner",
"DsRed_inhibitor": False,
"n_gen": 50,
"pop": False,
"obj_labels": ["ON_rel", "FI_rel"],
"max_time": 42,
"seeds": [seed, seed+1],
"folder_path": "./Amplifier_single_cell/"
}
with open(settings["folder_path"] + "settings.json", "w") as fid:
    json.dump(settings, fid)

def run(x:dict):
    '''Run the genetic algorithm for a given test case'''
    
    if settings["test_case"] == "Amplifier":
        test_case = Amplifier
    elif settings["test_case"] == "SignalConditioner":
        test_case = SignalConditioner
    elif settings["test_case"] == "PulseGenerator":
        test_case = PulseGenerator

    num_circuits = x["population_size"]*2
    num_one_part_circuits = int(
        num_circuits*x["population_ratio"]
    )
    num_two_part_circuits = int((num_circuits - 
                             num_one_part_circuits
    ))
    num_dict = {
        1: num_one_part_circuits,
        2: num_two_part_circuits
    }
    
    problem = test_case(
        promo_node="P1",
        dose_specs=[5, 75, 5],
        max_part=2,
        inhibitor=True,
        DsRed_inhibitor=settings["DsRed_inhibitor"],
        num_dict=num_dict,
        n_gen=settings["n_gen"],
        probability_crossover= x["crossover_rate"],
        probability_mutation= x["mutation_rate"],
        mutate_dose=True,
        pop=settings["pop"],
        obj_labels=settings["obj_labels"],
        max_time=settings["max_time"]
    )

    seeds = [seed, seed+1]
    fitness = []
    convergence = []
    for seed_val in seeds:
        # np.random.seed(seed_val)
        # population = sampling(
        #     problem.promo_node,
        #     problem.num_dict,
        #     problem.min_dose,
        #     problem.max_dose,
        #     problem.dose_interval,
        #     inhibitor=problem.inhibitor
        # )
        # num_circuits = len(population)
        # # calculate objective for each circuit 
        # # in initial population
        # obj = np.asarray(
        #     [problem.func(g[0]) for g in population])

        # run multi-objective GA if multiple objectives
        # (determined by type of first first objective 
        # in array- will be array if more than one)
        if len(problem.obj_labels) > 1:
            [pareto_objs, gen_converged] = multi_obj_GA(
                problem,
                # population,
                # num_circuits,
                # obj,
                seed_val
            )
            fitness.append(pareto_objs)
            convergence.append(gen_converged)
            
        else:
            [final_obj, gen_converged] = single_obj_GA(
                problem,
                # population,
                # num_circuits,
                # obj,
                seed_val
            )
            fitness.append(final_obj)
            convergence.append(gen_converged)

    average_fitness = np.mean(fitness)
    average_convergence = np.mean(convergence)

    return [average_fitness, average_convergence]

# Define fitness objective function as the first simulation output
def fitness(x,s): return s["GA"][0]

# Define convergence objective function as the second simulation output
def convergence(x,s): return s["GA"][1]


if __name__ == "__main__":
    np.random.seed(seed)

    # Create a libE_MOOP
    my_moop = libE_MOOP(LocalGPS)

    # Add design variables (the hyperparameters being optimized)

    # Population will be multiplied by 2 to guarantee even population
    # Must have a population of at least 3 for the chosen throuple in crossing over
    my_moop.addDesign({'name': "population_size",
                       'des_type': "integer", # Variable type
                       'lb': 2, # Lower bound
                       'ub': 50}) # Upper bound

    my_moop.addDesign({'name': "population_ratio",
                       'des_type': "continuous",
                       'lb': 0.0,
                       'ub': 1.0})

    my_moop.addDesign({'name': "mutation_rate",
                       'des_type': "continuous",
                       'lb': 0.0,
                       'ub': 1.0})

    my_moop.addDesign({'name': "crossover_rate",
                       'des_type': "continuous",
                       'lb': 0.0,
                       'ub': 1.0})

    # Add the simulation (note the budget of 20 sim evals during search phase)
    my_moop.addSimulation({'name': "GA",
                           'm': 2, # Number of outputs
                           'sim_func': run, # Variable for simulation function
                           'search': LatinHypercube, # Search method
                           'surrogate': GaussRBF, # Function used to build the surface
                           'hyperparams': {'search_budget': 5},  # number of experiments to get initial data
                           })

    # Add the objectives defined in the functions above
    my_moop.addObjective({'name': "Fitness", 'obj_func': fitness})
    my_moop.addObjective({'name': "Convergence", 'obj_func': convergence})


    # Add acquisition functions
    # The number of acquisition functions determines the number of simulations in each batch of solving
    for i in range(2):
        my_moop.addAcquisition({'acquisition': RandomConstraint, # Using default acquisition function
                                'hyperparams': {}})
        # hyperparams is where you would send hyperparameters to self-written and defined acquisition functions

    # Value in the parentheses determines the number of iterations of the optimizer after the initial search
    my_moop.solve(5)

    # Display the solution -- this "if" clause is needed when running with MPI
    if is_manager:
        # Get the pareto front as a pandas dataframe
        results = my_moop.getPF(format="pandas")

        # Save optimization data, which is the pareto front data
        with open(settings["folder_path"] + "pareto_front.pkl", "wb") as fid:
            pickle.dump(results, fid)

        # Save simulation data
        with open(settings["folder_path"] + "simulation_results.pkl", "wb") as fid:
            pickle.dump(my_moop.getSimulationData(), fid)

        # Save objective data, provides the hyperparameters for each simulation output
        with open(settings["folder_path"] + "all_hyperparameters.pkl", "wb") as fid:
            pickle.dump(my_moop.getObjectiveData(), fid)

        # Print the pareto front
        print(results)




    