from main_function import *
import pickle
import numpy as np
from parmoo.extras.libe import libE_MOOP
from parmoo.searches import LatinHypercube
from parmoo.surrogates import GaussRBF
from parmoo.acquisitions import RandomConstraint
from parmoo.optimizers import LocalGPS

# When running with MPI, we need to keep track of which thread is the manager
# using libensemble.tools.parse_args()
from libensemble.tools import parse_args
_, is_manager, _, _ = parse_args()


# Set the seed
seed = 20

# Define simulation function
def sim_func(x):

    # Define necessary problem variables
    promo_node = 'P1'
    min_dose = 5
    max_dose = 75
    dose_interval = 5
    max_part = 2
    inhibitor = True

    # Following four line calculate population ratio when size and ratio are varied
    total_population = x["population"]*2
    one_part = int(total_population*x["pop_ratio"])
    two_part = int(total_population - one_part)
    num_dict = {1: one_part, 2: two_part}

    # Uncomment this line to only vary population size instead of ratio
    # num_dict = {1: int(x["population"]), 2: int(x["population"])}

    # Run two simulations with the given seed and the seed plus 1
    # Two simulations are run because random chance allows a bad set of hyperparameters
    # To find the optimal solution quickly but not robustly
    value1 = full_sim(x["mut_rate"], x["crsovr_rate"], promo_node, num_dict, max_part, min_dose, max_dose, dose_interval, inhibitor, seed)
    value2 = full_sim(x["mut_rate"], x["crsovr_rate"], promo_node, num_dict, max_part, min_dose, max_dose, dose_interval, inhibitor, (seed + 1))

    # Calculate and return the average of the two simulations for both objectives
    average = [(value1[0] + value2[0])/2, (value1[1] + value2[1])/2]
    return average


# Define fitness objective function as the first simulation output
def fitness(x,s): return s["GA"][0]


# Define convergence objective function as the second simulation output
def convergence(x,s): return s["GA"][1]


# When using libEnsemble with Python MP, the "solve" command must be enclosed
# in an "if __name__ == '__main__':" block, as shown below
if __name__ == "__main__":
    # Fix the random seed for reproducibility
    np.random.seed(seed)

    # Create a libE_MOOP
    my_moop = libE_MOOP(LocalGPS)

    # Add design variables (the hyperparameters being optimized)

    # Population will be multiplied by 2 to guarantee even population
    # Must have a population of at least 3 for the chosen throuple in crossing over
    my_moop.addDesign({'name': "population",
                       'des_type': "integer", # Variable type
                       'lb': 2, # Lower bound
                       'ub': 50}) # Upper bound

    my_moop.addDesign({'name': "pop_ratio",
                       'des_type': "continuous",
                       'lb': 0.0,
                       'ub': 1.0})

    my_moop.addDesign({'name': "mut_rate",
                       'des_type': "continuous",
                       'lb': 0.0,
                       'ub': 1.0})

    my_moop.addDesign({'name': "crsovr_rate",
                       'des_type': "continuous",
                       'lb': 0.0,
                       'ub': 1.0})

    # Add the simulation (note the budget of 20 sim evals during search phase)
    my_moop.addSimulation({'name': "GA",
                           'm': 2, # Number of outputs
                           'sim_func': sim_func, # Variable for simulation function
                           'search': LatinHypercube, # Search method
                           'surrogate': GaussRBF, # Function used to build the surface
                           'hyperparams': {'search_budget': 200},  # number of experiments to get initial data
                           })

    # Add the objectives defined in the functions above
    my_moop.addObjective({'name': "Fitness", 'obj_func': fitness})
    my_moop.addObjective({'name': "Convergence", 'obj_func': convergence})


    # Add acquisition functions
    # The number of acquisition functions determines the number of simulations in each batch of solving
    for i in range(20):
        my_moop.addAcquisition({'acquisition': RandomConstraint, # Using default acquisition function
                                'hyperparams': {}})
        # hyperparams is where you would send hyperparameters to self-written and defined acquisition functions

    # Value in the parentheses determines the number of iterations of the optimizer after the initial search
    my_moop.solve(10)

    # Display the solution -- this "if" clause is needed when running with MPI
    if is_manager:
        # Get the pareto front as a pandas dataframe
        results = my_moop.getPF(format="pandas")

        # Save optimization data, which is the pareto front data
        with open("optimization_inhib_full.pkl", "wb") as fid:
            pickle.dump(results, fid)

        # Save simulation data
        with open("simulation_inhib_full.pkl", "wb") as fid:
            pickle.dump(my_moop.getSimulationData(), fid)

        # Save objective data, provides the hyperparameters for each simulation output
        with open("objective_inhib_full.pkl", "wb") as fid:
            pickle.dump(my_moop.getObjectiveData(), fid)

        # Print the pareto front
        print(results)