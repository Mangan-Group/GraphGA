from parmoo import MOOP
from parmoo.optimizers import LocalGPS

# Set the seed
seed = 20

# Create MOOP object using LocalGPS
my_moop = MOOP(LocalGPS)

# Add design variables (the hyperparameters being varied)
my_moop.addDesign({'name' : "mut_rate",
                   'des_type': "continuous", # Design variable type
                   'lb': 0.0, # Lower bound
                   'ub': 1.0}) # Upper bound

my_moop.addDesign({'name' : "cov_rate",
                   'des_type': "continuous",
                   'lb': 0.0,
                   'ub': 1.0})

# my_moop.addDesign({'name' : "population",
#                    'des_type': "integer",
#                    'lb': 2,
#                    'ub': 50})


from parmoo.searches import LatinHypercube
from parmoo.surrogates import GaussRBF
from main_function import *

# Define the simulation function
def sim_func(x):

    # Define the simulation setup variables
    promo_node = 'P1'
    num_dict = {1: 5, 2: 5}
    min_dose = 75
    max_dose = 75
    dose_interval = 5
    max_part = 2
    inhibitor = True

    # Uncomment to vary population size with a 50/50 ratio of 1-part to 2-part topologies
    # num_dict = {1: int(x['population']), 2: int(x['population'])}

    value = full_sim(x["mut_rate"], x["cov_rate"], promo_node, num_dict, max_part, min_dose, max_dose, dose_interval, inhibitor, seed)

    print("hi")
    return value

# hyperparams is where you would pass hyperparameters into the gaussRBF and latin hypercube if you defined your own methods
my_moop.addSimulation({'name': "TestSim",
                       'm': 2, # Number of simulation outputs
                       'sim_func': sim_func, # Variable for the function
                       'search': LatinHypercube, # Search method
                       'surrogate': GaussRBF, # Default surrogate function
                       'hyperparams': {'search_budget': 5}, # number of experiments to get initial data
                       })

# Define objective functions as their respective simulation ouptuts and add to MOOP object
def fitness(x,s): return s["TestSim"][0]
my_moop.addObjective({'name':"fitness", 'obj_func': fitness})
def convergence(x,s): return s["TestSim"][1]
my_moop.addObjective({'name':"convergence", 'obj_func': convergence})





from parmoo.acquisitions import RandomConstraint

# The range is the number of acquisition functions in each iteration after the initial search
for i in range(1):
    my_moop.addAcquisition({'acquisition': RandomConstraint, # Default acquisiton function
                            'hyperparams': {}})
    # Hyperparams is where you would send hyperparameters to self-defined acquisition functions


# Number in the solve method is the number of iterations
my_moop.solve(1)

# Results is the pareto front
results = my_moop.getPF(format = "pandas")

print(results)

# Save results
with open("optimization_consistent_pop1.pkl", "wb") as fid:
    pickle.dump(results, fid)

with open("simulation_consistent_pop1.pkl", "wb") as fid:
    pickle.dump(my_moop.getSimulationData(), fid)

with open("objective_consistent_pop1.pkl", "wb") as fid:
    pickle.dump(my_moop.getObjectiveData(), fid)




# Uncomment to open the pareto front on local device
# from parmoo.viz import scatter
# scatter(my_moop)