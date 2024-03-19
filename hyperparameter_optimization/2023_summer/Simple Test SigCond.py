from parmoo import MOOP
from parmoo.optimizers import LocalGPS

# Create MOOP object and set the seed
my_moop = MOOP(LocalGPS)
seed = 20

# Add design variables
my_moop.addDesign({'name' : "mut_rate",
                   'des_type': "continuous",
                   'lb': 0.0,
                   'ub': 1.0})

my_moop.addDesign({'name' : "cov_rate",
                   'des_type': "continuous",
                   'lb': 0.0,
                   'ub': 1.0})

my_moop.addDesign({'name' : "population",
                   'des_type': "integer",
                   'lb': 2,
                   'ub': 50})





from parmoo.searches import LatinHypercube
from parmoo.surrogates import GaussRBF
from main_function_sigcond import *
from pymoo.indicators.hv import HV

# Define hypervolume calculator object
ref_point = np.array([0, 0])
ind = HV(ref_point=ref_point)

# Define simulation function
def sim_func(x):

    # Define necessary problem variables
    promo_node = 'P1'
    min_dose = 5
    max_dose = 75
    dose_interval = 5
    max_part = 2
    inhibitor = True

    # Declare the population
    num_dict = {1: int(x["population"]), 2: int(x["population"])}

    values = full_sim(x["mut_rate"], x["cov_rate"], promo_node, num_dict, max_part, min_dose, max_dose, dose_interval, inhibitor, seed)

    # Output the hypervolume and convergence
    # ParMOO always minimizes, so return the negative value of the hypervolume
    return [-ind(values[0]), values[1]]

# Add simulation to MOOP object
my_moop.addSimulation({'name': "TestSim",
                       'm': 2,
                       'sim_func': sim_func,
                       'search': LatinHypercube,
                       'surrogate': GaussRBF,
                       'hyperparams': {'search_budget': 100}, # number of experiments to get initial data
                       })


# Define objective functions and add them to the MOOP object
def fitness(x,s):return s["TestSim"][0]
my_moop.addObjective({'name':"fitness", 'obj_func': fitness})
def convergence(x,s): return s["TestSim"][1]
my_moop.addObjective({'name':"convergence", 'obj_func': convergence})




# Define acquisition functions
from parmoo.acquisitions import RandomConstraint

for i in range(3):
    my_moop.addAcquisition({'acquisition': RandomConstraint,
                            'hyperparams': {}})


# Solve the MOOP
my_moop.solve(1)

# Output results
results = my_moop.getPF(format = "pandas")
print(results)

# Save MOOP results
with open("optimization_consistent_pop1.pkl", "wb") as fid:
    pickle.dump(results, fid)

with open("simulation_consistent_pop1.pkl", "wb") as fid:
    pickle.dump(my_moop.getSimulationData(), fid)

with open("objective_consistent_pop1.pkl", "wb") as fid:
    pickle.dump(my_moop.getObjectiveData(), fid)