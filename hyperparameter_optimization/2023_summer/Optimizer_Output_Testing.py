# Change the imported main function to change the output being tested
from main_function import *
import pandas as pd
from pymoo.indicators.hv import HV

# Define reference point and hypervolume calculator
ref_point = np.array([0, 0])
ind = HV(ref_point=ref_point)


# Define necessary simulation variables
cov_rate = 0.458389
mut_rate = 0.970409
promo_node = 'P1'
min_dose = 5
max_dose = 75
dose_interval = 5
max_part = 2
inhibitor = True
pop_ratio = 0.274159
total_population = 84 * 2

# Calculate the population ratio
one_part = int(total_population * pop_ratio)
two_part = int(total_population - one_part)

num_dict = {1: one_part, 2: two_part}


# Create an array of the fitnesses and convergences
objectives = []


# Iterate through 100 seeds as inputs to main_function
# Uses the same problem set up to explore the problem under random seeds
# Code currently set up for the signal conditioner test case
for i in range(88,101):

    # Run the simulation and store the results separately
    values = full_sim(mut_rate, cov_rate, promo_node, num_dict, max_part, min_dose, max_dose, dose_interval, inhibitor,i)
    fitness = values[0]
    convergence = values[1]

    objectives.append(values[0])

    # Print results (and calculate the hypervolume when working with a pareto front)
    print("Seed:", i, " HV:", -ind(fitness), " Conv:", convergence)

    # Uncomment the following code to plot pareto front with combinatorial pareto front
    # with open('SigCond_obj_pareto.pkl', 'rb') as f:
    #     circuit = pd.read_pickle(f)
    #
    # circuit = np.asarray(circuit)
    # plt.scatter(circuit[:, 0], circuit[:, 1], c = "Red")
    # plt.scatter(fitness[:, 0], fitness[:, 1], c = "Blue")
    # plt.savefig('SigCond Fixed/SigcCond_Pareto_%d' %i)

# Save the objective data
with open("optimizer_output_test.pkl", "wb") as fid:
    pickle.dump(objectives, fid)