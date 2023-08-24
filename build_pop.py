from GA import sampling
import pickle

# Define the variables for the sample you want to build
promo_node = 'P1'
num_dict = {1: 250, 2: 250}
min_dose = 5
max_dose = 75
dose_interval = 5
max_part = 2
inhibitor = True

# Create the population and store as a pickle file
population = sampling(promo_node, num_dict, min_dose, max_dose, dose_interval, inhibitor=True)

with open("init_pop_inhib_500.pkl", "wb") as fid:
    pickle.dump(population, fid)

print(population)