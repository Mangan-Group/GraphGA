from GA import *
from get_system_equations import system_equations
from define_problem import *


# Calculate the convergence by determining when the final GA output was first found
def first_seen(progression):
    gen_num = 0

    # Starting from the end, return the generation where the output first does not match the final output
    for gen in reversed(progression):
        if progression[-1] != gen:
            return len(progression) - gen_num
        gen_num += 1
    return 0


def simulate(topology, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_on = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('on', topology,))[-1, -1]
    return rep_on


def func(topology):
    return -simulate(topology) / Ref[topology.promo_node]['on']


# Define function by taking in necessary inputs to define the problem and set the seed
def full_sim(mut_rate, cov_rate, promo_node, num_dict, max_part, min_dose, max_dose, dose_interval, inhibitor, seed):

    # Set seed based on inputted number to allow for reproducibility
    np.random.seed(seed)

    problem = Problem(promo_node, max_part, min_dose, max_dose, dose_interval, inhibitor, func)

    # Generate a population using the commented code or use a set population
    population = sampling(problem.promo_node, num_dict, problem.min_dose, problem.max_dose, problem.dose_interval, inhibitor)

    # with open("init_pop_inhib_fixed.pkl", "rb") as fid:
    #      population = pickle.load(fid)

    # From this point, the code continues the same as the original main.py file
    num_circuits = len(population)
    n_gen = 50

    obj = np.asarray([problem.func(g[0]) for g in population])

    obj_min = np.zeros(n_gen + 1)
    circuit_min = []

    ind_min = np.argmin(obj)
    obj_min[0] = obj[ind_min]
    circuit_min.append(population[ind_min])

    for gen in range(n_gen):
        if np.random.uniform() < cov_rate:
            children = crossover(population, obj)
        else:
            children = deepcopy(population)
        mutate(problem, children, mut_rate, dose=True)
        obj_children = np.asarray([-simulate(g[0]) / Ref[g[0].promo_node]['on'] for g in children])
        obj = np.append(obj, obj_children)
        population = np.vstack((population, children))
        S = np.lexsort([obj])
        obj = obj[S[:num_circuits]]
        population = population[S[:num_circuits], :]

        ind_min = np.argmin(obj)

        obj_min[gen + 1] = obj[ind_min]
        circuit_min.append(population[ind_min])

    # Return the best fitness and the convergence found by the first_seen function
    return [obj_min[-1], first_seen(circuit_min)]