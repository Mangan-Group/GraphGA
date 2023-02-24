from GA import *
import sys
import matplotlib.pyplot as plt
from define_problem import *
# seed = int(sys.argv[1])
# np.random.seed(seed)
np.random.seed(20)


# class Problem:
#     def __init__(self, promo_node, max_part, min_dose, max_dose, dose_interval, inhibitor):
#         self.promo_node = promo_node
#         self.max_part = max_part
#         self.min_dose = min_dose
#         self.max_dose = max_dose
#         self.dose_interval = dose_interval
#         self.inhibitor = inhibitor
def simulate(topology, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_on = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('on', topology,))[-1, -1]
    return rep_on


def func(topology):
    return -simulate(topology) / Ref[topology.promo_node]['on']


promo_node = 'P1'
num_dict = {1: 10, 2: 30}
min_dose = 75
max_dose = 75
dose_interval = 5
max_part = 2
inhibitor = False

problem = Problem(promo_node, max_part, min_dose, max_dose, dose_interval, inhibitor, func)

# with open("init_population.pkl", "rb") as fid:
#    population = pickle.load(fid)

population = sampling(problem.promo_node, problem.num_dict, problem.min_dose, problem.max_dose, problem.dose_interval)

# with open("init_pop_%d.pkl" % seed, "wb") as fid:
#     pickle.dump(population, fid)

num_circuits = len(population)
n_gen = 30

obj = np.asarray([problem.func(g[0]) for g in population])

obj_min = np.zeros(n_gen + 1)
circuit_min = []

ind_min = np.argmin(obj)

ind_min = np.argmin(obj)
obj_min[0] = obj[ind_min]
circuit_min.append(population[ind_min])

for gen in range(n_gen):
    # if np.random.uniform() < prob:
    children = crossover(population, obj)
    # else:
    # children = deepcopy(population)
    mutate(problem, children, 1.)
    obj_children = np.asarray([-simulate(g[0]) / Ref[g[0].promo_node]['on'] for g in children])
    obj = np.append(obj, obj_children)
    population = np.vstack((population, children))
    S = np.lexsort([obj])
    obj = obj[S[:num_circuits]]
    population = population[S[:num_circuits], :]

    ind_min = np.argmin(obj)

    obj_min[gen + 1] = obj[ind_min]
    circuit_min.append(population[ind_min])

generations = np.arange(n_gen + 1)

#
# with open("min_topo_%d.pkl" % seed, "wb") as fid:
#     pickle.dump(circuit_min, fid)

a = 1
