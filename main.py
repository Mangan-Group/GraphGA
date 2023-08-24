from GA import *
import sys
from define_problem import *
import copy

seed = int(sys.argv[1])
np.random.seed(seed)

def first_seen(progression):

    looking = True
    gen_num = 0

    for gen in reversed(progression):
        if progression[-1] != gen:
            return len(progression) - gen_num
        gen_num += 1

    return 0


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
    # on_rel / Ref[topology.promo_node]['on']
    # should be in the range 30 - 69


promo_node = 'P1'
num_dict = {1: 10, 2: 30}
min_dose = 75
max_dose = 75
dose_interval = 5
max_part = 2
inhibitor = False

problem = Problem(promo_node, max_part, min_dose, max_dose, dose_interval, inhibitor, func)

# with open("init_pop_2.pkl", "rb") as fid:
#     population = pickle.load(fid)

population = sampling(problem.promo_node, num_dict, problem.min_dose, problem.max_dose, problem.dose_interval)

with open("init_pop_%d.pkl" % seed, "wb") as fid:
    pickle.dump(population, fid)

num_circuits = len(population)
n_gen = 35

needed_gens = []
results = []

for i in range(1,21):
    pop_i = copy.deepcopy(population)

    np.random.seed(i)

    obj = np.asarray([problem.func(g[0]) for g in pop_i])

    obj_min = np.zeros(n_gen + 1)
    circuit_min = []

    # ind_min = np.argmin(obj)

    ind_min = np.argmin(obj)
    obj_min[0] = obj[ind_min]
    circuit_min.append(pop_i[ind_min])

    for gen in range(n_gen):
        #if np.random.randint(2) == 1:
        if np.random.uniform() < .7:
            children = crossover(pop_i, obj)
        else:
            children = deepcopy(pop_i)
        if np.random.uniform() < 0.2:
            mutate(problem, children, 1.)
        obj_children = np.asarray([-simulate(g[0]) / Ref[g[0].promo_node]['on'] for g in children])
        obj = np.append(obj, obj_children)
        pop_i = np.vstack((pop_i, children))
        S = np.lexsort([obj])
        obj = obj[S[:num_circuits]]
        pop_i = pop_i[S[:num_circuits], :]

        ind_min = np.argmin(obj)

        obj_min[gen + 1] = obj[ind_min]
        circuit_min.append(pop_i[ind_min])

    needed_gens.append(first_seen(circuit_min))
    results.append(circuit_min[-1])

    print("Seed: ", i)

generations = np.arange(n_gen + 1)

with open("min_topo_%d.pkl" % seed, "wb") as fid:
    pickle.dump(circuit_min, fid)