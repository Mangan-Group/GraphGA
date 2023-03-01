from GA import *
import sys
import matplotlib.pyplot as plt
from define_problem import *
from rankcrowding import *

# seed = int(sys.argv[1])
# np.random.seed(seed)
np.random.seed(20)


def simulate(topology, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_off = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('off', topology,))[-1, -1]
    rep_on = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('on', topology,))[-1, -1]
    return rep_off, rep_on


def func(toplogy):
    rep_off, rep_on = simulate(toplogy)
    ON_rel = rep_on / Ref[toplogy.promo_node]['on']
    FI_rel = (rep_on / rep_off) / (Ref[toplogy.promo_node]['fi'])
    return [-ON_rel, -FI_rel]


promo_node = 'P1'
num_dict = {2: 30}
min_dose = 15
max_dose = 75
dose_interval = 5
max_part = 2
inhibitor = True

problem = Problem(promo_node, max_part, min_dose, max_dose, dose_interval, inhibitor, func)

population = sampling(problem.promo_node, num_dict, problem.min_dose, problem.max_dose, problem.dose_interval, problem.inhibitor)

# with open("init_pop_%d.pkl" % seed, "wb") as fid:
#     pickle.dump(population, fid)

nds = RankAndCrowding()
num_circuits = len(population)
n_gen = 10

all_obj = []
obj = [problem.func(g[0]) for g in population]
all_obj.append(obj)
obj = np.asarray(obj)


for gen in range(n_gen):
    _, rank_dict = nds.do(obj, num_circuits, return_rank=True)
    # if np.random.uniform() < prob:
    children = crossover(population, obj, rank_dict)
    # else:
    # children = deepcopy(population)
    mutate(problem, children, 1.)
    obj_children = [problem.func(g[0]) for g in children]
    all_obj.append(obj_children)
    obj_children = np.asarray(obj_children)

    obj = np.vstack((obj, obj_children))
    population = np.vstack((population, children))

    S = nds.do(obj, num_circuits)
    obj = obj[S]
    population = population[S, :]

    # ind_min = np.argmin(obj)
    #
    # obj_min[gen + 1] = obj[ind_min]
    # circuit_min.append(population[ind_min])
fronts = NonDominatedSorting().do(obj)
all_obj = np.asarray(all_obj).reshape(num_circuits*(1 + n_gen), 2)
a = 1
