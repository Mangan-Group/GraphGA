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

num_circuits = len(population)
n_gen = 30

obj = np.asarray([problem.func(g[0]) for g in population])

a = 1
