from GA import *
import sys
import matplotlib.pyplot as plt
from define_problem import *
seed = int(sys.argv[1])
np.random.seed(seed)

c_prob = int(sys.argv[2])
m_prob = int(sys.argv[3])

folder = "Amplifier/numdict_2_s_0/"
with open(folder + "init_pop.pkl", "rb") as fid:
	population = pickle.load(fid)

def simulate(topology, max_time=42):
	t = np.arange(0, max_time + 1, 1)
	rep_on = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('on', topology,))[-1, -1]
	return rep_on


def func(topology):
	return -simulate(topology) / Ref[topology.promo_node]['on']


promo_node = 'P1'
min_dose = 75
max_dose = 75
dose_interval = 5
max_part = 2
inhibitor = False

problem = Problem(promo_node, max_part, min_dose, max_dose, dose_interval, inhibitor, func)


num_circuits = len(population)
n_gen = 50

all_circuits = []
all_population = []
all_circuits.append(population)
all_population.append(population)

obj = np.asarray([problem.func(g[0]) for g in population])
obj_min = np.zeros(n_gen + 1)
circuit_min = []

ind_min = np.argmin(obj)
obj_min[0] = obj[ind_min]
circuit_min.append(population[ind_min])


print(c_prob*0.1)
print(m_prob*0.1)

for gen in range(n_gen):
	if np.random.uniform() < c_prob*0.1:
		children = crossover(population, obj)
	else:
		children = deepcopy(population)
	mutate(problem, children, m_prob*0.1)

	all_circuits.append(children)

	obj_children = np.asarray([-simulate(g[0]) / Ref[g[0].promo_node]['on'] for g in children])
	obj = np.append(obj, obj_children)
	population = np.vstack((population, children))
	S = np.lexsort([obj])
	obj = obj[S[:num_circuits]]
	population = population[S[:num_circuits], :]
	
	all_population.append(population)

	ind_min = np.argmin(obj)
	obj_min[gen + 1] = obj[ind_min]
	circuit_min.append(population[ind_min])


folder += "c%d_m%d/" % (c_prob, m_prob) 
with open(folder + "all_circuits_c%d_m%d_s%d.pkl" % (c_prob, m_prob, seed), "wb") as fid:
	pickle.dump(all_circuits, fid)

with open(folder + "all_population_c%d_m%d_s%d.pkl" % (c_prob, m_prob, seed), "wb") as fid:
        pickle.dump(all_population, fid)

with open(folder + "obj_min_c%d_m%d_s%d.pkl" % (c_prob, m_prob, seed), "wb") as fid:
        pickle.dump(obj_min, fid)

with open(folder + "circuit_min_c%d_m%d_s%d.pkl" % (c_prob, m_prob, seed), "wb") as fid:
        pickle.dump(circuit_min, fid)
