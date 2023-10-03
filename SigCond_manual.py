from GA import *
import sys
import matplotlib.pyplot as plt
from define_problem import *
from rankcrowding import *

# seed = int(sys.argv[1])
seed = 2
np.random.seed(seed)


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
num_dict = {1: 20, 2: 20}
min_dose = 5  # originally 15
max_dose = 75
dose_interval = 5
max_part = 2
inhibitor = True

problem = Problem(promo_node, max_part, min_dose, max_dose, dose_interval, inhibitor, func)

population = sampling(problem.promo_node, num_dict, problem.min_dose,problem.max_dose, problem.dose_interval, problem.inhibitor)

# with open("init_pop_inhib_100.pkl", "rb") as fid:
#      population = pickle.load(fid)


# Define a reference point and hypervolume calculator
from pymoo.indicators.hv import HV
ref_point = np.array([0, 0])
ind = HV(ref_point=ref_point)

# Store the progression of hypervolumes
hvs = []

nds = RankAndCrowding()
num_circuits = len(population)
n_gen = 50

all_obj = []
obj = [problem.func(g[0]) for g in population]
all_obj.append(obj)
obj = np.asarray(obj)

for gen in range(n_gen):
    _, rank_dict = nds.do(obj, num_circuits, return_rank=True)
    if np.random.uniform() < 1:
        children = crossover(population, obj, rank_dict)
    else:
        children = deepcopy(population)
    mutate(problem, children, 1, dose=True)
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
    hvs.append(ind(obj))

fronts = NonDominatedSorting().do(obj)
all_obj = np.asarray(all_obj).reshape(num_circuits*(1 + n_gen), 2)
a = 1

# Save the progression of hypervolumes over generations in a pickle file
with open("Hypervolumes.pkl", "wb") as fid:
    pickle.dump(hvs, fid)

# Create a plot of the population compared to the original pareto front
import pandas as pd
import seaborn as sns

# Create a vector based on whether the circuit does or does not have an inhibitor
types = []
for topo in population:
    inhib = "Activators"
    for part in topo[0].part_list:
        if part[0] == "I":
            inhib = "Inhibitors"
    types.append(inhib)

# Add types vector to the dataframe
df = pd.DataFrame(obj)
df['Type'] = types


# Load combinatorial pareto front for comparison
with open('SigCond_obj_pareto.pkl', 'rb') as f:
    circuit = pd.read_pickle(f)

circuit = pd.DataFrame(circuit)
circuit['color'] = 'Combo'

# Plot pareto front with the simulation output
sns.scatterplot(data = circuit, x = circuit[0], y = circuit[1], hue = 'color', palette = 'cividis')
sns.scatterplot(data=df, x= df[0], y= df[1], hue='Type')
plt.suptitle("Signal Conditioner Pareto Front")
plt.title("Population of 40, Varied Population Dose")
plt.ylabel("FI_rel")
plt.xlabel("ON_rel")
plt.savefig("VaryDose", bbox_inches="tight")
plt.show()


# Can adjust the weighted average below to pinpoint a specific location on the pareto front (change the 90)
# Both values are squared because you want the closest point to the origin

# values = np.array(obj)
# weighted_averages = values[:,0]**2 + (values[:,1]*90)**2
# print('Closest Point:', values[np.where(weighted_averages == min(weighted_averages))])