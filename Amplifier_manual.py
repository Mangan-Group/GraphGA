from GA4Graph import *
from scipy.integrate import odeint
from load_files import *
from get_system_equations import system_equations
from Test_amplifier import simulate, Amplifier
from pymoo.algorithms.soo.nonconvex.ga import GA, FitnessSurvival

import matplotlib.pyplot as plt
import scipy.stats as ss
from itertools import chain

problem = Amplifier(n_var=1, n_obj=1, n_ieq_constr=0)

def check_valid(g, num_parts):
    graph_parts = [i[0] for i in g.nodes]
    if ('P' not in graph_parts) or ('R' not in graph_parts) or ('Z' not in graph_parts):
        return 0
    elif (len(graph_parts) - 2) < num_parts:
        return 0

    for n in g.nodes:
        if n[0] == 'P':
            out_types = set([i[0] for i in list(g.successors(n))])
            if 'Z' not in out_types:
                return 0
        elif n[0] == 'R':
            in_types = set([i[0] for i in list(g.predecessors(n))])
            if (not in_types) or (('I' in in_types) and ('Z' not in in_types)):
                return 0
        else:
            if (n[0] == 'Z') and ('I' + n[1:] in g.nodes):
                if list(g.successors(n)) != list(g.successors('I' + n[1:])):
                    return 0
            in_nodes = list(g.predecessors(n))
            if not in_nodes:
                return 0
            elif in_nodes == [n]:
                return 0
            else:
                in_types = set([i[0] for i in in_nodes])
                if ('I' in in_types) and ('Z' not in in_types):
                    return 0
            if len(list(nx.all_simple_paths(g, n, 'Rep'))) == 0:
                return 0
    return 1


def is_equal(x1, x2):
    return compare_circuit(x1[0], x2[0])


def crossover(problem, X, **kwargs):
    X = X[np.newaxis, :].reshape(2, -1, 1)
    # The input of has the following shape (n_parents, n_matings, n_var)
    _, n_matings, n_var = X.shape
    Y = np.full_like(X, None, dtype=object)
    for k in range(n_matings):
        # get the first and the second parent
        parent1, parent2 = X[0, k, 0], X[1, k, 0]
        Y[0, k, 0], Y[1, k, 0] = crossover_structure(parent1, parent2)
    return Y.reshape(-1, 1)


def mutate(problem, X, prob, **kwargs):
    for i in range(len(X)):
        if np.random.uniform(0, 1) < prob:
            mutate_node_num(X[i, 0], problem.max_part, problem.min_dose, problem.max_dose, problem.dose_interval, problem.inhibitor)
        # if np.random.uniform(0, 1) < prob:
            mutate_node_type(X[i, 0], problem.min_dose, problem.max_dose, problem.dose_interval)
        # if np.random.uniform(0, 1) < prob:
            mutate_dose(X[i, 0], problem.min_dose, problem.max_dose, problem.dose_interval)
    # return X


def compare(X):
    dupes_dict = dict()
    for i, g in enumerate(X):
        if i not in list(chain.from_iterable(dupes_dict.values())):
            dupes_dict.update({i: [j for j in range(i+1, len(X)) if is_equal(g, X[j])]})
    return len(dupes_dict)


def grouping(F, delta=0.5):
    bins = np.arange(np.min(F), np.max(F)+delta, delta)
    n = np.full_like(F, 0)
    for i in range(len(F)):
        n[i] = np.argwhere(F[i] >= bins).flatten()[-1]
    return n, bins


def pheno_diversity(F, delta=0.5):
    n, bins = grouping(F, delta)
    return len(set(n))


def get_entropy(n):
    _, p = np.unique(n, return_counts=True)
    p = p/len(p)
    entropy = -np.sum(p*np.log(p))
    return entropy


def pseudo_iso(g):
    p_out = len(list(g.graph.successors('P1')))
    rep_in = len(list(g.graph.predecessors('Rep')))
    size = g.graph.size()
    return p_out, rep_in, size


def get_edit_distance(population, ind):
    ed1 = np.zeros(population.shape[0])
    for i in range(population.shape[0]):
        ed1[i] = nx.graph_edit_distance(population[i][0].graph, population[ind][0].graph)
    return ed1


promo_node = 'P1'
num_dict = {1: 10, 2: 20}
min_dose = 75
max_dose = 75
dose_interval = 5

num_circuits = sum(list(num_dict.values()))

np.random.seed(0)
n_gen = 20
prob = 0.9

population = sampling(promo_node, num_dict, min_dose, max_dose, dose_interval)
obj = np.asarray([-simulate(g[0])/Ref[g[0].promo_node]['on'] for g in population])

obj_min = np.zeros(n_gen+1)
circuit_min = []

geno = np.zeros(n_gen+1)
pheno = np.zeros_like(geno)
entropy = np.zeros_like(geno)
iso = np.zeros_like(geno)
ed_mean = np.zeros_like(geno)
ed_max = np.zeros_like(geno)
ed_median = np.zeros_like(geno)

geno[0] = compare(population)
n, bins = grouping(obj, 0.5)
entropy[0] = get_entropy(n)
pheno[0] = len(set(n))
iso[0] = len(np.unique(np.asarray([pseudo_iso(g[0]) for g in population]), axis=0))
ind_min = np.argmin(obj)
ed1 = get_edit_distance(population, ind_min)
ed_mean[0] = np.mean(ed1)
ed_max[0] = max(ed1)
ed_median[0] = np.median(ed1)

ind_min = np.argmin(obj)
obj_min[0] = obj[ind_min]
circuit_min.append(population[ind_min])

for gen in range(n_gen):
    if np.random.uniform() < 0.5:
        children = crossover(None, population)
        mutate(problem, children, 1.)
        obj_children = np.asarray([-simulate(g[0]) / Ref[g[0].promo_node]['on'] for g in children])
        obj = np.append(obj, obj_children)
        population = np.vstack((population, children))
        S = np.lexsort([obj])
        obj = obj[S[:num_circuits]]
        population = population[S[:num_circuits], :]

    geno[gen + 1] = compare(population)
    n, bins = grouping(obj, 0.5)
    entropy[gen + 1] = get_entropy(n)
    pheno[gen + 1] = len(set(n))
    iso[gen + 1] = len(np.unique(np.asarray([pseudo_iso(g[0]) for g in population]), axis=0))
    ind_min = np.argmin(obj)
    ed1 = get_edit_distance(population, ind_min)
    ed_mean[gen + 1] = np.mean(ed1)
    ed_max[gen + 1] = max(ed1)
    ed_median[gen + 1] = np.median(ed1)

    obj_min[gen + 1] = obj[ind_min]
    circuit_min.append(population[ind_min])

generations = np.arange(n_gen+1)

plt.figure(figsize=(15,10))
plt.plot(generations, geno, label='Genotype')
plt.plot(generations, pheno, label='Phenotype')
plt.plot(generations, iso, label='Pseudo-isomorphism')
plt.plot(generations, ed_mean, label='Mean edit distance')
plt.plot(generations, ed_max, label='Max edit distance')
plt.plot(generations, ed_median, label='Median edit distance')
plt.grid(True)
plt.legend(fontsize=15)

a = 1