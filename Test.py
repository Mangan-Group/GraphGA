import copy

from graph_operations import *
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.sampling import Sampling
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.duplicate import ElementwiseDuplicateElimination
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.optimize import minimize

def amplifier(g):
    rep_off, rep_on = g.simulate()
    ON_ratio = rep_on / Ref[g.promo_node]['on']
    FIn = rep_on / rep_off
    return np.array([-ON_ratio, -FIn])

def sigcond(g):
    rep_off, rep_on = g.simulate()
    ON_ratio = rep_on / Ref[g.promo_node]['on']
    OFF_ratio = rep_off / Ref[g.promo_node]['off']
    return np.array([-ON_ratio, OFF_ratio])

class MyProblem(ElementwiseProblem):
    def __init__(self, promo_node='P1', max_part=3, max_dose=75, min_dose=10, inhibitor=False, func=amplifier):
        super().__init__(n_var=1, n_obj=2, n_ieq_constr=0)
        self.promo_node = promo_node
        self.max_part = max_part
        self.max_dose = max_dose
        self.min_dose = min_dose
        self.inhibitor = inhibitor
        self.func = func
        self.X = []
        self.F = []

    def _evaluate(self, x, out, *args, **kwargs):
        out["F"] = self.func(x[0])
        self.X.append(x)
        self.F.append(out["F"])

# def amplifier(g):
#     rep_off, rep_on = g.simulate()
#     ON_ratio = rep_on/ Ref[g.promo_node]['on']
#     return -round(ON_ratio)
#
# class MyProblem(ElementwiseProblem):
#     def __init__(self, promo_node='P1', max_part=3, max_dose=75, min_dose=10, inhibitor=False, func=amplifier):
#         super().__init__(n_var=1, n_obj=1, n_ieq_constr=0)
#         self.promo_node = promo_node
#         self.max_part = max_part
#         self.max_dose = max_dose
#         self.min_dose = min_dose
#         self.inhibitor = inhibitor
#         self.func = func
#         self.X = []
#         self.F = []
#
#     def _evaluate(self, x, out, *args, **kwargs):
#         out["F"] = self.func(x[0])
#         self.X.append(x)
#         self.F.append(out["F"])

class MySampling(Sampling):
    def _do(self, problem, num_circuit, **kwargs):
        X = np.full((num_circuit, 1), None, dtype=object)
        X[:, 0] = sample_circuit(problem.promo_node, num_circuit, problem.max_part, problem.max_dose, problem.min_dose, problem.inhibitor)
        return X

class MyCrossover(Crossover):
    def __init__(self):
        super().__init__(n_parents=2, n_offsprings=2, prob=1.0)

    def _do(self, problem, X, **kwargs):
        # The input of has the following shape (n_parents, n_matings, n_var)
        _, n_matings, n_var = X.shape
        Y = np.full_like(X, None, dtype=object)
        for k in range(n_matings):
            # get the first and the second parent
            parent1, parent2 = X[0, k, 0], X[1, k, 0]
            Y[0, k, 0], Y[1, k, 0] = crossover_naive(parent1, parent2)
        # Y = copy.deepcopy(X)
        return Y

class MyMutation(Mutation):
    def __init__(self):
        super().__init__()

    def _do(self, problem, X, **kwargs):
        for i in range(len(X)):
            r = np.random.uniform(0, 1)
            if r < 0.3:
                mutate_node_num(X[i, 0], problem.max_part)
            elif r < 0.6:
                mutate_node_type(X[i, 0], problem.min_dose, problem.max_dose)
                # mutate_edge(X[i, 0])
            elif r < 0.9:
                mutate_dose(X[i, 0], problem.min_dose, problem.max_dose)
        return X

class MyDuplicateElimination(ElementwiseDuplicateElimination):
    def is_equal(self, x1, x2):
        return compare_circuit(x1.X[0], x2.X[0])

algorithm = NSGA2(pop_size=100,
                  sampling=MySampling(),
                  crossover=MyCrossover(),
                  mutation=MyMutation(),
                  eliminate_duplicates=MyDuplicateElimination())

# problem = MyProblem(promo_node='P1', max_part=3, max_dose=75, min_dose=10, inhibitor=True, func=sigcond)
problem=MyProblem()
res = minimize(problem,
               algorithm,
               ('n_gen', 10),
               seed=1,
               save_history=True)

X, F = res.opt.get("X", "F")
X = X[:, 0]
for i in range(len(X)):
    X[i].check_valid()
    if X[i].valid == 0:
        print(X[i].valid)

hist = res.history
n_evals = []
hist_F = []
hist_X = []
for algo in hist:
    # store the number of function evaluations
    n_evals.append(algo.evaluator.n_eval)
    # retrieve the optimum from the algorithm
    opt = algo.opt
    hist_F.extend(opt.get("F"))
    hist_X.append(opt.get('X'))
all_F = np.array(problem.F).reshape(-1, 50).T

plt.figure(figsize=(20,10))
plt.boxplot(all_F)
plt.xlabel('Iterations', fontsize=15)
plt.ylabel('Values of objective function', fontsize=15)

# from pymoo.indicators.hv import Hypervolume
# approx_ideal = F.min(axis=0)
# approx_nadir = F.max(axis=0)
# metric = Hypervolume(ref_point= np.array([1.1, 1.1]),
#                      norm_ref_point=False,
#                      zero_to_one=True,
#                      ideal=approx_ideal,
#                      nadir=approx_nadir)
# hv = [metric.do(_F) for _F in hist_F]
# plt.figure(figsize=(7, 5))
# plt.plot(n_evals, hv,  color='black', lw=0.7, label="Avg. CV of Pop")
# plt.scatter(n_evals, hv,  facecolor="none", edgecolor='black', marker="p")
# plt.title("Convergence")
# plt.xlabel("Function Evaluations")
# plt.ylabel("Hypervolume")
# #
# from pymoo.util.running_metric import RunningMetricAnimation
# running = RunningMetricAnimation(delta_gen=5, key_press=False,)
# for algorithm in res.history:
#     running.do(0, algorithm)

a = 1