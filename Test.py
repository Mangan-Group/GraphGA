from graph_operations import *
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.sampling import Sampling
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.optimize import minimize

def amplifier(g):
    rep_off, rep_on = g.simulate()
    ON_ratio = rep_on / Ref[g.promo_node]['on']
    FIn = rep_on / rep_off
    return -np.array([ON_ratio, FIn])

class MyProblem(ElementwiseProblem):
    def __init__(self, promo_node='P1', max_part=3, max_dose=75, min_dose=10, inhibitor=False,func=amplifier):
        super().__init__(n_var=1, n_obj=2, n_ieq_constr=0)
        self.promo_node = promo_node
        self.max_part = max_part
        self.max_dose = max_dose
        self.min_dose = min_dose
        self.inhibitor = inhibitor
        self.func = func

    def _evaluate(self, x, out, *args, **kwargs):
        out["F"] = self.func(x[0])

# def amplifier(g):
#     rep_off, rep_on = g.simulate()
#     ON_ratio = rep_on/ Ref[g.promo_node]['on']
#     return -ON_ratio
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
#     def _evaluate(self, x, out, *args, **kwargs):
#         out["F"] = self.func(x[0])

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
        # Offsprings
        Y = np.full_like(X, None, dtype=object)
        # for each mating provided
        for k in range(n_matings):
            # get the first and the second parent
            parent1, parent2 = X[0, k, 0], X[1, k, 0]
            Y[0, k, 0], Y[1, k, 0] = crossover_naive(parent1, parent2)
        return Y

class MyMutation(Mutation):
    def __init__(self):
        super().__init__()

    def _do(self, problem, X, **kwargs):
        return X

algorithm = NSGA2(pop_size=50,
                  sampling=MySampling(),
                  crossover=MyCrossover(),
                  mutation=MyMutation(),
                  eliminate_duplicates=False)


res = minimize(MyProblem(),
               algorithm,
               ('n_gen', 5),
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
for algo in hist:
    # store the number of function evaluations
    n_evals.append(algo.evaluator.n_eval)
    # retrieve the optimum from the algorithm
    opt = algo.opt
    hist_F.append(opt.get("F"))

plt.figure()
for i in range(len(hist_F)):
    plt.scatter(hist_F[i][:, 0], hist_F[i][:, 1])
plt.scatter(F[:, 0], F[:, 1], c='r')

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