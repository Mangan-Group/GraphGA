from GA4Graph import *
from scipy.integrate import odeint
from load_files import *
from get_system_equations import system_equations
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.optimize import minimize

def simulate(topology, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_on = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('on', topology,))[-1, -1]
    return rep_on

class Amplifier(ElementwiseProblem):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.promo_node = 'P1'
        self.max_part = 2
        self.min_dose = int(150/self.max_part)
        self.max_dose = int(150/self.max_part)
        self.dose_interval = 5
        self.inhibitor = False
        self.num_circuit = 3720
        self.n_gen = 1
        self.X = []
        self.F = []

    def _evaluate(self, x, out, *args, **kwargs):
        rep_on = simulate(x[0])
        out['F'] = - rep_on / Ref[x[0].promo_node]['on']
        self.X.append(x[0])
        self.F.append(out["F"])

problem = Amplifier(n_var=1, n_obj=1, n_ieq_constr=0)

algorithm = GA(pop_size=problem.num_circuit,
                  sampling=MySampling(),
                  crossover=MyCrossover(prob=1.),
                  mutation=MyMutation(prob=1.),
                  eliminate_duplicates=MyDuplicateElimination())
                  #   eliminate_duplicates=None)

res = minimize(problem,
               algorithm,
               ('n_gen', problem.n_gen),
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

# all_F = np.array(problem.F).reshape(-1, 50).T
#
# plt.figure(figsize=(20, 10))
# plt.boxplot(all_F)
# plt.xlabel('Iterations', fontsize=15)
# plt.ylabel('Values of objective function', fontsize=15)

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