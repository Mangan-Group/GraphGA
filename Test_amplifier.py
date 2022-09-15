from GA4Graph import *
from pymoo.algorithms.moo.nsga2 import NSGA2
# from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.optimize import minimize
from get_test_case import amplifier

test_case = amplifier()
algorithm = NSGA2(pop_size=test_case.num_circuit,
                  sampling=MySampling(),
                  crossover=MyCrossover(prob=1.0),
                  mutation=MyMutation(prob=0.9),
                  eliminate_duplicates=MyDuplicateElimination())
                  #   eliminate_duplicates=None)

problem = MyProblem(promo_node=test_case.promo_node, max_part=test_case.max_part,
                    min_dose=test_case.min_dose, max_dose=test_case.max_dose, dose_interval=test_case.dose_interval,
                    inhibitor=test_case.inhibitor, func=test_case.objective,
                    n_var=1, n_obj=test_case.n_obj, n_ieq_constr=test_case.n_ieq_constr)
res = minimize(problem,
               algorithm,
               ('n_gen', test_case.n_gen),
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