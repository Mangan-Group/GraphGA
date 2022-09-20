from GA4Graph import *
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from get_test_case import signal_conditioner

test_case = signal_conditioner()

algorithm = NSGA2(pop_size=test_case.num_circuit,
                  sampling=MySampling(),
                  crossover=MyCrossover(prob=1.0),
                  mutation=MyMutation(probs=0.9),
                  eliminate_duplicates=MyDuplicateElimination())
                  #   eliminate_duplicates=None)

problem = MyProblem(promo_node=test_case.promo_node, max_part=test_case.max_part,
                    min_dose=test_case.min_dose, max_dose=test_case.max_dose, dose_interval=test_case.dose_interval,
                    inhibitor=test_case.inhibitor, func=test_case.objective, constr=test_case.constr,
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

a=1