from GA4Graph import *
from scipy.integrate import odeint
from load_files import *
from get_system_equations import system_equations
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

def simulate(topology, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_off = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('off', topology,))[-1, -1]
    rep_on = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('on', topology,))[-1, -1]
    return rep_off, rep_on

class SignalConditioner(ElementwiseProblem):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.promo_node = 'P1'
        self.max_part = 2
        self.min_dose = 10
        self.max_dose = int(150/self.max_part)
        self.dose_interval = 5
        self.inhibitor = True
        self.num_circuit = 50
        self.n_gen = 1
        self.X = []
        self.F = []

    def _evaluate(self, x, out, *args, **kwargs):
        rep_off, rep_on = simulate(x[0])

        ON_rel = rep_on / Ref[x[0].promo_node]['on']
        FI_rel = (rep_on/rep_off)/(Ref[x[0].promo_node]['fi'])
        out['F'] = [- ON_rel, - FI_rel]
        self.X.append(x[0])
        self.F.append(out["F"])

problem = SignalConditioner(n_var=1, n_obj=2, n_ieq_constr=0)

algorithm = NSGA2(pop_size=problem.num_circuit,
                  sampling=MySampling(),
                  crossover=MyCrossover(prob=1.0),
                  mutation=MyMutation_full(prob=0.9),
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

a=1