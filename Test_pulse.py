from GA4Graph import *
from scipy.integrate import odeint
from scipy.signal import find_peaks
from load_files import *
from get_system_equations import system_equations_DsRed
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.optimize import minimize


def simulate(topology, max_time=48):
    t = np.arange(0, max_time + 1, 1)
    rep_on = odeint(system_equations_DsRed, np.zeros(topology.num_states * 2), t, args=('on', topology,))[:, -1]
    return rep_on

class Pulse(ElementwiseProblem):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.promo_node = 'P1'
        self.max_part = 2
        self.min_dose = 10
        self.max_dose = int(150/self.max_part)
        self.dose_interval = 5
        self.inhibitor = True
        self.num_circuit = 150
        self.n_gen = 2
        self.X = []
        self.F = []

    def _evaluate(self, x, out, *args, **kwargs):
        rep_on = simulate(x[0])
        out["F"] = -np.max(rep_on) / Ref[x[0].promo_node]['on']
        self.X.append(x[0])
        self.F.append(out["F"])
        out["G"] = 1 - len(find_peaks(rep_on, prominence=50.)[0])

problem = Pulse(n_var=1, n_obj=1, n_ieq_constr=1)

algorithm = GA(pop_size=problem.num_circuit,
                  sampling=MySampling(),
                  crossover=MyCrossover(prob=1.0),
                  mutation=MyMutation(prob=0.9),
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

a = 1