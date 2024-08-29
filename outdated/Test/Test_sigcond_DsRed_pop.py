from GA4Graph import *
from scipy.integrate import odeint
from multiprocessing import Pool
from load_files_pop import *
from get_system_equations_pop import system_equations_DsRed_pop
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

def simulate_cell(topology, Z_list, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_off = odeint(system_equations_DsRed_pop, np.zeros(topology.num_states * 2), t, args=('off', Z_list, topology,))[-1, -1]
    rep_on = odeint(system_equations_DsRed_pop, np.zeros(topology.num_states * 2), t, args=('on', Z_list, topology,))[-1, -1]
    return rep_off, rep_on

def simulate_pop(topology, Z, num_processes, max_time=42):
    nc = len(Z)
    zipped_args = list(zip([topology]*nc, Z, [max_time]*nc))
    with Pool(num_processes) as pool:
        results = pool.starmap(
            simulate_cell,
            zipped_args,
        )

    pool.close()
    pool.join()

    pop_rep_off, pop_rep_on = zip(*results)
    rep_off_mean = np.mean(pop_rep_off)
    rep_on_mean = np.mean(pop_rep_on)

    return rep_off_mean, rep_on_mean

class SignalConditioner_DsRed(ElementwiseProblem):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.promo_node = 'P1'
        self.max_part = 2
        self.min_dose = 5
        self.max_dose = int(150/self.max_part)
        self.dose_interval = 5
        self.inhibitor = True
        self.num_circuit = 50
        self.n_gen = 1
        self.X = []
        self.F = []

    def _evaluate(self, x, out, *args, **kwargs):
        rep_off, rep_on = simulate_pop(x[0], Z_20, 28)

        ON_rel = rep_on / Ref_pop20[x[0].promo_node]['on']
        FI_rel = (rep_on/rep_off)/(Ref_pop20[x[0].promo_node]['fi'])
        out['F'] = [- ON_rel, - FI_rel]
        self.X.append(x[0])
        self.F.append(out["F"])

problem = SignalConditioner_DsRed(n_var=1, n_obj=2, n_ieq_constr=0)

algorithm = NSGA2(pop_size=problem.num_circuit,
                  sampling=MySampling(),
                  crossover=MyCrossover(prob=1.0),
                  mutation=MyMutation_full(prob=0.9),
                  eliminate_duplicates=MyDuplicateElimination())
                  #   eliminate_duplicates=None)

if __name__ == '__main__':

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