from GA4Graph import *
from scipy.integrate import odeint
from multiprocessing import Pool
from scipy.signal import find_peaks, peak_prominences
from load_files_pop import *
from get_system_equations_pop import system_equations_DsRed_pop
from pymoo.algorithms.soo.nonconvex.ga import GA
from pymoo.optimize import minimize

def simulate_cell(topology, Z_list, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_on = odeint(system_equations_DsRed_pop, np.zeros(topology.num_states * 2), t, args=('on', Z_list, topology,))[:, -1]
    return rep_on

def simulate_pop(topology, Z, num_processes, max_time=42):
    nc = len(Z)
    zipped_args = list(zip([topology]*nc, Z, [max_time]*nc))
    with Pool(num_processes) as pool:
        pop_rep_on = pool.starmap(
            simulate_cell,
            zipped_args,
        )

    pool.close()
    pool.join() 

    rep_on_mean = [np.mean(k) for k in zip(*pop_rep_on)]

    rep_on_rel = [i/Ref_pop20[topology.promo_node]['on'] for i in rep_on_mean]

    return rep_on_rel 


class Pulse(ElementwiseProblem):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.promo_node = 'P1'
        self.max_part = 2
        self.min_dose = 5
        self.max_dose = int(150/self.max_part)
        self.dose_interval = 5
        self.inhibitor = True
        self.num_circuit = 700
        self.n_gen = 5
        self.X = []
        self.F = []

    def _evaluate(self, x, out, *args, **kwargs):
        rep_on_rel = simulate_pop(x[0])
        peak_rel = np.max(rep_on_rel)
        peaks, _ = find_peaks(rep_on_rel, prominence= 0.1*peak_rel)
        prominence_rel = peak_prominences(rep_on_rel, peaks)[0]
        out["F"] = [- peak_rel, - prominence_rel]
        self.X.append(x[0])
        self.F.append(out["F"])
        out["G"] = 1 - len(peaks)

problem = Pulse(n_var=1, n_obj=2, n_ieq_constr=1)

algorithm = GA(pop_size=problem.num_circuit,
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

    a = 1