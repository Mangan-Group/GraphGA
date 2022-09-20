from graph_operations import *
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.sampling import Sampling
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation
from pymoo.core.duplicate import ElementwiseDuplicateElimination


class MySampling(Sampling):
    def _do(self, problem, num_circuit, **kwargs):
        X = np.full((num_circuit, 1), None, dtype=object)
        X[:, 0] = sample_circuit(problem.promo_node, num_circuit, problem.max_part, problem.min_dose, problem.max_dose, problem.dose_interval, problem.inhibitor)
        return X

class MyCrossover(Crossover):
    def __init__(self, **kwargs):
        super(MyCrossover, self).__init__(n_parents=2, n_offsprings=2, **kwargs)

    def _do(self, problem, X, **kwargs):
        # The input of has the following shape (n_parents, n_matings, n_var)
        _, n_matings, n_var = X.shape
        Y = np.full_like(X, None, dtype=object)
        for k in range(n_matings):
            # get the first and the second parent
            parent1, parent2 = X[0, k, 0], X[1, k, 0]
            Y[0, k, 0], Y[1, k, 0] = crossover_structure(parent1, parent2)
        # Y = copy.deepcopy(X)
        return Y

class MyMutation(Mutation):
    def __init__(self, prob=0.9):
        super().__init__()
        self.prob = float(prob)/3

    def _do(self, problem, X, **kwargs):
        for i in range(len(X)):
            r = np.random.uniform(0, 1)
            if r < self.prob:
                mutate_node_num(X[i, 0], problem.max_part, problem.min_dose, problem.max_dose, problem.dose_interval, problem.inhibitor)
            elif r < self.prob*2:
                mutate_node_type(X[i, 0], problem.min_dose, problem.max_dose, problem.dose_interval)
                # mutate_edge(X[i, 0])
            elif r < self.prob*3:
                mutate_dose(X[i, 0], problem.min_dose, problem.max_dose, problem.dose_interval)
        return X

class MyDuplicateElimination(ElementwiseDuplicateElimination):
    def is_equal(self, x1, x2):
        return compare_circuit(x1.X[0], x2.X[0])

class MyProblem(ElementwiseProblem):
    def __init__(self, promo_node='P1', max_part=3, min_dose=10, max_dose=75, dose_interval=5, inhibitor=False, func=None, constr=None, **kwargs):
        super().__init__(**kwargs)
        self.promo_node = promo_node
        self.max_part = max_part
        self.min_dose = min_dose
        self.max_dose = max_dose
        self.dose_interval = dose_interval
        self.inhibitor = inhibitor
        self.func = func
        self.constr = constr
        self.X = []
        self.F = []

    def _evaluate(self, x, out, *args, **kwargs):
        out["F"] = self.func(x[0])
        self.X.append(x)
        self.F.append(out["F"])
        if self.n_ieq_constr > 0:
            out["G"] = self.constr(x[0])