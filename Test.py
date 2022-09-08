from GA import *
from pymoo.core.problem import ElementwiseProblem
from pymoo.core.sampling import Sampling
from pymoo.core.crossover import Crossover
from pymoo.core.mutation import Mutation

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize

class MyProblem(ElementwiseProblem):
    def __init__(self, promo_node='P1', max_part=3, max_dose=75, min_dose=10, inhibitor=False):
        super().__init__(n_var=1, n_obj=2, n_ieq_constr=0)
        self.promo_node = promo_node
        self.max_part = max_part
        self.max_dose = max_dose
        self.min_dose = min_dose
        self.inhibitor = inhibitor

    def _evaluate(self, x, out, *args, **kwargs):
        out["F"] = -amplifier_obj(x[0])

class MySampling(Sampling):
    def _do(self, problem, num_circuit, **kwargs):
        X = np.full((num_circuit, 1), None, dtype=object)
        X[:, 0] = sample_circuit(problem.promo_node, num_circuit, problem.max_part, problem.max_dose, problem.min_dose, problem.inhibitor)
        return X

class MyCrossover(Crossover):
    def __init__(self):
        # define the crossover: number of parents and number of offsprings
        super().__init__(n_parents=2, n_offsprings=2, prob=1.0)

    def _do(self, problem, X, **kwargs):
        # print(X.shape)
        # The input of has the following shape (n_parents, n_matings, n_var)
        _, n_matings, n_var = X.shape
        # Offsprings
        Y = np.full_like(X, None, dtype=object)
        # for each mating provided
        for k in range(n_matings):
            # get the first and the second parent
            parent1, parent2 = X[0, k, 0], X[1, k, 0]
            # print(type(parent1))
            Y[0, k, 0], Y[1, k, 0] = crossover(parent1, parent2, naive=True)

        return Y

class MyMutation(Mutation):
    def __init__(self):
        super().__init__()

    def _do(self, problem, X, **kwargs):
        return X
        # # for each individual
        # for i in range(len(X)):
        #
        #     r = np.random.uniform(0, 1)
        #
        #     # with a probabilty - change the node
        #     if r < 0.4:
        #         old_node = np.random.choice(X[i].part_list)
        #         if old_node[0] == 'Z':
        #             node_avail = list(set(tf_list).difference(set(X[i].part_list)))
        #         else:
        #             node_avail = list(set(in_list).difference(set(X[i].part_list)))
        #         new_node = np.random.choice(node_avail)
        #
        #         X[i] = Topo(switch_naive(X[i], old_node, new_node), X[i].dose, X[i].promo_node)
        #
        #     # also with a probabilty - remove (or add) a node
        #     elif r < 0.8:
        #         if len(X[i].part_list) > 1:
        #             old_node = np.random.choice(X[i].part_list)
        #             X[i].graph.remove_node(old_node)
        #             X[i].dose.pop(old_node)
        #             X[i] = Topo(list(X[i].graph.edges), X[i].dose, X[i].promo_node)
        #         else:
        #             pass
        #
        # return X

algorithm = NSGA2(pop_size=100,
                  sampling=MySampling(),
                  crossover=MyCrossover(),
                  mutation=MyMutation(),
                  eliminate_duplicates=False)

res = minimize(MyProblem(),
               algorithm,
               ('n_gen', 10),
               seed=1,
               save_history=True)

a = 1