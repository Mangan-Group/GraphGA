from GA import *
from define_problem import *
from rankcrowding import *


# Calculate the convergence by determining when the final hypervolume output was first found
def first_seen(progression):

    gen_num = 0

    # Starting from the end, return the generation where the output first does not match the final output
    for gen in reversed(progression):
        if progression[-1] != gen:
            return len(progression) - gen_num
        gen_num += 1

    return 0


def simulate(topology, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_off = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('off', topology,))[-1, -1]
    rep_on = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('on', topology,))[-1, -1]
    return rep_off, rep_on


def func(toplogy):
    rep_off, rep_on = simulate(toplogy)
    ON_rel = rep_on / Ref[toplogy.promo_node]['on']
    FI_rel = (rep_on / rep_off) / (Ref[toplogy.promo_node]['fi'])
    return [-ON_rel, -FI_rel]


# Define function by taking in necessary inputs to define the problem and set the seed
def full_sim(mut_rate, cov_rate, promo_node, num_dict, max_part, min_dose, max_dose, dose_interval, inhibitor, seed):

    # Inhibitor must be true in the signal conditioner case
    inhibitor = True

    # Set seed based on inputted number to allow for reproducibility
    np.random.seed(seed)

    problem = Problem(promo_node, max_part, min_dose, max_dose, dose_interval, inhibitor, func)

    # Generate a population using the code or use a set population
    population = sampling(problem.promo_node, num_dict, problem.min_dose, problem.max_dose, problem.dose_interval, problem.inhibitor)

    # with open("init_pop_inhib_fixed.pkl", "rb") as fid:
    #      population = pickle.load(fid)

    # Define a reference point and hypervolume calculator
    from pymoo.indicators.hv import HV
    ref_point = np.array([0, 0])
    ind = HV(ref_point=ref_point)

    # Store the progression of hypervolumes
    hvs = []

    # From this point, code continues exactly as the original SigCond_manual.py file
    nds = RankAndCrowding()
    num_circuits = len(population)
    n_gen = 50

    all_obj = []
    obj = [problem.func(g[0]) for g in population]
    all_obj.append(obj)
    obj = np.asarray(obj)

    for gen in range(n_gen):
        _, rank_dict = nds.do(obj, num_circuits, return_rank=True)
        if np.random.uniform() < cov_rate:
            children = crossover(population, obj, rank_dict)
        else:
            children = deepcopy(population)
        mutate(problem, children, mut_rate, dose=True)
        obj_children = [problem.func(g[0]) for g in children]
        all_obj.append(obj_children)
        obj_children = np.asarray(obj_children)

        obj = np.vstack((obj, obj_children))
        population = np.vstack((population, children))

        S = nds.do(obj, num_circuits)
        obj = obj[S]
        population = population[S, :]

        # Add hypervolume to the progression of hypervolumes
        hvs.append(ind(obj))

    # Return the final population and convergence
    return [obj, first_seen(hvs)]
