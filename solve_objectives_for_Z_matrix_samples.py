import numpy as np
import pickle
import pandas as pd
from multiprocessing import Pool
from scipy.stats import sem

def select_Z_matrix_sampling_topologies(
        results_path: str,
        obj_threshold: str=None
):
    """Collects the circuits and objective
    function values for selected topologies
    from population model GA run for simulation
    with 10 manifestations of the Z matrix. If
    obj_threshold is specified, it is assumed that
    the GA run was for single objective optimization
    (based on files loaded). Otherwise, it is assumed
    that multi-objective optimization was used, and the
    full pareto front is loaded."""

    # load the (unique) circuit objects and objectives from
    # the single objective GA run
    if obj_threshold:
        with open(results_path+"unique_circuits.pkl", "rb") as fid:
            unique_circuits = pickle.load(fid)
        with open(results_path+"unique_objectives.pkl", "rb") as fid:
            unique_objectives = pickle.load(fid)
        
        # subset objectives to only keep those with objective
        #above threshold (assuming maximizing objective, which
        # results in minization of -1*objective in GA).
        unique_objectives = unique_objectives*-1
        index_selected_objectives_threshold = np.argwhere(
            (unique_objectives[:, 0] >=obj_threshold)
        ).flatten()
        selected_objectives = unique_objectives[
            index_selected_objectives_threshold
        ]
        # subset circuits with the corresponding
        # objectives
        selected_circuits = unique_circuits[
            index_selected_objectives_threshold
        ]
    
    # load the pareto front circuit objects and objectives
    # from the multi-objective GA run
    else:
        with open(results_path+"final_population.pkl", "rb") as fid:
            pareto_front_circuits = pickle.load(fid)
        
        objectives_fname = "final_objectives_df.pkl"
        pareto_front_objectives = pd.read_pickle(
            results_path+objectives_fname
        )

        # drop duplicate objective function values and
        # corresponding circuits
        selected_objectives = (
            pareto_front_objectives.drop_duplicates()
        )
        selected_circuits = pareto_front_circuits[
            selected_objectives.index
        ]

    return selected_circuits, selected_objectives


def solve_single_objective_for_topology(
        problem: object,
        topology: np.ndarray, # this is an array because of formatting in GA: 1D array of length 1 with a single topology object
        Z_mat_list: list,
        Ref_list: list,
        idx: np.ndarray
):
    """***Single objective optimization***
    Simulates the ODEs for the given
    topology and calculates the objective
    function value for 10 manifestations of
    the Z matrix, and calculates metrics
    associated with the objective function 
    values."""
    
    obj_list = []
    # simulate the ODEs for the given topology
    # and calculate the objective function value
    # for each Z matrix in the list
    for i, Z_mat in enumerate(Z_mat_list):
        problem.Z = Z_mat
        problem.ref = Ref_list[i]
        obj = abs(problem.func(topology[0]))
        obj_list.append(obj)

    # save the objective function values and associated
    # metrics in a dictionary
    topology_dict = {}
    # get name of objective
    label = problem.obj_labels[0]
    # calculate range of objective function values across
    # simulations
    topology_dict[label+"_range"] = max(obj_list) - min(obj_list)
    # calculate average of objective function values across
    # simulations
    topology_dict[label+"_mean"] = np.mean(obj_list)
    # calculate minimum and maximum of objective function 
    # values across simulations
    topology_dict[label+"_min_max"] = [min(obj_list), max(obj_list)]
    # calculate standard error of objective function values
    # across simulations
    topology_dict[label+"_std_error"] = sem(obj_list, ddof=1)
    # save list of objective function values
    topology_dict[label+"_list"] = obj_list
    # save topology object
    topology_dict["topology"] = topology[0]
    print("topology "+str(idx[0])+" done")

    return topology_dict

def solve_multi_objective_for_topology(
        problem: object, 
        topology: np.ndarray, 
        Z_mat_list: list, 
        idx: np.ndarray
):
    """***Multi-objective optimization***
    Simulates the ODEs for the given
    topology and calculates the objective
    function values for 10 manifestations of
    the Z matrix, and calculates metrics
    associated with the objective function 
    values."""

    obj_list = []
    # simulate the ODEs for the given topology
    # and calculate the objective function values
    # for each Z matrix in the list
    for Z_mat in Z_mat_list:
        problem.Z = Z_mat
        obj_negative = problem.func(topology[0])
        obj = [abs(i) for i in obj_negative]
        obj_list.append(obj)

    topology_dict = {}
    obj_list = np.asarray(obj_list)
    # save the objective function values and
    # associated metrics in a dictionary
    # this is the same process as in 
    # select_Z_matrix_sampling_topologies, and
    # is done individually for each of the 
    # objectives
    for i in range(len(problem.obj_labels)):
        label = problem.obj_labels[i]
        topology_dict[label+"_range"] = max(obj_list[:, i]) - min(obj_list[:, i])
        topology_dict[label+"_mean"] = np.mean(obj_list[:, i])
        topology_dict[label+"_min_max"] = [min(obj_list[:, i]), max(obj_list[:, i])]
        topology_dict[label+"_std_error"] = sem(obj_list[:, i], ddof=1)
        topology_dict[label+"_list"] = obj_list[:, i]

    topology_dict["topology"] = topology[0]
    print("topology "+str(idx[0])+" done")

    return topology_dict


def solve_all_topology_objectives(
        problem: object, 
        topologies: np.ndarray, 
        Z_mat_list: list, 
        Ref_list: list
):
    """Solves objective(s) for all
    selected topologies and saves metrics."""    

    # get topology index list for printing
    # status
    num_topologies = len(topologies)
    idx = np.arange(0, num_topologies, 1).reshape(-1, 1)
    # choose appropriate solver function based
    # on number of objective labels in problem
    # object
    if len(problem.obj_labels) == 1:
        solver = solve_single_objective_for_topology
    else:
        solver = solve_multi_objective_for_topology
    
    # format function arguments for runnning with parallelization
    zipped_args = list(zip([problem]*num_topologies, topologies, [
        Z_mat_list]*num_topologies, [Ref_list]*num_topologies, idx
    ))
    with Pool(problem.num_processes) as pool:
        topology_dict_list = pool.starmap(solver, zipped_args)
    
    # save list of dictionaries for each topology as a df
    Z_matrix_sampling = pd.DataFrame(topology_dict_list)    

    return Z_matrix_sampling
    

def get_objective_errors(
        Z_matrix_sampling: pd.DataFrame,
        CI_metrics: list
):
    """Calculates the maximum value of the
    specified metrics (range and/or standard 
    error) across all simulated topologies for
    plotting CI."""
    
    CI_metric_maxes = []
    # calculate maximum value for each metric in the
    # given list
    for metric in CI_metrics:
        CI_metric_list = Z_matrix_sampling[metric].tolist()
        CI_metric_max = max(CI_metric_list)
        CI_metric_maxes.append(CI_metric_max)
    
    return CI_metric_maxes
