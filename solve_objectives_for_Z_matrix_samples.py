import numpy as np
import pickle
import pandas as pd
from multiprocessing import Pool
from scipy.stats import sem

def select_Z_matrix_sampling_topologies(results_path: str, obj_threshold: str=None):

    if obj_threshold:
        with open(results_path+"unique_circuits.pkl", "rb") as fid:
            unique_circuits = pickle.load(fid)
        with open(results_path+"unique_objectives.pkl", "rb") as fid:
            unique_objectives = pickle.load(fid)

        unique_objectives = unique_objectives*-1
        index_selected_objectives_threshold = np.argwhere(
            (unique_objectives[:, 0] >=obj_threshold)
        ).flatten()
        selected_objectives = unique_objectives[
            index_selected_objectives_threshold
        ]
        selected_circuits = unique_circuits[
            index_selected_objectives_threshold
        ]
    
    else:
        with open(results_path+"final_population.pkl", "rb") as fid:
            pareto_front_circuits = pickle.load(fid)
        
        objectives_fname = "final_objectives_df.pkl"
        pareto_front_objectives = pd.read_pickle(
            results_path+objectives_fname
        )
        selected_objectives = (
            pareto_front_objectives.drop_duplicates()
        )
        selected_circuits = pareto_front_circuits[
            selected_objectives.index
        ]


    return selected_circuits, selected_objectives

def solve_single_objective_for_topology(problem: object, topology: list, Z_mat_list: list, Ref_list: list, idx: np.ndarray):
    obj_list = []
    for i, Z_mat in enumerate(Z_mat_list):
        problem.Z = Z_mat
        problem.ref = Ref_list[i]
        obj = abs(problem.func(topology[0]))
        obj_list.append(obj)

    topology_dict = {}
    label = problem.obj_labels[0]
    topology_dict[label+"_range"] = max(obj_list) - min(obj_list)
    topology_dict[label+"_mean"] = np.mean(obj_list)
    topology_dict[label+"_min_max"] = [min(obj_list), max(obj_list)]
    topology_dict[label+"_std_error"] = sem(obj_list, ddof=1)
    topology_dict[label+"_list"] = obj_list
    topology_dict["topology"] = topology[0]
    print("topology "+str(idx[0])+" done")
    return topology_dict

def solve_multi_objective_for_topology(problem: object, topology: list, Z_mat_list: list, idx: np.ndarray):
    obj_list = []
    for Z_mat in Z_mat_list:
        problem.Z = Z_mat
        obj_negative = problem.func(topology[0])
        obj = [abs(i) for i in obj_negative]
        obj_list.append(obj)

    topology_dict = {}
    obj_list = np.asarray(obj_list)
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

def solve_all_topology_objectives(problem, topologies, Z_mat_list, Ref_list):
    num_topologies = len(topologies)
    idx = np.arange(0, num_topologies, 1).reshape(-1, 1)
    if len(problem.obj_labels) == 1:
        solver = solve_single_objective_for_topology
    else:
        solver = solve_multi_objective_for_topology
    zipped_args = list(zip([problem]*num_topologies, topologies, [Z_mat_list]*num_topologies, [Ref_list]*num_topologies, idx))
    with Pool(problem.num_processes) as pool:
        topology_dict_list = pool.starmap(solver, zipped_args)
    
    Z_matrix_sampling = pd.DataFrame(topology_dict_list)    

    return Z_matrix_sampling
    

def get_objective_errors(Z_matrix_sampling: pd.DataFrame, CI_metrics: list):
    
    CI_metric_maxes = []
    for metric in CI_metrics:
        CI_metric_list = Z_matrix_sampling[metric].tolist()
        CI_metric_max = max(CI_metric_list)
        CI_metric_maxes.append(CI_metric_max)
    
    return CI_metric_maxes
