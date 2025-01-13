import numpy as np
import pickle
import pandas as pd
from copy import deepcopy
from plot_search_results import (
    plot_hypervolumes_set,
    plot_hypervolumes_set_vs_combo,
    plot_pareto_front_set,
    plot_pareto_front_set_3D
)

def get_unique_in_a_seed_pareto_circuits(pareto_circuits):

    #add (part, dose) to circuit edges
    circuit_edge_lists = []
    for circuit in pareto_circuits:
        circuit_edges = deepcopy(circuit[0].edge_list)
        for key, val in circuit[0].dose.items():
            circuit_edges.append((key, str(val)))
        circuit_edge_lists.append(circuit_edges)

    # restructure edge_lists: combine first and
    # second edge into one string, e.g., 
    # ("Z1", "Z2") --> ("Z1Z2")
    combo_edges_lists = []
    for edges in circuit_edge_lists:
        edge_combos = [edge[0] + edge[1] for
            edge in edges]
        combo_edges_lists.append(edge_combos)

    # get unique circuits based on combo_edges_lists
    # and index for selecting them from pareto_circuits
    seed_unique_edge_combo = []
    index_list = []
    seen = set()
    for i, combo_list in enumerate(
        combo_edges_lists):
        combo_set = frozenset(combo_list)
        if combo_set not in seen:
            seen.add(combo_set)
            index_list.append(i)
            seed_unique_edge_combo.append(combo_list)

    seed_unique_circuits = pareto_circuits[index_list]
    # print(index_list)

    return seed_unique_edge_combo, seed_unique_circuits, seen

def get_unique_in_each_seed_pareto_circuits(results_path, seed_folder):
    
    seed_unique_edge_combo_list = []
    seed_unique_circuits_list = []
    seen_list = []
    for seed in range(0, 10):
        #import final population circuits
        full_path = results_path + seed_folder + str(seed) + "/"
        with open(full_path + "final_population.pkl", "rb") as fid:
            pareto_circuits = pickle.load(fid)
        (seed_unique_edge_combo, seed_unique_circuits,
         seen) = get_unique_in_a_seed_pareto_circuits(pareto_circuits)
        seed_unique_circuits = seed_unique_circuits.flatten()
        seed_unique_edge_combo_list.append(seed_unique_edge_combo)
        seed_unique_circuits_list.append(seed_unique_circuits)
        seen_list.append(seen)

    # print("unique in each seed: ", seed_unique_edge_combo_list)
    return seen_list, seed_unique_edge_combo_list, seed_unique_circuits_list

def get_common_pareto_circuits(seen_list, seed_unique_edge_combo_list, seed_unique_circuits_list):

    common_edge_combo_set = set.intersection(*seen_list)

    common_edge_combo_list = []
    common_circuits = []
    for i, edge_combo in enumerate(seed_unique_edge_combo_list[0]):
        combo_set = frozenset(edge_combo)
        if combo_set in common_edge_combo_set:
            common_circuits.append(seed_unique_circuits_list[0][i])
            common_edge_combo_list.append(edge_combo)

    return common_edge_combo_set, common_edge_combo_list, common_circuits

def get_uncommon_pareto_circuits(seed_unique_edge_combo_list, common_edge_combo_set, seed_unique_circuits_list):

    #with unique circuits in each seed and common circuits across all 
    # seeds, get circuits not common to all seeds (for each seed)
    all_uncommon_circuit_edge_combo_sets = []
    all_uncommon_circuits_lists = []
    all_uncommon_circuit_edge_combo_lists = []
    for i, edge_combo_list in enumerate(seed_unique_edge_combo_list):
        seed_unique_circuits = seed_unique_circuits_list[i]
        uncommon_circuits_list = []
        uncommon_circuit_edge_combo_list = []
        uncommon_circuits_edge_combo_set = set()
        for j, combo_edges in enumerate(edge_combo_list):
            combo_set = frozenset(combo_edges)
            if combo_set not in common_edge_combo_set:
                uncommon_circuit_edge_combo_list.append(combo_edges)
                uncommon_circuits_edge_combo_set.add(combo_set)
                uncommon_circuits_list.append(seed_unique_circuits[j])
        all_uncommon_circuit_edge_combo_lists.append(uncommon_circuit_edge_combo_list)
        all_uncommon_circuit_edge_combo_sets.append(uncommon_circuits_edge_combo_set)        
        all_uncommon_circuits_lists.append(uncommon_circuits_list)
    
    return (all_uncommon_circuit_edge_combo_sets, all_uncommon_circuit_edge_combo_lists, 
            all_uncommon_circuits_lists)

def get_circuits_unique_to_seed(all_uncommon_circuit_edge_combo_sets, 
                                all_uncommon_circuit_edge_combo_lists,
                                all_uncommon_circuits_lists
):

    # get circuits that only appear once across all seeds
    [l0, l1, l2, l3, l4, l5, l6, l7, l8, l9] = all_uncommon_circuit_edge_combo_sets
    circuits_freq1 = l0 ^ l1 ^ l2 ^ l3 ^ l4 ^ l5 ^ l6 ^ l7 ^ l8 ^ l9
    
    all_uncommon_circuits_edge_combo_lists_flat = [
        edge_lists 
        for seed_lists in all_uncommon_circuit_edge_combo_lists 
        for edge_lists in seed_lists
    ]
    # print(all_uncommon_circuits_lists)
    all_uncommon_circuits_lists_flat = [
        circuit_list
        for seed_lists in all_uncommon_circuits_lists
        for circuit_list in seed_lists
    ]
    # print(all_uncommon_circuits_lists_flat)
    # print(all_uncommon_circuits_lists)
    edge_combos_more_than_one_seed = []
    circuits_more_than_one_seed = []
    for i, edge_list in enumerate(all_uncommon_circuits_edge_combo_lists_flat):
        edge_set = frozenset(edge_list)
        if edge_set not in circuits_freq1:
            edge_combos_more_than_one_seed.append(edge_list)
            circuits_more_than_one_seed.append(all_uncommon_circuits_lists_flat[i])

    circuits_unique_to_seed_list = []
    edge_combos_unique_to_seed_list = []
    for i, uncommon_circuits_edge_combos in enumerate(all_uncommon_circuit_edge_combo_lists):
        edge_combos_unique_to_seed = []
        circuits_unique_to_seed = []
        uncommon_circuits_list = all_uncommon_circuits_lists[i]
        for j, circuit in enumerate(uncommon_circuits_edge_combos):
            circuit_set = frozenset(circuit)
            if circuit_set in circuits_freq1:
                edge_combos_unique_to_seed.append(circuit)
                circuits_unique_to_seed.append(uncommon_circuits_list[j])
        circuits_unique_to_seed_list.append(circuits_unique_to_seed)
        edge_combos_unique_to_seed_list.append(edge_combos_unique_to_seed)
    
    return (circuits_unique_to_seed_list, edge_combos_unique_to_seed_list, 
            edge_combos_more_than_one_seed, circuits_more_than_one_seed)


def get_dose_varied_common_circuits(
        circuits_unique_to_seed_list
):

    all_dose_varied_circuit_edge_combo_sets = []
    all_dose_varied_edge_combo_lists = []
    all_dose_varied_circuits_lists = []
    for i, circuit_list in enumerate(circuits_unique_to_seed_list):
        dose_varied_edge_combo_list = []
        dose_varied_edge_combo_set = set()
        dose_varied_circuit_list = []
        # for each unique to seed circuit
        for circuit in circuit_list:
            circuit_edges = deepcopy(circuit.edge_list)
            # edge combo list without doses
            edge_combo = [edge[0] + edge[1] for
            edge in circuit_edges]
            combo_set = frozenset(edge_combo)
            # only get circuits with unique edges, not
            # considering dose
            if combo_set not in dose_varied_edge_combo_set:
                dose_varied_edge_combo_set.add(combo_set)
                dose_varied_edge_combo_list.append(edge_combo)
                dose_varied_circuit_list.append(circuit)
        all_dose_varied_circuit_edge_combo_sets.append(dose_varied_edge_combo_set)
        all_dose_varied_edge_combo_lists.append(dose_varied_edge_combo_list)
        all_dose_varied_circuits_lists.append(dose_varied_circuit_list)

    (_,
     dose_varied_common_edge_combo_list,
     dose_varied_common_circuits) = get_common_pareto_circuits(
         all_dose_varied_circuit_edge_combo_sets,
         all_dose_varied_edge_combo_lists,
         all_dose_varied_circuits_lists
     )

    return dose_varied_common_edge_combo_list, dose_varied_common_circuits

def compare_pareto_circuits_across_seeds(results_path, seed_folder, selected_seed):

    (seen_list, seed_unique_edge_combo_list, 
     seed_unique_circuits_list) = get_unique_in_each_seed_pareto_circuits(results_path, seed_folder)
    
    (common_edge_combo_set,
     common_edge_combo_list, 
     common_circuits) = get_common_pareto_circuits(seen_list, seed_unique_edge_combo_list, seed_unique_circuits_list)

    (all_uncommon_circuits_sets, 
     all_uncommon_circuit_edge_combo_lists, all_uncommon_circuits_lists) = get_uncommon_pareto_circuits(
         seed_unique_edge_combo_list, common_edge_combo_set, seed_unique_circuits_list    
     )
    
    (circuits_unique_to_seed_list, edge_combos_unique_to_seed_list, 
    edge_combos_more_than_one_seed, circuits_more_than_one_seed) = get_circuits_unique_to_seed(
        all_uncommon_circuits_sets, all_uncommon_circuit_edge_combo_lists,
        all_uncommon_circuits_lists
    )

    (dose_varied_common_edge_combo_list, 
     dose_varied_common_circuits) = get_dose_varied_common_circuits(
         circuits_unique_to_seed_list
     )
    
    compare_hypervolumes(results_path, seed_folder, selected_seed)


    return (seed_unique_edge_combo_list, seed_unique_circuits_list,
            common_edge_combo_list, common_circuits, 
            edge_combos_unique_to_seed_list, circuits_unique_to_seed_list,
            dose_varied_common_edge_combo_list, dose_varied_common_circuits)

def compare_parteo_fronts(results_path=str, seed_folder=str, 
                          num_obj=int, obj_labels=list
):
    
    pareto_obj_dfs_list = []
    for seed in range(0, 10):
        #import final population circuits
        full_path = results_path + seed_folder + str(seed) + "/"
        with open(full_path + "final_objectives_df.pkl", "rb") as fid:
            pareto_front = pickle.load(fid)
        pareto_obj_dfs_list.append(pareto_front)

    pareto_plot = "pareto_front_set.svg"
    if num_obj == 2:
        plot_pareto_front_set(
            results_path+pareto_plot,
            pareto_obj_dfs_list,
            obj_labels
        )
    elif num_obj == 3:
        plot_pareto_front_set_3D(
            results_path+pareto_plot,
            pareto_obj_dfs_list,
            obj_labels
        )

def compare_hypervolumes(results_path, seed_folder, selected_seed):

    final_hypervolumes = []
    hypervolume_lists = []
    for seed in range(0, 10):
        #import final population circuits
        full_path = results_path + seed_folder + str(seed) + "/"
        with open(full_path + "hypervolumes.pkl", "rb") as fid:
            hypervolumes = pickle.load(fid)
        final_hypervolumes.append(hypervolumes[-1])
        hypervolume_lists.append(hypervolumes)

    n_gens = len(hypervolume_lists[0])
    hypervolumes_plot = "all_hypervolume_progressions.svg"
    plot_hypervolumes_set(results_path + hypervolumes_plot, 
                          n_gens, hypervolume_lists)
    
    y_lower_lim = min(final_hypervolumes) - min(final_hypervolumes)*0.05
    hypervolumes_plot_zoomed = "all_hypervolume_progressions_zoomed.svg"
    plot_hypervolumes_set(results_path + hypervolumes_plot_zoomed, 
                          n_gens, hypervolume_lists, y_lower_lim)
    hypervolumes_vs_combo = "all_hvs_vs_combo_paper.svg"
    plot_hypervolumes_set_vs_combo(results_path+hypervolumes_vs_combo,
                                  n_gens, hypervolume_lists, 45.89854218733082,
                                  selected_seed, y_lower_lim=0)
    
    return final_hypervolumes

def compile_common_results(results_path,
        common_edge_combo_list, common_circuits,
        dose_varied_common_edge_combo_list, dose_varied_common_circuits
):
    
    common_results_dict = {
        "edge_list": common_edge_combo_list + dose_varied_common_edge_combo_list,
        "circuit": common_circuits + dose_varied_common_circuits,
        "dose varied?": ["no"]*len(common_circuits) + ["yes"]*len(dose_varied_common_circuits)
    }

    common_pareto_circuits = pd.DataFrame(data=common_results_dict)
    filename = "common_pareto_circuits.csv"
    common_pareto_circuits.to_csv(results_path+filename)

def compile_seed_comparisons(results_path,
        seed_unique_edge_combo_list, seed_unique_circuits_list,
        edge_combos_unique_to_seed_list, circuits_unique_to_seed_list        
):
    seed_comparisons_dict = {
        "seed": np.arange(10),
        "unique in seed edge lists": seed_unique_edge_combo_list,
        "unique in seed circuits list": seed_unique_circuits_list,
        "number of unique in seed circuits": [len(i) for i in seed_unique_circuits_list],
        "unique to seed edge lists": edge_combos_unique_to_seed_list,
        "unique to seed circuits list": circuits_unique_to_seed_list,
        "number of unique to seed circuits": [len(i) for i in circuits_unique_to_seed_list],
    }

    seed_comparisons = pd.DataFrame(data=seed_comparisons_dict)
    filename = "seed_comparisons.csv"
    seed_comparisons.to_csv(results_path+filename)