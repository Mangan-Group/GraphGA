import numpy as np
import pickle
from get_comparisons_across_seeds import (
    compare_pareto_circuits_across_seeds, 
    compare_parteo_fronts,
    compile_common_results,
    compile_seed_comparisons
)
from plot_search_results import plot_obj_progression_set



def run_seed_results_comparisons(
        results_path,
        seed_folder,
        num_obj,
        obj_labels,
        ngens,
        selected_seed
):
    """Runs the seed comparison functions in 
    get_comparisons_across_seeds.py."""

    # for single objective optimization
    if num_obj == 1:
        objs_list = []
        # loop through seeds
        for seed in range(0, 10):
            full_path = results_path+seed_folder+str(seed)+"/minimum_obj_all_gens.pkl"
            # load the list of minimum objectives for each generation
            with open(full_path, "rb") as fid:
                objs = np.abs(pickle.load(fid))
            # append objective list to list for all seeds
            objs_list.append(objs)

        figure_path_zoomed = results_path+obj_labels[0]+"_progression_zoomed_paper.svg"
        # plot the maximum objective progression across
        # generations for each seed for comparison (multiplies
        # the negative of the objective value by -1 before
        # plotting)
        plot_obj_progression_set(figure_path_zoomed, ngens, 
                                objs_list, obj_labels[0],
                                selected_seed, 
                                opt_obj=max(objs_list[0])
        )

    # for multi-objective optimization
    else:
        # get common and unique pareto front
        # circuits across seeds
        (seed_unique_edge_combo_list, 
        seed_unique_circuits_list,
        common_edge_combo_list, common_circuits, 
        edge_combos_unique_to_seed_list, 
        circuits_unique_to_seed_list,
        dose_varied_common_edge_combo_list, 
        dose_varied_common_circuits) = compare_pareto_circuits_across_seeds(
            results_path, seed_folder, selected_seed
        )

        # plot pareto front for each seed for comparison
        compare_parteo_fronts(results_path, seed_folder,
                              num_obj, obj_labels)

        # compile and save common pareto front circuits
        # across seeds
        compile_common_results(results_path, common_edge_combo_list, 
                            common_circuits,
                            dose_varied_common_edge_combo_list,
                            dose_varied_common_circuits
        )

        # compile and save unique pareto front circuits
        # across seeds
        compile_seed_comparisons(results_path,
                                seed_unique_edge_combo_list,
                                seed_unique_circuits_list,
                                edge_combos_unique_to_seed_list, 
                                circuits_unique_to_seed_list
        )


### set up to run the function ###
path_test_case = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/GA_results_quest/run2_ngen10000/"
results_runs = "Pulse_single_seed_"
run_seed_results_comparisons(path_test_case, results_runs, 2, ["t_pulse", "prominence_rel"], 10000, 6)

