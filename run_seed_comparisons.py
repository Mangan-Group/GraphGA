import numpy as np
import pickle
from get_comparisons_across_seeds import (
    compare_pareto_circuits_across_seeds, 
    compare_parteo_fronts,
    compile_common_results,
    compile_seed_comparisons
)
from plot_search_results import plot_obj_progression_set

def run_seed_results_comparisons(results_path, seed_folder, num_obj, obj_labels=list, ngens=None):

    if num_obj == 1:
        objs_list = []
        for seed in range(0, 10):
            full_path = results_path+seed_folder+str(seed)+"/minimum_obj_all_gens.pkl"
            with open(full_path, "rb") as fid:
                objs = np.abs(pickle.load(fid))
            objs_list.append(objs)

        y_lower_lim = np.max(objs_list) - 0.1*np.max(objs_list)
        y_ticks = [0, 20, 40, 60]
        figure_path = results_path+obj_labels[0]+"_progression.svg"
        plot_obj_progression_set(figure_path, ngens, 
                                objs_list, "Maximum "+obj_labels[0],
                                y_ticks = y_ticks)
        
        figure_path_zoomed = results_path+obj_labels[0]+"_progression_zoomed.svg"
        plot_obj_progression_set(figure_path_zoomed, ngens, 
                                objs_list, "Maximum "+obj_labels[0], 
                                y_lower_lim=y_lower_lim)

    else:
        (seed_unique_edge_combo_list, 
        seed_unique_circuits_list,
        common_edge_combo_list, common_circuits, 
        edge_combos_unique_to_seed_list, 
        circuits_unique_to_seed_list,
        dose_varied_common_edge_combo_list, 
        dose_varied_common_circuits) = compare_pareto_circuits_across_seeds(
            results_path, seed_folder
        )

        compare_parteo_fronts(results_path, seed_folder,
                              num_obj, obj_labels)

        compile_common_results(results_path, common_edge_combo_list, 
                            common_circuits,
                            dose_varied_common_edge_combo_list,
                            dose_varied_common_circuits
        )

        compile_seed_comparisons(results_path,
                                seed_unique_edge_combo_list,
                                seed_unique_circuits_list,
                                edge_combos_unique_to_seed_list, 
                                circuits_unique_to_seed_list
        )

path_results = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
path_pulse = "SC_seed_single_cell_DsRED_inhibitor/Optimized_hyperparams_fixed_pop_ngen50_max_hv/run2_ngen80/"
results_runs = "2024-10-03_SignalCond_single_cell_DsRED_opt_hp_fixedpop_ngen50_seed_"

run_seed_results_comparisons(path_results+path_pulse, results_runs, 2, ["ON_rel", "FI_rel"])
