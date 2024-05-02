from get_comparisons_across_seeds import (
    compare_pareto_circuits_across_seeds, 
    compile_common_results,
    compile_seed_comparisons
)

def run_seed_results_comparisons(results_path, seed_folder):

    (seed_unique_edge_combo_list, 
     seed_unique_circuits_list,
     common_edge_combo_list, common_circuits, 
     edge_combos_unique_to_seed_list, 
     circuits_unique_to_seed_list,
     dose_varied_common_edge_combo_list, 
     dose_varied_common_circuits) = compare_pareto_circuits_across_seeds(
         results_path, seed_folder
     )

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
path_sc = "Pulse_seed_single_cell_DsRED_inhibitor/Optimized_hyperparams_vary_pop_t_pulse_opt_stdev_ngen50_nseed4/run1_ngen80/"
results_runs = "2024-04-29_Pulse_single_cell_DsRED_inhibitor_t_pulse_opt_hp_stdev_ngen80_nseed4_seed_"

run_seed_results_comparisons(path_results+path_sc, results_runs)
