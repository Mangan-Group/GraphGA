import numpy as np
import pandas as pd
from typing import Union
from plot_search_results import plot_graph
from define_circuit import Topo
from saving import make_results_analysis_directory
from get_selected_results import (
    import_single_obj_GA_files, 
    import_multi_obj_GA_files, 
    get_single_obj_selected_results,
    get_multi_obj_selected_results,
    get_selected_all_cell_metrics,
    plot_all_cell_objs
)

def run_results_analysis(results_analysis_settings: dict):
    
    folder_path = make_results_analysis_directory(
        results_analysis_settings["repository_path"] + 
        results_analysis_settings["results_path"],
        results_analysis_settings
    )

    if not results_analysis_settings["multi_obj"]:
        (
            settings, 
            unique_obj, 
            unique_circuits
        ) = import_single_obj_GA_files(
            results_analysis_settings["repository_path"],
            results_analysis_settings["results_path"]
        )
        selected_results_df = get_single_obj_selected_results(
            settings, unique_obj, unique_circuits, 
            results_analysis_settings["obj_range"]
        )
    else:
        (
            settings,
            pareto_unique_obj, 
            pareto_unique_circuits
        ) = import_multi_obj_GA_files(
            results_analysis_settings["repository_path"],
            results_analysis_settings["results_path"]
        )
        selected_results_df = get_multi_obj_selected_results(
            settings, pareto_unique_obj,
            pareto_unique_circuits,
            results_analysis_settings["obj_range"]
        )

    all_cell_results_df, all_cell_metrics_df = get_selected_all_cell_metrics(settings, selected_results_df)

    if results_analysis_settings["plot_topologies"]:
        for i, circuit in enumerate(selected_results_df["Topology"].tolist()):
            graph_file_name = ("selected_circuit_" + 
                               results_analysis_settings["selected_results_name"] 
                               + "_" + str(i) + ".svg")
            plot_graph(folder_path + graph_file_name, circuit)

    if results_analysis_settings["plot_all_cell_results"]:
        plot_all_cell_objs(
            folder_path, settings,
            all_cell_results_df
        )

    selected_results_file_name = (
        "selected_results_" + 
        results_analysis_settings["selected_results_name"] +
        ".csv"
    )
    # selected_results_df.to_csv(folder_path + selected_results_file_name)
    all_cell_results_file_name = "all_cell_" + selected_results_file_name
    # all_cell_results_df.to_csv(folder_path + all_cell_results_file_name)
    all_cell_metrics_file_name = "all_cell_metrics_" + results_analysis_settings["selected_results_name"]
    all_cell_metrics_df.to_csv(folder_path + all_cell_metrics_file_name + ".csv")



if __name__ == "__main__":
    repo_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/"

    #amplifier vary dose
    results_path_amp_vary_dose = "Amp_seed_pop_vary_dose/2024-03-08_Amplifier_pop_vary_dose_new_dose_terms_seed_0/"
    amp_vary_dose_results_name = "high_ON_rel"
    amp_vary_dose_obj_range = [63.4]

    #signal conditioner all ZFs
    results_path_sc_all_ZFs = "SC_seed_pop_DsRED_inhibitor/2024-03-06_Signal_Cond_pop_DsRED_inhibitor_ngen80_new_dose_terms_seed_0/"
    sc_all_ZFs_results_name = "high_FI_rel"
    sc_all_ZFs_obj_range = {"FI_rel": [1.0]}

    #pulse all ZFs t_pulse
    results_path_pulse_all_ZFs_t_pulse = "Pulse_seed_pop_DsRED_inhibitor/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ngen80_new_dose_terms_seed_0/"
    pulse_all_ZFs_t_pulse_results_name = "low_t_pulse"
    pulse_all_ZFs_t_pulse_obj_range = {"t_pulse": [0.0, 24.0]}

    #pulse ZF1 and ZF2 t_pulse
    results_path_pulse_ZF1_ZF2_t_pulse = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_new_dose_terms_seed_0/"
    pulse_ZF1_ZF2_t_pulse_results_name = "low_t_pulse"
    pulse_ZF1_ZF2_t_pulse_obj_range = {"t_pulse": [0.0, 24.0]}

    results_paths = [results_path_amp_vary_dose, results_path_sc_all_ZFs, results_path_pulse_all_ZFs_t_pulse, results_path_pulse_ZF1_ZF2_t_pulse]
    results_names = [amp_vary_dose_results_name, sc_all_ZFs_results_name, pulse_all_ZFs_t_pulse_results_name, pulse_ZF1_ZF2_t_pulse_results_name]
    obj_ranges = [amp_vary_dose_obj_range, sc_all_ZFs_obj_range, pulse_all_ZFs_t_pulse_obj_range, pulse_ZF1_ZF2_t_pulse_obj_range]
    multi_objs = [False, True, True, True]

    # results_analysis_settings = {
    #     "repository_path": "/Users/kdreyer/Documents/Github/GraphGA/GA_results/",
    #     "results_path": "",
    #     "selected_results_name": "",
    #     "obj_range": "",
    #     "multi_obj": False,
    #     "plot_topologies": True,
    #     "plot_all_cell_results": True
    # }

    # for i in range(len(results_paths)):
    #     results_analysis_settings["results_path"] = results_paths[i]
    #     results_analysis_settings["selected_results_name"] = results_names[i]
    #     results_analysis_settings["obj_range"] = obj_ranges[i]
    #     results_analysis_settings["multi_obj"] = multi_objs[i]

    #     run_results_analysis(results_analysis_settings)

    results_analysis_settings = {
        "repository_path": "/Users/kdreyer/Documents/Github/GraphGA/GA_results/",
        "results_path": "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_3obj_126h_ZF1_ZF2_new_dose_terms_seed_0/",
        "selected_results_name": "low_t_pulse",
        "obj_range": {"t_pulse": [0.0, 24.0], "prominence_rel": [0.2]},
        "multi_obj": True,
        "plot_topologies": True,
        "plot_all_cell_results": True
    }

    run_results_analysis(results_analysis_settings)
