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
    selected_results_df.to_csv(folder_path + selected_results_file_name)
    all_cell_results_file_name = "all_cell_" + selected_results_file_name
    all_cell_results_df.to_csv(folder_path + all_cell_results_file_name)
    if all_cell_metrics_df:
        all_cell_metrics_file_name = "all_cell_metrics_" + results_analysis_settings["selected_results_name"]
        all_cell_metrics_df.to_csv(folder_path + all_cell_metrics_file_name + ".csv")



if __name__ == "__main__":
    repo_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/"

    #amplifier vary dose
    results_path_amp_vary_dose = "Amp_seed_pop_vary_dose/Original_hyperparams_worked_well/2024-04-23_Amplifier_pop_vary_dose_original_hp_seed_0/"
    amp_vary_dose_results_name = "high_ON_rel"
    amp_vary_dose_obj_range = [63.11786016]



    results_analysis_settings = {
        "repository_path": "/Users/kdreyer/Documents/Github/GraphGA/GA_results/",
        "results_path": results_path_amp_vary_dose,
        "selected_results_name": amp_vary_dose_results_name,
        "obj_range": amp_vary_dose_obj_range,
        "multi_obj": False,
        "plot_topologies": True,
        "plot_all_cell_results": True
    }

    run_results_analysis(results_analysis_settings)
