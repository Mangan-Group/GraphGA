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
        if "objs_to_analyze" in results_analysis_settings:
            objs_ = results_analysis_settings["objs_to_analyze"]
        else:
            objs_ = "pareto"
        (
            settings,
            pareto_unique_obj, 
            pareto_unique_circuits
        ) = import_multi_obj_GA_files(
            results_analysis_settings["repository_path"],
            results_analysis_settings["results_path"],
            objs=objs_
        )
        selected_results_df = get_multi_obj_selected_results(
            settings, pareto_unique_obj,
            pareto_unique_circuits,
            results_analysis_settings["obj_range"]
        )

    
    if settings["pop"]:
        all_cell_results_df, all_cell_metrics_df = get_selected_all_cell_metrics(settings, selected_results_df)
    else:
        all_cell_results_df = None
        all_cell_metrics_df = None

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
        results_analysis_settings["selected_results_name"]
    )
    selected_results_df.to_csv(folder_path + selected_results_file_name + ".csv")
    selected_results_df.to_pickle(folder_path + selected_results_file_name + ".pkl")
    if all_cell_results_df is not None:
        all_cell_results_file_name = "all_cell_" + selected_results_file_name
        all_cell_results_df.to_csv(folder_path + all_cell_results_file_name + ".csv")
        all_cell_results_df.to_pickle(folder_path + all_cell_results_file_name + ".pkl")
    if all_cell_metrics_df is not None:
        all_cell_metrics_file_name = "all_cell_metrics_" + results_analysis_settings["selected_results_name"]
        all_cell_metrics_df.to_csv(folder_path + all_cell_metrics_file_name + ".csv")
        all_cell_metrics_df.to_pickle(folder_path + all_cell_metrics_file_name + ".pkl")




if __name__ == "__main__":

    results_path = "2024-10-10_Pulse_pop_DsRED_inhibitor_t_frac_pulse_pop200/"
    results_name = "sub_opt"
    results_analysis_settings = {
        # "repository_path": "/Users/kdreyer/Library/CloudStorage/OneDrive-NorthwesternUniversity/KatieD_LL/GCAD_Collab/Selected_GA_results_paper/",
        "repository_path": "/Users/kdreyer/Documents/Github/GraphGA/GA_results/",
        "results_path": results_path,
        "selected_results_name": results_name,
        "obj_range": {"frac_pulse": [0.7, 1.0]},
        "multi_obj": True,
        "objs_to_analyze": "sub_opt",
        "plot_topologies": True,
        "plot_all_cell_results": True

    }

    run_results_analysis(results_analysis_settings)
