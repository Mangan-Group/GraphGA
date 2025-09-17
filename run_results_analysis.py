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

def run_results_analysis(
    results_analysis_settings: dict
):
    """Run the results analysis."""
    
    # make the directory to store analysis outputs
    folder_path = make_results_analysis_directory(
        results_analysis_settings["repository_path"] + 
        results_analysis_settings["results_path"],
        results_analysis_settings
    )

    # if not multi-objective, run single objecitive
    # results analysis
    if not results_analysis_settings["multi_obj"]:
        # import files
        (
            settings, 
            unique_obj, 
            unique_circuits
        ) = import_single_obj_GA_files(
            results_analysis_settings["repository_path"],
            results_analysis_settings["results_path"]
        )
        # subset results files for selected results
        selected_results_df = get_single_obj_selected_results(
            settings, unique_obj, unique_circuits, 
            results_analysis_settings["obj_range"]
        )

    # otherwise, run multi-objective results analysis
    else:
        # if specifying objectives other than pareto
        # to analyze, set objs_ as in settings
        if "objs_to_analyze" in results_analysis_settings:
            objs_ = results_analysis_settings["objs_to_analyze"]
        else:
            objs_ = "pareto"
        # import files
        (
            settings,
            pareto_unique_obj, 
            pareto_unique_circuits
        ) = import_multi_obj_GA_files(
            results_analysis_settings["repository_path"],
            results_analysis_settings["results_path"],
            objs=objs_
        )
        # subset results files for selected results
        selected_results_df = get_multi_obj_selected_results(
            settings, pareto_unique_obj,
            pareto_unique_circuits,
            results_analysis_settings["obj_range"]
        )

    # if the population model was used, get metrics
    # for indivudual cells
    if settings["pop"]:
        all_cell_results_df, all_cell_metrics_df = (
            get_selected_all_cell_metrics(
                settings, 
                selected_results_df,
                folder_path
            )
        )
    else:
        all_cell_results_df = None
        all_cell_metrics_df = None

    # if plot_topologies is True, plot the topology 
    # network graphs
    if results_analysis_settings["plot_topologies"]:
        for i, circuit in enumerate(selected_results_df["Topology"].tolist()):
            graph_file_name = ("selected_circuit_" + 
                               results_analysis_settings["selected_results_name"] 
                               + "_" + str(i) + ".svg")
            plot_graph(folder_path + graph_file_name, circuit)

    # if plot_all_cell_results is True, plot the objectives
    # for individual cells (only relevant for population 
    # model)
    if results_analysis_settings["plot_all_cell_results"]:
        plot_all_cell_objs(
            folder_path, settings,
            all_cell_results_df
        )

    # save results files
    selected_results_file_name = (
        "selected_results_" + 
        results_analysis_settings["selected_results_name"]
    )
    selected_results_df.to_csv(
        folder_path + selected_results_file_name + ".csv"
    )
    selected_results_df.to_pickle(
        folder_path + selected_results_file_name + ".pkl"
    )
    if all_cell_results_df is not None:
        all_cell_results_file_name = "all_cell_" + selected_results_file_name
        all_cell_results_df.to_csv(
            folder_path + all_cell_results_file_name + ".csv"
        )
        all_cell_results_df.to_pickle(
            folder_path + all_cell_results_file_name + ".pkl"
        )
    if all_cell_metrics_df is not None:
        all_cell_metrics_file_name = (
            "all_cell_metrics_" + results_analysis_settings["selected_results_name"]
        )
        all_cell_metrics_df.to_csv(
            folder_path + all_cell_metrics_file_name + ".csv"
        )
        all_cell_metrics_df.to_pickle(
            folder_path + all_cell_metrics_file_name + ".pkl"
        )




if __name__ == "__main__":
    # define run settings and run analysis

    results_path = ("Pulse_single_cell/Optimized_hyperparams/"
    "t_pulse_fixed_pop_max_hv/run1_ngen80/"
    "2024-10-09_Pulse_single_DsRED_t_pulse_opt_hps_ngen80_seed_2/")
    results_name = "full_pareto"
    results_analysis_settings = {
        "repository_path": ("/Users/kdreyer/Library/CloudStorage/"
                "OneDrive-NorthwesternUniversity/KatieD_LL/GCAD_Collab/"
                "Selected_GA_results_paper/"),
        "results_path": results_path,
        "selected_results_name": results_name,
        "obj_range": {"t_pulse": [0, 70.0]},
        "multi_obj": True,
        "plot_topologies": True,
        "plot_all_cell_results": False
    }

    run_results_analysis(results_analysis_settings)


