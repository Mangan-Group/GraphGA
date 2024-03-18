import numpy as np
import pickle
import pandas as pd
from saving import make_main_directory
from load_Z_mat_samples import Z_mat_list
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator
from solve_objectives_for_Z_matrix_samples import (
    select_Z_matrix_sampling_topologies,
    solve_all_topology_objectives,
    get_objective_errors
)
from plot_search_results import (
    plot_1D_obj_confidence_interval, 
    plot_2D_obj_confidence_interval, 
    plot_3D_obj_confidence_interval
)

def run_Z_matrix_sampling(
        testcase: object,
        settings: dict,
        Z_mat_list: list,
):
    folder_path = make_main_directory(settings)

    topologies, objectives = select_Z_matrix_sampling_topologies(
        settings["repository_path"]+settings["results_path"],
        obj_threshold=settings["objective_threshold"]
    )
    selected_topologies_file_name = "Z_matrix_sampling_topologies.pkl"
    with open(
        folder_path + "/" + selected_topologies_file_name, "wb"
    ) as fid:
        pickle.dump(topologies, fid)
    
    problem = testcase(
        promo_node=settings["promo_node"],
        dose_specs=settings["dose_specs"],
        max_part=settings["max_part"],
        inhibitor=settings["inhibitor"],
        DsRed_inhibitor=settings["DsRed_inhibitor"],
        num_dict=settings["num_dict"],
        n_gen=settings["n_gen"],
        probability_crossover=settings["probability_crossover"],
        probability_mutation=settings["probability_mutation"],
        pop=True,
        num_processes=settings["num_processes"],
        obj_labels=settings["obj_labels"],
        max_time=settings["max_time"],
        single_cell_tracking=False
    )

    Z_matrix_sampling = solve_all_topology_objectives(
        problem, topologies, Z_mat_list
    )

    Z_matrix_sampling_file_name = "Z_matrix_sampling_for_CI.pkl"
    Z_matrix_sampling.to_pickle(
        folder_path + "/" + Z_matrix_sampling_file_name
    )

    for i in range(len(settings["CI_metrics"])):
        CI_metric_maxes = get_objective_errors(Z_matrix_sampling, settings["CI_metrics"][i])
        if len(CI_metric_maxes) == 1:
            figure_path = folder_path + "/" + settings["CI_metrics"][i][0]+"_CI.svg"
            plot_1D_obj_confidence_interval(
                settings["repository_path"]+settings["results_path"],
                figure_path, CI_metric_maxes[0], settings["obj_labels"], settings["CI_ylim"]
            )
        elif len(CI_metric_maxes) == 2:
            figure_path = folder_path + "/" + settings["CI_metrics"][i][0]+"_"+settings["CI_metrics"][i][1]+"_CI.svg"
            plot_2D_obj_confidence_interval(
                objectives, settings["repository_path"]+settings["results_path"],
                figure_path, CI_metric_maxes, settings["obj_labels"]
            )

        elif len(CI_metric_maxes) == 3:
            figure_path = folder_path + "/" + settings["CI_metrics"][i][0]+"_"+settings["CI_metrics"][i][1]+"_"+settings["CI_metrics"][i][2]+"_CI.svg"
            plot_3D_obj_confidence_interval(
                objectives, settings["repository_path"]+settings["results_path"],
                figure_path, CI_metric_maxes, settings["obj_labels"]
            )

    return Z_matrix_sampling

settings = {
    "test_case": "PulseGenerator",
    "promo_node": "P1",
    "dose_specs": [5, 75, 5],
    "max_part": None,
    "inhibitor": True,
    "DsRed_inhibitor": True,
    "num_dict": None,
    "n_gen": None,
    "probability_crossover": None,
    "probability_mutation": None,
    "mutate_dose": True,
    "pop": True,
    "num_processes": 8,
    "obj_labels": ["t_pulse", "prominence_rel"],
    "objective_threshold": None,
    "max_time": 126,
    "CI_metrics": [["t_pulse_range", "prominence_rel_range"], ["t_pulse_std_error", "prominence_rel_std_error"]],
    "CI_ylim": False,
    "repository_path": "/Users/kdreyer/Documents/Github/GraphGA/",
    "results_path": "GA_results/Pulse_seed_pop_DsRED_inhibitor/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ngen80_new_dose_terms_seed_0/",
    "folder_name": "Pulse_pop_DsRED_inhibitor_Z_matrix_sampling"
}

if __name__ == "__main__":
    if settings["test_case"] == "Amplifier":
            test_case = Amplifier
    elif settings["test_case"] == "SignalConditioner":
        test_case = SignalConditioner
    elif settings["test_case"] == "PulseGenerator":
        test_case = PulseGenerator
    else:
        raise Exception("Error: test case not defined")
    

    Z_matrix_sampling = run_Z_matrix_sampling(test_case, settings, Z_mat_list)