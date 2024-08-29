import numpy as np
import pandas as pd
import json
import pickle
from multiprocessing import Pool
from typing import Union
from plot_search_results import (
    plot_1D_all_cell_obj,
    plot_2D_all_cell_obj,
    plot_all_cell_time_series
)
from define_circuit import Topo
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator


def import_single_obj_GA_files(repository_path, results_path):

    #import settings
    with open(repository_path + results_path + 
              "settings.json", "rb") as fid:
        settings = json.load(fid)

    #import unique obj
    with open(repository_path + results_path + 
              "unique_objectives.pkl", "rb") as fid:
        unique_obj = pickle.load(fid)
        unique_obj = np.abs(unique_obj)
    #import unique circuits
    with open(repository_path + results_path + 
              "unique_circuits.pkl", "rb") as fid:
        unique_circuits = pickle.load(fid)

    return settings, unique_obj, unique_circuits

def import_multi_obj_GA_files(repository_path, results_path):

    #import settings
    with open(repository_path + results_path + 
              "settings.json", "rb") as fid:
        settings = json.load(fid)

    #import final pareto front obj
    pareto_obj = pd.read_pickle(repository_path + results_path +
                                "final_objectives_df.pkl")
    #get unique pareto front obj
    pareto_unique_obj = pareto_obj.drop_duplicates()
    pareto_unique_obj = pareto_unique_obj.copy()
    #convert obj column values to abs values
    for column in pareto_unique_obj.columns:
        if column != "type":
            pareto_unique_obj[column] = pareto_unique_obj[column].abs()
    #get indices with unique objs
    pareto_unique_obj_idx = pareto_unique_obj.index.tolist()
    # print(pareto_unique_obj_idx)
    #import final population circuits
    with open(repository_path + results_path + 
              "final_population.pkl", "rb") as fid:
        pareto_circuits = pickle.load(fid)
    #get unique pareto front circuits
    pareto_unique_circuits = pareto_circuits[pareto_unique_obj_idx]

    #reset index for future code selecting by index in obj and unique circuits
    pareto_unique_obj.reset_index(inplace=True)

    return settings, pareto_unique_obj, pareto_unique_circuits

def get_single_obj_selected_results(
        settings, unique_obj, unique_circuits, obj_range
):

    # sort descending- need to get idx to sort circuits too
    if len(obj_range) == 1:
        selected_obj_idx = np.where(unique_obj > obj_range[0])[0]
        selected_unique_obj = unique_obj[selected_obj_idx].flatten()
        selected_unique_circuits = unique_circuits[selected_obj_idx].flatten()

    elif len(obj_range) == 2:
        selected_obj_idx = np.where(
            np.logical_and(
                unique_obj > obj_range[0], 
                unique_obj < obj_range[1]))[0]
        selected_unique_obj = unique_obj[selected_obj_idx].flatten()
        selected_unique_circuits = unique_circuits[selected_obj_idx].flatten()

    descending_obj_idx = np.argsort(selected_unique_obj)[::-1]
    sorted_selected_unique_obj = selected_unique_obj[descending_obj_idx]
    sorted_selected_unique_circuits = selected_unique_circuits[descending_obj_idx]

    edge_lists = []
    dose_dict_list = []
    for circuit in sorted_selected_unique_circuits:
        edge_lists.append(circuit.edge_list)
        dose_dict_list.append(circuit.dose)

    selected_results_dict = {
        "Topology": sorted_selected_unique_circuits,
        "Edge list": edge_lists,
        "Doses": dose_dict_list,
        settings["obj_labels"][0]: sorted_selected_unique_obj
    }
    selected_results_df = pd.DataFrame.from_dict(selected_results_dict)

    return selected_results_df

def get_multi_obj_selected_results(
        settings, pareto_unique_obj, 
        pareto_unique_circuits, obj_range:dict
):

    for key, val_list in obj_range.items():
        obj_name = key
        range_list = val_list
        if len(range_list) == 1:
            selected_pareto_obj = pareto_unique_obj[
                pareto_unique_obj[obj_name] > range_list[0]
            ]
        elif len(range_list) == 2:
            selected_pareto_obj = pareto_unique_obj[
                (pareto_unique_obj[obj_name] > range_list[0]) &
                (pareto_unique_obj[obj_name] < range_list[1])
            ]
        pareto_unique_obj = selected_pareto_obj.copy()

    selected_obj_idx = selected_pareto_obj.index.tolist()
    selected_pareto_circuits = pareto_unique_circuits[selected_obj_idx]
    selected_pareto_obj.reset_index(inplace=True)
    sort_by = list(obj_range.keys())[0]
    sorted_selected_pareto_obj = selected_pareto_obj.sort_values(sort_by, ascending=False)
    sorted_obj_idx = sorted_selected_pareto_obj.index.tolist()
    sorted_selected_pareto_circuits = selected_pareto_circuits[sorted_obj_idx].flatten()

    edge_lists = []
    dose_dict_list = []
    for circuit in sorted_selected_pareto_circuits:
        edge_lists.append(circuit.edge_list)
        dose_dict_list.append(circuit.dose)

    selected_results_dict = {
        "Topology": sorted_selected_pareto_circuits,
        "Edge list": edge_lists,
        "Doses": dose_dict_list
    }
    for obj in settings["obj_labels"]:
        selected_results_dict[obj] = sorted_selected_pareto_obj[obj].tolist()
    
    selected_results_df = pd.DataFrame.from_dict(selected_results_dict)

    return selected_results_df


def get_selected_all_cell_metrics(settings, selected_results_df):

    if settings["test_case"] == "Amplifier":
        test_case = Amplifier
    elif settings["test_case"] == "SignalConditioner":
        test_case = SignalConditioner
    elif settings["test_case"] == "PulseGenerator":
        test_case = PulseGenerator
    else:
        raise Exception("Error: test case not defined")   
    
    problem = test_case(
        promo_node=settings["promo_node"],
        dose_specs=settings["dose_specs"],
        max_part=settings["max_part"],
        inhibitor=settings["inhibitor"],
        DsRed_inhibitor=settings["DsRed_inhibitor"],
        num_dict=settings["num_dict"],
        n_gen=settings["n_gen"],
        probability_crossover=settings["probability_crossover"],
        probability_mutation=settings["probability_mutation"],
        mutate_dose=settings["mutate_dose"],
        pop=settings["pop"],
        num_processes=settings["num_processes"],
        obj_labels=settings["obj_labels"],
        max_time=settings["max_time"],
        single_cell_tracking=True
    )

    with Pool(settings["num_processes"]) as pool:
        obj_all_cells_dict_list = pool.imap(
            problem.func, selected_results_df["Topology"].tolist()
        )
        pool.close()
        pool.join()
    obj_all_cells_dict_list = list(obj_all_cells_dict_list)
    obj_all_cells_dict_list = np.asarray(obj_all_cells_dict_list, dtype=object)
    # extract list of all_cells_dfs from output list
    all_cells_dict_list = obj_all_cells_dict_list[:,1].tolist()
    all_cell_results_df = pd.DataFrame.from_dict(all_cells_dict_list)

    for label in settings["obj_labels"]:
        all_cell_results_df[label + "_mean"] = selected_results_df[label].tolist()

    if isinstance(problem, PulseGenerator):
        all_cell_results_df["single_cell_peaks"] = 0
        all_cell_results_df["single_cell_peaks"] = all_cell_results_df["single_cell_peaks"].astype(object)
        all_cell_results_df["single_cell_prominence"] = 0
        all_cell_results_df["single_cell_prominence"] = all_cell_results_df["single_cell_prominence"].astype(object)
        for index, row in all_cell_results_df.iterrows():
            peak_cell_list = []
            prom_cell_list = []
            for i in range(20):
                peak_cell = max(row["Rep_rel time series for each cell"][i])
                peak_cell_list.append(peak_cell)
                prom_cell =  test_case.calc_prominence_rel(row["Rep_rel time series for each cell"][i], peak_cell)
                prom_cell_list.append(prom_cell)
            all_cell_results_df.at[index, "single_cell_peaks"] = peak_cell_list
            all_cell_results_df.at[index, "single_cell_prominence"] = prom_cell_list

        all_cell_metrics_df = all_cell_results_df.copy().drop(['Rep_rel time series for each cell', 'Rep_rel time series mean'], axis=1)

    else:
        all_cell_metrics_df = None
        
    return all_cell_results_df, all_cell_metrics_df

def plot_all_cell_objs(base_path, settings, all_cell_results_df):

    if settings["test_case"] == "Amplifier":
        plot_function = plot_1D_all_cell_obj
    elif settings["test_case"] == "SignalConditioner":
        plot_function = plot_2D_all_cell_obj
    elif settings["test_case"] == "PulseGenerator":
        plot_function = plot_all_cell_time_series

    for index, row in all_cell_results_df.iterrows():
        figure_path = base_path + "all_cells_labeled" + str(index) + ".svg"
        plot_function(figure_path, settings, row)