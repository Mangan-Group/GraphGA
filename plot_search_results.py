import numpy as np
import networkx as nx
import pickle
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('/Users/kdreyer/Documents/Github/GraphGA/paper.mplstyle.py')

def plot_graph(
        figure_path: str,
        topology: object
):
    plt.figure()
    plt.tight_layout()
    nx.draw_networkx(topology.graph,
        arrows=True,
        arrowsize=15, 
        node_size=600,
        node_shape='s'
    )
    plt.savefig(figure_path, )

def plot_obj_distribution(
        figure_path: str,
        obj_list: list, 
        x_label: str,
        text: str, 
        n_bins: int =30,
):
    plt.figure(figsize=(3,3))
    plt.hist(obj_list, bins=n_bins)
    plt.xlabel(x_label)
    plt.ylabel("count")
    plt.text(0.01, 100, text)
    
    # plt.show()
    plt.savefig(figure_path)
    

def plot_metric(
        figure_path: str,
        metric_list: np.ndarray,
        label: str
):
    generations = np.arange(len(metric_list))

    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    ax.plot(generations, metric_list)
    ax.set_xlabel("generation")
    ax.set_ylabel(label)
    plt.savefig(figure_path)

def plot_hypervolume(
        figure_path: str,
        hypervolumes_list: np.ndarray,
        gens_to_plot: int,
        y_lower_lim: float=None
):
    generations = np.arange(50)[-gens_to_plot:]
    hvs_to_plot = hypervolumes_list[-gens_to_plot:]

    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    for hvs in hypervolumes_list:
        ax.plot(generations, hvs)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Hypervolume")
    if y_lower_lim:
        ax.set_ylim(bottom=y_lower_lim)
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_pareto_front(
        figure_path: str, 
        obj_df: pd.DataFrame,
        obj_labels: list
):
        if obj_df[obj_labels[0]].to_list()[0] < 0:
            obj_df[obj_labels[0]] = obj_df[
                obj_labels[0]]*-1
            obj_df[obj_labels[1]] = obj_df[
                obj_labels[1]]*-1
        
        palette = ["gray", "black"]
        fig, ax = plt.subplots(1, 1, figsize= (4, 4))
        sns.scatterplot(data=obj_df, x= obj_df[obj_labels[0]],
                        y= obj_df[obj_labels[1]], hue='type', 
                        palette=palette, ax=ax)
        plt.xlabel(obj_labels[0])
        plt.ylabel(obj_labels[1])
        # plt.show()
        plt.savefig(figure_path, bbox_inches="tight")

def plot_1D_obj_scatter(
        figure_path: str,
        obj_vals: np.ndarray,
        obj_label: str,
        y_lower_lim: int=None
):
    if obj_vals.flatten()[0] < 0:
        obj_vals = obj_vals*-1
    x_vals = [1]*len(obj_vals)
    jittered_x = x_vals + 0.1*np.random.rand(
        len(x_vals))
    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    ax.plot(jittered_x, obj_vals, linestyle="None",
             marker="o", color="gray")
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_ylabel(obj_label)
    if y_lower_lim:
        ax.set_ylim(lower = y_lower_lim)
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

# def plot_2D_obj_scatter(
#         figure_path: str,
#         obj_vals: np.ndarray,
#         obj_labels: list,
#         y_lower_lim: int=None
# ):
#     if obj_df[obj_labels[0]].to_list()[0] < 0:
#         obj_df[obj_labels[0]] = obj_df[
#         obj_labels[0]]*-1
#     obj_df[obj_labels[1]] = obj_df[
#         obj_labels[1]]*-1
#     # if obj_vals.flatten()[0] < 0:
#     #     obj_vals = obj_vals*-1
#     # x_vals = [1]*len(obj_vals)
#     # jittered_x = x_vals + 0.1*np.random.rand(
#     #     len(x_vals))
#     # fig, ax = plt.subplots(1, 1, figsize= (4, 4))
#     # ax.plot(jittered_x, obj_vals, linestyle="None",
#     #          marker="o", color="gray")
#     # ax.set_xticklabels([])
#     # ax.set_xticks([])
#     # ax.set_ylabel(obj_labels)
#     if y_lower_lim:
#         ax.set_ylim(lower = y_lower_lim)
#     # plt.show()
#     plt.savefig(figure_path, bbox_inches="tight")

# Amplifier single cell GA with constant dose, Z_mat sampling
# with population model
# path = ("/Users/kdreyer/Documents/Github/GraphGA/GA_results/" + 
#         "2023-10-20_Amplifier_Z_matrix_samples_obj/" + 
#         "Z_mat_sampling_all_GA_circuits.pkl")
# amp_const_dose = pd.read_pickle(path)

# ON_rel_range = amp_const_dose["objectives_range"].tolist()
# ON_rel_range_mean = round(np.mean(ON_rel_range), 4)
# fig_text = "mean = " + str(ON_rel_range_mean)
# fig_path = ("/Users/kdreyer/Documents/Github/GraphGA/GA_results/" + 
#         "2023-10-20_Amplifier_Z_matrix_samples_obj/" + 
#         "ON_rel_range_distribution.svg")
# plot_obj_distribution(fig_path, ON_rel_range, 
#                       "ON_rel range", fig_text)

# Amplifier single cell GA with varied dose, Z_mat sampling
# with population model
# path = ("/Users/kdreyer/Documents/Github/GraphGA/GA_results/" + 
#         "2023-10-20_Amplifier_Z_matrix_samples_obj_vary_dose/" + 
#         "Z_mat_sampling_all_GA_circuits.pkl")
# amp_varied_dose = pd.read_pickle(path)

# ON_rel_range_varied_dose = amp_varied_dose["objectives_range"].tolist()
# ON_rel_range_vd_mean = round(np.mean(ON_rel_range_varied_dose), 4)
# fig_text_vd = "mean = " + str(ON_rel_range_vd_mean)
# fig_path = ("/Users/kdreyer/Documents/Github/GraphGA/GA_results/" + 
#         "2023-10-20_Amplifier_Z_matrix_samples_obj_vary_dose/" + 
#         "ON_rel_range_distribution.svg")
# plot_obj_distribution(fig_path, ON_rel_range_varied_dose,
#                       "ON_rel range", fig_text_vd)

# Signal Conditioner single cell GA with DsRED inhibitor, Z_mat sampling
# with population model
# path = ("/Users/kdreyer/Documents/Github/GraphGA/GA_results/" +
#         "2023-10-24_SignalConditioner_DsRED_Z_matrix_samples_obj/" +
#         "Z_mat_sampling_all_GA_circuits.pkl")
# SC_varied_dose = pd.read_pickle(path)
# FI_rel_range = SC_varied_dose["objectives[1]_range"].tolist()
# FI_rel_range_mean = round(np.mean(FI_rel_range), 4)
# fig_text = "mean = " + str(FI_rel_range_mean)
# fig_path = ("/Users/kdreyer/Documents/Github/GraphGA/GA_results/" +
#             "2023-10-24_SignalConditioner_DsRED_Z_matrix_samples_obj/" +
#             "FI_rel_range_distribution.svg")
# plot_obj_distribution(fig_path, FI_rel_range, 
#                       "FI_rel range", fig_text)

# Signal Conditioner single cell GA with inhibitor, Z_mat sampling
# with population model
# path = ("/Users/kdreyer/Documents/Github/GraphGA/GA_results/" +
#         "2023-10-24_SignalConditioner_Z_matrix_samples_obj/" +
#         "Z_mat_sampling_all_GA_circuits.pkl")
# SC_varied_dose = pd.read_pickle(path)
# FI_rel_range = SC_varied_dose["objectives[1]_range"].tolist()
# FI_rel_range_mean = round(np.mean(FI_rel_range), 4)
# fig_text = "mean = " + str(FI_rel_range_mean)
# fig_path = ("/Users/kdreyer/Documents/Github/GraphGA/GA_results/" +
#             "2023-10-24_SignalConditioner_Z_matrix_samples_obj/" +
#             "FI_rel_range_distribution.svg")
# plot_obj_distribution(fig_path, FI_rel_range, 
#                       "FI_rel range", fig_text) 

# Signal Conditioner single cell GA hypervolumes for all seeds
# hv_list = []
# for seed in range(0, 20):
#     path = ("/Users/kdreyer/Documents/Github/GraphGA/GA_results/" +
#             "SC_seed_single_cell_DsRED_inhibitor/2023-10-24_Signal_cond_single_cell_DsRED_inhibitor_seed_" +
#             str(seed) + "/hypervolumes.pkl")
#     with open(path, "rb") as fid:
#         hv = pickle.load(fid)
#     hv_list.append(hv[-40:])

# print(hv_list)

# plot_hypervolumes(hv_list, "signal_conditioner", 40, 47.3, "/Users/kdreyer/Documents/Github/GraphGA/GA_results SC_seed_single_cell_DsRED_inhibitor/")
    
