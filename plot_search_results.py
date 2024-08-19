import numpy as np
import networkx as nx
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import seaborn as sns
from scipy.interpolate import interp2d
from math import ceil, floor

plt.style.use('/Users/kdreyer/Documents/Github/GraphGA/paper.mplstyle.py')
sky_blue = [i/255 for i in [86, 180, 233]]

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
    plt.savefig(figure_path)

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

    fig, ax = plt.subplots(1, 1, figsize= (2.5, 2.5))
    ax.plot(generations, metric_list)
    ax.set_xlabel("generation")
    ax.set_ylabel(label)
    plt.savefig(figure_path)

def plot_hypervolume(
        figure_path: str,
        n_gens: int,
        hypervolumes_list: np.ndarray,
        y_lower_lim: float=None
):
    generations = np.arange(n_gens)

    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    ax.plot(generations, hypervolumes_list)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Hypervolume")
    if y_lower_lim:
        ax.set_ylim(bottom=y_lower_lim)
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_hypervolumes_set(
        figure_path: str,
        n_gens: int,
        hypervolumes_lists: np.ndarray,
        y_lower_lim: float=None
):
    generations = np.arange(n_gens)

    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    for hv_list in hypervolumes_lists:
        ax.plot(generations, hv_list)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Hypervolume")
    if y_lower_lim:
        ax.set_ylim(bottom=y_lower_lim)
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_pareto_front(
        figure_path: str, 
        obj_df: pd.DataFrame,
        obj_labels: list,
        types: bool
):
        if np.any(np.array(obj_df[obj_labels[0]].to_list()) < 0):
            obj_df[obj_labels[0]] = obj_df[
                obj_labels[0]]*-1
        if np.any(np.array(obj_df[obj_labels[1]].to_list()) < 0):
            obj_df[obj_labels[1]] = obj_df[
                obj_labels[1]]*-1
            
        if types:
            palette = ["gray", "black"]
            fig, ax = plt.subplots(1, 1, figsize= (2.25, 2.25))
            sns.scatterplot(data=obj_df, x= obj_df[obj_labels[0]],
                            y= obj_df[obj_labels[1]], hue='type', 
                            palette=palette, ax=ax, s=8)
            plt.legend(fontsize="8")
            
        else:
            fig, ax = plt.subplots(1, 1, figsize= (2.25, 2.25))
            sns.scatterplot(data=obj_df, x= obj_df[obj_labels[0]],
                            y= obj_df[obj_labels[1]], 
                            color="black", ax=ax, s=8)

        plt.xlabel(obj_labels[0])
        plt.ylabel(obj_labels[1])
        # plt.xticks([0, 20, 40, 60])
        # plt.yticks([0, 1, 2, 3])
        # plt.show()
        plt.savefig(figure_path, bbox_inches="tight")

def plot_pareto_front3D(
        figure_path: str, 
        obj_df: pd.DataFrame,
        obj_labels: list,
        types: bool
):
    if np.any(np.array(obj_df[obj_labels[0]].to_list()) < 0):
        obj_df[obj_labels[0]] = obj_df[
            obj_labels[0]]*-1
    if np.any(np.array(obj_df[obj_labels[1]].to_list()) < 0):
        obj_df[obj_labels[1]] = obj_df[
            obj_labels[1]]*-1
    if np.any(np.array(obj_df[obj_labels[2]].to_list()) < 0):
        obj_df[obj_labels[2]] = obj_df[
            obj_labels[2]]*-1
            
    # print(obj_df.tail(n=50))
    fig = plt.figure(figsize= (4, 4))
    ax = fig.add_subplot(projection='3d')

    ax.scatter(
        xs=obj_df[obj_labels[0]], ys=obj_df[obj_labels[1]],
        zs=obj_df[obj_labels[2]], color="k", 
    )
    ax.view_init(elev=10, azim=-115)
    ax.set_xlabel(obj_labels[0])
    ax.set_ylabel(obj_labels[1])
    ax.set_zlabel(obj_labels[2])
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_1D_obj_scatter(
        figure_path: str,
        obj_vals: np.ndarray,
        obj_labels: list,
        y_lower_lim: int=None
):
    if obj_vals.flatten()[0] < 0:
        obj_vals = obj_vals*-1
    x_vals = [1]*len(obj_vals)
    jittered_x = x_vals + 0.1*np.random.rand(
        len(x_vals))
    fig, ax = plt.subplots(1, 1, figsize= (2.25, 2))
    ax.plot(jittered_x, obj_vals, linestyle="None",
             marker="o", markersize=1, color="gray")
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_ylabel(obj_labels[0])
    if y_lower_lim:
        ax.set_ylim(lower = y_lower_lim)
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_1D_obj_confidence_interval(
        results_path: str,
        figure_path: str,
        CI_metric_max: float,
        obj_labels: list,
        y_lim: int=False
):
    unique_objectives = pd.read_pickle(results_path+"unique_objectives.pkl")
    unique_objectives = unique_objectives.flatten()*-1
    max_objective = max(unique_objectives)

    x_vals = [1]*len(unique_objectives)
    jittered_x = x_vals + 0.1*np.random.rand(
        len(x_vals))
    lower_bound = [max_objective-CI_metric_max]*len(unique_objectives)
    upper_bound = [max_objective]*len(unique_objectives)
    fig, ax = plt.subplots(1, 1, figsize= (2.25, 2))
    ax.plot(jittered_x, unique_objectives, linestyle="None",
             marker="o", markersize=1, color="black", zorder=1) #markersize=1
    jittered_x.sort()
    ax.fill_between(jittered_x, lower_bound, upper_bound, alpha=0.4, color=sky_blue, zorder=2)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_ylabel(obj_labels[0])
    if y_lim:
        # if CI_metric_max <= 0.5:
        #     y_lim = [floor(max_objective)-0.5, ceil(max_objective)]
        # else:
        y_lim = [floor(max_objective)-CI_metric_max, ceil(max_objective)]
        ax.set_ylim(y_lim)
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_2D_obj_confidence_interval(
        objectives: pd.DataFrame,
        results_path: str,
        figure_path: str,
        CI_metric_maxes: list,
        obj_labels: list,
):
    all_objectives = np.abs(pd.read_pickle(results_path+"all_objectives.pkl"))
    
    obj1_CI = CI_metric_maxes[0]
    obj2_CI = CI_metric_maxes[1]

    # upper_obj1 = np.array([i+obj1_CI_list[0] for i in (objectives[obj_labels[0]]*-1)])
    if "t_pulse" in '\t'.join(obj_labels):
        obj1 = np.abs(np.array(objectives[obj_labels[0]]))
        upper_obj1 = np.array([i+obj1_CI for i in obj1])
        sorted_upper_idx = np.argsort(upper_obj1)
        sorted_upper_obj1 = upper_obj1[sorted_upper_idx]
        lower_obj1 = obj1
        sorted_lower_idx = np.argsort(lower_obj1)
        sorted_lower_obj1 = lower_obj1[sorted_lower_idx]
    else:
        upper_obj1 = np.abs(np.array(objectives[obj_labels[0]]))
        sorted_upper_idx = np.argsort(upper_obj1)
        sorted_upper_obj1 = upper_obj1[sorted_upper_idx]
        lower_obj1 = np.array([i-obj1_CI for i in upper_obj1])
        sorted_lower_idx = np.argsort(lower_obj1)
        sorted_lower_obj1 = lower_obj1[sorted_lower_idx]

    # upper_obj2 = np.array([i+obj2_CI_list[0] for i in (objectives[obj_labels[1]]*-1)])
    upper_obj2 = np.abs(np.array(objectives[obj_labels[1]]))
    sorted_upper_obj2 = upper_obj2[sorted_upper_idx]
    lower_obj2 = np.array([i-obj2_CI for i in upper_obj2])
    sorted_lower_obj2 = lower_obj2[sorted_lower_idx]

    xfill = np.sort(np.concatenate([upper_obj1, lower_obj1]))
    if "t_pulse" in '\t'.join(obj_labels):
        y1fill = np.interp(xfill, sorted_upper_obj1, sorted_lower_obj2)
        y2fill = np.interp(xfill, sorted_lower_obj1, sorted_upper_obj2)
    else:
        y1fill = np.interp(xfill, sorted_upper_obj1, sorted_upper_obj2)
        y2fill = np.interp(xfill, sorted_lower_obj1, sorted_lower_obj2)
    fig, ax = plt.subplots(1, 1, figsize= (2.25, 2))
    ax.fill_between(xfill, y1fill, y2fill, alpha=0.4, color=sky_blue, zorder=2)
    ax.plot(all_objectives[:, 0], all_objectives[:, 1], linestyle="none", 
            marker="o", markersize=1, color="black", zorder=1)

    ax.set_label(obj_labels[0])
    ax.set_ylabel(obj_labels[1])
    if np.any(y1fill < 0) or np.any(y2fill < 0):
        ax.set_ylim(bottom=0)
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_3D_obj_confidence_interval(
        objectives: pd.DataFrame,
        results_path: str,
        figure_path: str,
        CI_metric_maxes: list,
        obj_labels: list,
):
    all_objectives = pd.read_pickle(results_path+"all_objectives.pkl")
    
    obj1_CI = CI_metric_maxes[0]
    obj2_CI = CI_metric_maxes[1]
    obj3_CI = CI_metric_maxes[2]

    upper_obj1 = np.array(objectives[obj_labels[0]])
    sorted_upper_idx = np.argsort(upper_obj1)
    sorted_upper_obj1 = upper_obj1[sorted_upper_idx]
    lower_obj1 = np.array([i-obj1_CI for i in (objectives[obj_labels[0]])])
    sorted_lower_idx = np.argsort(lower_obj1)
    sorted_lower_obj1 = lower_obj1[sorted_lower_idx]

    upper_obj2 = np.array(objectives[obj_labels[1]]*-1)
    sorted_upper_obj2 = upper_obj2[sorted_upper_idx]
    lower_obj2 = np.array([i-obj2_CI for i in (objectives[obj_labels[1]]*-1)])
    sorted_lower_obj2 = lower_obj2[sorted_lower_idx]

    upper_obj3 = np.array(objectives[obj_labels[2]]*-1)
    sorted_upper_obj3 = upper_obj3[sorted_upper_idx]
    lower_obj3 = np.array([i-obj3_CI for i in (objectives[obj_labels[2]]*-1)])
    sorted_lower_obj3 = lower_obj3[sorted_lower_idx]

    # obj1_obj1_upper, obj2_obj2_upper = np.meshgrid(sorted_upper_obj1, sorted_upper_obj2)
    # obj3_obj3_upper = np.array(sorted_upper_obj3*len(sorted_upper_obj3)).reshape(len(sorted_upper_obj3), len(sorted_upper_obj3))
    # obj_3_upper_interpolation = interp2d(obj1_obj1_upper, obj2_obj2_upper, obj3_obj3_upper)
    all_obj1_vals = np.sort(np.concatenate([upper_obj1, lower_obj1]))
    all_obj2_vals = np.sort(np.concatenate([upper_obj2, lower_obj2]))

    # obj3_obj3_upper = np.repeat([sorted_upper_obj3], len(sorted_upper_obj3), axis=0)
    obj_3_upper_function = interp2d(sorted_upper_obj1, sorted_upper_obj2, sorted_upper_obj3)
    obj_3_upper_interpolation = obj_3_upper_function(all_obj1_vals, all_obj2_vals)[0, :]

    # obj3_obj3_lower = np.repeat([sorted_lower_obj3], len(sorted_lower_obj3), axis=0)
    obj_3_lower_function = interp2d(sorted_lower_obj1, sorted_lower_obj2, sorted_lower_obj3)
    obj_3_lower_interpolation = obj_3_lower_function(all_obj1_vals, all_obj2_vals)[0, :]

    upper_vertices = [[obj1_i, obj2_i, obj3_i] for obj1_i, obj2_i, obj3_i in zip(all_obj1_vals, all_obj2_vals, obj_3_upper_interpolation)]
    lower_vertices = [[obj1_i, obj2_i, obj3_i] for obj1_i, obj2_i, obj3_i in zip(all_obj1_vals, all_obj2_vals, obj_3_lower_interpolation)]
    all_vertices = [upper_vertices]+[lower_vertices]
    print(all_vertices)
    fig = plt.figure(figsize= (2.25, 2.25))
    ax = fig.add_subplot(projection='3d', computed_zorder=False)
    ax.scatter(
        xs=all_objectives[:, 0], ys=all_objectives[:, 1]*-1,
        zs=all_objectives[:, 2]*-1, color="black", zorder=1#marker="o", markersize=1,
        # color="black", zorder=1
    )
    ax.add_collection3d(Poly3DCollection(all_vertices, alpha=0.4, facecolors=sky_blue, zorder=2))
    ax.view_init(elev=10, azim=-115)
    ax.set_xlabel(obj_labels[0])
    ax.set_ylabel(obj_labels[1])
    ax.set_zlabel(obj_labels[2])
    plt.show()
    # plt.savefig(figure_path, bbox_inches="tight")

def plot_1D_all_cell_obj(figure_path, settings, all_cell_results_df_row):
    obj_all_cells = all_cell_results_df_row[
        settings["obj_labels"][0] + " for each cell"][0]
    obj_mean = all_cell_results_df_row[
        settings["obj_labels"][0] + "_mean"]
    x_vals = [1]*len(obj_all_cells)
    jittered_x = x_vals + 0.1*np.random.rand(
        len(x_vals))
    fig, ax = plt.subplots(1, 1, figsize= (2.25, 2))
    for i in range(len(obj_all_cells)):
        ax.plot(jittered_x[i], obj_all_cells[i], linestyle="None",
                marker="o", markersize=2)
    ax.plot(1, obj_mean, linestyle="None",
             marker = "o",markersize=2, color="k", 
             label="population mean"
    )
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_ylabel(settings["obj_labels"][0])
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_2D_all_cell_obj(figure_path, settings, all_cell_results_df_row):
    obj1_all_cells = all_cell_results_df_row[
        settings["obj_labels"][0] + " for each cell"][0]
    obj2_all_cells = all_cell_results_df_row[
        settings["obj_labels"][1] + " for each cell"][0]
    obj1_mean = all_cell_results_df_row[
        settings["obj_labels"][0] + "_mean"]
    obj2_mean = all_cell_results_df_row[
        settings["obj_labels"][1] + "_mean"]
    fig, ax = plt.subplots(1, 1, figsize= (2.25, 2))
    for i in range(len(obj1_all_cells)):
        plt.plot(obj1_all_cells[i], obj2_all_cells[i],
                 linestyle="None", marker = "o",
                 markersize=1
        )
    plt.plot(obj1_mean, obj2_mean, linestyle="None",
             marker = "o",markersize=1, color="k", 
             label="population mean"
    )
    # plt.legend()
    plt.xlabel(settings["obj_labels"][0])
    plt.ylabel(settings["obj_labels"][1])
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_all_cell_time_series(figure_path, settings, all_cell_results_df_row):

    all_cells_rep_rel = all_cell_results_df_row[
        "Rep_rel time series for each cell"]
    rep_rel_max = np.amax(np.asarray(all_cells_rep_rel))
    # print(rep_rel_max)
    rep_rel_mean = all_cell_results_df_row[
        "Rep_rel time series mean"]
    all_cells_peaks = all_cell_results_df_row["single_cell_peaks"]
    t = np.arange(0, settings["max_time"] + 1, 1)

    fig, axs = plt.subplots(1, 2, figsize= (5, 2.5))
    for i in range(len(all_cells_rep_rel)):
        axs[0].plot(t[:43], all_cells_rep_rel[i][:43],
                 lw="1", label=("cell " + str(i) + " peak " + str(round(all_cells_peaks[i], 1))))
        axs[1].plot(t[:43], all_cells_rep_rel[i][:43],
                 lw="1")
    axs[0].plot(t[:43], rep_rel_mean[:43], color="k",
             label="population mean", lw="2"
    )
    axs[1].plot(t[:43], rep_rel_mean[:43], color="k",
             label="population mean", lw="2"
    )
    axs[0].legend()
    axs[0].set_xlabel("time (hours)")
    axs[1].set_xlabel("time (hours)")
    axs[0].set_ylabel("Reporter_rel")
    axs[1].set_ylabel("Reporter_rel")
    axs[1].set_ylim(top=(rep_rel_max/2))
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_pulse_ensemble_time_series(figure_path, all_cell_results_df_row, reference_all_cell):

    all_cells_rep_rel = all_cell_results_df_row[
        "Rep_rel time series for each cell"
    ]
    time_points = [14, 18, 22, 26, 38, 42, 46]
    for i in range(len(all_cells_rep_rel)[0:1]):
        print(all_cells_rep_rel)
        for j, time in enumerate(time_points):
            time_list = [time]*len(all_cells_rep_rel)
            rep_rel_t_exp = all_cells_rep_rel[time_points]
            print(rep_rel_t_exp)


def plot_obj_progression_set(
        figure_path: str,
        n_gens: int,
        objectives_lists: np.ndarray,
        obj_label: str,
        y_lower_lim: float=None
):
    generations = np.arange(n_gens+1)

    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    for obj_list in objectives_lists:
        ax.plot(generations, obj_list)
    ax.set_xlabel("Generation")
    ax.set_ylabel(obj_label)
    if y_lower_lim:
        ax.set_ylim(bottom=y_lower_lim)
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_pareto_front_set(
        figure_path: str, 
        obj_df_list: list[pd.DataFrame],
        obj_labels: list,
):
        
        fig, ax = plt.subplots(1, 1, figsize= (3, 3))
        for obj_df in obj_df_list:
            if np.any(np.array(obj_df[obj_labels[0]].to_list()) < 0):
                obj_df[obj_labels[0]] = obj_df[
                    obj_labels[0]]*-1
            if np.any(np.array(obj_df[obj_labels[1]].to_list()) < 0):
                obj_df[obj_labels[1]] = obj_df[
                    obj_labels[1]]*-1
                
            sns.scatterplot(data=obj_df, x= obj_df[obj_labels[0]],
                            y= obj_df[obj_labels[1]], ax=ax, s=30)
            
        for dots in ax.collections:
            facecolors = dots.get_facecolors()
            dots.set_edgecolors(facecolors.copy())
            dots.set_facecolors('none')
            dots.set_linewidth(0.75)

        plt.xlabel(obj_labels[0])
        plt.ylabel(obj_labels[1])
        plt.savefig(figure_path, bbox_inches="tight")

def plot_pareto_front_set_3D(
        figure_path: str, 
        obj_df_list: list[pd.DataFrame],
        obj_labels: list,
):
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    fig = plt.figure(figsize= (4, 4))
    ax = fig.add_subplot(projection='3d')
    for i, obj_df in enumerate(obj_df_list):
        if np.any(np.array(obj_df[obj_labels[0]].to_list()) < 0):
            obj_df[obj_labels[0]] = obj_df[
                obj_labels[0]]*-1
        if np.any(np.array(obj_df[obj_labels[1]].to_list()) < 0):
            obj_df[obj_labels[1]] = obj_df[
                obj_labels[1]]*-1
        if np.any(np.array(obj_df[obj_labels[2]].to_list()) < 0):
            obj_df[obj_labels[2]] = obj_df[
                obj_labels[2]]*-1

        ax.scatter(
            xs=obj_df[obj_labels[0]], ys=obj_df[obj_labels[1]],
            zs=obj_df[obj_labels[2]], fc="none", color=colors[i],
            marker="o", s=30, linewidth=0.75
        )
    ax.view_init(elev=10, azim=-115)
    ax.set_xlabel(obj_labels[0])
    ax.set_ylabel(obj_labels[1])
    ax.set_zlabel(obj_labels[2])
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")



# path_results = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# path_sc = "Pulse_seed_pop_DsRED_inhibitor/3_obj/Optimized_hyperparams_vary_pop_3obj_opt_stdev_ngen70/"
# results_runs = "2024-07-05_Pulse_pop_DsRED_inhibitor_3obj_vary_pop_opt_hp_stdev_ngen70_seed_"

# pareto_obj_dfs_list = compare_parteo_fronts(path_results+path_sc, results_runs, ["t_pulse", "peak_rel", "prominence_rel"])
# plot_pareto_front_set_3D(path_results+path_sc+"pareto_front_set.svg", pareto_obj_dfs_list, ["t_pulse", "peak_rel", "prominence_rel"])
# plot_pareto_front_set(path_results+path_sc+"pareto_front_set.svg", pareto_obj_dfs_list, ["t_pulse", "peak_rel", "prominence_rel"])

#2024-06-07_Pulse_single_cell_DsRED_inhibitor_3_obj_vary_pop_opt_hp_stdev_ngen80_nseed4_run2_seed_


repo_path = "/Users/kdreyer/Library/CloudStorage/OneDrive-NorthwesternUniversity/KatieD_LL/GCAD_Collab/GA_results/"
results_path = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/Pulse_pop_DsRED_inhibitor_3obj_126h_ZF1_ZF2_seed_0/Results_analysis/"
file_name = "all_cell_selected_results_low_t_pulse.csv"
all_cell_results = pd.read_csv(repo_path+results_path+file_name)
# print(type(all_cell_results.loc[0, "Rep_rel time series for each cell"]))
all_cell_time_series_opt = all_cell_results.copy().iloc[:1]
# print(all_cell_time_series_opt)
all_cell_time_series_opt["Rep_rel time series for each cell"] = all_cell_time_series_opt["Rep_rel time series for each cell"].astype(object)
for index, row in all_cell_time_series_opt.iterrows():
#     for cell_time_series in row["Rep_rel time series for each cell"].split(", ["):
#         # cell_time_series = cell_time_series.strip(' []')
#         print(cell_time_series)
        # cell_time_series = [float(val) for val in cell_time_series]
        all_cell_time_series_opt.at[index, "Rep_rel time series for each cell"] = eval(row["Rep_rel time series for each cell"])
    # row["Rep_rel time series for each cell"] = [cell_time_series.strip(' []') for cell_time_series in row["Rep_rel time series for each cell"].split(",")]
    #                                           #  for val in cell_time_series]
print(type(all_cell_time_series_opt.at[0, "Rep_rel time series for each cell"]))
# print(type(all_cell_time_series_opt[0]))
# for all_cell_time_series in all_cell_time_series_opt:
    # for time_series in all_cell_time_series:
    #     print(time_series)
# figure_path = 
# plot_pulse_ensemble_time_series()