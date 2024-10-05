import numpy as np
import networkx as nx
import pandas as pd
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import seaborn as sns
from scipy.interpolate import interp2d
from math import ceil, floor
from load_files_pop import Z_20, Ref_pop20
# from Reference_pop import simulate_reference_time_series

plt.style.use('/Users/kdreyer/Documents/Github/GraphGA/paper.mplstyle.py')
orange_ = [i/255 for i in [230, 159, 0]]
sky_blue = [i/255 for i in [86, 180, 233]]
pink_ = [i/255 for i in [204, 121, 167]]
bluish_green = [i/255 for i in [0, 158, 115]]
vermillion = [i/255 for i in [213, 94, 0]]
yellow_ = [i/255 for i in [240, 228, 66]]

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

    fig, ax = plt.subplots(1, 1, figsize= (2, 2))
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
    fig, ax = plt.subplots(1, 1, figsize= (1.9, 1.75))
    ax.plot(jittered_x, obj_vals, linestyle="None",
             marker="o", markersize=1, color="k",zorder=1)
    ax.plot(1.05, max(obj_vals), linestyle="none", marker="o",
            markersize=1.25, color=yellow_, zorder=2)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_ylabel(obj_labels[0])
    if y_lower_lim:
        ax.set_ylim(lower = y_lower_lim)
    # ax.set_yticks([0, 20, 40, 60])
    ax.set_box_aspect(1)
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


### not currently executable 
# def plot_3D_obj_confidence_interval(
#         objectives: pd.DataFrame,
#         results_path: str,
#         figure_path: str,
#         CI_metric_maxes: list,
#         obj_labels: list,
# ):
#     all_objectives = pd.read_pickle(results_path+"all_objectives.pkl")
    
#     obj1_CI = CI_metric_maxes[0]
#     obj2_CI = CI_metric_maxes[1]
#     obj3_CI = CI_metric_maxes[2]

#     upper_obj1 = np.array(objectives[obj_labels[0]])
#     sorted_upper_idx = np.argsort(upper_obj1)
#     sorted_upper_obj1 = upper_obj1[sorted_upper_idx]
#     lower_obj1 = np.array([i-obj1_CI for i in (objectives[obj_labels[0]])])
#     sorted_lower_idx = np.argsort(lower_obj1)
#     sorted_lower_obj1 = lower_obj1[sorted_lower_idx]

#     upper_obj2 = np.array(objectives[obj_labels[1]]*-1)
#     sorted_upper_obj2 = upper_obj2[sorted_upper_idx]
#     lower_obj2 = np.array([i-obj2_CI for i in (objectives[obj_labels[1]]*-1)])
#     sorted_lower_obj2 = lower_obj2[sorted_lower_idx]

#     upper_obj3 = np.array(objectives[obj_labels[2]]*-1)
#     sorted_upper_obj3 = upper_obj3[sorted_upper_idx]
#     lower_obj3 = np.array([i-obj3_CI for i in (objectives[obj_labels[2]]*-1)])
#     sorted_lower_obj3 = lower_obj3[sorted_lower_idx]

#     # obj1_obj1_upper, obj2_obj2_upper = np.meshgrid(sorted_upper_obj1, sorted_upper_obj2)
#     # obj3_obj3_upper = np.array(sorted_upper_obj3*len(sorted_upper_obj3)).reshape(len(sorted_upper_obj3), len(sorted_upper_obj3))
#     # obj_3_upper_interpolation = interp2d(obj1_obj1_upper, obj2_obj2_upper, obj3_obj3_upper)
#     all_obj1_vals = np.sort(np.concatenate([upper_obj1, lower_obj1]))
#     all_obj2_vals = np.sort(np.concatenate([upper_obj2, lower_obj2]))

#     # obj3_obj3_upper = np.repeat([sorted_upper_obj3], len(sorted_upper_obj3), axis=0)
#     obj_3_upper_function = interp2d(sorted_upper_obj1, sorted_upper_obj2, sorted_upper_obj3)
#     obj_3_upper_interpolation = obj_3_upper_function(all_obj1_vals, all_obj2_vals)[0, :]

#     # obj3_obj3_lower = np.repeat([sorted_lower_obj3], len(sorted_lower_obj3), axis=0)
#     obj_3_lower_function = interp2d(sorted_lower_obj1, sorted_lower_obj2, sorted_lower_obj3)
#     obj_3_lower_interpolation = obj_3_lower_function(all_obj1_vals, all_obj2_vals)[0, :]

#     upper_vertices = [[obj1_i, obj2_i, obj3_i] for obj1_i, obj2_i, obj3_i in zip(all_obj1_vals, all_obj2_vals, obj_3_upper_interpolation)]
#     lower_vertices = [[obj1_i, obj2_i, obj3_i] for obj1_i, obj2_i, obj3_i in zip(all_obj1_vals, all_obj2_vals, obj_3_lower_interpolation)]
#     all_vertices = [upper_vertices]+[lower_vertices]
#     print(all_vertices)
#     fig = plt.figure(figsize= (2.25, 2.25))
#     ax = fig.add_subplot(projection='3d', computed_zorder=False)
#     ax.scatter(
#         xs=all_objectives[:, 0], ys=all_objectives[:, 1]*-1,
#         zs=all_objectives[:, 2]*-1, color="black", zorder=1#marker="o", markersize=1,
#         # color="black", zorder=1
#     )
#     ax.add_collection3d(Poly3DCollection(all_vertices, alpha=0.4, facecolors=sky_blue, zorder=2))
#     ax.view_init(elev=10, azim=-115)
#     ax.set_xlabel(obj_labels[0])
#     ax.set_ylabel(obj_labels[1])
#     ax.set_zlabel(obj_labels[2])
#     plt.show()
    # plt.savefig(figure_path, bbox_inches="tight")


def plot_1D_all_cell_obj(
        figure_path:str,
        settings:dict,
        all_cell_results_df_row:pd.Series
):
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


def plot_2D_all_cell_obj(
        figure_path:str, 
        settings:dict,
        all_cell_results_df_row:pd.Series
):
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


def plot_all_cell_time_series(
        figure_path:str,
        settings:dict,
        all_cell_results_df_row:pd.Series
):
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


def plot_pulse_ensemble_time_series(
        figure_path:str,
        all_cell_results_df_row:pd.Series,
        condition:int
):

    all_cells_rep_rel = all_cell_results_df_row[
        "Rep_rel time series for each cell"
    ]
    time_points = [14, 18, 22, 26, 38, 42, 46]
    fig, ax = plt.subplots(1, 1, figsize= (3, 3))
    for i in range(len(all_cells_rep_rel)):
        rep_rel_t_exp = [all_cells_rep_rel[i][t] for t in time_points]
        ax.plot(time_points, rep_rel_t_exp, linestyle="none", marker="o", markersize=4, color="k")
        ax.set_xlabel("Time(h)")
        ax.set_ylabel("Rep_rel")
        ax.set_title("Pulse exp " + str(condition))
        ax.set_box_aspect(1)
    plt.savefig(figure_path+"pulse_"+str(condition)+"_time_series.svg")


def plot_pulse_ensemble_violin(
        figure_path:str,
        all_cell_results_df_row:pd.Series,
        condition:int
):
    
    all_cells_rep_rel = all_cell_results_df_row[
        "Rep_rel time series for each cell"
    ]
    time_points = [14, 18, 22, 26, 38, 42, 46]
    pulse_rel_all_cells = []
    fig, ax = plt.subplots(1, 1, figsize= (3, 3))
    for i in range(len(all_cells_rep_rel)):
        pulse_rep_rel = all_cells_rep_rel[i]
        rep_rel_t_exp = [pulse_rep_rel[t] for t in time_points]
        pulse_rel_all_cells.append(rep_rel_t_exp)
    pulse_rel_t_exp_df = pd.DataFrame(pulse_rel_all_cells, columns = [str(i) + "h" for i in time_points])
    pulse_rel_t_exp_df_T = pulse_rel_t_exp_df.transpose().copy()
    pulse_rel_t_exp_df_T["Time (h)"] = pulse_rel_t_exp_df_T.index
    pulse_rel_t_exp_df_T_plot = pd.melt(frame=pulse_rel_t_exp_df_T,
                                      id_vars="Time (h)",
                                      var_name="column_name",
                                      value_name="Rep_rel")
    plot = sns.violinplot(data=pulse_rel_t_exp_df_T_plot, x="Time (h)", y="Rep_rel", ax=ax)
    plot.set_xticks(range(len(pulse_rel_t_exp_df_T.index)))
    plot.set_xticklabels(time_points)
    ax.set_xlabel("Time(h)")
    ax.set_ylabel("Rep_rel")
    ax.set_title("Pulse exp " + str(condition))
    ax.set_box_aspect(1)
    plt.savefig(figure_path+"pulse_"+str(condition)+"_violin_plot.svg")


def plot_ref_ensemble_time_series(
        figure_path:str,
        ref_all_cell_time_series:np.ndarray
):
    ref_on = Ref_pop20["P1"]["on"]
    time_points = [14, 18, 22, 26, 38, 42, 46]
    ref_rel_all_cells = []
    fig, ax = plt.subplots(1, 1, figsize= (3, 3))
    for i in range(len(ref_all_cell_time_series)):
        ref_rep_rel = ref_all_cell_time_series[i]/ref_on
        rep_rel_t_exp = [ref_rep_rel[t] for t in time_points]
        ref_rel_all_cells.append(rep_rel_t_exp)
        ax.plot(time_points, rep_rel_t_exp, linestyle="none", marker="o", markersize=4, color="k")
        ax.set_xlabel("Time(h)")
        ax.set_ylabel("Rep_rel")
        ax.set_title("Reference")
        ax.set_box_aspect(1)
    plt.savefig(figure_path+"reference_time_series.svg")


def plot_ref_ensemble_violin(
        figure_path:str,
        ref_all_cell_time_series:np.ndarray
):
    ref_on = Ref_pop20["P1"]["on"]
    time_points = [14, 18, 22, 26, 38, 42, 46]
    ref_rel_all_cells = []
    fig, ax = plt.subplots(1, 1, figsize= (3, 3))
    for i in range(len(ref_all_cell_time_series)):
        ref_rep_rel = ref_all_cell_time_series[i]/ref_on
        rep_rel_t_exp = [ref_rep_rel[t] for t in time_points]
        ref_rel_all_cells.append(rep_rel_t_exp)
    ref_rel_t_exp_df = pd.DataFrame(ref_rel_all_cells, columns = [str(i) + "h" for i in time_points])
    ref_rel_t_exp_df_T = ref_rel_t_exp_df.transpose().copy()
    ref_rel_t_exp_df_T["Time (h)"] = ref_rel_t_exp_df_T.index
    ref_rel_t_exp_df_T_plot = pd.melt(frame=ref_rel_t_exp_df_T,
                                      id_vars="Time (h)",
                                      var_name="column_name",
                                      value_name="Rep_rel")
    plot = sns.violinplot(data=ref_rel_t_exp_df_T_plot, x="Time (h)", y="Rep_rel", ax=ax)
    plot.set_xticks(range(len(ref_rel_t_exp_df_T.index)))
    plot.set_xticklabels(time_points)
    ax.set_xlabel("Time(h)")
    ax.set_ylabel("Rep_rel")
    ax.set_title("Reference")
    ax.set_box_aspect(1)
    plt.savefig(figure_path+"reference_violin_plot.svg")


def plot_obj_progression_set(
        figure_path: str,
        n_gens: int,
        objectives_lists: np.ndarray,
        obj_label: str,
        y_ticks: list=None,
        y_lower_lim: float=None
):
    generations = np.arange(n_gens+1)

    cmap = plt.get_cmap("plasma_r", len(objectives_lists))
    fig, ax = plt.subplots(1, 1, figsize= (1.9, 1.75))
    for i, obj_list in enumerate(objectives_lists):
        ax.plot(generations, obj_list, color=cmap(i), label="seed"+str(i))
    ax.set_xlabel("Generation")
    ax.set_xticks(np.arange(0, n_gens+1, 10))
    ax.set_xlim(left=0)
    ax.set_ylabel(obj_label)
    if y_lower_lim:
        ax.set_ylim(bottom=y_lower_lim)
    else:
        ax.set_ylim(bottom=0)
    # plt.show()
    plt.legend()
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
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd',
              '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
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