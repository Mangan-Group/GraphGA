import numpy as np
import networkx as nx
import pandas as pd
import pickle
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import seaborn as sns
from scipy.interpolate import interp2d
from math import ceil, floor
from brokenaxes import brokenaxes
from load_files_pop import Z_20, Ref_pop20
# from Reference_pop import simulate_reference_time_series

plt.style.use('/Users/kdreyer/Documents/Github/GraphGA/paper.mplstyle.py')
orange_ = [i/255 for i in [230, 159, 0]]
sky_blue = [i/255 for i in [86, 180, 233]]
pink_ = [i/255 for i in [204, 121, 167]]
bluish_green = [i/255 for i in [0, 158, 115]]
vermillion = [i/255 for i in [213, 94, 0]]
yellow_ = [i/255 for i in [240, 228, 66]]
grey_ = [(i/255) for i in [150, 150, 150]]

combinatorial_yellow = [(i/255) for i in [221, 204, 119]]

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
    max_hv = max(np.array(hypervolumes_lists).flatten())
    fig, ax = plt.subplots(1, 1, figsize= (2, 2))
    for hv_list in hypervolumes_lists:
        ax.plot(generations, hv_list)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Hypervolume")
    ax.axhline(max_hv, xmin=0, xmax=generations[-1], linestyle="dashed", color="grey", label="max hv="+str(round(max_hv, 3)))
    if y_lower_lim:
        ax.set_ylim(bottom=y_lower_lim)
    # plt.show()
    ax.legend()
    plt.savefig(figure_path, bbox_inches="tight")

def plot_hypervolumes_set_vs_combo(
        figure_path: str,
        n_gens: int,
        hypervolumes_lists: np.ndarray,
        opt_combo_hv:float,
        selected_seed: int,
        y_lower_lim: float=None
):
    generations = np.arange(n_gens)
    mpl.rcParams["figure.autolayout"] = False
    # fig = plt.figure(figsize= (1.585, 1.6))
    fig, ax = plt.subplots(1, 1, figsize= (1.585, 1.6))
    # bax = brokenaxes(ylims=((0, 1), (35, 46)), hspace=0.1) #signal conditioner
    for i, hv_list in enumerate(hypervolumes_lists):
        if i == selected_seed:
            color_="k"
            zorder_=(i+1)*100
            linewidth_="0.75"
        else:
            color_=grey_
            zorder_=i
            linewidth_="0.5"

        ax.plot(generations, hv_list, linewidth=linewidth_,
                 color=color_, zorder=zorder_)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Hypervolume")
    ax.axhline(opt_combo_hv, xmin=0, xmax=generations[-1], 
               linestyle="dashed", color="k", label="opt hv="+str(round(opt_combo_hv, 3)))
    if y_lower_lim:
        ax.set_ylim(bottom=y_lower_lim)
    # plt.show()
    ax.legend()
    ax.set_xticks(np.arange(0, n_gens+1, 25))
    ax.set_ylim(bottom=0)
    ax.set_xlim(left=0)
    plt.savefig(figure_path, bbox_inches="tight")
        # bax.plot(generations, hv_list, linewidth=linewidth_,
        #          color=color_, zorder=zorder_)

    # bax.axhline(opt_combo_hv, xmin=0, xmax=generations[-1]+5,
    #             linestyle="dashed", linewidth="0.75", color="k",
    #             label="opt hv="+str(round(opt_combo_hv, 3)),
    #             zorder=20)
    # bax.set_xlim([0, generations[-1]+5])
    # bax.axs[1].set_yticks([0])
    # bax.axs[0].set_yticks([36, 38, 40, 42, 44])
    # bax.set_xlabel("Generation in GA")
    # bax.set_xticks(np.arange(0, n_gens+1, 2000))
    # bax.set_ylabel("Hypervolume")
    # if y_lower_lim:
    #     ax.set_ylim(bottom=y_lower_lim)
    # bax.legend()
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
            fig, ax = plt.subplots(1, 1, figsize= (1.955, 1.955))
            sns.scatterplot(data=obj_df, x= obj_df[obj_labels[0]],
                            y= obj_df[obj_labels[1]], 
                            color="black", ax=ax, s=6)

        plt.xlabel(obj_labels[0])
        plt.ylabel(obj_labels[1])
        plt.xticks([0, 20, 40, 60])
        # plt.yticks([0, 1, 2, 3])
        plt.xlim(left=0)
        plt.ylim(bottom=0, top=1.75)
        plt.yticks([0, 0.5, 1.0, 1.5])
        ax.set_box_aspect(1)
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
            
    prom_rel_exp = [2.345337964907094, 2.0821875761716178, 2.681593101318251, 3.053945976792732]

    obj_df.drop(obj_df[obj_df["prominence_rel"] == prom_rel_exp[0]].index, inplace=True, axis=0)
    obj_df.drop(obj_df[obj_df["prominence_rel"] == prom_rel_exp[1]].index, inplace=True, axis=0)
    obj_df.drop(obj_df[obj_df["prominence_rel"] == prom_rel_exp[2]].index, inplace=True, axis=0)
    obj_df.drop(obj_df[obj_df["prominence_rel"] == prom_rel_exp[3]].index, inplace=True, axis=0)

    fig = plt.figure(figsize= (1.955, 1.955))
    ax = fig.add_subplot(projection='3d')

    ax.scatter(
        xs=obj_df[obj_labels[0]], ys=obj_df[obj_labels[1]],
        zs=obj_df[obj_labels[2]], color="grey", s=2,
    )

    pulse_blue = [(i/255) for i in [51, 34, 136]]
    frac_p_exp = [0.45, 0.45, 0.4, 0.35]
    t_p_exp = [13, 11, 13, 13]
    ax.scatter(
        xs=frac_p_exp, ys=t_p_exp,
        zs=prom_rel_exp, color=pulse_blue, s=10,
    )

    # ax.view_init(elev=15, azim=-60)
    ax.view_init(elev=10, azim=-115)
    ax.set_xlabel(obj_labels[0])
    ax.set_ylabel(obj_labels[1])
    ax.set_zlabel(obj_labels[2])
    ax.set_xticks([0, 0.2, 0.4, 0.6])
    ax.set_yticks([0, 20, 40, 60])
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
    np.random.seed(0)
    jittered_x = x_vals + 0.1*np.random.rand(
        len(x_vals))
    fig, ax = plt.subplots(1, 1, figsize= (2, 1.75))
    ax.plot(jittered_x, obj_vals, linestyle="None",
             marker="o", markersize=1, color="k",zorder=1)
    ax.plot(1.05, max(obj_vals), linestyle="none", marker="o",
            markersize=2.5, color=combinatorial_yellow, zorder=2)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_ylabel(obj_labels[0])
    # if y_lower_lim:
    #     ax.set_ylim(bottom = y_lower_lim)
    ax.set_ylim(bottom = 0)
    ax.set_yticks([0, 20, 40, 60])
    ax.set_box_aspect(1)
    # plt.show()
    plt.savefig(figure_path, bbox_inches="tight")


def plot_1D_obj_confidence_interval(
        results_path: str,
        figure_path: str,
        CI_metric_max: float,
        obj_labels: list,
        y_lim_bottom: float=None 
):
    unique_objectives = pd.read_pickle(results_path+"unique_objectives.pkl")
    unique_objectives = unique_objectives.flatten()*-1
    max_objective = max(unique_objectives)
    np.random.seed(0)
    x_vals = [1]*len(unique_objectives)
    jittered_x = x_vals + 0.1*np.random.rand(
        len(x_vals))
    lower_bound = [max_objective-CI_metric_max]*(len(unique_objectives)+2)
    upper_bound = [max_objective]*(len(unique_objectives)+2)
    fig, ax = plt.subplots(1, 1, figsize= (2, 1.75))
    ax.plot(jittered_x, unique_objectives, linestyle="None",
             marker="o", markersize=1, color="black", zorder=1)
    x_CI_fill = np.append(jittered_x, [0.995, 1.105])
    x_CI_fill.sort()
    ax.fill_between(x_CI_fill, lower_bound, upper_bound, alpha=0.75, color=grey_, zorder=2,
                    linewidth=0.25)
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_xlim([0.995, 1.105])
    ax.set_ylabel(obj_labels[0])
    if y_lim_bottom is not None:
        y_lim = [y_lim_bottom, round(max_objective+CI_metric_max, 1)]
    else:
        y_lim = [round(max_objective-CI_metric_max, 1), round(max_objective+0.05, 2)]
        print(y_lim)
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
    # axs[0].legend()
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
    fig, ax = plt.subplots(1, 1, figsize= (2, 2))
    for i in range(len(all_cells_rep_rel)):
        rep_rel_t_exp = [all_cells_rep_rel[i][t] for t in time_points]
        ax.plot(time_points, rep_rel_t_exp, linestyle="none", marker="o", markersize=4, color="k")
        ax.set_xlabel("Time(h)")
        ax.set_ylabel("Rep_rel")
        ax.set_title("Pulse exp " + str(condition))
        ax.set_box_aspect(1)
    plt.savefig(figure_path+"pulse_"+str(condition)+"_time_series_paper.svg")


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


def plot_obj_progression_set(
        figure_path: str,
        n_gens: int,
        objectives_lists: np.ndarray,
        obj_label: str,
        selected_seed: int,
        opt_obj: float,
        y_lower_lim: float,
        y_ticks: list=None,
        # y_lower_lim: float=None
):
    generations = np.arange(n_gens+1)

    mpl.rcParams["figure.autolayout"] = False
    fig = plt.figure(figsize= (1.585, 1.6))
    bax = brokenaxes(ylims=((0, 1), (47, 64)), hspace=0.1)
    for i, obj_list in enumerate(objectives_lists):
        if i == selected_seed:
            color_="k"
            zorder_=(i+1)*100
            linewidth_="0.75"
        else:
            color_=grey_
            zorder_=i
            linewidth_="0.5"
        bax.plot(generations, obj_list, linewidth=linewidth_,
                 color=color_, zorder=zorder_)
    bax.axhline(opt_obj, xmin=0, xmax=generations[-1]+5,
                linestyle="dashed", linewidth="0.75", color="black", 
                label="opt obj="+str(round(opt_obj, 3)),
                zorder=20)
    # axt.set_ylim([55, 64])
    # axb.set_ylim([0, 10])
    bax.set_xlim([0, 55])
    bax.axs[1].set_yticks([0])
    bax.axs[0].set_yticks([50, 55, 60, 65])
    # bax.draw_diags()
    bax.set_xlabel("Generation in GA")

    bax.set_xticks(np.arange(0, n_gens+1, 10))
    # ax.set_xlim(left=0)
    bax.set_ylabel("objective ("+obj_label+")")
    # bax.axs[0].set_box_aspect(1)
    # if y_lower_lim:
    #     ax.set_ylim(bottom=y_lower_lim)
    # else:
    #     ax.set_ylim(bottom=0)
    # plt.show()
    bax.legend()
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