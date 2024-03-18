import numpy as np
import networkx as nx
import pickle
import pandas as pd
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
            fig, ax = plt.subplots(1, 1, figsize= (4, 4))
            sns.scatterplot(data=obj_df, x= obj_df[obj_labels[0]],
                            y= obj_df[obj_labels[1]], hue='type', 
                            palette=palette, ax=ax)
            
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
    fig = plt.figure(figsize= (2.25, 2))
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


#plot amplifier single cell 1D objective scatter plot
# amp_all_obj_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/all_objectives.pkl"

# with open(amp_all_obj_path, "rb") as fid:
#     amp_all_obj = pickle.load(fid)

# amp_unique_obj = np.unique(amp_all_obj)

# plot_1D_obj_scatter("/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/obj_scatter_plot_small.svg", amp_all_obj, "ON_rel")

# sc_obj_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/final_objectives_df_with_type.pkl"
# sc_pareto = pd.read_pickle(sc_obj_path)
# print(sc_pareto)
# plot_pareto_front("/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/final_population_pareto_front_small.svg", sc_pareto, ["ON_rel", "FI_rel"], False)
# sc_pareto = sc_pareto[["ON_rel", "FI_rel"]]
# print(sc_pareto)

# sc_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/"
# unique_obj_df = pd.read_pickle(sc_path + "unique_objectives_df.pkl")*-1
# graph_file_name = "unique_obj_scatter_plot.svg"
# plot_pareto_front(
#     sc_path + "/" + graph_file_name,
#     unique_obj_df,
#     ["ON_rel", "FI_rel"],
#     types=False
# )
# print(unique_obj_df[(unique_obj_df["FI_rel"] > 2.5) & (unique_obj_df["FI_rel"] <= 3.0)])
# print(unique_obj_df)
