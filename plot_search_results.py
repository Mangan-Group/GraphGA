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
            fig, ax = plt.subplots(1, 1, figsize= (4, 4))
            sns.scatterplot(data=obj_df, x= obj_df[obj_labels[0]],
                            y= obj_df[obj_labels[1]], 
                            color="black", ax=ax)

        plt.xlabel(obj_labels[0])
        plt.ylabel(obj_labels[1])
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
            
    print(obj_df.tail(n=50))
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
    # plt.savefig(figure_path, bbox_inches="tight")


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

#############################################################
### Signal Conditioner Population Model ###
#############################################################

######## SC pop GA hypervolumes with DsRED inhibitor, 50 gens ########
# results_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_pop_DsRED_inhibitor/2023-10-31_Signal_Cond_pop_DsRED_inhibitor_seed_0/"
# hv_fname = "hypervolumes.pkl"
# with open(results_path + hv_fname, "rb") as fid:
#     hv = pickle.load(fid)
# print(hv)
# hv_plot_fname = "hypervolume_progression_full.svg"
# plot_hypervolume(results_path + hv_plot_fname, 50, hv)

######## SC pop GA hypervolumes/pareto with DsRED inhibitor, 60 gens ########
# results_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/"
# results_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_pop_DsRED_inhibitor/2023-12-01_Signal_Cond_pop_DsRED_inhibitor_ngen70_seed_0/"
# results_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_pop_DsRED_inhibitor/2023-12-02_Signal_Cond_pop_DsRED_inhibitor_ngen100_seed_0/"
# hv_fname = "hypervolumes.pkl"
# # final_obj_fname = "final_objectives_df_with_type.pkl"

# with open(results_path + hv_fname, "rb") as fid:
#     hv = pickle.load(fid)
# print(hv)
# hv_plot_fname = "hypervolume_progression.svg"
# plot_hypervolume(results_path+hv_plot_fname, 100, hv)
# obj_df = pd.read_pickle(results_path + final_obj_fname)
# graph_file_name = "final_population_pareto_front.svg"
# plot_pareto_front(
#         results_path + graph_file_name,
#         obj_df,
#         ["ON_rel", "FI_rel"],
#         types=False
#     )


# plot_hypervolumes(hv_list, "signal_conditioner", 40, 47.3, "/Users/kdreyer/Documents/Github/GraphGA/GA_results SC_seed_single_cell_DsRED_inhibitor/")
    

### Pulse hypervolumes ###
# repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"

# # prob_c = 0.32, prob_m = 0.57
# path_pulse1 = repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-10-31_Pulse_pop_DsRED_inhibitor_seed_0/hypervolumes.pkl"
# with open(path_pulse1, "rb") as fid:
#     hvs1 = pickle.load(fid)
# # print(hvs1)
# # plot_hypervolume(repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-10-31_Pulse_pop_DsRED_inhibitor_seed_0/hypervolume_progression.svg",
# #                  [hvs1], 50)

# # prob_c = 0.5, prob_m = 0.57
# path_pulse2 = repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-11-07_Pulse_pop_DsRED_inhibitor_c0_5_seed_0/hypervolumes.pkl"
# with open(path_pulse2, "rb") as fid:
#     hvs2 = pickle.load(fid)
# # print(hvs2)
# # plot_hypervolume(repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-11-07_Pulse_pop_DsRED_inhibitor_c0_5_seed_0/hypervolume_progression.svg",
# #                  [hvs2], 50)

# # prob_c = 0.75, prob_m = 0.57
# path_pulse3 = repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-11-07_Pulse_pop_DsRED_inhibitor_c0_75_seed_0/hypervolumes.pkl"
# with open(path_pulse3, "rb") as fid:
#     hvs3 = pickle.load(fid)
# # print(hvs3)
# # plot_hypervolume(repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-11-07_Pulse_pop_DsRED_inhibitor_c0_75_seed_0/hypervolume_progression.svg",
# #                  [hvs3], 50)

# # prob_c = 0.32, prob_m = 0.75
# path_pulse4 = repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-11-07_Pulse_pop_DsRED_inhibitor_m0_75_seed_0/hypervolumes.pkl"
# with open(path_pulse4, "rb") as fid:
#     hvs4 = pickle.load(fid)
# # print(hvs4)
# # plot_hypervolume(repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-11-07_Pulse_pop_DsRED_inhibitor_m0_75_seed_0/hypervolume_progression.svg",
# #                  [hvs4], 50)

# # prob_c = 0.32, prob_m = 1.0
# path_pulse5 = repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-11-07_Pulse_pop_DsRED_inhibitor_m1_seed_0/hypervolumes.pkl"
# with open(path_pulse5, "rb") as fid:
#     hvs5 = pickle.load(fid)
# # print(hvs5)
# # plot_hypervolume(repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-11-07_Pulse_pop_DsRED_inhibitor_m1_seed_0/hypervolume_progression.svg",
# #                  [hvs5], 50)

# # prob_c = 0.5, prob_m = 0.75
# path_pulse6 = repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-11-08_Pulse_pop_DsRED_inhibitor_c0_5_m0_75_seed_0/hypervolumes.pkl"
# with open(path_pulse6, "rb") as fid:
#     hvs6 = pickle.load(fid)
# # print(hvs5)
# plot_hypervolume(repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-11-08_Pulse_pop_DsRED_inhibitor_c0_5_m0_75_seed_0/hypervolume_progression.svg",
#                  [hvs6], 50)


#plot amplifier single cell 1D objective scatter plot
# amp_all_obj_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_single_cell_const_dose/2023-10-31_Amplifier_single_cell_seed_0/all_objectives.pkl"

# with open(amp_all_obj_path, "rb") as fid:
#     amp_all_obj = pickle.load(fid)

# amp_unique_obj = np.unique(amp_all_obj)

# plot_1D_obj_scatter("/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_single_cell_const_dose/2023-10-31_Amplifier_single_cell_seed_0/unique_obj_scatter.svg", amp_unique_obj, "ON_rel")

# obj_df = pd.DataFrame(data=[[1, 2, 3], [4, 5, 6], [7, 8, 9]], columns=["obj1", "obj2", "obj3"])
# # print(obj_df)

# plot_pareto_front3D("/Users/kdreyer/Desktop/3d_pareto.svg", obj_df, ["obj1", "obj2", "obj3"], False)
    
# repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# path_final_obj = "2024-02-05_Pulse_3obj_42h_seed_0/final_objectives_df.pkl"
# df_obj = pd.read_pickle(repo_path + path_final_obj)
# plot_pareto_front3D(repo_path+"2024-02-05_Pulse_3obj_42h_seed_0/final_popluation_pareto_front2.svg", df_obj, ["t_pulse", "peak_rel", "prominence_rel"], False)