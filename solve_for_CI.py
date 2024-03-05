import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator
# from run_test_case import run_combinitorial_pop_samples
from load_Z_mat_samples import Z_mat_list

seed = 0
np.random.seed(seed)
def plot_1D_obj_scatter_for_ci(
        figure_path: str,
        obj_vals: np.ndarray,
        obj_label: str,
        error: list=None,
        y_lower_lim: int=None
):
    # if obj_vals.flatten()[0] < 0:
    #     obj_vals = obj_vals*-1
    x_vals = [1]*len(obj_vals[0])
    jittered_x = x_vals + 0.1*np.random.rand(
        len(x_vals))
    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    ax.plot(jittered_x, obj_vals[0], linestyle="None",
             marker="o", color="gray")
    if error:
        ax.errorbar(jittered_x[:len(obj_vals[1])], obj_vals[1], yerr=error, 
                    linestyle="None", fillstyle="none", marker="o", 
                    color="black", capsize=2)
    # else:
    #     ax.plot(jittered_x[:len(obj_vals[1])], obj_vals[1], linestyle="None",
    #             fillstyle="none", marker="o", color="black")
    ax.set_xticklabels([])
    ax.set_xticks([])
    ax.set_ylabel(obj_label)
    if y_lower_lim:
        ax.set_ylim(bottom = y_lower_lim, top=71)
    plt.show()
    # plt.savefig(figure_path, bbox_inches="tight")

# import amplifier const dose results
# path_amp_const_obj = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_const_dose/2023-10-31_Amplifier_pop_const_dose_seed_0/unique_objectives.pkl"
# path_amp_const_circuits = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_const_dose/2023-10-31_Amplifier_pop_const_dose_seed_0/unique_circuits.pkl"
# with open(path_amp_const_obj, "rb") as fid:
#     unique_obj = pickle.load(fid)

# with open(path_amp_const_circuits, "rb") as fid:
#     unique_circuits = pickle.load(fid)

# unique_obj = unique_obj*-1
# unique_obj = unique_obj.flatten()

# unique_circuits = unique_circuits.flatten()

# unique_obj_ascending_idx = unique_obj.argsort()
# unique_obj_descending_idx = unique_obj_ascending_idx[::-1]
# unique_obj = unique_obj[unique_obj_descending_idx]

# unique_circuits = unique_circuits[unique_obj_descending_idx]
# # print(unique_obj)
# unique_circuits[0].plot_graph()

# # index_unique_obj_high = np.argwhere(unique_obj >= 68.8)
# # index_unique_obj_high = np.argwhere(unique_obj >= 67.5)
# index_unique_obj_high = np.argwhere(unique_obj >= 67)
# # print(index_unique_obj_high)
# unique_obj_high = unique_obj[index_unique_obj_high.flatten()]
# unique_circuits_high = unique_circuits[index_unique_obj_high.flatten()]
# # print(unique_circuits_high)
# unique_circuits_high = np.asarray(unique_circuits_high).reshape(len(unique_circuits_high), 1)
# # print(unique_circuits_high)
# unique_obj = unique_obj.tolist()
# unique_obj_high = unique_obj_high.tolist()
# # print(len(unique_obj_high))
# # print(unique_obj_high)
# all_obj = [unique_obj, unique_obj_high]
# # print(len(all_obj[1]))

# path_before_ci = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_const_dose/2023-10-31_Amplifier_pop_const_dose_seed_0/high_obj_scatter_67.svg"
# # plot_1D_obj_scatter_for_ci(path_before_ci, all_obj, "ON_rel", y_lower_lim=66)

# folder_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_const_dose/2023-10-31_Amplifier_pop_const_dose_seed_0/"
# file_name = "high_unique_objs_67.pkl"
# with open(folder_path+file_name, "rb") as fid:
#     high_obj = pickle.load(fid)
# # with open(folder_path+file_name, "wb") as fid:
# #     pickle.dump(unique_obj_high, fid)
# print(high_obj)
# file_name = "high_unique_circuits_67.pkl"
# with open(folder_path+file_name, "rb") as fid:
#     high_circuit = pickle.load(fid)
# high_circuit[0][0].plot_graph()
# # with open(folder_path+file_name, "wb") as fid:
# #     pickle.dump(unique_circuits_high, fid)
          
# settings = {
#     "promo_node":"P1",
#     "dose_specs": [75, 75, 5],
#     "max_part": 2,
#     "inhibitor": True,
#     "DsRed_inhibitor": False,
#     "num_dict": {1: 26, 2: 26},
#     "n_gen": 40,
#     "probability_crossover": 0.55,
#     "probability_mutation": 1.0,
#     "mutate_dose": False,
#     "pop": True,
#     "CI": None,
#     "num_processes": 1,
#     "get_unique": False,
#     "plot": False,
#     "seed": 0,
#     "repository_path": "/Users/kdreyer/Documents/Github/GraphGA/",
#     "folder_name": "Amplifier_pop_const_dose_sampling_high_circuits_67",
#     "results_path": folder_path+file_name
# }

# # run_combinitorial_pop_samples(Amplifier, settings, Z_mat_list)

# z_mat_sampling_file_name = "2023-11-01_Amplifier_pop_const_dose_sampling_high_circuits_67/Z_mat_sampling_high_obj_circuits.pkl"
# with open(folder_path+z_mat_sampling_file_name, "rb") as fid:
#     z_mat_sampling = pickle.load(fid)

# ON_rel_ranges = z_mat_sampling["objectives_range"].tolist()
# # # print(ON_rel_ranges)
# fig_path = folder_path+"2023-11-01_Amplifier_pop_const_dose_sampling_high_circuits_67/high_obj_scatter_ci.svg"
# plot_1D_obj_scatter_for_ci(fig_path, all_obj, "ON_rel", error = ON_rel_ranges, y_lower_lim=66)

##################################################################################################################
##################################################################################################################
##################################################################################################################
# import amplifier vary dose results
# path_amp_vary_obj = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/unique_objectives.pkl"
# path_amp_vary_circuits = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/unique_circuits.pkl"
# with open(path_amp_vary_obj, "rb") as fid:
#     unique_obj = pickle.load(fid)

# with open(path_amp_vary_circuits, "rb") as fid:
#     unique_circuits = pickle.load(fid)

# unique_obj = unique_obj*-1
# # unique_obj = unique_obj.flatten()

# # unique_circuits = unique_circuits.flatten()

# # unique_obj_ascending_idx = unique_obj.argsort()
# # unique_obj_descending_idx = unique_obj_ascending_idx[::-1]
# # unique_obj = unique_obj[unique_obj_descending_idx]

# # unique_circuits = unique_circuits[unique_obj_descending_idx]
# # # print(unique_obj[:40])

# index_unique_obj_high = np.argwhere((unique_obj[:, 0] >=68)).flatten()
# # print(index_unique_obj_high)
# # print(index_unique_obj_high)
# unique_obj_high = unique_obj[index_unique_obj_high]
# # print(len(unique_obj_high))
# unique_circuits_high = unique_circuits[index_unique_obj_high]
# print(unique_circuits_high)
# unique_circuits_high = np.asarray(unique_circuits_high).reshape(len(unique_circuits_high), 1)
# for i, obj in enumerate(unique_obj_high):
#     # p1 -> Z2 auto -> Z6 auto -> rep (idx 19), 75ng each
#     print("circuit: ", str(i))
#     print("edges, doses: ", unique_circuits_high[i][0].edge_list)
#     print("obj:", obj)
#     unique_circuits_high[i][0].plot_graph()

# print("edges, doses: ", unique_circuits_high[19][0].edge_list)
# print("obj:", unique_obj_high[19])
# unique_circuits_high[19][0].plot_graph()


# unique_obj = unique_obj.tolist()
# unique_obj_high = unique_obj_high.tolist()
# # print(len(unique_obj_high))
# # print(unique_obj_high)
# all_obj = [unique_obj, unique_obj_high]
# print(len(all_obj[1]))

# path_before_ci = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/high_obj_scatter_68.svg"
# plot_1D_obj_scatter_for_ci(" ", all_obj, "ON_rel", y_lower_lim=39)

folder_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/"
# file_name = "high_unique_objs_68.pkl"
# with open(folder_path+file_name, "rb") as fid:
#     high_obj = pickle.load(fid)
# # # with open(folder_path+file_name, "wb") as fid:
# # #     pickle.dump(unique_obj_high, fid)

# file_name = "high_unique_circuits_68.pkl"
# with open(folder_path+file_name, "rb") as fid:
#     high_circuit = pickle.load(fid)
# high_circuit[0][0].plot_graph()
# # with open(folder_path+file_name, "wb") as fid:
# #     pickle.dump(unique_circuits_high, fid)
          
# settings = {
#     "promo_node":"P1",
#     "dose_specs": [5, 75, 5],
#     "max_part": 2,
#     "inhibitor": True,
#     "DsRed_inhibitor": False,
#     "num_dict": {1: 26, 2: 26},
#     "n_gen": 50,
#     "probability_crossover": 0.55,
#     "probability_mutation": 1.0,
#     "mutate_dose": False,
#     "pop": True,
#     "CI": None,
#     "num_processes": 1,
#     "get_unique": False,
#     "plot": False,
#     "seed": 0,
#     "repository_path": "/Users/kdreyer/Documents/Github/GraphGA/",
#     "folder_name": "Amplifier_pop_vary_dose_sampling_high_circuits_68",
#     "results_path": folder_path+file_name
# }

# # run_combinitorial_pop_samples(Amplifier, settings, Z_mat_list)

# z_mat_sampling_file_name = "2023-11-01_Amplifier_pop_vary_dose_sampling_high_circuits_68/Z_mat_sampling_high_obj_circuits.pkl"
# with open(folder_path+z_mat_sampling_file_name, "rb") as fid:
#     z_mat_sampling = pickle.load(fid)
# print(z_mat_sampling)

# print(z_mat_sampling["objectives_range"].max)
# ON_rel_ranges = z_mat_sampling["objectives_range"].tolist()
# # print(ON_rel_ranges)
# fig_path = folder_path+"2023-11-01_Amplifier_pop_vary_dose_sampling_high_circuits_68/high_obj_scatter_ci.svg"
# plot_1D_obj_scatter_for_ci(fig_path, all_obj, "ON_rel", error = ON_rel_ranges, y_lower_lim=66)
