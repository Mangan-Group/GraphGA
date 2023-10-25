import numpy as np
import networkx as nx
import pickle
import pandas as pd
import matplotlib.pyplot as plt

plt.style.use('/Users/kdreyer/Documents/Github/GraphGA/paper.mplstyle.py')



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
    

def plot_metric(metric_list: np.ndarray, label: str):

    generations = np.arange(len(metric_list))

    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    ax.plot(generations, metric_list)
    ax.set_xlabel("Generation")
    ax.set_ylabel(label)
    plt.savefig("Amplifier_"+label+"_"+str(len(generations)-1)+".svg")

def plot_hypervolumes(hypervolumes_list: np.ndarray,
                      testcase: str, gens_to_plot: int,
                      y_lower: float, path: str):

    generations = np.arange(50)[-gens_to_plot:]

    # gens_to_plot = generations[-gens_to_plot:]
    # hvs_to_plot = [l[-gens_to_plot:] for l in hypervolumes_list]

    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    for hvs in hypervolumes_list:
        ax.plot(generations, hvs)
    ax.set_xlabel("Generation")
    ax.set_ylabel("Hypervolume")
    # ax.set_ylim(bottom=y_lower)
    plt.show()
    # plt.savefig(path + testcase + "_hypervolumes.svg")

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
hv_list = []
for seed in range(0, 20):
    path = ("/Users/kdreyer/Documents/Github/GraphGA/GA_results/" +
            "SC_seed_single_cell_DsRED_inhibitor/2023-10-24_Signal_cond_single_cell_DsRED_inhibitor_seed_" +
            str(seed) + "/hypervolumes.pkl")
    with open(path, "rb") as fid:
        hv = pickle.load(fid)
    hv_list.append(hv[-40:])

# print(hv_list)

plot_hypervolumes(hv_list, "signal_conditioner", 40, 47.3, "/Users/kdreyer/Documents/Github/GraphGA/GA_results SC_seed_single_cell_DsRED_inhibitor/")
    
