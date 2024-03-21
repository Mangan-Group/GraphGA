import pickle
import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from pulse_generator_problem import PulseGenerator
from define_circuit import Topo
from amplifier_problem import Amplifier
from load_files_pop import Ref_pop20
from scipy.stats import sem
from copy import deepcopy
from diversity_metrics import first_seen

plt.style.use('/Users/kdreyer/Documents/Github/GraphGA/paper.mplstyle.py')

#############################################################
### Pulse population model and single cell model ###
#############################################################
# pulse = PulseGenerator(
#         promo_node="P1", dose_specs=[5, 75, 5], max_part=2,
#         inhibitor=True, DsRed_inhibitor=True, num_dict={1: 1, 2: 2}, 
#         n_gen=1, probability_crossover=1, probability_mutation=1,
#         mutate_dose=True, pop=True
# )

# all_cell_results = pulse.all_cells

# load files- population model
# repo_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# file_path_all_cells = "Pulse_seed_pop_DsRED_inhibitor/2023-10-31_Pulse_pop_DsRED_inhibitor_seed_0/All_circuits_all_cell_results.pkl"
# file_path_all_obj = "Pulse_seed_pop_DsRED_inhibitor/2023-10-31_Pulse_pop_DsRED_inhibitor_seed_0/all_objectives.pkl"
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/2023-10-31_Pulse_pop_DsRED_inhibitor_seed_0/final_objectives_df_with_type.pkl"
# file_path_final_circuits = "Pulse_seed_pop_DsRED_inhibitor/2023-10-31_Pulse_pop_DsRED_inhibitor_seed_0/final_population.pkl"

#load files- population model 126h
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/2023-11-27_Pulse_126h_pop_DsRED_inhibitor_seed_0/final_objectives_df_with_type.pkl"
# file_path_final_circuits = "Pulse_seed_pop_DsRED_inhibitor/2023-11-27_Pulse_126h_pop_DsRED_inhibitor_seed_0/final_population.pkl"

# load files- populstion model t_pulse
# file_path_all_cells = "Pulse_seed_pop_DsRED_inhibitor/2023-11-28_Pulse_pop_DsRED_inhibitor_t_pulse_seed_0/All_circuits_all_cell_results.pkl"
# file_path_all_obj = "Pulse_seed_pop_DsRED_inhibitor/2023-11-28_Pulse_pop_DsRED_inhibitor_t_pulse_seed_0/all_objectives.pkl"
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/2023-11-28_Pulse_pop_DsRED_inhibitor_t_pulse_seed_0/final_objectives_df_with_type.pkl"

# load files- population model t_pulse 126h
# file_path_all_cells = "Pulse_seed_pop_DsRED_inhibitor/2023-11-28_Pulse_126h_pop_DsRED_inhibitor_t_pulse_seed_0/All_circuits_all_cell_results.pkl"
# file_path_all_obj = "Pulse_seed_pop_DsRED_inhibitor/2023-11-28_Pulse_126h_pop_DsRED_inhibitor_t_pulse_seed_0/all_objectives.pkl"
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/2023-11-28_Pulse_126h_pop_DsRED_inhibitor_t_pulse_seed_0/final_objectives_df_with_type.pkl"

# #load files- single cell
# # file_path_final_obj = "Pulse_seed_single_cell_DsRED_inhibitor/2023-10-31_Pulse_single_cell_DsRED_inhibitor_seed_0/final_objectives_df_with_type.pkl"
# # file_path_final_circuits = "Pulse_seed_single_cell_DsRED_inhibitor/2023-10-31_Pulse_single_cell_DsRED_inhibitor_seed_0/final_population.pkl"


# df_final_obj = pd.read_pickle(repo_path + file_path_final_obj)
# # print(df_final_obj)
# obj_labels: list=["t_pulse (hr)",
#     "prominence_rel"
# ]
# # print(df_final_obj[obj_labels[1]].to_list()[0])

# if np.any(np.array(df_final_obj[obj_labels[1]].to_list()) < 0):
#     print("contains negative values")

# if df_final_obj[obj_labels[1]].to_list()[0] < 0:
#     df_final_obj[obj_labels[1]] = df_final_obj[
#         obj_labels[1]]*-1
# print(df_final_obj)

# with open(repo_path + file_path_all_obj, "rb") as fid:
#             all_obj = pickle.load(fid)

# with open(repo_path + file_path_final_circuits, "rb") as fid:
#             final_circuits = pickle.load(fid)
# final_circuits[3, 0].plot_graph()
# print(final_circuits[3, 0].dose)

# df_all_cells = pd.read_pickle(repo_path + file_path_all_cells)
# print(df_all_cells)

# get unique values and indices in all_obj
# unique_all_obj, unique_all_obj_idx = np.unique(all_obj, return_index=True, axis=0)
# print(len(unique_all_obj))
# print(unique_all_obj_idx)

# # get rows of df_all_cells that yielded unique_all_obj
# unique_df_all_cells = df_all_cells.loc[unique_all_obj_idx]
# unique_df_all_cells = unique_df_all_cells.reset_index()

# # get unique objective values on pareto front
# final_objs_list = df_final_obj.values.tolist()
# unique_final_objs, unique_final_obj_idx = np.unique(final_objs_list, return_index=True, axis=0)
# print(unique_final_objs)
# unique_final_circuits = final_circuits[unique_final_obj_idx]
# # print(unique_final_circuits)
# mean_ts_lists = []
# t = np.arange(127)
# plt.figure()
# for circuit in unique_final_circuits[1:]:
#     t, mean_ts = pulse.simulate(circuit[0])
#     mean_ts_rel = pulse.calc_rep_rel(circuit[0], mean_ts)
#     mean_ts_lists.append(mean_ts_rel)
#     metrics = pulse.func(circuit[0])
#     print("t_pulse = ", metrics[0])
#     print("prominence = ", metrics[1]*-1)
#     plt.plot(t, mean_ts_rel)
# plt.show() 
# pareto_df_all_cells = pulse.all_cells
# # print(all_obj)

# # get indices in unique_all_obj where values = a set
# # of unique values on pareto front
# pareto_all_obj_indices = []
# i = 0
# for row1 in unique_all_obj:
#     for row2 in unique_final_objs:
#         if (row1[0] == row2[0]) & (row1[1] == row2[1]):
#             # print(row1, row2)
#             pareto_all_obj_indices.append(i)
#             break
#     i += 1
# print(pareto_all_obj_indices)

# # get subset of all_cells_df with unique objectives,
# # add objs from population mean, and add dose info
# pareto_df_all_cells = unique_df_all_cells.loc[pareto_all_obj_indices]
# pareto_df_all_cells["Mean Objs"] = list(unique_all_obj[pareto_all_obj_indices])
# pareto_df_all_cells = pareto_df_all_cells.reset_index()
# # print(pareto_df_all_cells)
# pareto_df_all_cells["Doses"] = [0]*len(pareto_all_obj_indices)
# for index, row in pareto_df_all_cells.iterrows():
#     # print([[row["Topology"].dose]])
#     pareto_df_all_cells.at[index, "Doses"] = [row["Topology"].dose]
# print(pareto_df_all_cells[["Doses", "Mean Objs"]])
# pareto_df_all_cells.to_pickle(repo_path + "Pulse_seed_pop_DsRED_inhibitor/2023-10-31_Pulse_pop_DsRED_inhibitor_seed_0/pareto_df_all_cells.pkl")

# # plot single cell time series for each topology on pareto
# # front
# # ***choose subset from t_pulse results?***
# for index, row in pareto_df_all_cells.iloc[1:].iterrows():
#     fig, ax = plt.subplots(1, 1, figsize= (4, 4))
#     t = np.arange(0, 43, 1)
#     topology = pareto_df_all_cells.loc[index, "Topology"]
#     # topology.plot_graph()
#     all_cell_ts = pareto_df_all_cells.loc[index, "Rep ON state time series for each cell"][0]
#     mean_cell_ts = [np.mean(k) for k in zip(*all_cell_ts)]
#     mean_cell_ts_rel = pulse.calc_rep_rel(topology, mean_cell_ts)
#     for cell_ts in all_cell_ts:
#         cell_ts_rel = pulse.calc_rep_rel(topology, cell_ts)
#         plt.plot(t, cell_ts_rel)
#     plt.plot(mean_cell_ts_rel, color="k", label="population mean")
#     plt.legend()
#     plt.xlabel("Time (hr)")
#     plt.ylabel("Relative reporter expression")
#     plt.ylim([0, max(mean_cell_ts_rel)])
#     # plt.show()
#     plt.savefig(repo_path+"Pulse_seed_pop_DsRED_inhibitor/2023-10-31_Pulse_pop_DsRED_inhibitor_seed_0/All_cell_time_series_zoomed_"+str(index)+".svg")

# # print(unique_all_obj[pareto_all_obj_indices])

#############################################################
### Signal Conditioner population model ###
#############################################################
# # load files
# repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# file_path_all_cells = "SC_seed_pop_DsRED_inhibitor/Original_dose_terms/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/All_circuits_all_cell_results.pkl"
# file_path_all_obj = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/all_objectives.pkl"
# file_path_final_obj = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/final_objectives_df_with_type.pkl"
# file_path_final_circuits = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/final_population.pkl"

# df_final_obj = pd.read_pickle(repo_path + file_path_final_obj)

# with open(repo_path + file_path_all_obj, "rb") as fid:
#             all_obj = pickle.load(fid)

# df_all_cells = pd.read_pickle(repo_path + file_path_all_cells)
# print(df_all_cells)
# with open(repo_path + file_path_final_circuits, "rb") as fid:
#             final_circuits = pickle.load(fid)

# # get unique values and indices in all_obj
# unique_all_obj, unique_all_obj_idx = np.unique(all_obj, return_index=True, axis=0)
# # print(len(unique_all_obj))
# unique_all_obj = unique_all_obj*-1

# # get rows of df_all_cells that yielded unique_all_obj
# unique_df_all_cells = df_all_cells.loc[unique_all_obj_idx]
# unique_df_all_cells = unique_df_all_cells.reset_index()

# get unique objective values on pareto front
# final_objs_list = df_final_obj[["ON_rel", "FI_rel"]].values.tolist()
# unique_final_objs, unique_final_obj_idx = np.unique(final_objs_list, return_index=True, axis=0)
# unique_final_objs = np.asarray(unique_final_objs*-1)
# unique_final_circuits = final_circuits[unique_final_obj_idx]
# # print(unique_final_objs)

# # get selected values from unique_final_objs
# selected_unique_final_objs = unique_final_objs[np.where(unique_final_objs[:, 1] > 2.0)]
# selected_unique_final_circuits = unique_final_circuits[np.where(unique_final_objs[:, 1] > 2.0)]
# # print(selected_unique_final_circuits)
# for i, obj in enumerate(selected_unique_final_objs):
#     print("edges: ", selected_unique_final_circuits[i][0].edge_list)
#     print("doses: ", selected_unique_final_circuits[i][0].dose)
#     print("obj:", obj)
#     selected_unique_final_circuits[i][0].plot_graph()
# print(selected_unique_final_objs)

# # get indices in unique_all_obj where values = a set
# # of unique values on pareto front
# pareto_all_obj_indices = []
# i = 0
# for row1 in unique_all_obj:
#     for row2 in selected_unique_final_objs:
#         if (row1[0] == row2[0]) & (row1[1] == row2[1]):
#             pareto_all_obj_indices.append(i)
#             break
#     i += 1

# # get subset of all_cells_df with unique objectives,
# # add objs from population mean, and add dose info
# pareto_df_all_cells = unique_df_all_cells.loc[pareto_all_obj_indices]
# pareto_df_all_cells["Mean Objs"] = list(unique_all_obj[pareto_all_obj_indices])
# pareto_df_all_cells = pareto_df_all_cells.reset_index()
# pareto_df_all_cells["Doses"] = [0]*len(pareto_df_all_cells)
# for index, row in pareto_df_all_cells.iterrows():
#     pareto_df_all_cells.at[index, "Doses"] = [row["Topology"].dose]
# # print(pareto_df_all_cells)

# # calculate ON_rel and FI_rel for each cell in 
# # population for each topology
# def calc_ON_rel(rep_on):
#     # print(rep_on)
#     reference_on = Ref_pop20["P1"]['on']
#     # print(reference_on)
#     ON_rel = rep_on/reference_on
#     return ON_rel

# def calc_FI(off, on):
#     FI = on/off
#     return FI

# def calc_FI_rel(FI):
#     FI_ref = Ref_pop20["P1"]['fi']
#     FI_rel = FI/FI_ref
#     return FI_rel

# pareto_df_all_cells["ON rel for each cell"] = [[0]]*len(pareto_df_all_cells)
# pareto_df_all_cells["FI rel for each cell"] = [[0]]*len(pareto_df_all_cells)
# for index, row in pareto_df_all_cells.iterrows():
#     # rep_off and rep_on for each cell in population
#     rep_off_all_cells = row["Rep OFF state for each cell"][0]
#     rep_on_all_cells = row["Rep ON state for each cell"][0]
#     # print("rep off all: ", rep_off_all_cells)
#     # print("rep on all: ", rep_on_all_cells)
    
#     # for each cell, calculate ON_rel
#     ON_rel_all_cells = []
#     FI_rel_all_cells = []
#     for i in range(20):
#         ON_rel = calc_ON_rel(rep_on_all_cells[i])
#         ON_rel_all_cells.append(ON_rel)
#         FI = calc_FI(rep_off_all_cells[i], rep_on_all_cells[i])
#         FI_rel = calc_FI_rel(FI)
#         FI_rel_all_cells.append(FI_rel)
#     # print("FI rel all: ", FI_rel_all_cells)
#     pareto_df_all_cells.at[index, "ON rel for each cell"] = [ON_rel_all_cells]
#     pareto_df_all_cells.at[index, "FI rel for each cell"] = [FI_rel_all_cells]

# print(pareto_df_all_cells["Topology"])

# for i, topology in enumerate(pareto_df_all_cells["Topology"].values.tolist()):
#      topology.plot_graph()
#      print(topology.edge_list)
#      print(topology.dose)
#      print(pareto_df_all_cells["Mean Objs"].loc[i])
     
# pareto_df_all_cells.to_pickle(repo_path + "SC_seed_pop_DsRED_inhibitor/2023-10-31_Signal_Cond_pop_DsRED_inhibitor_seed_0/pareto_df_all_cells.pkl")

# for index, row in pareto_df_all_cells.iterrows():
#     ON_rel_all_cells = pareto_df_all_cells.at[index, "ON rel for each cell"][0] #CHECK all for first set of rows
#     # print("ON_rel: ", ON_rel_all_cells)
#     FI_rel_all_cells = pareto_df_all_cells.at[index, "FI rel for each cell"][0]
#     # print("FI_rel: ", FI_rel_all_cells)
#     metrics_all_cells = list(zip(ON_rel_all_cells, FI_rel_all_cells))
#     # print("metrics: ", metrics_all_cells)
#     # print("mean_metrics: ", pareto_df_all_cells["Mean Objs"][0][0]) #also test [0][0]
#     fig, ax = plt.subplots(1, 1, figsize= (4, 4))
#     for cell_metrics in metrics_all_cells:
#         plt.plot(cell_metrics[0], cell_metrics[1],
#                  linestyle="None", marker = "o")
#     plt.plot(pareto_df_all_cells["Mean Objs"][0][0],
#                 pareto_df_all_cells["Mean Objs"][0][1],
#                 linestyle="None", marker = "o",
#                 color="k", label="population mean")
#     plt.legend()
#     plt.xlabel("ON_rel")
#     plt.ylabel("FI_rel")
#     if index == 0 or index == 4:
#         plt.yticks([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
#     # plt.show()
#     plt.savefig(repo_path + "SC_seed_pop_DsRED_inhibitor/2023-10-31_Signal_Cond_pop_DsRED_inhibitor_seed_0/All_cell_metrics"+str(index)+".svg")




# repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# path_final_obj = "2024-02-05_Pulse_3obj_42h_seed_0/final_objectives_df.pkl"
# path_final_circuits = "2024-02-05_Pulse_3obj_42h_seed_0/final_population.pkl"
# df_obj = pd.read_pickle(repo_path + path_final_obj)
# df_obj_unique = df_obj.drop_duplicates()
# print(df_obj_unique)
# print(df_obj["prominence_rel"].to_list())
# print(df_obj["peak_rel"].to_list())
# for circuit in unique_final_circuits[1:]:
#     t, mean_ts = pulse.simulate(circuit[0])
#     mean_ts_rel = pulse.calc_rep_rel(circuit[0], mean_ts)
#     mean_ts_lists.append(mean_ts_rel)
#     metrics = pulse.func(circuit[0])
#     print("t_pulse = ", metrics[0])
#     print("prominence = ", metrics[1]*-1)
#     plt.plot(t, mean_ts_rel)
# plt.show()

# obj_list = np.array([[1, 2], [3, 4], [5, 6]])
# obj_list = np.array(obj_list)
# obj_list = [1, 2, 3, 4, 5]
# obj_labels = ["ON_rel", "FI_rel"]
# # print(obj_list)
# # topology_dict = {}
# # for i in range(len(obj_labels)):
# #     label = obj_labels[i]
# #     topology_dict[label+"_range"] = max(obj_list[:, i]) - min(obj_list[:, i])
# #     topology_dict[label+"_mean"] = np.mean(obj_list[:, i])
# # print(topology_dict)
# topology_dict_list = []
# for topology in range(3):
#     topology_dict = {}
#     topology_dict["objectives_range"] = max(obj_list) - min(obj_list)
#     topology_dict["objectives_mean"] = np.mean(obj_list)
#     topology_dict["objectives [min, max]"] = [min(obj_list), max(obj_list)]
#     topology_dict["objectives_std_error"] = sem(obj_list, ddof=1)
#     topology_dict["objectives_list"] = obj_list
#     topology_dict_list.append(topology_dict)
#     obj_list = [i*2 for i in obj_list]

# Z_mat_sampling = pd.DataFrame(topology_dict_list)
# print(Z_mat_sampling)



### 28 Feb- load pulse t_pulse and 3obj for all ZFs and ZF1/2 only
# load files- population model
repo_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/2024-02-05_Pulse_pop_DsRED_inhibitor_t_pulse_126h_seed_0/final_objectives_df.pkl"
# file_path_final_circuits = "Pulse_seed_pop_DsRED_inhibitor/2024-02-05_Pulse_pop_DsRED_inhibitor_t_pulse_126h_seed_0/final_population.pkl"
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/2024-02-05_Pulse_pop_DsRED_inhibitor_3obj_126h_seed_0/final_objectives_df.pkl"
# file_path_final_circuits = "Pulse_seed_pop_DsRED_inhibitor/2024-02-05_Pulse_pop_DsRED_inhibitor_3obj_126h_seed_0/final_population.pkl"
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/2024-02-26_Pulse_pop_DsRED_inhibitor_ZF1_ZF2_t_pulse_126h_seed_0/final_objectives_df.pkl"
# file_path_final_circuits = "Pulse_seed_pop_DsRED_inhibitor/2024-02-26_Pulse_pop_DsRED_inhibitor_ZF1_ZF2_t_pulse_126h_seed_0/final_population.pkl"
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/2024-02-26_Pulse_pop_DsRED_inhibitor_ZF1_ZF2_3_obj_126h_seed_0/final_objectives_df.pkl"
# file_path_final_circuits = "Pulse_seed_pop_DsRED_inhibitor/2024-02-26_Pulse_pop_DsRED_inhibitor_ZF1_ZF2_3_obj_126h_seed_0/final_population.pkl"

# pulse = PulseGenerator("P1", [5, 75, 5], 2, True, True, {1: 46, 2: 122}, 2, 0.32, 0.57, mutate_dose=False, pop=True, max_time=126, obj_labels=["t_pulse (hr)", "peak_rel", "prominence_rel"])

# final_obj_unique = pd.read_pickle(repo_path+file_path_final_obj).drop_duplicates()
# print(final_obj_unique)
# unique_obj_idx = final_obj_unique.index.to_list()
# unique_obj_idx.remove(1)
# final_circuits_unique = pd.read_pickle(repo_path+file_path_final_circuits)[unique_obj_idx]
# for circuit in final_circuits_unique:
#     # print(circuit[0].in_dict, circuit[0].dose)

#     doses = deepcopy(circuit[0].dose)
#     edges = deepcopy(circuit[0].edge_list)
#     del doses["Rep"]
#     circuit_exp = Topo(edges, doses, "P1")
#     # print(circuit_exp.in_dict)
#     # print({'Z2': circuit[0].dose["Z2"], 'I1': circuit[0].dose["I1"]})
#     # obj_topo = pulse.func(circuit[0])
#     # objs = pulse.func(circuit_exp)
#     # print(obj_topo)
#     print(circuit_exp.in_dict, circuit_exp.dose)
    # print(objs)

amp_min_obj_path = "Amp_seed_single_cell_vary_dose/2024-03-12_Amplifier_single_cell_vary_dose_new_dose_terms_seed_0/minimum_obj_all_gens.pkl"
amp_obj_path = "Amp_seed_single_cell_vary_dose/2024-03-12_Amplifier_single_cell_vary_dose_new_dose_terms_seed_0/all_objectives.pkl"
amp_min_circuit_path = "Amp_seed_single_cell_vary_dose/2024-03-12_Amplifier_single_cell_vary_dose_new_dose_terms_seed_0/min_obj_circuit_all_gens.pkl"

with open(repo_path+amp_min_circuit_path, "rb") as fid:
    amp_min_circuit = pickle.load(fid)

with open(repo_path+amp_min_obj_path, "rb") as fid:
    amp_min_obj = pickle.load(fid)
# print(first_seen(amp_min_obj))

with open(repo_path+amp_obj_path, "rb") as fid:
    amp_all_obj = pickle.load(fid)
amp_all_obj = np.sort(amp_all_obj.flatten())
amp_all_obj = np.unique(amp_all_obj)
# print(amp_all_obj[:20])

print(amp_min_circuit[-1][0].dose)