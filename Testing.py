# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 14:16:07 2022

@author: Katie_Dreyer
"""
import os
import json
import numpy as np
import networkx as nx
import time
import pickle
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

import pandas as pd
import seaborn as sns
from define_circuit import Topo
from itertools import product, combinations, permutations, chain
from scipy.integrate import odeint
# from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
# from rankcrowding import RankAndCrowding
from plot_search_results import *
from GA import sampling, check_valid
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator
from pymoo.indicators.hv import HV
from plot_search_results import plot_graph, plot_hypervolumes_set
from math import exp
from scipy.interpolate import interp2d
from get_selected_results import get_selected_all_cell_metrics, plot_all_cell_objs
from diversity_metrics import first_seen

plt.style.use('/Users/kdreyer/Documents/Github/GraphGA/paper.mplstyle.py')
sky_blue = [i/255 for i in [86, 180, 233]]

# from load_files_pop import (
#     Z_20,  
#     Ref_pop20, 
#     Ref, 
#     tf_list, 
#     inhibitor_list
# )
# # print(tf_list + inhibitor_list + ['Rep'])
# edge_list = [('P1', 'Z1'), ('P1', 'I6'), ('I6', 'Z2'), ('Z1', 'Z2'), ('Z1', 'Rep')]
# dose_list = {'Z2': 25, 'Z1': 75, 'I6': 75}
# promo_node = 'P1'
# topology = Topo(
#     edge_list,
#     dose_list, 
#     promo_node
# )
# print(topology.in_dict["Z2"]["Z"][0])
# edge_list = [('P1', 'Z1'), ('Z1', 'Z1'), ('Z1', 'Rep')]
# dose_list = {'Z1': 75}
# promo_node = 'P1'
# topology = Topo(
#     edge_list,
#     dose_list, 
#     promo_node
# )

# print(type(topology.part_list[0]))
# print(list(topology.graph.predecessors('Z1')))
# print(list(topology.graph.successors('Z1')))


# print(topology.pool)
# print(topology.in_dict)
# print((topology.in_dict["Z1"]['P'] != []) + (topology.in_dict["Z2"]['P'] != []))

# with open("Amplifier/Amplifier_objective_1.pkl", "rb") as fid:
#     amp_opt1 = pickle.load(fid)

# with open("Amplifier/Amplifier_objective_2.pkl", "rb") as fid:
#     amp_opt2 = pickle.load(fid)

# with open("Amplifier/Amplifier_combo_1.pkl", "rb") as fid:
#     amp_topos1 = pickle.load(fid)

# with open("Amplifier/Amplifier_combo_2.pkl", "rb") as fid:
#     amp_topos2 = pickle.load(fid)

# amp_on_rel1 = [-i/Ref['P1']['on'] for i in amp_opt1]
# sort_idx1 = np.argsort(amp_on_rel1)
# obj_sorted_sc1 = -1*(np.asarray(amp_on_rel1)[sort_idx1])
# topo_sorted_sc1 = np.asarray(amp_topos1)[sort_idx1]
# # print(obj_sorted_s1)
# amp_on_rel2 = [-i/Ref['P1']['on'] for i in amp_opt2]
# sort_idx2 = np.argsort(amp_on_rel2)
# obj_sorted_sc2 = -1*(np.asarray(amp_on_rel2)[sort_idx2])
# topo_sorted_sc2 = np.asarray(amp_topos2)[sort_idx2]

# amp_topos_sorted_all = np.concatenate((topo_sorted_sc1, topo_sorted_sc2))
# # with open("Amplifier/Amplifier_topos_all.pkl", "wb") as fid:
# #     pickle.dump(amp_topos_sorted_all, fid)

# print(amp_topos_sorted_all[0])

# with open("Amplifier/Amplifier_ON_rel_1.pkl", 'rb') as fid:
#     obj_sorted_sc1 = pickle.load(fid)

# with open("Amplifier/Amplifier_ON_rel_2.pkl", 'rb') as fid:
#     obj_sorted_sc2 = pickle.load(fid)

# obj_sorted_sc = np.concatenate((obj_sorted_sc1, obj_sorted_sc2))
# # print(obj_sorted_sc)

# # with open("Amplifier/Amplifier_ON_rel_all.pkl", "wb") as fid:
# #     pickle.dump(obj_sorted_sc, fid)

# # hist, bins, _ = plt.hist(obj_sorted_sc, bins= 30)
# # plt.xlabel("ON_rel")
# # plt.ylabel("Count")
# # plt.show()ß
# # plt.savefig("230322_amp_all_sc_ON_rel_dist.svg")
# # print(bins[-5])

# # idx_high = np.argwhere(obj_sorted_sc >= bins[-5])
# # obj_sc_high = obj_sorted_sc[idx_high]
# # print(obj_sc_high)
# # print(idx_high)


# with open("Results/Combinatorial_results/230329_Amplifier_pop_sampling.pkl", "rb") as fid:
#     pop_sampling = pickle.load(fid)

# num_topos = 432
# obj_range_list = []
# for i in range(num_topos):
#     obj_info = pop_sampling["topologies_list["+str(i)+"]"]
#     obj_range = obj_info["objectives_range"]
#     obj_range_list.append(obj_range)

# with open("Results/Combinatorial_results/230329_Amplifier_pop_sampling2.pkl", "rb") as fid:
#     pop_sampling2 = pickle.load(fid)

# with open("Results/Combinatorial_results/230403_Amplifier_pop_sampling3.pkl", "rb") as fid:
#     pop_sampling3 = pickle.load(fid)

# num_topos = 432
# obj_range_list3 = []
# for i in range(num_topos):
#     obj_info = pop_sampling3["topologies_list["+str(i)+"]"]
#     obj_range = obj_info["objectives_range"]
#     obj_range_list3.append(obj_range)


# # print(obj_range_list)
# plt.hist(obj_range_list3, bins= 30)
# plt.xlabel("ON_rel range")
# plt.ylabel("Count")
# plt.show()

# list1 = [[1], [2], [3]]
# arr1 = np.asarray(list1)
# list2 = np.array([0, 2])
# print(list1, arr1)

# print(arr1[list2])

# list1 = ['Z1', 'Z2', 'Z3']
# list2 = ['I4', 'I5', 'I6']
# list1 = [5, 10, 15]

# list3 = combinations(list1, 1)
# list4 = combinations(list1, 1)
# # print(list(list3))
# prod = list(product(list1, list1))
# print(prod)
# combo_two = [(i[0] + i[1]) for i in list(product(list3, list4))]
# print(combo_two)

# l = [2]*5
# print(l)

# print(np.random.choice(10, 1))

# def add_to_list(l1):
#     l1.append(2)


# l1 = [1, 2, 3]
# add_to_list(l1)
# print(l1)

# with open("init_pop_2.pkl", "rb") as fid:
#     population = pickle.load(fid)

# print(population)

# with open("/Users/kdreyer/Documents/Github/GraphGA/SigCond/SigCond_combo_2.pkl", "rb") as fid:
#     sig_cond = pickle.load(fid)

# # print(sig_cond)
# print(len(sig_cond))

# with open("Amplifier/Amplifier_combo_1.pkl", "rb") as fid:
#     amp1 = pickle.load(fid)

# print(len(amp1))

# with open("Amplifier/Amplifier_combo_2.pkl", "rb") as fid:
#     amp2 = pickle.load(fid)

# print(len(amp2))

# with open("/Users/kdreyer/Desktop/SigCond_combo_65.pkl", "rb") as fid:
#     sig_cond = pickle.load(fid)

# print(len(sig_cond))

# with open('SigCond_combo_pareto.pkl', 'rb') as f:
#     sig_cond_obj = pd.read_pickle(f)

# print(len(sig_cond_obj))
# # print(sig_cond_obj)

# print(np.arange(5, 75+1, 5))

# print(np.random.randint(1))
# path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/2023-10-19_Amplifier_single_cell_test/initial_population.pkl"
# with open(path, "rb") as fid:
#     pop = pickle.load(fid)

# # print(type(pop))

# path = path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/2023-10-19_Amplifier_single_cell_test/top_num_circuit_circuits_each_gen.pkl"
# with open(path, "rb") as fid:
#     top_circuits = pickle.load(fid)

# print(type(top_circuits))

# l1 = [[1, 2], [3, 4], [1, 2]]

# l2 = np.asarray(l1)

# unique_rows, unique_index = np.unique(l2, axis=0, return_index=True)

# print(unique_rows, unique_index)


# arr = np.array([[1], [2], [3], [4], [5]])

# num_vals = 3

# top_num_vals_obj = arr[-num_vals:, 0]
# print(top_num_vals_obj)

# path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/2023-10-20_Amplifier_pop_const_dose_seed_0/all_objectives.pkl"
# with open(path, "rb") as fid:
#     all_obj = pickle.load(fid)



# with open("/Users/kdreyer/Documents/Github/GraphGA/GA_results/2023-10-24_Signal_cond_pop_inhibitor_seed_0/all_objectives.pkl", "rb") as fid:
#     all_obj = pickle.load(fid)

# # print(all_obj)

# with open("/Users/kdreyer/Documents/Github/GraphGA/GA_results/2023-10-24_Signal_cond_pop_inhibitor_seed_0/all_circuits.pkl", "rb") as fid:
#     all_circuits = pickle.load(fid)

# with open("/Users/kdreyer/Documents/Github/GraphGA/GA_results/2023-10-24_Signal_cond_pop_inhibitor_seed_0/final_objectives_df_with_type.pkl", 'rb') as f:
#     obj_df_type = pd.read_pickle(f)
# obj_df_type["ON_rel"] = obj_df_type["ON_rel"]*-1
# obj_df_type["FI_rel"] = obj_df_type["FI_rel"]*-1
# print(obj_df_type)

# plot_pareto_front("", obj_df_type, ["ON_rel", "FI_rel"])
# path1 = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/outdated/Amp_seed_pop_const_dose/2023-10-23_Amplifier_pop_const_dose_seed_0/all_circuits.pkl"
# with open(path1, "rb") as fid:
#     all_circuits = pickle.load(fid)

# # path2 = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/outdated/Amp_seed_pop_const_dose/2023-10-23_Amplifier_pop_const_dose_seed_0/all_objectives.pkl"
# # with open(path2, "rb") as fid:
# #     all_obj = pickle.load(fid)
# circuit_edge_lists = []
# for circuit in all_circuits:
#     circuit_edges = circuit[0].edge_list
#     for key, val in circuit[0].dose.items():
#         circuit_edges.append((key, str(val)))
#     circuit_edge_lists.append(circuit_edges)
#     print("finished circuit " + str(len(circuit_edge_lists)))
# print(circuit_edge_lists[-1])
# combo_edges_lists = []
# for edges in circuit_edge_lists:
#     edge_combos = list([edge[0] + edge[1] for edge in edges])
#     combo_edges_lists.append(edge_combos)
# print(type(combo_edges_lists))

# unique_edge_combo = []
# index_list = []
# seen = set()
# for i, combo_list in enumerate(combo_edges_lists):
#     combo_set = frozenset(combo_list)
#     if combo_set not in seen:
#         seen.add(combo_set)
#         index_list.append(i)
#         unique_edge_combo.append(combo_list)
#     print("finished circuit " + str(i))
# # print(index_list, unique_edge_combo)
# print(len(index_list))
# # print(seen)
# unique_obj = all_obj[index_list]
# unique_obj_pos = unique_obj*-1
# print(len(np.argwhere(unique_obj_pos >= 67.5)))

 


#################################################################
########### Amplifier experimental circuit new predictions ######
#################################################################
# repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# amp_results_path = "Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/all_circuits.pkl"
# amp_unique_obj_path = "Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/unique_objectives.pkl"
# amp_unique_circuits_path = "Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/unique_circuits.pkl"
# amp_obj_path = "Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/all_objectives.pkl"
# amp_results = pd.read_pickle(repo_path+amp_results_path)
# amp_obj = pd.read_pickle(repo_path+amp_obj_path)*-1
# amp_obj = amp_obj.flatten()
# print(max(amp_obj))
# print(amp_obj)
# amp = Amplifier("P1", [5, 75, 5], 2, True, False, {1: 46, 2: 122}, 2, 0.32, 0.57, False, True)

# amp_unique_obj = pd.read_pickle(repo_path+amp_unique_obj_path).flatten()*-1
# amp_unique_circuits = pd.read_pickle(repo_path+amp_unique_circuits_path)
# amp_unique_obj_high = np.where(amp_unique_obj > 69.5)[0]
# unique_circuits_high = amp_unique_circuits[amp_unique_obj_high]
# unique_obj_high = amp_unique_obj[amp_unique_obj_high]
# print(unique_obj_high[31])
# print(unique_circuits_high)


# print(len(amp_unique_obj_high[0]))
### Circuit 1 from experiments
# print(np.where(np.logical_and(amp_obj>69.982, amp_obj<69.98215)))
# print(amp_obj[[1874, 1940, 1951, 1965, 2027, 2074, 2108, 2112, 2165, 2175, 2226,
#        2236, 2253, 2260, 2278, 2300, 2311, 2324, 2331, 2342, 2497, 2516,
#        2573]])
# amp_results[1874][0].plot_graph()
# circuit_1 = amp_results[1874][0]
# print(circuit_1.dose, circuit_1.in_dict)
# circuit_1_exp = Topo([('P1', 'Z9'), ('Z9', 'Z9'), ('Z9', 'Z6'), ('Z6', 'Z6'), ('Z6', 'Rep'), ('Z6', 'Z9')], {'Z6': 75, 'Z9': 75}, "P1")
# obj = amp.func(circuit_1_exp)
# print(obj)
# edge_list = circuit_1_exp.edge_list
# flattened_edges = list(chain.from_iterable(circuit_1_exp.edge_list))
# if "Z6" in flattened_edges:
#     print(edge_list)
# ON_rel_neg = amp.func(circuit_1_exp)
# print(ON_rel_neg)

### Circuit 2 from experiments
# print(np.where(np.logical_and(amp_obj>69.97675, amp_obj<69.9768)))
# print(amp_obj[[1692, 1712, 1889, 2009, 2038, 2350, 2543, 2621, 2624, 2647]])
# amp_results[1692][0].plot_graph()
# circuit_2 = amp_results[1692][0]
# print(circuit_2.dose, circuit_2.in_dict)
# print(circuit_2.edge_list)
# circuit_2_exp = Topo([('P1', 'Z2'), ('Z2', 'Z2'), ('Z2', 'Z6'), ('Z6', 'Z6'), ('Z6', 'Rep')], {'Z6': 75, 'Z2': 75}, "P1")
# ON_rel_neg = amp.func(circuit_2_exp)
# print(ON_rel_neg)

### Circuit 3 from experiments
# print(np.where(np.logical_and(amp_obj>39.8308, amp_obj<39.83095)))
# print(amp_obj[1251])
# amp_results[1251][0].plot_graph()
# circuit_3 = amp_results[1251][0]
# print(circuit_3.dose, circuit_3.in_dict)
# print(circuit_3.edge_list)
# circuit_3_exp = Topo([('P1', 'Z2'), ('Z2', 'Z2'), ('Z2', 'Z9'), ('Z9', 'Z9'), ('Z9', 'Rep')], {'Z2': 75, 'Z9': 75}, "P1")
# ON_rel_neg = amp.func(circuit_3_exp)
# print(ON_rel_neg)

### Circuit 4 from experiments
# print(np.where(np.logical_and(amp_obj>39.8018, amp_obj<39.802)))
# print("GCAD design scenario: ", amp_obj[22])
# print(amp_results[22][0].plot_graph())
# circuit_4 = amp_results[22][0]
# print(circuit_4.dose, circuit_4.in_dict)
# circuit_4_exp = Topo([('P_exp_amp', 'Z9'), ('Z9', 'Rep'), ('Z9', 'Z9')], {'Z9': 75}, "P_exp_amp")
# print(circuit_4.dose)
# ON_rel_neg = amp.func(circuit_4_exp)
# print(ON_rel_neg)


#############################################################################
########### Signal Conditioner experimental circuit new predictions #########
#############################################################################
# repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# sc_results_path = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/final_population.pkl"
# sc_obj_path = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/final_objectives_df.pkl"
# sc_all_obj_path = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/all_objectives.pkl"
# sc_all_circuits_path = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/all_circuits.pkl"
# sc_all_obj = pd.read_pickle(repo_path+sc_all_obj_path)*-1
# sc_pareto_obj = pd.read_pickle(repo_path+sc_obj_path)
# print(sc_all_obj)
# print(sc_pareto_obj)
# unique_sc_pareto_obj = (
#     sc_pareto_obj.drop_duplicates()

# print(unique_sc_pareto_obj.index.tolist())
# sc_pareto_circuits = pd.read_pickle(repo_path+sc_results_path)
# unique_sc_pareto_circuits = sc_pareto_circuits[unique_sc_pareto_obj.index]
# print(unique_sc_pareto_circuits, len(unique_sc_pareto_circuits))
# sc_all_circuits = pd.read_pickle(repo_path+sc_all_circuits_path)
# unique_all_obj, unique_all_obj_idx = np.unique(sc_all_obj, return_index=True, axis=0)
# unique_all_circuits = sc_all_circuits[unique_all_obj_idx]
# high_obj = unique_all_obj[np.where(unique_all_obj[:, 1] > 1.0)] 
# print(high_obj[49])
# high_circuits = unique_all_circuits[np.where(unique_all_obj[:, 1] > 1.0)]
# for i, circuit in enumerate(high_circuits):
#     # print(circuit[0].part_list)
#     if ("Z1" in circuit[0].part_list) and ("I1" in circuit[0].part_list):
#         print(i, circuit[0].edge_list)
#         print(circuit[0].dose)
#     elif ("Z1" in circuit[0].part_list) and ("I2" in circuit[0].part_list):
#         print(i, circuit[0].edge_list)
#         print(circuit[0].dose)
#     elif ("Z2" in circuit[0].part_list) and ("I2" in circuit[0].part_list):
#         print(i, circuit[0].edge_list)
#         print(circuit[0].dose)
#     elif ("Z2" in circuit[0].part_list) and ("I1" in circuit[0].part_list):
#         print(i, circuit[0].edge_list)
#         print(circuit[0].dose)
#     if circuit[0].part_list[0][1:] == circuit[0].part_list[1][1:]:
#         print(i, circuit[0].edge_list)
#         print(circuit[0].dose)

# sc_results = pd.read_pickle(repo_path+sc_results_path)
# print(len(sc_results))
# sc_obj = pd.read_pickle(repo_path+sc_obj_path)[["ON_rel", "FI_rel"]]*-1
# print(sc_obj[sc_obj["FI_rel"] > 2].drop_duplicates().sort_values(by="FI_rel", ascending=False))
# sc = SignalConditioner('P1', [5, 75, 5], 2, True, True, {1: 46, 2: 122}, 2, 0.32, 0.57, False, True)


# not_equal = []
# for i, circuit in enumerate(sc_results[:1000]):
#     in_dict_list = list(circuit[0].in_dict.keys())
#     in_dict_list.remove("Rep") 
#     dose_dict_list = list(circuit[0].dose.keys())
#     dose_dict_list.remove("Rep")
#     # print(in_dict_list)
#     if in_dict_list != dose_dict_list:
#         not_equal.append(circuit)
#         print(i)
# print(len(not_equal))
        # print("not equal")
### Circuit 1 from experiments
# sc_results[1][0].plot_graph()
# circuit_1 = sc_results[1][0]

# print(circuit_1.in_dict, circuit_1.dose)
# circuit_1_exp = Topo([('P1', 'Z12'), ('Z12', 'Rep'), ('Z12', 'I1'), ('P1', 'I1'), ('I1', 'Rep')], {'Z12': 60, 'I1': 60}, 'P1')
# objs = sc.func(circuit_1_exp)
# print("c1: ",objs)
# print(circuit_1_exp.part_list)
# if ("Z12" in circuit_1_exp.part_list) and ("I1" in circuit_1_exp.part_list):
#     print(circuit_1_exp.edge_list)
# objs_neg, FI_sc = sc.func(circuit_1_exp)
# print(objs_neg, FI_sc)

### Circuit 2 from experiments
# sc_results[7][0].plot_graph()
# circuit_2 = sc_results[7][0]
# print(circuit_2.dose, circuit_2.in_dict)
# circuit_2_exp = Topo([('P1', 'Z12'), ('P1', 'I11'), ('Z12', 'Rep'), ('Z12', 'I11'), ('I11', 'Rep')], {'Z12': 55, 'I11': 75}, "P1")
# objs_neg = sc.func(circuit_2_exp)
# print("c2: ", objs_neg)

### Circuit 3 from experiments
# sc_results[12][0].plot_graph()
# circuit_3 = sc_results[12][0]
# print(circuit_3.dose, circuit_3.in_dict)
# circuit_3_exp = Topo([('P1', 'Z9'), ('P1', 'I11'), ('Z9', 'Rep'), ('Z9', 'I11'), ('I11', 'Rep')], {'I11': 35, 'Z9': 75}, 'P1')
# objs_neg = sc.func(circuit_3_exp)
# print("c3: ", objs_neg)

### Circuit 4 from experiments
# sc_results[24][0].plot_graph()
# circuit_4 = sc_results[24][0]
# print(circuit_4.dose, circuit_4.in_dict)
# print(circuit_4.edge_list)
# circuit_4_exp = Topo([('P1', 'Z12'), ('Z12', 'Rep'), ('Z12', 'I1'), ('I1', 'Rep')], {'Z12': 55, 'I1': 5}, "P1")
# objs_neg = sc.func(circuit_4_exp)
# print("c4: ", objs_neg)

### Circuit 5 from experiments
# sc_results[21][0].plot_graph()
# circuit_5 = sc_results[21][0]
# print(circuit_5.dose, circuit_5.in_dict)
# print(circuit_5.edge_list, circuit_5.dose)
# circuit_5_exp = Topo([('P1', 'Z12'), ('Z12', 'Rep'), ('Z12', 'I15'), ('I15', 'Rep')], {'Z12': 55, 'I15': 20}, "P1")
# objs_neg = sc.func(circuit_5_exp)
# print("c5: ", objs_neg)

### Circuit 6 from experiments
# sc_results[11][0].plot_graph()
# circuit_6 = sc_results[11][0]
# print(circuit_6.dose, circuit_6.in_dict)
# print(circuit_6.edge_list, circuit_6.dose)
# circuit_6_exp = Topo([('P1', 'Z12'), ('Z12', 'Rep'), ('Z12', 'I11'), ('I11', 'Rep')], {'I11': 5, "Z12": 50}, "P1")
# objs_neg = sc.func(circuit_6_exp)
# print("c6: ", objs_neg)


### Original GA run for signal conditioner ZF1 and ZF2 only
# repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# sc_results_path = "SC_seed_pop_DsRED_inhibitor/2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/final_population.pkl"
# sc_obj_path = "SC_seed_pop_DsRED_inhibitor/2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/final_objectives_df.pkl"
# sc_all_obj_path = "SC_seed_pop_DsRED_inhibitor/2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/all_objectives.pkl"
# sc_all_circuits_path = "SC_seed_pop_DsRED_inhibitor/2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/all_circuits.pkl"
# sc_all_cell_results_path = "SC_seed_pop_DsRED_inhibitor/Original_dose_terms/2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/All_circuits_all_cell_results.pkl"
# sc_all_obj = pd.read_pickle(repo_path+sc_all_obj_path)*-1
# # print(sc_all_obj)
# sc_pareto_obj = pd.read_pickle(repo_path+sc_obj_path)
# sc_pareto_circuits = pd.read_pickle(repo_path+sc_results_path)
# sc_pareto_obj = sc_pareto_obj.drop_duplicates()
# # print(sc_pareto_obj)
# sc_pareto_obj["ON_rel"] = sc_pareto_obj["ON_rel"]*-1
# sc_pareto_obj["FI_rel"] = sc_pareto_obj["FI_rel"]*-1
# sc_pareto_obj_high = sc_pareto_obj[sc_pareto_obj["FI_rel"] > 1.3]
# print(sc_pareto_obj_high)
# for i, circuit in enumerate(sc_pareto_circuits[sc_pareto_obj_high.index.tolist()]):
#     print(circuit[0].edge_list, circuit[0].dose)
#     plot_graph(repo_path +"2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/topology_" + str(i) + ".svg", circuit[0])
# sc_all_circuits = pd.read_pickle(repo_path+sc_all_circuits_path)
# unique_all_obj, unique_all_obj_idx = np.unique(sc_all_obj, return_index=True, axis=0)
# unique_all_circuits = sc_all_circuits[unique_all_obj_idx]
# high_obj = unique_all_obj[np.where(unique_all_obj[:, 1] > 1.9)]
# # print(high_obj[np.argsort(high_obj[:, -1])])
# high_circuits = unique_all_circuits[np.where(unique_all_obj[:, 1] > 1.9)]
# circuit_idx = []
# for i, circuit in enumerate(high_circuits):
#     # if len(circuit[0].edge_list) == 4:
#     #     circuit_idx.append(i)
#     print(circuit[0].edge_list)
#     print(circuit[0].dose)
# print(high_obj)
# all_cells = pd.read_pickle(repo_path+sc_all_cell_results_path)
# print(all_cells)

####################################################################
########## NEW GA RUNS WITH UPDATED DOSE TERMS ######################
####################################################################
### SC all ZFs, 80gen
# repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# sc_results_path = "SC_seed_pop_DsRED_inhibitor/2024-03-06_Signal_Cond_pop_DsRED_inhibitor_ngen80_new_dose_terms_seed_0/final_population.pkl"
# sc_obj_path = "SC_seed_pop_DsRED_inhibitor/2024-03-06_Signal_Cond_pop_DsRED_inhibitor_ngen80_new_dose_terms_seed_0/final_objectives_df.pkl"
# sc_all_obj_path = "SC_seed_pop_DsRED_inhibitor/2024-03-06_Signal_Cond_pop_DsRED_inhibitor_ngen80_new_dose_terms_seed_0/all_objectives.pkl"
# sc_all_circuits_path = "SC_seed_pop_DsRED_inhibitor/2024-03-06_Signal_Cond_pop_DsRED_inhibitor_ngen80_new_dose_terms_seed_0/all_circuits.pkl"
# sc_all_cells_path = "SC_seed_pop_DsRED_inhibitor/2024-03-06_Signal_Cond_pop_DsRED_inhibitor_ngen80_new_dose_terms_seed_0/All_circuits_all_cell_results.pkl"
# # sc_all_obj = pd.read_pickle(repo_path+sc_all_obj_path)*-1
# # print(sc_all_obj)
# sc_pareto_obj = pd.read_pickle(repo_path+sc_obj_path)
# sc_pareto_circuits = pd.read_pickle(repo_path+sc_results_path)
# sc_pareto_obj = sc_pareto_obj.drop_duplicates()
# print(sc_pareto_obj)
# sc_pareto_obj["ON_rel"] = sc_pareto_obj["ON_rel"]*-1
# sc_pareto_obj["FI_rel"] = sc_pareto_obj["FI_rel"]*-1
# sc_pareto_obj_high = sc_pareto_obj[sc_pareto_obj["FI_rel"] >= 1.0]
# print(sc_pareto_obj_high, len(sc_pareto_obj_high))
# for i, circuit in enumerate(sc_pareto_circuits[sc_pareto_obj_high.index.tolist()]):
#     print(circuit[0].edge_list, circuit[0].dose)
#     print(((float(circuit[0].dose["I13"]) / circuit[0].pool["I13"])/200)**0.5)
# print((45/200)**0.5)
#     circuit[0].plot_graph()
#     plot_graph(repo_path +"2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/topology_" + str(i) + ".svg", circuit[0])

# all_cell_results = pd.read_pickle(repo_path+sc_all_cells_path)
# print(all_cell_results)

### SC ZF1 and ZF2 only, 60gen
# repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# sc_results_path = "SC_seed_pop_DsRED_inhibitor/2024-03-05_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_new_ZF_reg_dose_term_seed_0/final_population.pkl"
# sc_obj_path = "SC_seed_pop_DsRED_inhibitor/2024-03-05_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_new_ZF_reg_dose_term_seed_0/final_objectives_df.pkl"
# sc_all_obj_path = "SC_seed_pop_DsRED_inhibitor/2024-03-05_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_new_ZF_reg_dose_term_seed_0/all_objectives.pkl"
# sc_all_circuits_path = "SC_seed_pop_DsRED_inhibitor/2024-03-05_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_new_ZF_reg_dose_term_seed_0/all_circuits.pkl"
# # sc_all_obj = pd.read_pickle(repo_path+sc_all_obj_path)*-1
# # print(sc_all_obj)
# sc_pareto_obj = pd.read_pickle(repo_path+sc_obj_path)
# sc_pareto_circuits = pd.read_pickle(repo_path+sc_results_path)
# sc_pareto_obj = sc_pareto_obj.drop_duplicates()
# # # print(sc_pareto_obj)
# sc_pareto_obj["ON_rel"] = sc_pareto_obj["ON_rel"]*-1
# sc_pareto_obj["FI_rel"] = sc_pareto_obj["FI_rel"]*-1
# sc_pareto_obj_high = sc_pareto_obj[sc_pareto_obj["FI_rel"] >= 2.0]
# # print(sc_pareto_obj_high, len(sc_pareto_obj_high))
# for i, circuit in enumerate(sc_pareto_circuits[sc_pareto_obj_high.index.tolist()]):
#     print(circuit[0].edge_list, circuit[0].dose)
#     circuit[0].plot_graph()
#     plot_graph(repo_path +"2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/topology_" + str(i) + ".svg", circuit[0])

# ### Pulse all ZFs, t_pulse 126h (80 gen)
# repo_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ngen80_new_dose_terms_seed_0/final_objectives_df.pkl"
# file_path_final_circuits = "Pulse_seed_pop_DsRED_inhibitor/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ngen80_new_dose_terms_seed_0/final_population.pkl"

# final_obj_unique = pd.read_pickle(repo_path+file_path_final_obj).drop_duplicates()
# final_obj_unique["prominence_rel"] = final_obj_unique["prominence_rel"]*-1
# final_obj_unique = final_obj_unique[final_obj_unique["prominence_rel"] > 1.0]
# final_obj_unique_reset_idx = final_obj_unique.reset_index()
# # print(final_obj_unique_reset_idx)
# unique_obj_idx = final_obj_unique.index.to_list()
# # unique_obj_idx.remove(0)
# # print(unique_obj_idx)
# final_circuits_unique = pd.read_pickle(repo_path+file_path_final_circuits)[unique_obj_idx]
# # pulse = PulseGenerator("P1", [5, 75, 5], 2, True, True, {1: 46, 2: 122}, 2, 0.32, 0.57, mutate_dose=False, pop=True, max_time=126, obj_labels=["t_pulse (hr)", "prominence_rel"])
# for i, circuit in enumerate(final_circuits_unique):
# #     print(circuit[0].edge_list)
#     print(circuit[0].dose)
#     print([final_obj_unique_reset_idx.loc[i, "t_pulse"], final_obj_unique_reset_idx.loc[i, "prominence_rel"]])
# #         circuit[0].dose.pop("Rep")
# #         circuit_topo = Topo(circuit[0].edge_list, circuit[0].dose, "P1")
# #         objs = pulse.func(circuit_topo)
# #         print(objs)
#     circuit[0].plot_graph()


# ### Pulse ZF1 and ZF2 only, t_pulse 126h (50 gen)
# repo_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_new_dose_terms_seed_0/final_objectives_df.pkl"
# file_path_final_circuits = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_new_dose_terms_seed_0/final_population.pkl"
# file_path_unique_obj = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_new_dose_terms_seed_0/unique_objectives_df.pkl"
# file_path_unique_circuits = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_new_dose_terms_seed_0/unique_circuits.pkl"

# final_obj_unique = pd.read_pickle(repo_path+file_path_final_obj).drop_duplicates()
# final_obj_unique["prominence_rel"] = final_obj_unique["prominence_rel"]*-1
# final_obj_unique = final_obj_unique[final_obj_unique["prominence_rel"] > 0]

# unique_obj = pd.read_pickle(repo_path+file_path_unique_obj)
# unique_circuits = pd.read_pickle(repo_path+file_path_unique_circuits)
# unique_obj_high_prom = unique_obj[unique_obj["prominence_rel"] < -0.2]
# unique_obj_high_prom_abs = unique_obj_high_prom.copy()
# unique_obj_high_prom_abs["prominence_rel"] = abs(unique_obj_high_prom["prominence_rel"])
# for index, row in unique_obj_high_prom_abs.iterrows():
#     if row["prominence_rel"] in final_obj_unique["prominence_rel"].tolist():
#         unique_obj_high_prom_abs.drop(index, inplace=True)

# unique_circuits_high_prom = unique_circuits[unique_obj_high_prom_abs.index.tolist()]
# unique_obj_high_prom_abs.reset_index(inplace=True)
# sorted_selected_pareto_obj = unique_obj_high_prom_abs.sort_values("t_pulse", ascending=False)
# sorted_obj_idx = sorted_selected_pareto_obj.index.tolist()
# sorted_selected_pareto_circuits = unique_circuits_high_prom[sorted_obj_idx].flatten()
# sorted_selected_pareto_obj.reset_index(inplace=True)
# # print(sorted_selected_pareto_obj)
# edge_lists = []
# dose_dict_list = []
# # objs_list_p1 = []
# # all_cells_dict_list_p1 = []
# pulse_p1 = PulseGenerator("P1", [5, 75, 5], 2, True, True, {1: 46, 2: 122}, 2, 0.32, 0.57, mutate_dose=False, pop=True, max_time=126, obj_labels=["t_pulse (hr)", "prominence_rel"], single_cell_tracking=False)

# # pulse_idx_p0 = []
# # objs_list_p0 = []
# # all_cells_dict_list_p0 = []
# # pulse_p0 = PulseGenerator("P0", [5, 75, 5], 2, True, True, {1: 46, 2: 122}, 2, 0.32, 0.57, mutate_dose=False, pop=True, max_time=126, obj_labels=["t_pulse (hr)", "prominence_rel"], single_cell_tracking=True)

# # pulse_idx_p2 = []
# # objs_list_p2 = []
# # all_cells_dict_list_p2 = []
# # pulse_p2 = PulseGenerator("P2", [5, 75, 5], 2, True, True, {1: 46, 2: 122}, 2, 0.32, 0.57, mutate_dose=False, pop=True, max_time=126, obj_labels=["t_pulse (hr)", "prominence_rel"], single_cell_tracking=True)

# for i, circuit in enumerate(sorted_selected_pareto_circuits):
#       edge_lists.append(circuit.edge_list)
#       dose_dict_list.append(circuit.dose)
#       keys_order = list(circuit.in_dict.keys())
#       keys_order.remove("Rep")
#       circuit.dose.pop("Rep")
#       doses_ordered = {key: circuit.dose[key] for key in keys_order}

#       circuit_topo_p1 = Topo(circuit.edge_list, doses_ordered, "P1")
#       [objs_p1, all_cells_dict_p1] = pulse_p1.func(circuit_topo_p1)
# #       objs_list.append(objs_p1)
# #       all_cells_dict_list.append(all_cells_dict_p1)
#       circuit_topo_p0 = Topo(circuit.edge_list, doses_ordered, "P0")
#       [objs_p0, all_cells_dict_p0] = pulse_p0.func(circuit_topo_p0)
#       objs_list_p0.append(objs_p0)
#       all_cells_dict_list_p0.append(all_cells_dict_p0)
#       if objs_p0[0] != 0.0:
#             pulse_idx_p0.append(i)

#       circuit_topo_p2 = Topo(circuit.edge_list, doses_ordered, "P2")
#       [objs_p2, all_cells_dict_p2] = pulse_p2.func(circuit_topo_p2)
#       objs_list_p2.append(objs_p2)
#       all_cells_dict_list_p2.append(all_cells_dict_p2)
#       if objs_p2[0] != 0.0:
#             pulse_idx_p2.append(i)
      # print(objs)
      # print(circuit.edge_list)
# print(np.abs(objs_list_p0))
# print(np.abs(objs_list_p2))

# selected_results_dict = {
#       "Topology": sorted_selected_pareto_circuits,
#       "Edge list": edge_lists,
#       "Doses": dose_dict_list
# }
# for obj in ["t_pulse", "prominence_rel"]:
#         selected_results_dict[obj] = sorted_selected_pareto_obj[obj].tolist()
# selected_results_df = pd.DataFrame.from_dict(selected_results_dict)
# # # selected_results_df.to_csv("/Users/kdreyer/Documents/Github/GraphGA/GA_results/Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_new_dose_terms_seed_0/2024-03-19_results_analysis_sub_opt/" + "selected_results_sub_opt.csv")# print(unique_circuits_high_prom[-1])
# settings_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_new_dose_terms_seed_0/settings.json"
# with open(settings_path, "rb") as fid:
#       settings = json.load(fid)
# if __name__ == "__main__":
#         all_cell_results_df = get_selected_all_cell_metrics(settings, selected_results_df)
#         # all_cell_results_tp20 = all_cell_results_df.copy()[all_cell_results_df["t_pulse_mean"] >= 20]
#         all_cell_results_df["single_cell_peaks"] = 0
#         all_cell_results_df["single_cell_peaks"] = all_cell_results_df["single_cell_peaks"].astype(object)
#         all_cell_results_df["single_cell_prominence"] = 0
#         all_cell_results_df["single_cell_prominence"] = all_cell_results_df["single_cell_prominence"].astype(object)
#         for index, row in all_cell_results_df.iterrows():
#                 peak_cell_list = []
#                 prom_cell_list = []
#                 # print("t_pulse = ", row["t_pulse_mean"], " h topology: ")
#                 for i in range(20):
#                 # print(row["t_pulse_mean"])
#                         peak_cell = max(row["Rep_rel time series for each cell"][i])
#                         peak_cell_list.append(peak_cell)
#                         prom_cell =  pulse_p1.calc_prominence_rel(row["Rep_rel time series for each cell"][i], peak_cell)
#                         prom_cell_list.append(prom_cell)
#                         # print(prom_cell)
#                 all_cell_results_df.at[index, "single_cell_peaks"] = peak_cell_list
#                 all_cell_results_df.at[index, "single_cell_prominence"] = prom_cell_list
#         # print(all_cell_results_df.columns)
#         all_cell_results_df_save = all_cell_results_df.copy().drop(['Rep_rel time series for each cell', 'Rep_rel time series mean'], axis=1)
#         # print(all_cell_results_df_save)
#         # all_cell_results_df_save.to_csv("/Users/kdreyer/Documents/Github/GraphGA/GA_results/Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_new_dose_terms_seed_0/2024-03-19_results_analysis_sub_opt/" + "all_cell_metrics_sub_opt.csv")
#         # all_cell_results_df = pd.read_csv("/Users/kdreyer/Documents/Github/GraphGA/GA_results/Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_new_dose_terms_seed_0/2024-03-19_results_analysis_sub_opt/all_cellselected_results_sub_opt.csv")
#         # print(all_cell_results_df.iloc[0:3])
#         # print(type(all_cell_results_df["Rep_rel time series for each cell"][0]))
#         # df_plot = all_cell_results_df.copy()[all_cell_results_df["t_pulse_mean"] >= 20]
#         base_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_new_dose_terms_seed_0/2024-03-19_results_analysis_sub_opt/"
#         plot_all_cell_objs(base_path, settings, all_cell_results_df)
# for circuit in unique_circuits_high_prom:
#     print(circuit[0].edge_list)
#     print(circuit[0].dose)
#     circuit[0].plot_graph()

# print(final_obj_unique)
# final_obj_unique_reset_idx = final_obj_unique.reset_index()
# # print(final_obj_unique_reset_idx)
# unique_obj_idx = final_obj_unique.index.to_list()
# # unique_obj_idx.remove(0)
# # print(unique_obj_idx)
# final_circuits_unique = pd.read_pickle(repo_path+file_path_final_circuits)[unique_obj_idx]
# print(final_circuits_unique)
# # pulse = PulseGenerator("P1", [5, 75, 5], 2, True, True, {1: 46, 2: 122}, 2, 0.32, 0.57, mutate_dose=False, pop=True, max_time=126, obj_labels=["t_pulse (hr)", "prominence_rel"])

# for i, circuit in enumerate(final_circuits_unique):
#     print(circuit[0].edge_list)
#     print(circuit[0].dose)
#     print([final_obj_unique_reset_idx.loc[i, "t_pulse"], final_obj_unique_reset_idx.loc[i, "prominence_rel"]])
# #         circuit[0].dose.pop("Rep")
# #         circuit_topo = Topo(circuit[0].edge_list, circuit[0].dose, "P1")
# #         objs = pulse.func(circuit_topo)
# #         print(objs)
#     circuit[0].plot_graph()


### Pulse ZF1 and ZF2 only, 3 obj 126h (50 gen)
# repo_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# file_path_final_obj = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_3obj_126h_ZF1_ZF2_new_dose_terms_seed_0/final_objectives_df.pkl"
# file_path_final_circuits = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/2024-03-07_Pulse_pop_DsRED_inhibitor_3obj_126h_ZF1_ZF2_new_dose_terms_seed_0/final_population.pkl"

# final_obj_unique = pd.read_pickle(repo_path+file_path_final_obj).drop_duplicates()
# final_obj_unique["prominence_rel"] = final_obj_unique["prominence_rel"]*-1
# final_obj_unique = final_obj_unique[final_obj_unique["prominence_rel"] > 0]
# final_obj_unique_reset_idx = final_obj_unique.reset_index()
# # print(final_obj_unique_reset_idx)
# unique_obj_idx = final_obj_unique.index.to_list()
# # unique_obj_idx.remove(0)
# # print(unique_obj_idx)
# final_circuits_unique = pd.read_pickle(repo_path+file_path_final_circuits)[unique_obj_idx]
# # pulse = PulseGenerator("P1", [5, 75, 5], 2, True, True, {1: 46, 2: 122}, 2, 0.32, 0.57, mutate_dose=False, pop=True, max_time=126, obj_labels=["t_pulse (hr)", "prominence_rel"])
# for i, circuit in enumerate(final_circuits_unique):
# #     print(circuit[0].edge_list)
# #     print(circuit[0].dose)
#     print([final_obj_unique_reset_idx.loc[i, "t_pulse"], final_obj_unique_reset_idx.loc[i, "prominence_rel"]])
# #         circuit[0].dose.pop("Rep")
# #         circuit_topo = Topo(circuit[0].edge_list, circuit[0].dose, "P1")
# #         objs = pulse.func(circuit_topo)
# #         print(objs)
#     circuit[0].plot_graph()

# from get_system_equations_pop import system_equations_pop

# topology = Topo([('P1', 'Z2'), ('Z2', 'Z2'), ('Z2', 'Z6'), ('Z6', 'Z6'), ('Z6', 'Rep')], {'Z6': 75, 'Z2': 75}, "P1")


# def simulate(topology, max_time=42):
#     t = np.arange(0, max_time + 1, 1)
#     rep_on = odeint(system_equations_pop, np.zeros(topology.num_states * 2), t, args=('on', np.ones(5), topology,))[:, -1]
#     return rep_on

# rep_on = simulate(topology)
# print(rep_on)

# from statistics import mean

# vals = [10, 20]
# mean1 = np.mean(vals)
# mean2 = (vals[0] + vals[1])/2
# print(mean1, mean2)

# with open("libE_history_at_abort_0.npy", 'rb') as fid:
#     libe = np.load(fid)

# print(libe)

# with open("libE_persis_info_at_abort_0.pickle", 'rb') as fid:
#     libe = pickle.load(fid)

# print(libe)



#########################################################
### Loading results for comparison to opt hyperparams ###
#########################################################

### AMPLIFIER single cell###
# path_amp_const = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_single_cell_const_dose/2024-03-12_Amplifier_single_cell_new_dose_terms_seed_0/minimum_obj_all_gens.pkl"
# with open(path_amp_const, "rb") as fid:
#     min_objs_const = pickle.load(fid)
# print(min_objs_const, len(min_objs_const))

# path_amp_vary = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_single_cell_const_dose/2024-03-12_Amplifier_single_cell_new_dose_terms_seed_0/minimum_obj_all_gens.pkl"
# with open(path_amp_vary, "rb") as fid:
#     min_objs_vary = pickle.load(fid)
# print(min_objs_vary)

### AMPLIFIER population###

### SIGNAL CONDITIONER single cell###
# path_sc = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_single_cell_inhibitor/2024-03-12_Signal_Cond_single_cell_inhibitor_new_dose_terms_seed_0/hypervolumes.pkl"
# with open(path_sc, "rb") as fid:
#     hvs1 = pickle.load(fid)
# print("SC: ", first_seen(hvs1))

# path_sc_dsr = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_single_cell_DsRED_inhibitor/2024-03-11_Signal_Cond_single_cell_DsRED_inhibitor_new_dose_terms_seed_0/hypervolumes.pkl"
# with open(path_sc_dsr, "rb") as fid:
#     hvs2 = pickle.load(fid)
# print("SC DsR: ", first_seen(hvs2))

### SIGNAL CONDITIONER population###
# path_sc_pop = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_pop_inhibitor/2024-03-07_Signal_Cond_pop_inhibitor_new_dose_terms_seed_0/hypervolumes.pkl"
# with open(path_sc_pop, "rb") as fid:
#     hvs1 = pickle.load(fid)
# print("SC: ", hvs1)

# path_sc_pop = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_pop_inhibitor/2024-03-07_Signal_Cond_pop_inhibitor_ngen80_new_dose_terms_seed_0/hypervolumes.pkl"
# with open(path_sc_pop, "rb") as fid:
#     hvs1 = pickle.load(fid)
# print("SC: ", first_seen(hvs1))

# path_sc_dsr = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_single_cell_DsRED_inhibitor/2024-03-11_Signal_Cond_single_cell_DsRED_inhibitor_new_dose_terms_seed_0/hypervolumes.pkl"
# with open(path_sc_dsr, "rb") as fid:
#     hvs2 = pickle.load(fid)
# print("SC DsR: ", first_seen(hvs2))


## Pulse single cell###
# path3 = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Pulse_seed_single_cell_DsRED_inhibitor/2024-03-11_Pulse_single_cell_DsRED_inhibitor_126h_new_dose_terms_seed_0/hypervolumes.pkl"
# with open(path3, "rb") as fid:
#     hvs3 = pickle.load(fid)
# print(hvs3, first_seen(hvs3))

# path4 = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Pulse_seed_single_cell_DsRED_inhibitor/2024-03-11_Pulse_single_cell_DsRED_inhibitor_t_pulse_126h_new_dose_terms_seed_0/hypervolumes.pkl"
# with open(path4, "rb") as fid:
#     hvs4 = pickle.load(fid)
# print(hvs4)

# path5 = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Pulse_seed_single_cell_DsRED_inhibitor/2024-03-11_Pulse_single_cell_DsRED_inhibitor_3obj_126h_new_dose_terms_seed_0/hypervolumes.pkl"
# with open(path5, "rb") as fid:
#     hvs5 = pickle.load(fid)
# print(hvs5, first_seen(hvs5))


#########################################################
### Loading results for opt hyperparams 10 seed runs  ###
#########################################################
# path_amp_const = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/2024-03-25_Amplifier_single_cell_const_dose_opt_hp_seed_"
# min_obj_path = "minimum_obj_all_gens.pkl"
# min_circuit_path = "min_obj_circuit_all_gens.pkl"
# for seed in range(0,10):
#     full_path = path_amp_const + str(seed) + "/"
#     with open(full_path + min_obj_path, "rb") as fid:
#         min_objs_const = pickle.load(fid)
#     with open(full_path + min_circuit_path, "rb") as fid:
#         min_circuit_const = pickle.load(fid)
#     print(min_objs_const[-1])
#     print(min_circuit_const[-1][0].edge_list)
#     print(min_circuit_const[-1][0].dose)

# path_amp_vary = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/2024-03-25_Amplifier_single_cell_vary_dose_opt_hp_seed_"
# min_obj_path = "minimum_obj_all_gens.pkl"
# min_circuit_path = "min_obj_circuit_all_gens.pkl"
# for seed in range(0,10):
#     full_path = path_amp_vary + str(seed) + "/"
#     with open(full_path + min_obj_path, "rb") as fid:
#         min_objs_vary = pickle.load(fid)
#     with open(full_path + min_circuit_path, "rb") as fid:
#         min_circuit_vary = pickle.load(fid)
#     print(min_objs_vary[-1], first_seen(min_objs_vary))
#     print(min_circuit_vary[-1][0].edge_list)
#     print(min_circuit_vary[-1][0].dose)

# path_sc = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_single_cell_inhibitor/Optimized_hyperparams_fixed_pop_ngen50/"
# results_runs = "2024-03-26_Signal_cond_single_cell_inhibitor_opt_hp_seed_"
# hv_path = "hypervolumes.pkl"
# hv_lists = []
# for seed in range(0,10):
#     full_path = path_sc + results_runs + str(seed) + "/"
#     with open(full_path + hv_path, "rb") as fid:
#         hypervolumes = pickle.load(fid)
#     print(hypervolumes[-1], first_seen(hypervolumes))
#     hv_lists.append(hypervolumes)
# plot_hypervolumes_set(path_sc+"all_hypervolume_progressions_zoomed.svg", 60, hv_lists, 42)

# path_sc_vary = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_single_cell_inhibitor/Optimized_hyperparams_vary_pop_ngen50/"
# results_runs = "2024-03-28_Signal_cond_single_cell_inhibitor_opt_hp_seed_"
# hv_path = "hypervolumes.pkl"
# hv_lists = []
# for seed in range(0,10):
#     full_path = path_sc_vary + results_runs + str(seed) + "/"
#     with open(full_path + hv_path, "rb") as fid:
#         hypervolumes = pickle.load(fid)
#     print(hypervolumes[-1], first_seen(hypervolumes))
#     hv_lists.append(hypervolumes)
# plot_hypervolumes_set(path_sc_vary+"all_hypervolume_progressions_zoomed.svg", 50, hv_lists, 42.0)

# path_sc_dsr = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_single_cell_DsRED_inhibitor/Optimized_hyperparams_fixed_pop_ngen50/"
# results_runs = "2024-03-28_Signal_cond_single_cell_DsRED_inhibitor_opt_hp_seed_"
# hv_path = "hypervolumes.pkl"
# hv_lists = []
# for seed in range(0,10):
#     full_path = path_sc_dsr + results_runs + str(seed) + "/"
#     with open(full_path + hv_path, "rb") as fid:
#         hypervolumes = pickle.load(fid)
#     print(hypervolumes[-1], first_seen(hypervolumes))
#     hv_lists.append(hypervolumes)
# plot_hypervolumes_set(path_sc_dsr+"all_hypervolume_progressions.svg", 50, hv_lists)

# path_sc_dsr_ngen80 = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/SC_seed_single_cell_DsRED_inhibitor/Optimized_hyperparams_fixed_pop_ngen80/"
# results_runs = "2024-03-28_Signal_cond_single_cell_DsRED_inhibitor_opt_hp_seed_"
# hv_path = "hypervolumes.pkl"
# hv_lists = []
# for seed in range(0,10):
#     full_path = path_sc_dsr_ngen80 + results_runs + str(seed) + "/"
#     with open(full_path + hv_path, "rb") as fid:
#         hypervolumes = pickle.load(fid)
#     print(hypervolumes[-1], first_seen(hypervolumes))
#     hv_lists.append(hypervolumes)
# plot_hypervolumes_set(path_sc_dsr_ngen80+"all_hypervolume_progressions.svg", 80, hv_lists)


# path_pulse = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/2024-03-25_Pulse_single_cell_opt_hp_seed_"
# hv_path = "hypervolumes.pkl"
# for seed in range(0,2):
#     full_path = path_pulse + str(seed) + "/"
#     with open(full_path + hv_path, "rb") as fid:
#         hypervolumes = pickle.load(fid)
#     print(hypervolumes[-1])

# path_pulse = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Pulse_seed_single_cell_DsRED_inhibitor/Optimized_hyperparams_fixed_pop_opt_stderr/"
# results_runs = "2024-04-12_Pulse_single_cell_opt_stderr_hp_seed_"
# hv_path = "hypervolumes.pkl"
# hv_lists = []
# for seed in range(0,10):
#     full_path = path_pulse + results_runs + str(seed) + "/"
#     with open(full_path + hv_path, "rb") as fid:
#         hypervolumes = pickle.load(fid)
#     print(hypervolumes[-1], first_seen(hypervolumes))
#     hv_lists.append(hypervolumes)
# plot_hypervolumes_set(path_pulse+"all_hypervolume_progressions.svg", 50, hv_lists)

# repository_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# results_base_path = "Amp_seed_pop_vary_dose/Single_cell_model_opt_hyperparams/2024-04-02_Amplifier_pop_vary_dose_single_cell_opt_hp"
# for seed in range(0, 10):
#     results_full_path = results_base_path.removesuffix("_seed_" + str(seed-1))
#     results_full_path = results_full_path + "_seed_" + str(seed) + "/"
#     with open(repository_path + results_full_path + 
#         "all_objectives.pkl", "rb") as fid:
#         all_obj = pickle.load(fid)
#     all_obj = np.abs(all_obj)
#     max_obj_idx = np.argmax(all_obj)
#     max_obj = all_obj[max_obj_idx]
#     #import unique circuits
#     with open(repository_path + results_full_path + 
#             "all_circuits.pkl", "rb") as fid:
#         all_circuits = pickle.load(fid)
#     max_circuit = all_circuits[max_obj_idx]
#     print(max_circuit[0].dose)
#     max_circuit[0].plot_graph()
    
repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
results_path = "Amp_seed_pop_vary_dose/Original_hyperparams_worked_well/2024-04-23_Amplifier_pop_vary_dose_original_hp_seed_2/"
amp = Amplifier("P1", [5, 75, 5], 2, True, False, {1: 46, 2: 122}, 2, 0.32, 0.57, False, True)

with open(repo_path + results_path + 
            "all_circuits.pkl", "rb") as fid:
        all_circuits = pickle.load(fid)

for circuit in all_circuits:
      edge_list = circuit[0].edge_list
      flattened_edges = list(chain.from_iterable(edge_list))
      if "Z9" in flattened_edges and "Z2" in flattened_edges and len(flattened_edges) == 10:
             print(edge_list)
             print(amp.func(circuit[0]))

        # if edge_list == [('P1', 'Z9'), ('Z9', 'Z9'), ('Z9', 'Z6'), ('Z6', 'Z6'), ('Z6', 'Z9'), ('Z6', 'Rep')]:
        #     if edge_list == [('P1', 'Z9'), ('Z9', 'Rep'), ('Z9', 'Z9')]:
        # print(amp.func(circuit[0]))



# obj = amp.func(circuit_1_exp)
# print(obj)
# edge_list = circuit_1_exp.edge_list
# flattened_edges = list(chain.from_iterable(circuit_1_exp.edge_list))
# if "Z6" in flattened_edges:
#     print(edge_list)
# ON_rel_neg = amp.func(circuit_1_exp)
# print(ON_rel_neg)