# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 14:16:07 2022

@author: Katie_Dreyer
"""
import os
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
# from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
# from rankcrowding import RankAndCrowding
from plot_search_results import *
from GA import sampling, check_valid
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator
from pymoo.indicators.hv import HV
from plot_search_results import plot_graph
from math import exp
from scipy.interpolate import interp2d

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
# # amp_results_path = "Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/all_circuits.pkl"
# amp_unique_obj_path = "Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/unique_objectives.pkl"
# # amp_unique_circuits_path = "Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/unique_circuits.pkl"
# amp_obj_path = "Amp_seed_pop_vary_dose/2023-10-31_Amplifier_pop_vary_dose_seed_0/all_objectives.pkl"
# # amp_results = pd.read_pickle(repo_path+amp_results_path)
# amp_obj = pd.read_pickle(repo_path+amp_obj_path)*-1
# amp_obj = amp_obj.flatten()
# print(max(amp_obj))
# amp_unique_obj = pd.read_pickle(repo_path+amp_unique_obj_path)
# amp_unique_obj = amp_unique_obj.flatten()*-1
# print(max(amp_unique_obj))
# amp_obj = amp_obj.flatten()
# print(amp_obj)
# pulse = PulseGenerator("P1", [5, 75, 5], 2, True, False, {1: 46, 2: 122}, 2, 0.32, 0.57, False, True)
# circuit = Topo([('P1', 'Z6'), ('Z6', 'Rep')], {'Z6': 75}, "P1")
# # print(circuit.dose)
# pulse.simulate(circuit)
# amp_unique_obj = pd.read_pickle(repo_path+amp_unique_obj_path).flatten()*-1
# amp_unique_circuits = pd.read_pickle(repo_path+amp_unique_circuits_path)
# amp_unique_obj_high = np.where(amp_unique_obj > 69.5)[0]
# unique_circuits_high = amp_unique_circuits[amp_unique_obj_high]
# unique_obj_high = amp_unique_obj[amp_unique_obj_high]
# print(unique_obj_high[31])
# print(unique_circuits_high)
# edges_no_ZF6 = []
# edges_no_ZF6_idx = []
# for i, circuit in enumerate(unique_circuits_high):
#     edge_list = circuit[0].edge_list
#     edge_list_flattened = list(chain.from_iterable(edge_list))
#     if "Z6" not in edge_list_flattened:
#         edges_no_ZF6.append(edge_list)
#         edges_no_ZF6_idx.append(i)
# print(edges_no_ZF6)
# for i, circuit in enumerate(unique_circuits_high):
#     edge_list = circuit[0].edge_list
#     edge_list_flattened = list(chain.from_iterable(edge_list))
#     if len(edge_list_flattened) == 10:
#         print(edge_list, i)

# print(len(amp_unique_obj_high[0]))
### Circuit 1 from experiments
# print(np.where(np.logical_and(amp_obj>69.982, amp_obj<69.98215)))
# print(amp_obj[[1874, 1940, 1951, 1965, 2027, 2074, 2108, 2112, 2165, 2175, 2226,
#        2236, 2253, 2260, 2278, 2300, 2311, 2324, 2331, 2342, 2497, 2516,
#        2573]])
# amp_results[1874][0].plot_graph()
# circuit_1 = amp_results[1874][0]
# print(circuit_1.dose)
# circuit_1_exp = Topo([('P_exp_amp', 'Z9'), ('Z9', 'Z9'), ('Z9', 'Z6'), ('Z6', 'Z6'), ('Z6', 'Rep'), ('Z6', 'Z9')], {'Z6': 75, 'Z9': 75}, "P_exp_amp")
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
# print(circuit_2.edge_list)
# circuit_2_exp = Topo([('P_exp_amp', 'Z2'), ('Z2', 'Z2'), ('Z2', 'Z6'), ('Z6', 'Z6'), ('Z6', 'Rep')], {'Z6': 75, 'Z2': 75}, "P_exp_amp")
# ON_rel_neg = amp.func(circuit_2_exp)
# print(ON_rel_neg)

### Circuit 3 from experiments
# print(np.where(np.logical_and(amp_obj>39.8308, amp_obj<39.83095)))
# print(amp_obj[1251])
# amp_results[1251][0].plot_graph()
# circuit_3 = amp_results[1251][0]
# print(circuit_3.edge_list)
# circuit_3_exp = Topo([('P_exp_amp', 'Z2'), ('Z2', 'Z2'), ('Z2', 'Z9'), ('Z9', 'Z9'), ('Z9', 'Rep')], {'Z2': 75, 'Z9': 75}, "P_exp_amp")
# ON_rel_neg = amp.func(circuit_3_exp)
# print(ON_rel_neg)

### Circuit 4 from experiments
# print(np.where(np.logical_and(amp_obj>39.8018, amp_obj<39.802)))
# print(amp_obj[22])
# print(amp_results[22][0].plot_graph())
# circuit_4 = amp_results[22][0]
# circuit_4_exp = Topo([('P_exp_amp', 'Z9'), ('Z9', 'Rep'), ('Z9', 'Z9')], {'Z9': 75}, "P_exp_amp")
# print(circuit_4.dose)
# ON_rel_neg = amp.func(circuit_4_exp)
# print(ON_rel_neg)


#############################################################################
########### Signal Conditioner experimental circuit new predictions #########
#############################################################################
repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
sc_results_path = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/final_population.pkl"
sc_obj_path = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/final_objectives_df.pkl"
sc_all_obj_path = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/all_objectives.pkl"
sc_all_circuits_path = "SC_seed_pop_DsRED_inhibitor/2023-11-30_Signal_Cond_pop_DsRED_inhibitor_ngen60_seed_0/all_circuits.pkl"
sc_all_obj = pd.read_pickle(repo_path+sc_all_obj_path)*-1
sc_pareto_obj = pd.read_pickle(repo_path+sc_obj_path)
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

sc_results = pd.read_pickle(repo_path+sc_all_circuits_path)
print(len(sc_results))
sc_obj = pd.read_pickle(repo_path+sc_obj_path)[["ON_rel", "FI_rel"]]*-1
# print(sc_obj[sc_obj["FI_rel"] > 2].drop_duplicates().sort_values(by="FI_rel", ascending=False))

sc = SignalConditioner("P1", [5, 75, 5], 2, True, True, {1: 46, 2: 122}, 2, 0.32, 0.57, False, True)
not_equal = []
for i, circuit in enumerate(sc_results[:1000]):
    in_dict_list = list(circuit[0].in_dict.keys())
    in_dict_list.remove("Rep") 
    dose_dict_list = list(circuit[0].dose.keys())
    dose_dict_list.remove("Rep")
    # print(in_dict_list)
    if in_dict_list != dose_dict_list:
        not_equal.append(circuit)
        print(i)
print(len(not_equal))
        # print("not equal")
### Circuit 1 from experiments
# sc_results[1][0].plot_graph()
# circuit_1 = sc_results[1][0]
# print(circuit_1.in_dict, circuit_1.dose)
# circuit_1_exp = Topo([('P1', 'Z12'), ('Z12', 'Rep'), ('Z12', 'I1'), ('P1', 'I1'), ('I1', 'Rep')], {'Z12': 60, 'I1': 60}, "P1")
# print(circuit_1_exp.part_list)
# if ("Z12" in circuit_1_exp.part_list) and ("I1" in circuit_1_exp.part_list):
#     print(circuit_1_exp.edge_list)
# objs_neg, FI_sc = sc.func(circuit_1_exp)
# print(objs_neg, FI_sc)

### Circuit 2 from experiments
# sc_results[7][0].plot_graph()
circuit_2 = sc_results[7][0]
print(circuit_2.dose, circuit_2.in_dict)
# circuit_2_exp = Topo([('P_exp_sc', 'Z12'), ('P_exp_sc', 'I11'), ('Z12', 'Rep'), ('Z12', 'I11'), ('I11', 'Rep')], {'Z12': 55, 'I11': 75}, "P_exp_sc")
# objs_neg, FI_sc = sc.func(circuit_2_exp)
# print(objs_neg, FI_sc)

### Circuit 3 from experiments
# sc_results[12][0].plot_graph()
# circuit_3 = sc_results[12][0]
# # print(circuit_3.dose, circuit_3.in_dict)
# circuit_3_exp = Topo([('P1', 'Z9'), ('P1', 'I11'), ('Z9', 'Rep'), ('Z9', 'I11'), ('I11', 'Rep')], {'I11': 35, 'Z9': 75}, "P1")
# objs_neg = sc.func(circuit_3_exp)
# print(objs_neg)

### Circuit 4 from experiments
# sc_results[24][0].plot_graph()
# circuit_4 = sc_results[24][0]
# print(circuit_4.edge_list)
# circuit_4_exp = Topo([('P_exp_sc', 'Z12'), ('Z12', 'Rep'), ('Z12', 'I1'), ('I1', 'Rep')], {'Z12': 55, 'I1': 5}, "P_exp_sc")
# objs_neg, FI_sc = sc.func(circuit_4_exp)
# print(objs_neg, FI_sc)

### Circuit 5 from experiments
# sc_results[21][0].plot_graph()
# circuit_5 = sc_results[21][0]
# print(circuit_5.edge_list, circuit_5.dose)
# circuit_5_exp = Topo([('P1', 'Z12'), ('Z12', 'Rep'), ('Z12', 'I15'), ('I15', 'Rep')], {'Z12': 55, 'I15': 20}, "P_exp_sc")
# objs_neg, FI_sc = sc.func(circuit_5_exp)
# print(objs_neg, FI_sc)

### Circuit 6 from experiments
# sc_results[11][0].plot_graph()
# circuit_6 = sc_results[11][0]
# print(circuit_6.edge_list, circuit_6.dose)
# circuit_6_exp = Topo([('P_exp_sc', 'Z12'), ('Z12', 'Rep'), ('Z12', 'I11'), ('I11', 'Rep')], {'I11': 5, "Z12": 50}, "P_exp_sc")
# objs_neg, FI_sc = sc.func(circuit_6_exp)
# print(objs_neg, FI_sc)

# repo_path ="/Users/kdreyer/Documents/Github/GraphGA/GA_results/"
# sc_results_path = "2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/final_population.pkl"
# sc_obj_path = "2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/final_objectives_df.pkl"
# sc_all_obj_path = "2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/all_objectives.pkl"
# sc_all_circuits_path = "2024-02-15_Signal_Cond_pop_DsRED_inhibitor_ZF1_ZF2_seed_0/all_circuits.pkl"
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





#### test plotting for CI with pareto front ####
# fig, ax = plt.subplots(1, 3, figsize=(12, 4))
# # x = [1, 1.5, 2, 3, 4, 5]
# # y = [exp(-i) for i in x]
# ON_rel_err = [2.0]*len(sc_pareto_obj["ON_rel"].tolist())
# FI_rel_err = [0.25]*len(sc_pareto_obj["FI_rel"].tolist())
# # upper_ON_rel = np.array([i+ON_rel_err[0] for i in (sc_pareto_obj["ON_rel"]*-1)])
# upper_ON_rel = np.array(sc_pareto_obj["ON_rel"]*-1)
# sorted_upper_idx = np.argsort(upper_ON_rel)
# sorted_upper_ON_rel = upper_ON_rel[sorted_upper_idx]
# # print(sorted_upper_ON_rel)
# lower_ON_rel = np.array([i-ON_rel_err[0] for i in (sc_pareto_obj["ON_rel"]*-1)])
# sorted_lower_idx = np.argsort(lower_ON_rel)
# sorted_lower_ON_rel = lower_ON_rel[sorted_lower_idx]
# # upper_FI_rel = np.array([i+FI_rel_err[0] for i in (sc_pareto_obj["FI_rel"]*-1)])
# upper_FI_rel = np.array(sc_pareto_obj["FI_rel"]*-1)
# sorted_upper_FI_rel = upper_FI_rel[sorted_upper_idx]
# lower_FI_rel = np.array([i-FI_rel_err[0] for i in (sc_pareto_obj["FI_rel"]*-1)])
# sorted_lower_FI_rel = lower_FI_rel[sorted_lower_idx]
# # print(sorted_lower_FI_rel)
# xfill = np.sort(np.concatenate([upper_ON_rel, lower_ON_rel]))
# y1fill = np.interp(xfill, sorted_upper_ON_rel, sorted_upper_FI_rel)
# y2fill = np.interp(xfill, sorted_lower_ON_rel, sorted_lower_FI_rel)
# # print(y2fill)
# ax[0].errorbar(sc_pareto_obj["ON_rel"]*-1, sc_pareto_obj["FI_rel"]*-1, xerr=ON_rel_err,
#                yerr=FI_rel_err, linestyle="None", marker="o", capsize=4)
# ax[1].fill_between(xfill, y1fill, y2fill, alpha=0.3)
# ax[1].plot(sc_pareto_obj["ON_rel"]*-1, sc_pareto_obj["FI_rel"]*-1, linestyle="None", marker="o")
# ax[1].plot(lower_ON_rel, lower_FI_rel, linestyle="None", marker="o")
# ax[1].plot(upper_ON_rel, upper_FI_rel, linestyle="None", marker="o")
# ax[2].plot(sc_all_obj[:, 0], sc_all_obj[:, 1], linestyle="None", marker="o", markersize=4)
# ax[2].fill_between(xfill, y1fill, y2fill, alpha=0.3)

# plt.show()

# all_vertices = [[[-31.0, -6.53452205539751, 0.0], [-29.0, -6.483134823574519, 0.0], [-28.0, -5.9889566791622055, 0.0], [-27.0, -5.9554042087325865, 0.0], [-26.0, -5.6419117140686605, 0.0], [-25.0, -5.4046441221955455, 0.0], [-24.0, -5.15955583568886, 0.0], [-23.0, -4.904103702363367, 0.0], [-21.0, -4.3487887373756, 0.0], [-19.0, -4.186657512107043, 0.0], [-17.0, -3.253241160149568, 0.0], [-13.0, -1.4178456939579807, 0.0], [0.0, -0.8621307599677381, 0.0], [1.0, 0.04792253338347758, 0.019736509901236654], [2.0, 0.09930976520646771, 0.03947301980247331], [3.0, 0.5934879096187822, 0.08691635738331432], [4.0, 0.6270403800484003, 0.3954655909939231], [5.0, 0.9405328747123269, 0.4277589480926045], [6.0, 1.1778004665854414, 0.7262357246427164], [6.0, 1.4228887530921275, 0.7262357246427164], [7.0, 1.6783408864176201, 0.9517388696511238], [8.0, 2.2336558514053877, 1.1847211015661638], [10.0, 2.3957870766739435, 1.4254010544752462], [12.0, 3.3292034286314194, 1.9376563950272714], [14.0, 3.63168294042376, 2.0606017603094613], [18.0, 5.164598894823007, 2.6219693891107267], [32.0, 5.720313828813249, 4.089996151753775], [34.0, 10.214127529204747, 4.142305897492051], [37.0, 90.74101026536327, 4.220770516099464], [65.0, 97.32345485414426, 1.4813962181584142]], [[-31.0, -6.53452205539751, -1.210437574896098], [-29.0, -6.483134823574519, -1.1709645550936247], [-28.0, -5.9889566791622055, -1.1235212175127838], [-27.0, -5.9554042087325865, -0.814971983902175], [-26.0, -5.6419117140686605, -0.7826786268034935], [-25.0, -5.4046441221955455, -0.4842018502533817], [-24.0, -5.15955583568886, -0.2586987052449743], [-23.0, -4.904103702363367, -0.02571647332993421], [-21.0, -4.3487887373756, 0.2149634795791482], [-19.0, -4.186657512107043, 0.7272188201311733], [-17.0, -3.253241160149568, 0.8501641854133632], [-13.0, -1.4178456939579807, 1.4115318142146287], [0.0, -0.8621307599677381, 2.7746995223831736], [1.0, 0.04792253338347758, 2.879558576857677], [2.0, 0.09930976520646771, 2.905713449726815], [3.0, 0.5934879096187822, 2.931868322595953], [4.0, 0.6270403800484003, 2.958023195465091], [5.0, 0.9405328747123269, 2.9841780683342285], [6.0, 1.1778004665854414, 3.0103329412033664], [6.0, 1.4228887530921275, 3.0103329412033664], [7.0, 1.6783408864176201, 2.9124981448483287], [8.0, 2.2336558514053877, 2.814663348493291], [10.0, 2.3957870766739435, 2.618993755783216], [12.0, 3.3292034286314194, 2.4233241630731412], [14.0, 3.63168294042376, 2.227654570363066], [18.0, 5.164598894823007, 1.8363153849429164], [32.0, 5.720313828813249, 0.4666282359723911], [34.0, 10.214127529204747, 0.2709586432623161], [37.0, 90.74101026536327, 0.2709586432623161], [65.0, 97.32345485414426, 0.2709586432623161]]]
# all_vertices = [[[-31.0, -6.53452205539751, 0.0], [-29.0, -6.483134823574519, 0.0],
#                 [-28.0, -5.9889566791622055, 0.0], [-27.0, -5.9554042087325865, 0.0]],
#                 [[-31.0, -6.53452205539751, -1.210437574896098], [-29.0, -6.483134823574519, -1.1709645550936247],
#                 [-28.0, -5.9889566791622055, -1.1235212175127838], [-27.0, -5.9554042087325865, -0.814971983902175]]]
# x1, y1, z1 = list(zip(*all_vertices[0]))
# x2, y2, z2 = list(zip(*all_vertices[1]))

# fig = plt.figure(figsize= (2.25, 2))
# ax = fig.add_subplot(projection='3d')
# ax.scatter(xs=x1, ys=y1, zs=z1, color="k")
# ax.scatter(xs=x2, ys=y2, zs=z2, color="k")
# ax.add_collection3d(Poly3DCollection(all_vertices, facecolors=sky_blue))
# ax.view_init(elev=10, azim=-115)
# plt.show()

# dose = {'I1': 5, 'Z2': 55, 'Rep': 1}
# print(list(dose.keys()))

# in_dict = {'I1': {'P': [], 'Z': ['Z2'], 'I': []}, 'Z2': {'P': ['P1'], 'Z': ['Z2'], 'I': ['I1']}, 'Rep': {'P': [], 'Z': ['Z2'], 'I': ['I1']}}
# print(list(in_dict.keys()))

# if list(in_dict.keys()) != list(dose.keys()):
#     print("not equal")