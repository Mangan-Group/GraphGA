# -*- coding: utf-8 -*-
"""
Created on Wed Aug 31 14:16:07 2022

@author: Katie_Dreyer
"""
import numpy as np
import networkx as nx
import time
import pickle
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from define_circuit import Topo
from itertools import product, combinations, permutations
# from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
# from rankcrowding import RankAndCrowding
from plot_search_results import *


# from load_files_pop import (
#     Z_20,  
#     Ref_pop20, 
#     Ref, 
#     tf_list, 
#     inhibitor_list
# )
# # print(tf_list + inhibitor_list + ['Rep'])
edge_list = [('P1', 'Z1'), ('P1', 'I6'), ('I6', 'Z2'), ('Z1', 'Z2'), ('Z1', 'Rep')]
dose_list = {'Z2': 25, 'Z1': 75, 'I6': 75}
promo_node = 'P1'
topology = Topo(
    edge_list,
    dose_list, 
    promo_node
)
print(topology.in_dict["Z2"]["Z"][0])
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

# path2 = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/outdated/Amp_seed_pop_const_dose/2023-10-23_Amplifier_pop_const_dose_seed_0/all_objectives.pkl"
# with open(path2, "rb") as fid:
#     all_obj = pickle.load(fid)
# circuit_edge_lists = []
# for circuit in all_circuits:
#     circuit_edges = circuit[0].edge_list
#     for key, val in circuit[0].dose.items():
#         circuit_edges.append((key, str(val)))
#     circuit_edge_lists.append(circuit_edges)
#     print("finished circuit " + str(len(circuit_edge_lists)))
# # print(circuit_edge_lists[-1])
# combo_edges_lists = []
# for edges in circuit_edge_lists:
#     edge_combos = list([edge[0] + edge[1] for edge in edges])
#     combo_edges_lists.append(edge_combos)
# # print(len(combo_edges_lists))

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


##################################################################
### outdated converting to tuples first (not necessary)###
# tuple_combo = []
# for tuple_ in tuple_list:
#     tuple_new = tuple([sub_tuple[0] + sub_tuple[1] for sub_tuple in tuple_])
#     tuple_combo.append(tuple_new)

# tuple_list2 = []
# idx_list = []
# seen = set()
# for idx, tuple_ in enumerate(tuple_combo):
#     tuple_s = frozenset(tuple_)
#     # tuple_
#     if tuple_s not in seen:
#         seen.add(tuple_s)
#         idx_list.append(idx)
#         tuple_list2.append(tuple_)