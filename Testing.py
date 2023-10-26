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
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from rankcrowding import RankAndCrowding


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

# if pop:
    # use obj and population bc they are the final
    # generation (outside of for loop)
    # obj = obj*-1;
    # final_min_obj = obj_min[-1]*-1
    # index_obj_within_CI = np.argwhere(obj >= final_min_obj-CI)
    # final_obj_within_CI = obj[index_obj_within_CI]
    # final_circuits_within_CI = population[index_obj_within_CI]
    # save final_obj_within_CI and final_circuits_within_CI
    # also save all_objectives and all_circuits
    
    # for circuit in final_circuits_within_CI:
        # plot graph
# else:
    # save min_obj_all_gens
    # save min_obj_circuit_all_gens
    # plot circuit_min graph
# save all_objectives and all_circuits in both cases
# if get_unique:
#     do unique stuff




# CI = 0.4727
# num_circuit = 52
# final_obj = all_obj[-num_circuit:, 0]*-1
# # sorted_index = np.lexsort(final_obj)
# # sorted_final_obj = final_obj[sorted_index]
# final_obj_within_CI = final_obj[np.argwhere(final_obj >= final_obj[-1]-CI)]

# print(final_obj_within_CI)

# sorted_final_obj = sorted(final_obj)
# print(sorted_final_obj)

# for i in range(1, 20):
#     print(i)


#mask = (10 < a[:, 0]) & (a[:, 0] < 15)
# idx = np.flatnonzero(mask)
arr1 = np.array([[3, 11], [1, 9], [5, 10], [0, 8], [-1, 11]])

arr2 = np.array([[0, 7], [0, 11], [6, 9], [2, 11]])

# arr4 = np.array([[0, 7, 5], [0, 11, 5], [6, 9, 5], [2, 11, 5]])

# df = pd.DataFrame(arr4, columns=["col1", "col2", "col3"])
# print(df)

# print(df[["col1", "col2"]].to_numpy())
# print(arr1)
# all_results = []


# i = 0
# idx = []
# for row1 in arr1:

#     for row2 in arr2:
#         if ((row1[0] <= row2[0]) & (row1[1] <= row2[1])):
#             print(row1)
#             idx.append(i)
#             print(i)
#             break
#     i += 1

# print(arr1[idx])

# arr3 = np.argwhere(arr1[:, 0] > 1)

# arr3 = np.argwhere((arr1[:, 0] <= arr2[:, 0]) & (arr1[:, 1] <= arr2[:, 1]))#.flatten()
# print(arr3)
# print(idx, arr1[idx])

with open("/Users/kdreyer/Documents/Github/GraphGA/GA_results/2023-10-24_Signal_cond_pop_inhibitor_seed_0/all_objectives.pkl", "rb") as fid:
    all_obj = pickle.load(fid)

# print(all_obj)

with open("/Users/kdreyer/Documents/Github/GraphGA/GA_results/2023-10-24_Signal_cond_pop_inhibitor_seed_0/all_circuits.pkl", "rb") as fid:
    all_circuits = pickle.load(fid)

with open("/Users/kdreyer/Documents/Github/GraphGA/GA_results/2023-10-24_Signal_cond_pop_inhibitor_seed_0/final_objectives_df_with_type.pkl", 'rb') as f:
    obj_df_type = pd.read_pickle(f)
# print(obj_df_type)

obj = obj_df_type[["ON_rel", "FI_rel"]].to_numpy()
# print(obj)
types = list(obj_df_type[["type"]].to_numpy().flatten())
# print(types)


CI = [1.2891, 0.0042]
nds = RankAndCrowding()

S_all = nds.do(all_obj, len(all_obj))
sorted_all_obj = all_obj[S_all, :]*-1
sorted_all_circuits = all_circuits[S_all]

            # (all_obj[:, 0] >= obj[:, 0] - problem.CI[0]) &
            # (all_obj[:, 1] >= obj[:, 1] - problem.CI[1])

i = 0
index_obj_within_CI = []
for row1 in sorted_all_obj:
    for row2 in obj:
        if ((row1[0] >= row2[0] - CI[0]) & (row1[1] >= row2[1] - CI[1])):
            # print(row1, row2)
            index_obj_within_CI.append(i)
            # print(i)
            break
    i += 1

final_objs_within_CI = sorted_all_obj[index_obj_within_CI]
final_circuits_within_CI = sorted_all_circuits[index_obj_within_CI]

edges_circuits_within_CI = []
for circuit in final_circuits_within_CI:
    circuit_edges = circuit[0].edge_list
    for key, val in circuit[0].dose.items():
        circuit_edges.append((key, val))
    edges_circuits_within_CI.append(circuit_edges)
unique_edges_set = set(map(frozenset, edges_circuits_within_CI))


types_CI = []
for topo in final_circuits_within_CI:
    inhib = "Activators"
    for part in topo[0].part_list:
        if part[0] == "I":
            inhib = "Inhibitors"
    types_CI.append(inhib)
types.extend(types_CI)

obj_within_CI = np.vstack((obj, final_objs_within_CI))
obj_within_CI_df = pd.DataFrame(obj_within_CI, columns=["ON_rel", "FI_rel"])
obj_within_CI_df["type"] = types

# print(obj_within_CI_df)



# print(len(sorted_all_obj), len(index_obj_within_CI))
# print(len(unique_edges_set))
# print(arr1[idx])

# folder_path = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/2023-10-24_Signal_cond_pop_inhibitor_seed_0/"

# file_name = "all_obj_within_CI_df_with_type.pkl"
# obj_within_CI_df.to_pickle(folder_path + "/" + file_name)

# file_name = "final_objs_within_CI.pkl"
# with open(folder_path + "/" + file_name, "wb") as fid:
#     pickle.dump(final_objs_within_CI, fid)

# file_name = "final_circuits_within_CI.pkl"
# with open(folder_path + "/" + file_name, "wb") as fid:
#     pickle.dump(final_circuits_within_CI, fid)

# file_name = "unique_edges_set.pkl"
# with open(folder_path + "/" + file_name, "wb") as fid:
#     pickle.dump(unique_edges_set, fid)

# fig, ax = plt.subplots(1, 1, figsize= (4, 4))
# sns.scatterplot(data=obj_within_CI_df, x= obj_within_CI_df["ON_rel"],
#                 y= obj_within_CI_df["FI_rel"], hue='type', 
#                 palette="colorblind", ax=ax)
# plt.title("Signal Conditioner Population CI Pareto Front")
# plt.ylabel("FI_rel")
# plt.xlabel("ON_rel")
# plt.savefig(folder_path +"population_CI_pareto_front.svg", bbox_inches="tight")

# fig, ax = plt.subplots(1, 1, figsize= (4, 4))
# sns.scatterplot(data=obj_df_type, x= obj_df_type["ON_rel"],
#                 y= obj_df_type["FI_rel"], hue='type', 
#                 palette="colorblind", ax=ax)
# plt.title("Signal Conditioner Population Pareto Front")
# plt.ylabel("FI_rel")
# plt.xlabel("ON_rel")
# plt.savefig(folder_path +"population_pareto_front.svg", bbox_inches="tight")



# path1 = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Amp_seed_pop_const_dose/2023-10-23_Amplifier_pop_const_dose_seed_0/final_circuits_within_CI.pkl"
# with open(path1, "rb") as fid:
#     CI = pickle.load(fid)

# path2 = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/outdated/Amp_seed_pop_const_dose/2023-10-23_Amplifier_pop_const_dose_seed_0/unique_edges_set.pkl"
# with open(path2, "rb") as fid:
#     circuits = pickle.load(fid)

# print(type(circuits))
# print(circuits[:3])
# circuits_set_list = list(circuits)
# circuits_list = [list(edges) for edges in circuits_set_list]
# circuits_list_short = circuits_list[:3]
# print(circuits_list_short)
# circuits_tuple = [tuple(i) for i in circuits_list_short]
# print(circuits_tuple)
# print(circuits_list_short)
# print([list(i[0]) for i in circuits_list_short])
# circuits_list_new = [l[:-1] for l in circuits_list_short]
# print(circuits_list_new)

tuple_list = [(('Z9', 'Z6'),('Z15', 'Rep'), ('P1', 'Z9')), (('Z15', 'Rep'),('Z9', 'Z6'), ('P1', 'Z9'))]
# tuple_combo = [sub_tuple[0] + sub_tuple[1] for tuple_ in tuple_list for sub_tuple in tuple_]

# print(tuple_combo)
tuple_combo = []
for tuple_ in tuple_list:
    tuple_new = tuple([sub_tuple[0] + sub_tuple[1] for sub_tuple in tuple_])
    tuple_combo.append(tuple_new)

# print(tuple_combo)
# tuple_ = ("hi", "helo")
# print(tuple_[0] + tuple_[1])

tuple_list2 = []
idx_list = []
seen = set()
for idx, tuple_ in enumerate(tuple_combo):
    tuple_s = frozenset(tuple_)
    # tuple_
    if tuple_s not in seen:
        seen.add(tuple_s)
        idx_list.append(idx)
        tuple_list2.append(tuple_)

print(tuple_list2)
print(idx_list)
print(seen)

# edges_unique_circuits = []
# for circuit in sorted_all_circuits:
#     circuit_edges = circuit[0].edge_list
#     for key_val in circuit[0].dose.items():
#         circuit_edges.append((key, val))
#     edges_unique_circuits


