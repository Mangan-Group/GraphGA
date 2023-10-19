import numpy as np
from copy import deepcopy
from GA import *

### this check_valid function looks out of date, so now using one from
### GA file
# def check_valid(g, num_parts):
#     graph_parts = [i[0] for i in g.nodes]
#     if ('P' not in graph_parts) | ('R' not in graph_parts) | ('Z' not in graph_parts):
#         return 0
#     elif (len(graph_parts) - 2) < num_parts:
#         return 0

#     for n in g.nodes:
#         if n[0] == 'P':
#             out_types = set([i[0] for i in list(g.successors(n))])
#             if 'Z' not in out_types:
#                 return 0
#         elif n[0] == 'R':
#             in_types = set([i[0] for i in list(g.predecessors(n))])
#             if (not in_types) | (('I' in in_types) & ('Z' not in in_types)):
#                 return 0
#         else:
#             in_nodes = list(g.predecessors(n))
#             if not in_nodes:
#                 return 0
#             elif in_nodes == [n]:
#                 return 0
#             else:
#                 in_types = set([i[0] for i in in_nodes])
#                 if ('I' in in_types) & ('Z' not in in_types):
#                     return 0
#             if len(list(nx.all_simple_paths(g, n, 'Rep'))) == 0:
#                 return 0
#     return 1


def get_combinatorics(promo_node, num_part, min_dose, max_dose, dose_interval, inhibitor):
    combo = []

    if not inhibitor:
        #for i in range(1, max_part + 1): (outdated)
        # get all combinations of tfs (order doesn't matter- (Z6, Z2) = (Z2, Z6))
        combo.extend(list(combinations(tf_list, num_part)))
    else:
        #for num_part in range(2, max_part + 1): (outdated)
        # get all combinations of tf and inhibitor (1 of each for 2 part)
        for num_tf in range(1, num_part):
            num_in = num_part - num_tf
            list1 = combinations(tf_list, num_tf)
            list2 = combinations(inhibitor_list, num_in)
            combo_two = [(i[0] + i[1]) for i in list(product(list1, list2))]
            combo.extend(combo_two)

    circuits = []

    for i in range(len(combo)):
        part_list = combo[i]
        # all nodes have P1->node->rep, node1->node2, node2->node1,
        # node1->node1, node2->node2 if 2 parts or
        # P1->node->rep, node1->node1 if 1 part
        # and all permutations from part_list with number of parts
        full_edge_list = get_full_connected(part_list, promo_node)
        # set doses to max_dose for each part
        dose_list = dict(zip(part_list, [max_dose]*len(part_list)))
        # add fully connected circuit to list 
        circuits.append(Topo(full_edge_list, dose_list, promo_node))
        # max number of edges leaves at least 1 edge (validity
        # checked in check_valid)
        max_remove = len(full_edge_list) - 1 #2*len(part_list) (outdated)

        for j in range(1, max_remove+1):
            # get all possible combinations of all allowable amounts of
            # edges to remove
            combo_remove = list(combinations(full_edge_list, j))
            for k in combo_remove:
                #remove each edge in combo_remove
                edges = [e for e in full_edge_list if e not in k]
                #check that circuit is still valid
                if check_valid(nx.DiGraph(edges), promo_node, part_list) == 1:
                    circuits.append(Topo(edges, dose_list, promo_node))
    return circuits

def add_dose_variation(circuits, num_part, min_dose, max_dose, dose_interval):

    # get list of doses based on specifications
    doses = np.arange(min_dose, max_dose+1, dose_interval)

    dose_varied_circuits = []
    # if 1 part, use each possible dose for the part
    if num_part == 1:
        for circuit in circuits:
            part_list = circuit.part_list
            for dose in doses:
                circuit_copy = deepcopy(circuit)
                circuit_copy.dose[part_list[0]] = dose
                dose_varied_circuits.append(circuit_copy)

    if num_part == 2:
        dose_combos_2part = list(product(doses, repeat= 2))
        for circuit in circuits:
            part_list = circuit.part_list
            for dose in dose_combos_2part:
                circuit_copy = deepcopy(circuit)
                circuit_copy.dose[part_list[0]] = dose[0]
                circuit_copy.dose[part_list[1]] = dose[1]
                dose_varied_circuits.append(circuit_copy)

    return dose_varied_circuits


circuits_1part = get_combinatorics('P1', 1, 0, 75, 0, False)
# print(len(circuits))

dose_varied_circuits_1part = add_dose_variation(circuits_1part, 1, 5, 75, 5)
# for circuit in dose_varied_circuits_1part:
#     print(circuit.dose)
print(len(dose_varied_circuits_1part))

circuits_2part = get_combinatorics('P1', 2, 0, 75, 0, False)
# print(len(circuits_2part))
dose_varied_circuits_2part = add_dose_variation(circuits_2part, 2, 5, 75, 5)
print(len(dose_varied_circuits_2part))

# circuits = get_combinatorics('P1', 2, 0, 75, 0, False)
# print(len(circuits))

#sampling = get_combinatorics('P1', 3, 0, 75, 0, False)
#with open("Amplifier_3.pkl", "wb") as fid:
#	pickle.dump(sampling, fid)

# sampling = get_combinatorics('P1', 2, 0, 75, 0, True)
# print(len(sampling))
#with open("SigCond.pkl", "wb") as fid:
#	pickle.dump(sampling, fid)
