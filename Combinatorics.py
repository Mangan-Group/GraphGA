from graph_operations import *


def check_valid(g, num_parts):
    graph_parts = [i[0] for i in g.nodes]
    if ('P' not in graph_parts) | ('R' not in graph_parts) | ('Z' not in graph_parts):
        return 0
    elif (len(graph_parts) - 2) < num_parts:
        return 0

    for n in g.nodes:
        if n[0] == 'P':
            out_types = set([i[0] for i in list(g.successors(n))])
            if 'Z' not in out_types:
                return 0
        elif n[0] == 'R':
            in_types = set([i[0] for i in list(g.predecessors(n))])
            if (not in_types) | (('I' in in_types) & ('Z' not in in_types)):
                return 0
        else:
            in_nodes = list(g.predecessors(n))
            if not in_nodes:
                return 0
            elif in_nodes == [n]:
                return 0
            else:
                in_types = set([i[0] for i in in_nodes])
                if ('I' in in_types) & ('Z' not in in_types):
                    return 0
            if len(list(nx.all_simple_paths(g, n, 'Rep'))) == 0:
                return 0
    return 1


def get_combinatorics(promo_node, num_part, min_dose, max_dose, dose_interval, inhibitor):
    combo = []

    if not inhibitor:
        #for i in range(1, max_part + 1):
        combo.extend(list(combinations(tf_list, num_part)))
    else:
        #for num_part in range(2, max_part + 1):
        for num_tf in range(1, num_part):
            num_in = num_part - num_tf
            list1 = combinations(tf_list, num_tf)
            list2 = combinations(inhibitor_list, num_in)
            combo_two = [(i[0] + i[1]) for i in list(product(list1, list2))]
            combo.extend(combo_two)

    circuits = []

    for i in range(len(combo)):
        part_list = combo[i]
        full_edge_list = get_full_connected(part_list, promo_node)
        dose_list = dict(zip(part_list, [max_dose]*len(part_list)))
        circuits.append(Topo(full_edge_list, dose_list, promo_node))
        max_remove = len(full_edge_list) - 2*len(part_list)

        for j in range(1, max_remove+1):
            combo_remove = list(combinations(full_edge_list, j))
            for k in combo_remove:
                edges = [e for e in full_edge_list if e not in k]
                if check_valid(nx.DiGraph(edges), len(part_list)) == 1:
                    circuits.append(Topo(edges, dose_list, promo_node))
    return circuits




#sampling = get_combinatorics('P1', 3, 0, 75, 0, False)
#with open("Amplifier_3.pkl", "wb") as fid:
#	pickle.dump(sampling, fid)

#sampling = get_combinatorics('P1', 2, 0, 75, 0, True)
#with open("SigCond.pkl", "wb") as fid:
#	pickle.dump(sampling, fid)
