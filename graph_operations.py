import networkx as nx

from define_circuit import *
from itertools import combinations_with_replacement
from copy import deepcopy

def get_out_path(n, part_list):
    out_list = [k for k in part_list if k != n]
    out_path = [n]
    if len(out_list) > 0:
        num_connect = np.random.randint(len(out_list))
        out_path.extend(np.random.choice(out_list, num_connect, replace=False))
    out_path.append('Rep')
    edges = [(i, j) for i, j in zip(out_path[:-1], out_path[1:])]

    return edges

def get_in_path(n, promo_node, circuit_tf_list):
    # if promo_node is not None:
    #     in_node = np.random.choice(circuit_tf_list + [promo_node])
    # else:
    #     in_node = np.random.choice(circuit_tf_list)
    # edges = [(in_node, n)]
    # if in_node == n:
    #     new_tf_list = [k for k in circuit_tf_list if k != n]
    #     # in_node = np.random.choice(new_tf_list + [promo_node])
    #     if promo_node is not None:
    #         in_node = np.random.choice(new_tf_list + [promo_node])
    #     else:
    #         in_node = np.random.choice(new_tf_list)
    # edges.append((in_node, n))
    if promo_node is not None:
        in_node = np.random.choice(circuit_tf_list + [promo_node])
        edges = [(in_node, n)]
        if in_node != promo_node:
            in_path = [promo_node]
            num_connect = np.random.randint(len(circuit_tf_list))
            in_path.extend(np.random.choice(circuit_tf_list, num_connect, replace=False))
            in_path.append(n)
            edges.extend([(i, j) for i, j in zip(in_path[:-1], in_path[1:])])
    else:
        in_node = np.random.choice(circuit_tf_list)
        edges = [(in_node, n)]
        if in_node == n:
            in_node = np.random.choice([k for k in circuit_tf_list if k != n])
            edges.append((in_node, n))

    return edges

def get_edges(promo_node, part_list):
    circuit_tf_list = [k for k in part_list if k[0] == 'Z']
    edge_list = set([(promo_node, np.random.choice(circuit_tf_list)), (np.random.choice(circuit_tf_list), 'Rep')])
    for n in part_list:
        # edge_list.update([(np.random.choice(circuit_tf_list + [promo_node]), n)])
        in_edges = get_in_path(n, promo_node, circuit_tf_list)
        edge_list.update(in_edges)
        out_edges = get_out_path(n, part_list)
        edge_list.update(out_edges)

    return list(edge_list)

def get_dose(min_dose=10, max_dose=75, dose_interval=5, num_part=1):
    return np.random.choice(np.arange(min_dose, max_dose+1, dose_interval), size=num_part, replace=True)

def sample_circuit(promo_node, num_circuit, max_part, min_dose, max_dose, dose_interval, inhibitor):
    circuits = []
    if not inhibitor:
        for i in range(num_circuit):
            num_part = np.random.randint(1, max_part+1)
            part_list = np.random.choice(tf_list, num_part, replace=False)
            dose_list = dict(zip(part_list, get_dose(min_dose, max_dose, dose_interval, num_part)))
            edge_list = get_edges(promo_node, part_list)
            circuits.append(Topo(edge_list, dose_list, promo_node))
    else:
        for i in range(num_circuit):
            num_tf = np.random.randint(1, max_part)
            num_in = np.random.randint(1, max_part-num_tf+1)
            part_list = np.append(np.random.choice(tf_list, num_tf), np.random.choice(inhibitor_list, num_in))
            dose_list = dict(zip(part_list, get_dose(min_dose, max_dose, dose_interval, num_tf+num_in)))
            edge_list = get_edges(promo_node, part_list)
            circuits.append(Topo(edge_list, dose_list, promo_node))

    return circuits

def validate(g):
    # new_edges = set([k for k in g.graph.edges])
    circuit_tf_list = [k for k in g.part_list if k[0] == 'Z']
    if len(circuit_tf_list) == 0:
        raise Exception("Something's wrong. No TFs in the circuit.")

    if ('Rep' not in g.graph.nodes) | (len([k for k in g.graph.predecessors('Rep') if k[0] == 'Z']) == 0):
        g.graph.add_edges_from(get_in_path('Rep', None, circuit_tf_list))
        # g.graph.add_edges_from([(np.random.choice(circuit_tf_list), 'Rep')])

    for n in g.part_list:
        if n not in g.graph.nodes:
            g.graph.add_edges_from(get_in_path(n, g.promo_node, circuit_tf_list))
            g.graph.add_edges_from(get_out_path(n, g.part_list))
        else:
            viable_type = [k[-2][0] for k in nx.all_simple_paths(g.graph, g.promo_node, n)]
            if len(viable_type) == 0:
                g.graph.add_edges_from(get_in_path(n, g.promo_node, circuit_tf_list))
                # g.graph.add_edges_from([(g.promo_node, n)])
            else:
                if ('I' in viable_type) & ('Z' not in viable_type):
                    g.graph.add_edges_from(get_in_path(n, None, circuit_tf_list))

            if len(list(nx.all_simple_paths(g.graph, n, 'Rep'))) == 0:
                g.graph.add_edges_from(get_out_path(n, g.part_list))
            #     node_avail = set(g.graph.successors(n))
            #     node_avail.update([n])
            #     new_edges.update([np.random.choice(list(node_avail)), 'Rep'])
            #     # new_edges.update(get_out_path(n, g.part_list))

    g.update(list(g.graph.edges))
    # return list(new_edges)

def compare_circuit(g1, g2):
        ind = (set(g1.edge_list) == set(g2.edge_list)) & (g1.dose == g2.dose)

        return ind

def get_crosspt(list1, list2):
    same = set(list1).intersection(set(list2))
    same = list(same)
    if len(same) > 0:
        pt1 = np.random.choice(same)
        pt2 = pt1
    else:
        pt1 = np.random.choice(list1)
        if pt1[0] == 'Z':
            pt2 = np.random.choice([k for k in list2 if k[0] == 'Z'])
        else:
            pt2 = np.random.choice([k for k in list2 if k[0] == 'I'])

    return pt1, pt2

def switch_node(g, old_node, new_node):
    child_edge = []
    for edge in list(g.graph.edges):
        source, target = edge
        if source == old_node:
            source = new_node
        if target == old_node:
            target = new_node
        edge = (source, target)
        child_edge.append(tuple(edge))

    return child_edge

def crossover_node(g1, g2):
    pt1, pt2 = get_crosspt(g1.part_list, g2.part_list)

    child1_edge = switch_node(g1, pt1, pt2)
    child2_edge = switch_node(g2, pt2, pt1)
    child1_dose = {k: g1.dose[k] for k in g1.part_list if k != pt1}
    child1_dose.update({pt2: g2.dose[pt2]})
    child2_dose = {k: g2.dose[k] for k in g2.part_list if k != pt2}
    child2_dose.update({pt1: g1.dose[pt1]})

    child1 = Topo(child1_edge, child1_dose, g1.promo_node)
    child2 = Topo(child2_edge, child2_dose, g2.promo_node)

    return child1, child2

def match_node(new_node, part_list, promo_node, circuit_tf_list, circuit_in_list, pt2, node_list2):
     # = []
    for n in node_list2:
        if n == pt2:
            new_node.append(pt2)
        elif n in (part_list + [promo_node, 'Rep']):
            new_node.append(n)
        elif n[0] == 'Z':
            node_avail = set(circuit_tf_list).difference(set(new_node))
            if len(node_avail) > 0:
                n_new = np.random.choice(list(node_avail))
                new_node.append(n_new)
            # else:
            #     n_new = promo_node
            #     new_node.append(n_new)
        elif n[0] == 'I':
            node_avail = set(circuit_in_list).difference(set(new_node))
            if len(node_avail) > 0:
                n_new = np.random.choice(list(node_avail))
                new_node.append(n_new)

    # return

def switch_edge(g1, pt1, pt2, in_list2, out_list2, dose2):
    child = deepcopy(g1)

    child.part_list.remove(pt1)
    child.dose.pop(pt1)
    child.graph.remove_node(pt1)

    circuit_tf_list = [k for k in child.part_list if k[0] == 'Z']
    circuit_in_list = [k for k in child.part_list if k[0] == 'I']

    in_node = []
    common_list2 = [k for k in in_list2 if k in out_list2]
    match_node(in_node, child.part_list, child.promo_node, circuit_tf_list, circuit_in_list, pt2, common_list2)
    out_node = [k for k in in_node if k[0] != 'P']

    in_list2 = [k for k in in_list2 if k not in common_list2]
    match_node(in_node, child.part_list, child.promo_node, circuit_tf_list, circuit_in_list, pt2, in_list2)
    out_list2 = [k for k in out_list2 if k not in common_list2]
    match_node(out_node, child.part_list, child.promo_node, circuit_tf_list, circuit_in_list, pt2, out_list2)

    new_edges = set([(k, pt2) for k in in_node])
    new_edges.update([(pt2, k) for k in out_node])

    child.part_list.append(pt2)
    child.dose.update({pt2: dose2})
    child.graph.add_edges_from(new_edges)

    # new_edges = validate(child)
    # child.update(new_edges)
    validate(child)

    return child

def crossover_structure(g1, g2):
    pt1, pt2 = get_crosspt(g1.part_list, g2.part_list)

    child1 = switch_edge(g1, pt1, pt2, list(g2.graph.predecessors(pt2)), list(g2.graph.successors(pt2)), g2.dose[pt2])
    child2 = switch_edge(g2, pt2, pt1, list(g1.graph.predecessors(pt1)), list(g1.graph.successors(pt1)), g1.dose[pt1])

    # child1.check_valid()
    # child2.check_valid()
    # if child1.valid == 0 | child2.valid == 0:
    #     print(pt1)
    #     print(pt2)
    #     print(g1.edge_list)
    #     print(g2.edge_list)
    #     print(child1.edge_list)
    #     print(child2.edge_list)

    return child1, child2

def mutate_dose(g, min_dose=10, max_dose=75, dose_interval=5):
    n = np.random.choice(g.part_list)
    g.dose.update({n: get_dose(min_dose, max_dose, dose_interval, 1)[0]})

def mutate_node_type(g, min_dose=10, max_dose=75, dose_interval=5):
    old_node = np.random.choice(g.part_list)
    if old_node[0] == 'Z':
        node_avail = list(set(tf_list).difference(set(g.part_list)))
    else:
        node_avail = list(set(inhibitor_list).difference(set(g.part_list)))
    new_node = np.random.choice(node_avail)
    g.dose.pop(old_node)
    g.dose.update({new_node: get_dose(min_dose, max_dose, dose_interval, 1)[0]})
    new_edges = switch_node(g, old_node, new_node)
    g.graph.remove_node(old_node)
    g.update(new_edges)

def add_node(g, circuit_tf_list, min_dose=10, max_dose=75, dose_interval=5, inhibitor=False):
    if not inhibitor:
        node_avail = [k for k in tf_list if k not in circuit_tf_list]
    else:
        node_avail = [k for k in parts.keys() if k not in g.part_list]
    new_node = np.random.choice(node_avail)
    g.dose.update({new_node: get_dose(min_dose, max_dose, dose_interval, 1)[0]})
    new_edges = set([k for k in g.edge_list])
    new_edges.update(get_in_path(new_node, g.promo_node, circuit_tf_list))
    new_edges.update(get_out_path(new_node, g.part_list))
    g.update(list(new_edges))

def remove_node(g, circuit_tf_list):
    circuit_in_list = [k for k in g.part_list if k[0] == 'I']

    if len(circuit_tf_list) > 1 | len(circuit_in_list) > 1:
        if len(circuit_in_list) <= 1:
            old_node = np.random.choice(circuit_tf_list)
        elif len(circuit_tf_list) <= 1:
            old_node = np.random.choice(circuit_in_list)
        else:
            old_node = np.random.choice(g.part_list)

        g.part_list.remove(old_node)
        g.dose.pop(old_node)
        g.graph.remove_node(old_node)
        # new_edges = validate(g)
        # g.update(new_edges)
        validate(g)

def mutate_node_num(g, max_part, min_dose=10, max_dose=75, dose_interval=5, inhibitor=False):
    circuit_tf_list = [k for k in g.part_list if k[0] == 'Z']
    if max_part > 1:
        if len(g.part_list) == 1:
            add_node(g, circuit_tf_list, min_dose, max_dose, dose_interval, inhibitor)
        elif len(g.part_list) < max_part:
            if np.random.uniform() < 0.5:
                add_node(g, circuit_tf_list, min_dose, max_dose, dose_interval, inhibitor)
            else:
                remove_node(g, circuit_tf_list)
        else:
            remove_node(g, circuit_tf_list)

def get_full_connected(part_list, promo_node):
    edge_list = [(promo_node, k) for k in part_list]
    edge_list.extend([(k, 'Rep') for k in part_list])
    edge_list.extend(list(combinations_with_replacement(part_list, 2)))
    return edge_list

def mutate_edge(g):
    edge_full = get_full_connected(g.part_list, g.promo_node)
    if g.graph.size() < len(edge_full):
        edge_avail = [k for k in edge_full if k not in g.graph.edges]
        ind = np.random.choice(len(edge_avail))
        g.edge_list.append(edge_avail[ind])
        g.update(g.edge_list)
    else:
        pass

a=1