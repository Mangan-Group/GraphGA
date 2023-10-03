import networkx as nx
from define_circuit import *
from copy import deepcopy, copy

from scipy.integrate import odeint
from load_files import *
from get_system_equations import *

from itertools import combinations, permutations, product


# Randomly generate a path from a part n to the reporter
def get_out_path(n, part_list):
    # Build a list of potential parts and declare path
    out_list = [k for k in part_list if k != n]
    out_path = [n]

    # Add a random circuit part to the path
    if len(out_list) > 0:
        num_connect = np.random.randint(len(out_list))
        out_path.extend(np.random.choice(out_list, num_connect, replace=False))
    out_path.append('Rep')

    # Create edges using path
    edges = [(i, j) for i, j in zip(out_path[:-1], out_path[1:])]

    return edges


# Randomly generate a parth from the promoter to a part n
def get_in_path(n, promo_node, circuit_tf_list):

    # Randomly choose a tf or the promoter to direct toward n
    if promo_node is not None:
        in_node = np.random.choice(circuit_tf_list + [promo_node])
        edges = [(in_node, n)]

        # Add promoter to the start of the path and choose a random set of parts for the path
        if in_node != promo_node:
            in_path = [promo_node]
            num_connect = np.random.randint(len(circuit_tf_list))
            in_path.extend(np.random.choice(circuit_tf_list, num_connect, replace=False))
            in_path.append(n)
            edges.extend([(i, j) for i, j in zip(in_path[:-1], in_path[1:])])

    # Choose a tf to direct toward n
    else:
        in_node = np.random.choice(circuit_tf_list)
        edges = [(in_node, n)]

        # Choose another tf to direct toward n if added edge is a self-regulation
        if in_node == n:
            in_node = np.random.choice([k for k in circuit_tf_list if k != n])
            edges.append((in_node, n))

    return edges



def get_edges(promo_node, part_list):
    circuit_tf_list = []
    same_list = []
    for k in part_list:
        if k[0] == 'Z':
            circuit_tf_list.append(k)
            if ('I' + k[1:]) in part_list:
                same_list.append([k, ('I' + k[1:])])

    # Declare edge_list as a list instead of a set
    edge_list = [(promo_node, np.random.choice(circuit_tf_list)), (np.random.choice(circuit_tf_list),'Rep')]

    for n in part_list:
        if not any(n in sublist for sublist in same_list):
            in_edges = get_in_path(n, promo_node, circuit_tf_list)

            # Add from in_edges to edge_list if not already in edge_list
            for edge in in_edges:
                if edge not in edge_list:
                    edge_list.append(edge)

            out_edges = get_out_path(n, part_list)

            # Add from out_edges to edge_list if not already in edge_list
            for edge in out_edges:
                if edge not in edge_list:
                    edge_list.append(edge)

    for z, i in same_list:
        in_edges_z = get_in_path(z, promo_node, circuit_tf_list)

        # Add from in_edges_z to edge_list if not already in edge_list
        for edge in in_edges_z:
            if edge not in edge_list:
                edge_list.append(edge)

        in_edges_i = get_in_path(i, promo_node, circuit_tf_list)

        # Add from in_edges_i to edge_list if not already in edge_list
        for edge in in_edges_i:
            if edge not in edge_list:
                edge_list.append(edge)

        out_edges_z = get_out_path(z, part_list)

        # Add from out_edges_z to edge_list if not already in edge_list
        for edge in out_edges_z:
            if edge not in edge_list:
                edge_list.append(edge)

        all_out_edges_z = [edge for edge in (in_edges_z + in_edges_i) if edge[0] == z] + out_edges_z

        # Add from all_out_edges_z to edge_list if not already in edge_list
        for edge in all_out_edges_z:
            if edge not in edge_list:
                edge_list.append(edge)

    return edge_list


# Randomly generate a dose amount within the given constraints
def get_dose(min_dose=10, max_dose=75, dose_interval=5, num_part=1):
    return np.random.choice(np.arange(min_dose, max_dose + 1, dose_interval), size=num_part, replace=True)


# def sample_circuit(promo_node, num_circuit, max_part, min_dose, max_dose, dose_interval, inhibitor):
#     circuits = []
#     if not inhibitor:
#         for i in range(num_circuit):
#             num_part = np.random.randint(1, max_part + 1)
#             part_list = np.random.choice(tf_list, num_part, replace=False)
#             dose_list = dict(zip(part_list, get_dose(min_dose, max_dose, dose_interval, num_part)))
#             edge_list = get_edges(promo_node, part_list)
#             circuits.append(Topo(edge_list, dose_list, promo_node))
#     else:
#         for i in range(num_circuit):
#             num_tf = np.random.randint(1, max_part)
#             num_in = np.random.randint(1, max_part - num_tf + 1)
#             part_list = np.append(np.random.choice(tf_list, num_tf), np.random.choice(inhibitor_list, num_in))
#             dose_list = dict(zip(part_list, get_dose(min_dose, max_dose, dose_interval, num_tf + num_in)))
#             edge_list = get_edges(promo_node, part_list)
#             circuits.append(Topo(edge_list, dose_list, promo_node))

#     return circuits


def sampling(promo_node, num_dict, min_dose, max_dose, dose_interval, inhibitor=False):
    combo = []
    for num_part, num_circuit in num_dict.items():
        if not inhibitor:
            part_combo = list(combinations(tf_list, num_part))
            ind = np.random.choice(len(part_combo), num_circuit)
            combo.extend([part_combo[i] for i in ind])
        else:

            for i in range(num_circuit):
                # Select a transcription factor because circuit must have at least 1
                guaranteed_tf_index = np.random.choice(len(tf_list), 1)
                guaranteed_tf = tf_list[guaranteed_tf_index[0]]

                # Make a list of remaining transcription factors combined with inhibitors
                remaining_tfs = np.delete(tf_list, guaranteed_tf_index[0])
                remaining_options = np.concatenate((remaining_tfs, inhibitor_list))

                # Select a random subset of the remaining parts to add to the guaranteed tf
                choices = np.random.choice(remaining_options, size=num_part - 1)

                # Append set of choices to the list of combinations
                combo.extend([np.append(choices, guaranteed_tf)])

    circuits = []

    for i in range(len(combo)):
        part_list = combo[i]
        # full_edge_list = get_full_connected(part_list, promo_node)
        edge_list = get_edges(promo_node, list(part_list))
        dose_list = dict(zip(part_list, get_dose(min_dose, max_dose, dose_interval, len(part_list))))
        circuits.append([Topo(edge_list, dose_list, promo_node)])

    return np.asarray(circuits)


# Make the inputted graph valid
def validate(g):
    circuit_tf_list = []
    same_list = []

    # Fill circuit_tf_list with tfs and same_list with their paired inhibitor if it exists
    for k in g.part_list:
        if k[0] == 'Z':
            circuit_tf_list.append(k)
            if ('I' + k[1:]) in g.part_list:
                same_list.append([k, ('I' + k[1:])])

    # Exception if no transcription factors
    if len(circuit_tf_list) == 0:
        raise Exception("Something's wrong. No TFs in the circuit.")

    # Must have a reporter and at least one ZF activator directed toward the reporter
    if ('Rep' not in g.graph.nodes) | (len([k for k in g.graph.predecessors('Rep') if k[0] == 'Z']) == 0):
        g.graph.add_edges_from(get_in_path('Rep', None, circuit_tf_list))

    for n in g.part_list:
        # Ensure all parts are in the circuit graph
        if n not in g.graph.nodes:
            g.graph.add_edges_from(get_in_path(n, g.promo_node, circuit_tf_list))
            g.graph.add_edges_from(get_out_path(n, g.part_list))

        else:
            # Check for paths from promoter to the part of the circuit
            viable_type = [k[-2][0] for k in nx.all_simple_paths(g.graph, g.promo_node, n)]

            # Add edge from promoter to part
            if len(viable_type) == 0:
                g.graph.add_edges_from(get_in_path(n, g.promo_node, circuit_tf_list))
            # Add transcription factor for each inhibitor because they must have the same path (IDK)
            else:
                if ('I' in viable_type) and (('Z' not in viable_type) and ('P' not in viable_type)):
                    # g.graph.add_edges_from(get_in_path(n, None, circuit_tf_list))
                    g.graph.add_edges_from(get_in_path(n, g.promo_node, circuit_tf_list))

            # Add path from every part to the reporter
            if len(list(nx.all_simple_paths(g.graph, n, 'Rep'))) == 0:
                g.graph.add_edges_from(get_out_path(n, g.part_list))

    # All edges in list must be in the circuit graph
    if set(g.graph.edges) != set(g.edge_list):
        g.update(list(g.graph.edges))

    # ZFs and their inhibitors must have the same successors, so add edges if not true
    for z, i in same_list:
        z_succ = set(g.graph.successors(z))
        i_succ = set(g.graph.successors(i))
        if z_succ != i_succ:
            z_out = list(g.graph.out_edges(z))
            i_out = list(g.graph.out_edges(i))
            i_out_new = [(i, k[1]) for k in z_out]
            g.graph.remove_edges_from(i_out)
            g.graph.add_edges_from(i_out_new)

    # Check that the edge list and graph edges are still the same
    if set(g.graph.edges) != set(g.edge_list):
        g.update(list(g.graph.edges))


# Check if the inputted graph is valid
def check_valid(g, promo_node, part_list):
    # Graph must have a promoter and reporter
    if (promo_node not in g.nodes) or ('Rep' not in g.nodes):
        return 0

    # Fill circuit_tf_list with tfs and same_list with their paired inhibitor if it exists
    circuit_tf_list = []
    same_list = []
    for k in part_list:
        if k[0] == 'Z':
            circuit_tf_list.append(k)
            if ('I' + k[1:]) in part_list:
                same_list.append([k, ('I' + k[1:])])

    # Graph must have a transcription factor
    if len(circuit_tf_list) == 0:
        return 0

    # Graph must have a reporter with a transcription factor directed towards it
    if ('Rep' not in g.nodes) | (len([k for k in g.predecessors('Rep') if k[0] == 'Z']) == 0):
        return 0

    # All parts must be nodes in the circuit graph
    for n in part_list:
        if n not in g.nodes:
            return 0
        else:
            viable_type = [k[-2][0] for k in nx.all_simple_paths(g, promo_node, n)]
            if len(viable_type) == 0:
                return 0
            else:
                if ('I' in viable_type) and (('Z' not in viable_type) and ('P' not in viable_type)):
                    return 0
            if len(list(nx.all_simple_paths(g, n, 'Rep'))) == 0:
                return 0

    for z, i in same_list:
        z_succ = set(g.successors(z))
        i_succ = set(g.successors(i))
        if z_succ != i_succ:
            return 0
    return 1


# def check_valid(g, num_parts):
#     graph_parts = [i[0] for i in g.nodes]
#     if ('P' not in graph_parts) or ('R' not in graph_parts) or ('Z' not in graph_parts):
#         return 0
#     elif (len(graph_parts) - 2) < num_parts:
#         return 0
#
#     for n in g.nodes:
#         if n[0] == 'P':
#             out_types = set([i[0] for i in list(g.successors(n))])
#             if 'Z' not in out_types:
#                 return 0
#         elif n[0] == 'R':
#             in_types = set([i[0] for i in list(g.predecessors(n))])
#             if (not in_types) or (('I' in in_types) and ('Z' not in in_types)):
#                 return 0
#         else:
#             if (n[0] == 'Z') and ('I' + n[1:] in g.nodes):
#                 if list(g.successors(n)) != list(g.successors('I' + n[1:])):
#                     return 0
#             in_nodes = list(g.predecessors(n))
#             if not in_nodes:
#                 return 0
#             elif in_nodes == [n]:
#                 return 0
#             else:
#                 in_types = set([i[0] for i in in_nodes])
#                 if ('I' in in_types) and ('Z' not in in_types):
#                     return 0
#             if len(list(nx.all_simple_paths(g, n, 'Rep'))) == 0:
#                 return 0
#     return 1


def compare_circuit(g1, g2):
    ind = (set(g1.edge_list) == set(g2.edge_list)) & (g1.dose == g2.dose)

    return ind


def get_crosspt(list1, list2):

    # Create a list of shared circuit parts by iterating through both lists
    same = []
    for item in list1:
        if item in list2:
            same.append(item)

    if len(same) > 0:
        pt1 = np.random.choice(same)
        pt2 = pt1
    else:
        pt1 = np.random.choice(list1)

        if pt1[0] == 'Z':
            pt2 = np.random.choice([k for k in list2 if k[0] == 'Z'])
        else:
            list2_inhibitors = [k for k in list2 if k[0] == 'I']

            # Follow the same process if both circuits have inhibitors
            if len(list2_inhibitors) > 0:
                pt2 = np.random.choice(list2_inhibitors)

            # New method for if the second circuit doesn't have an inhibitor
            else:

                # Switch inhibitor for tf if second circuit has more than 1 part
                if len(list2) > 1:
                    pt2 = np.random.choice(list2)

                # Have to switch tf for tf if second circuit only has 1 part
                else:
                    pt1 = np.random.choice([k for k in list1 if k[0] == 'Z'])
                    pt2 = np.random.choice([k for k in list2 if k[0] == 'Z'])

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
    for n in node_list2:
        if n == pt2:
            new_node.append(pt2)
        elif n in (part_list + [promo_node, 'Rep']):
            new_node.append(n)
        elif n[0] == 'Z':
            # Determine available nodes by iterating through list of transcription factors
            node_avail = []
            for item in circuit_tf_list:
                if item not in new_node:
                    node_avail.append(item)

            if len(node_avail) > 0:
                n_new = np.random.choice(list(node_avail))
                new_node.append(n_new)
        elif n[0] == 'I':
            # Determine available nodes by iterating through list of transcription factors
            node_avail = []
            for item in circuit_in_list:
                if item not in new_node:
                    node_avail.append(item)

            if len(node_avail) > 0:
                n_new = np.random.choice(list(node_avail))
                new_node.append(n_new)


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

    # Build an array of new edges by appending all unique edges from in_node and out_node
    new_edges = []
    for k in in_node:
        new_edges.append((k, pt2))
    for k in out_node:
        out_edge = (pt2,k)
        if out_edge not in new_edges:
            new_edges.append(out_edge)

    child.part_list.append(pt2)
    child.dose.update({pt2: dose2})
    child.graph.add_edges_from(new_edges)

    validate(child)

    return child


def crossover_structure(g1, g2):
    pt1, pt2 = get_crosspt(g1.part_list, g2.part_list)

    child1 = switch_edge(g1, pt1, pt2, list(g2.graph.predecessors(pt2)), list(g2.graph.successors(pt2)), g2.dose[pt2])
    child2 = switch_edge(g2, pt2, pt1, list(g1.graph.predecessors(pt1)), list(g1.graph.successors(pt1)), g1.dose[pt1])

    return child1, child2


def mutate_dose(g, min_dose=10, max_dose=75, dose_interval=5):
    n = np.random.choice(g.part_list)
    g.dose.update({n: get_dose(min_dose, max_dose, dose_interval, 1)[0]})


def mutate_node_type(g, min_dose=10, max_dose=75, dose_interval=5):
    old_node = np.random.choice(g.part_list)
    if old_node[0] == 'Z':
        same_type = 'I'
        node_avail = []
        for item in tf_list:
            if item not in g.part_list:
                node_avail.append(item)
    else:
        same_type = 'Z'
        node_avail = []
        for item in tf_list:
            if item not in g.part_list:
                node_avail.append(item)

    new_node = np.random.choice(node_avail)
    same = same_type + new_node[1:]
    g.dose.update({new_node: g.dose[old_node]})
    g.dose.pop(old_node)
    new_edges = switch_node(g, old_node, new_node)
    g.graph.remove_node(old_node)
    g.update(new_edges)
    if same in g.part_list:
        same_out = list(g.graph.out_edges(same))
        new_node_out = list(g.graph.out_edges(new_node))
        if same_type == 'Z':
            new_node_out_new = [(new_node, k[1]) for k in same_out]
            g.graph.remove_edges_from(new_node_out)
            g.graph.add_edges_from(new_node_out_new)
        else:
            same_out_new = [(same, k[1]) for k in new_node_out]
            g.graph.remove_edges_from(same_out)
            g.graph.add_edges_from(same_out_new)
        g.update(list(g.graph.edges))


def add_node(g, circuit_tf_list, min_dose=10, max_dose=75, dose_interval=5, inhibitor=False):
    if not inhibitor:
        node_avail = [k for k in tf_list if k not in circuit_tf_list]
    else:
        node_avail = [k for k in parts.keys() if k not in g.part_list]
    new_node = np.random.choice(node_avail)
    new_node_type = new_node[0]
    g.dose.update({new_node: get_dose(min_dose, max_dose, dose_interval, 1)[0]})

    # Create list of new edges from the graph edge list
    new_edges = []
    for k in g.edge_list:
        new_edges.append(k)

    # Add new edges from new_node_in to new_edges by iterating through the list
    new_node_in = get_in_path(new_node, g.promo_node, circuit_tf_list)
    for edge in new_node_in:
        if edge not in new_edges:
            new_edges.append(edge)

    if new_node_type == 'I':
        same_type = 'Z'
        same = same_type + new_node[1:]
        if same in g.part_list:
            same_out = list(g.graph.out_edges(same))
            all_out_edges_same = [edge for edge in new_node_in if edge[0] == same] + same_out

            # Add new edges from new_node_out to new_edges by iterating through the list
            for k in all_out_edges_same:
                possible_edge = (new_node, k[1])
                if possible_edge not in new_edges:
                    new_edges.append(possible_edge)
        else:
            # Add new edges from new_node_out to new_edges by iterating through the list
            new_node_out = get_out_path(new_node, g.part_list)
            for edge in new_node_out:
                if edge not in new_edges:
                    new_edges.append(edge)
    else:
        same_type = 'I'
        same = same_type + new_node[1:]
        if same in g.part_list:
            same_out = list(g.graph.out_edges(same))

            # Remove shared edges between new_edges and same_out from the new_edges array
            for i in range(len(new_edges)):
                if new_edges[i] in same_out:
                    new_edges.pop(i)

            new_node_out = get_out_path(new_node, g.part_list)

            # Add new edges from new_node_out to new_edges by iterating through the list
            for edge in new_node_out:
                if edge not in new_edges:
                    new_edges.append(edge)

            all_out_edges_new = [edge for edge in new_node_in if edge[0] == new_node] + new_node_out

            # Add new edges from the out edges to new_edges by iterating through the list
            for k in all_out_edges_new:
                possible_edge = (same, k[1])
                if possible_edge not in new_edges:
                    new_edges.append(possible_edge)
        else:
            # Add new edges from the out path to new_edges by iterating through the list
            new_node_out = get_out_path(new_node, g.part_list)
            for edge in new_node_out:
                if edge not in new_edges:
                    new_edges.append(edge)

    g.update(new_edges)


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
    edge_list.extend([(k, k) for k in part_list])
    edge_list.extend(list(permutations(part_list, 2)))
    return edge_list


def mutate_edge(g, inhibitor=False):
    if len(g.part_list) == 1:
        if g.graph.size() == 3:
            g.edge_list.remove((g.part_list[0], g.part_list[0]))
        else:
            g.edge_list.append((g.part_list[0], g.part_list[0]))
        g.update(g.edge_list)
    else:
        edge_full = get_full_connected(g.part_list, g.promo_node)
        if g.graph.size() == len(edge_full):
            if inhibitor:
                num_parts = len(g.part_list)

                # edge_avail = [k for k in g.edge_list]
                valid = 0
                while (valid == 0) and (len(edge_full) > 0):
                    g_graph = deepcopy(g.graph)
                    ind = np.random.choice(len(edge_full))
                    edge_removed = edge_full[ind]
                    g_graph.remove_edge(edge_removed[0], edge_removed[1])
                    valid = check_valid(g_graph, g.promo_node, g.part_list)
                    edge_full.remove(edge_removed)
                if valid == 1:
                    g.update(list(g_graph.edges))
            else:
                ind = np.random.choice(g.graph.size())
                g.edge_list.remove(edge_full[ind])
                g.update(g.edge_list)

        else:
            if np.random.uniform() < 0.5:
                edge_avail = [k for k in edge_full if k not in g.graph.edges]
                ind = np.random.choice(len(edge_avail))
                g.edge_list.append(edge_avail[ind])
                g.update(g.edge_list)
            else:
                num_parts = len(g.part_list)

                edge_avail = [k for k in g.edge_list]
                valid = 0
                while (valid == 0) and (len(edge_avail) > 0):
                    g_graph = deepcopy(g.graph)
                    ind = np.random.choice(len(edge_avail))
                    edge_removed = edge_avail[ind]
                    g_graph.remove_edge(edge_removed[0], edge_removed[1])
                    valid = check_valid(g_graph, g.promo_node, g.part_list)
                    edge_avail.remove(edge_removed)
                if valid == 1:
                    g.update(list(g_graph.edges))


def is_equal(x1, x2):
    return compare_circuit(x1[0], x2[0])


def tournament(X, obj, pressure=3):
    throuple = np.random.choice(range(len(X)), size=pressure, replace=False)
    return throuple[obj[throuple].argsort()][1:]


def crossover(X, obj, rank_dict=None, **kwargs):
    n_matings = int(len(X) / 2)
    Y = np.full_like(X, None, dtype=object)
    for k in range(n_matings):
        throuple = np.random.choice(range(len(X)), 3, replace=False)
        if rank_dict is None:
            parents = throuple[obj[throuple].argsort()][:-1]
        else:
            throuple_rank = np.asarray([rank_dict[i]['rank'] for i in throuple])
            parents = throuple[throuple_rank.argsort()][:-1]
        # if strategy == "tournament":
        #     parents = tournament(X, obj)
        # elif strategy == "random":
        #     parents = np.random.choice(len(X), size=2, replace=False)
        Y[2 * k, 0], Y[2 * k + 1, 0] = crossover_structure(X[parents[0], 0], X[parents[1], 0])
    return Y


def mutate(problem, X, prob, dose=False, **kwargs):
    for i in range(len(X)):
        if np.random.uniform(0, 1) < prob:
            if dose:
                r = np.random.choice(4)
            else:
                r = np.random.choice(3)

            if r == 0:
                mutate_node_num(X[i, 0], problem.max_part, problem.min_dose, problem.max_dose, problem.dose_interval,
                                problem.inhibitor)
            elif r == 1:
                mutate_node_type(X[i, 0], problem.min_dose, problem.max_dose, problem.dose_interval)
            elif r == 2:
                mutate_edge(X[i, 0], problem.inhibitor)
            else:
                mutate_dose(X[i, 0], problem.min_dose, problem.max_dose, problem.dose_interval)


            # r = np.random.uniform()
            #
            # if r < 0.52:
            #     mutate_node_num(X[i, 0], problem.max_part, problem.min_dose, problem.max_dose, problem.dose_interval,
            #                             problem.inhibitor)
            # elif r < 0.68:
            #     mutate_node_type(X[i, 0], problem.min_dose, problem.max_dose, problem.dose_interval)
            # elif r < 0.74:
            #     mutate_edge(X[i, 0], problem.inhibitor)
            # else:
            #     mutate_dose(X[i, 0], problem.min_dose, problem.max_dose, problem.dose_interval)

            # r = np.random.uniform()
            #
            # if r < 0.4:
            #     mutate_node_num(X[i, 0], problem.max_part, problem.min_dose, problem.max_dose, problem.dose_interval,
            #                             problem.inhibitor)
            # elif r < 0.6:
            #     mutate_node_type(X[i, 0], problem.min_dose, problem.max_dose, problem.dose_interval)
            # elif r < 0.8:
            #     mutate_edge(X[i, 0], problem.inhibitor)
            # else:
            #     mutate_dose(X[i, 0], problem.min_dose, problem.max_dose, problem.dose_interval)