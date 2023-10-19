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
    # if no other parts in circuit, add reporter to path
    if len(out_list) > 0:
        # random selection of number of connections from
        # number of possible parts to connect to (1 fewer
        # than total length)
        num_connect = np.random.randint(len(out_list))
        # add random parts to out_path of length num_connect
        out_path.extend(np.random.choice(out_list, num_connect, replace=False))
    # add reporter
    out_path.append('Rep')

    # Create edges using out_path
    edges = [(i, j) for i, j in zip(out_path[:-1], out_path[1:])]

    return edges


# Randomly generate a path from the promoter to a part n
def get_in_path(n, promo_node, circuit_tf_list):

    # Randomly choose a tf or the promoter to direct toward n
    if promo_node is not None:
        in_node = np.random.choice(circuit_tf_list + [promo_node])
        # add to edge list
        edges = [(in_node, n)]

        # Add promoter to the start of the path and choose a
        # random set of parts for the path
        if in_node != promo_node:
            in_path = [promo_node]
            # random selection of number of connections from
            # number of possible parts be connected to (1 fewer
            # than total length)
            num_connect = np.random.randint(len(circuit_tf_list))
            # add random parts to in_path of length num_connect
            # not using inhibitors because direct connection must 
            # be a tf or promoter
            in_path.extend(np.random.choice(circuit_tf_list, num_connect, replace=False))
            in_path.append(n)
            # create edges using in_path
            edges.extend([(i, j) for i, j in zip(in_path[:-1], in_path[1:])])

    # Choose a tf to direct toward n if no promoter
    else:
        in_node = np.random.choice(circuit_tf_list)
        edges = [(in_node, n)]

        # Choose another tf to direct toward n if added edge is 
        # a self-regulation
        if in_node == n:
            in_node = np.random.choice([k for k in circuit_tf_list if k != n])
            edges.append((in_node, n))

    return edges



def get_edges(promo_node, part_list):
    circuit_tf_list = []
    same_list = []
    # add tf and inhibitor to same list if same type
    for k in part_list:
        if k[0] == 'Z':
            circuit_tf_list.append(k)
            if ('I' + k[1:]) in part_list:
                same_list.append([k, ('I' + k[1:])])

    # Declare edge_list as a list instead of a set
    # so order is the same each time if part_list
    # is the same
    # add path P1 -> random tf -> reporter so it is
    # valid circuit
    edge_list = [(promo_node, np.random.choice(circuit_tf_list)),
                 (np.random.choice(circuit_tf_list),'Rep')
    ]

    for n in part_list:
        # for parts that are not in same_list (edges 
        # added later)
        if not any(n in sublist for sublist in same_list):
            in_edges = get_in_path(n, promo_node, circuit_tf_list)

            # Add from in_edges to edge_list if not
            # already in edge_list
            for edge in in_edges:
                if edge not in edge_list:
                    edge_list.append(edge)

            out_edges = get_out_path(n, part_list)

            # Add from out_edges to edge_list if not
            # already in edge_list
            for edge in out_edges:
                if edge not in edge_list:
                    edge_list.append(edge)

    for z, i in same_list:
        in_edges_z = get_in_path(z, promo_node, circuit_tf_list)

        # Add from in_edges_z to edge_list if not 
        # already in edge_list
        for edge in in_edges_z:
            if edge not in edge_list:
                edge_list.append(edge)

        in_edges_i = get_in_path(i, promo_node, circuit_tf_list)

        # Add from in_edges_i to edge_list if not 
        # already in edge_list
        for edge in in_edges_i:
            if edge not in edge_list:
                edge_list.append(edge)

        out_edges_z = get_out_path(z, part_list)

        # Add from out_edges_z to edge_list if not 
        # already in edge_list
        for edge in out_edges_z:
            if edge not in edge_list:
                edge_list.append(edge)

        # make list of all edges that will need to
        # be regulated by same type inhibitor:
        # self regulated tf edges, tf regulated
        # inhibitor (not in out_edges_z), out edges 
        # for tf
        all_out_edges_z = [edge for edge in (in_edges_z + in_edges_i) 
                           if edge[0] == z] + out_edges_z
        # update from 
        # edge_list.update([(i, k[1]) for k in all_out_edges_z])
        # change edges to inhibitor regulation
        inhibitor_edges = [(i, k[1]) for k in all_out_edges_z]

        # Add from inhibitor_edges to edge_list if not 
        # already in edge_list
        for edge in inhibitor_edges:
            if edge not in edge_list:
                edge_list.append(edge)

    return edge_list


# Randomly generate a dose amount within the given constraints
def get_dose(min_dose=10, max_dose=75, dose_interval=5, num_part=1):
    return (np.random.choice(np.arange(min_dose, max_dose + 1, dose_interval),
                             size=num_part, replace=True))


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

# generate population of circuits with specifications 
# as args
def sampling(promo_node, num_dict, min_dose, max_dose, dose_interval, inhibitor=False):
    combo = []
    # for each specified number of circuits with
    # specified number of parts
    for num_part, num_circuit in num_dict.items():
        if not inhibitor:
            # all combinations of parts for all tf
            # for specified number of parts
            part_combo = list(combinations(tf_list, num_part))
            # randomly choose from combinations num_circuit
            # indices
            ind = np.random.choice(len(part_combo), num_circuit)
            combo.extend([part_combo[i] for i in ind])
        else:

            for i in range(num_circuit):
                # Select a transcription factor because circuit 
                # must have at least 1
                guaranteed_tf_index = np.random.choice(len(tf_list), 1)
                guaranteed_tf = tf_list[guaranteed_tf_index[0]]

                # Make a list of remaining transcription factors 
                # combined with inhibitors
                remaining_tfs = np.delete(tf_list, guaranteed_tf_index[0])
                remaining_options = np.concatenate((remaining_tfs, inhibitor_list))

                # Select a random subset of the remaining parts to 
                # add to the guaranteed tf
                choices = np.random.choice(remaining_options, size=num_part - 1)

                # Append set of choices and guaranteed tf to the 
                # list of combinations
                combo.extend([np.append(choices, guaranteed_tf)])

    circuits = []

    for i in range(len(combo)):
        part_list = combo[i]
        # full_edge_list = get_full_connected(part_list, promo_node)
        edge_list = get_edges(promo_node, list(part_list))
        # get dose list in dict format
        dose_list = dict(zip(part_list, 
                             get_dose(min_dose, max_dose, dose_interval, len(part_list))))
        # append to circuits class instance based on specs
        circuits.append([Topo(edge_list, dose_list, promo_node)])

    return np.asarray(circuits)


# Make the inputted graph valid
def validate(g):
    circuit_tf_list = []
    same_list = []

    # Fill circuit_tf_list with tfs and same_list 
    # with their paired inhibitor if it exists
    for k in g.part_list:
        if k[0] == 'Z':
            circuit_tf_list.append(k)
            if ('I' + k[1:]) in g.part_list:
                same_list.append([k, ('I' + k[1:])])

    # Exception if no transcription factors
    if len(circuit_tf_list) == 0:
        raise Exception("Something's wrong. No TFs in the circuit.")

    # Must have a reporter and at least one ZF activator 
    # directed toward the reporter
    if (('Rep' not in g.graph.nodes) | 
        (len([k for k in g.graph.predecessors('Rep') if k[0] == 'Z']) == 0)):
        g.graph.add_edges_from(get_in_path('Rep', None, circuit_tf_list))

    for n in g.part_list:
        # Ensure all parts are in the circuit graph
        if n not in g.graph.nodes:
            g.graph.add_edges_from(get_in_path(n, g.promo_node, circuit_tf_list))
            g.graph.add_edges_from(get_out_path(n, g.part_list))

        else:
            # Check for paths from promoter to the part of the circuit
            # to see what direct connections to part are 
            viable_type = [k[-2][0] 
                           for k in nx.all_simple_paths(g.graph, g.promo_node, n)]

            # Add in_path edges that will include tf or promoter
            # if none exist
            if len(viable_type) == 0:
                g.graph.add_edges_from(get_in_path(n, g.promo_node, circuit_tf_list))
            # If direct connnections to part do not contain tf or promoter
            # add in_path edges that will include tf or promoter
            else:
                if (('I' in viable_type) and 
                    (('Z' not in viable_type) and ('P' not in viable_type))):
                    # g.graph.add_edges_from(get_in_path(n, None, circuit_tf_list))
                    g.graph.add_edges_from(get_in_path(n, g.promo_node, circuit_tf_list))

            # Add path from every part to the reporter
            if len(list(nx.all_simple_paths(g.graph, n, 'Rep'))) == 0:
                g.graph.add_edges_from(get_out_path(n, g.part_list))

    # All edges in list must be in the circuit graph
    # after making above updates to graph only
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


# Check if the inputted graph is valid (doesn't 
# update if not valid)
# based on validate function above so operations
# are mostly the same
def check_valid(g, promo_node, part_list):
    # Graph must have a promoter and reporter
    if (promo_node not in g.nodes) or ('Rep' not in g.nodes):
        return 0

    # Fill circuit_tf_list with tfs and same_list with their 
    # paired inhibitor if it exists
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

    # Graph must have a reporter with a transcription
    # factor directed towards it
    if (('Rep' not in g.nodes) |
        (len([k for k in g.predecessors('Rep') if k[0] == 'Z']) == 0)):
        return 0

    # All parts must be nodes in the circuit graph
    for n in part_list:
        if n not in g.nodes:
            return 0
        else:
            viable_type = [k[-2][0] 
                           for k in nx.all_simple_paths(g, promo_node, n)]
            if len(viable_type) == 0:
                return 0
            else:
                if (('I' in viable_type) and 
                    (('Z' not in viable_type) and ('P' not in viable_type))):
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

# not used
def compare_circuit(g1, g2):
    ind = (set(g1.edge_list) == set(g2.edge_list)) & (g1.dose == g2.dose)

    return ind

# gets node for crossover point between two circuits
def get_crosspt(list1, list2):

    # Create a list of shared circuit parts 
    # by iterating through both lists
    same = []
    for item in list1:
        if item in list2:
            same.append(item)

    # crossover point is same part if there are any
    if len(same) > 0:
        pt1 = np.random.choice(same)
        pt2 = pt1
    
    # otherwise, pick random part from 1st
    # cicuit as crossover point
    else:
        pt1 = np.random.choice(list1)

        # if TF, select crossover point for 2nd 
        # circuit as another TF
        if pt1[0] == 'Z':
            pt2 = np.random.choice([k for k in list2 if k[0] == 'Z'])

        # if inhibitor, select crossover point for
        # 2nd circuit as another inhibitor if possible
        else:
            list2_inhibitors = [k for k in list2 if k[0] == 'I']

            # Follow the same process if both circuits
            # have inhibitors
            if len(list2_inhibitors) > 0:
                pt2 = np.random.choice(list2_inhibitors)

            # New method for if the second circuit 
            # doesn't have an inhibitor
            else:

                # Switch inhibitor for tf if second 
                # circuit has more than 1 part
                if len(list2) > 1:
                    pt2 = np.random.choice(list2)

                # Have to switch tf for tf if second 
                # circuit only has 1 part
                else:
                    pt1 = np.random.choice([k for k in list1 if k[0] == 'Z'])
                    pt2 = np.random.choice([k for k in list2 if k[0] == 'Z'])

    return pt1, pt2

# switch two nodes and associated edges
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

# not used
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

# match node from pt2 predecessors or successors
# to node in g1 for switching edges
def match_node(new_node, part_list, promo_node, circuit_tf_list, 
               circuit_in_list, pt2, node_list2):
    for n in node_list2:
        # if node is crossover point 2,
        # use this node
        if n == pt2:
            new_node.append(pt2)
        # if node is in list of parts from child
        # use node
        elif n in (part_list + [promo_node, 'Rep']):
            new_node.append(n) #p1
        # if node is tf not in child
        elif n[0] == 'Z':
            # Determine available nodes by iterating 
            # through list of transcription factors
            # that aren't in existing new_node list    
            node_avail = []
            for item in circuit_tf_list:
                if item not in new_node:
                    node_avail.append(item)

            # randomly choose new_node from available
            # nodes
            if len(node_avail) > 0:
                n_new = np.random.choice(list(node_avail))
                new_node.append(n_new)
        elif n[0] == 'I':
            # same process if inhibitor but use
            # inhibitor list
            node_avail = []
            for item in circuit_in_list:
                if item not in new_node:
                    node_avail.append(item)

            if len(node_avail) > 0:
                n_new = np.random.choice(list(node_avail))
                new_node.append(n_new)

# switch node pt1 in g1 for node pt2 and edges from 
# pt2 to maintain structure as much as possible given 
# parts in g1
def switch_edge(g1, pt1, pt2, in_list2, out_list2, dose2):
    child = deepcopy(g1)

    # remove crossover pt1 (and all adjacent edges)
    # and pt1 dose from child graph
    child.part_list.remove(pt1)
    child.dose.pop(pt1)
    child.graph.remove_node(pt1)

    # get list of remaining tfs and inhibitors
    circuit_tf_list = [k for k in child.part_list if k[0] == 'Z']
    circuit_in_list = [k for k in child.part_list if k[0] == 'I'] 

    in_node = []
    # get list of parts in other graph that are both
    # predecessors and successors of pt2 (e.g., self 
    # regulation)
    common_list2 = [k for k in in_list2 if k in out_list2]
    # get new in and out node- either same part as in
    # common_list2 or same type 
    match_node(in_node, child.part_list, child.promo_node, circuit_tf_list, circuit_in_list, pt2, common_list2)
    out_node = [k for k in in_node if k[0] != 'P']

    # for parts not in common_list2- get new in and out
    # nodes and add to respective lists; use same part
    # or same type
    in_list2 = [k for k in in_list2 if k not in common_list2]
    match_node(in_node, child.part_list, child.promo_node, circuit_tf_list, circuit_in_list, pt2, in_list2)
    out_list2 = [k for k in out_list2 if k not in common_list2]
    match_node(out_node, child.part_list, child.promo_node, circuit_tf_list, circuit_in_list, pt2, out_list2)

    # Build an array of new edges by appending all unique 
    # edges from in_node and out_node
    new_edges = []
    for k in in_node:
        new_edges.append((k, pt2))
    for k in out_node:
        out_edge = (pt2,k) 
        if out_edge not in new_edges:
            new_edges.append(out_edge)

    # update part list, dose, and edges for child
    child.part_list.append(pt2)
    child.dose.update({pt2: dose2})
    child.graph.add_edges_from(new_edges)

    # ensure child is valid circuit by updating if
    # needed
    validate(child)

    return child

# switch crossover node and associated edges
def crossover_structure(g1, g2):
    pt1, pt2 = get_crosspt(g1.part_list, g2.part_list)

    child1 = switch_edge(g1, pt1, pt2, list(g2.graph.predecessors(pt2)),
                          list(g2.graph.successors(pt2)), g2.dose[pt2])
    child2 = switch_edge(g2, pt2, pt1, list(g1.graph.predecessors(pt1)), 
                         list(g1.graph.successors(pt1)), g1.dose[pt1])

    return child1, child2

# choose a random part from circuit and 
# select a new random dose for part
def mutate_dose(g, min_dose=10, max_dose=75, dose_interval=5):
    n = np.random.choice(g.part_list)
    g.dose.update({n: get_dose(min_dose, max_dose, dose_interval, 1)[0]})

# change tf or inhibitor variant of selected
# part (e.g., Z1 -> Z6) and maintain edges 
def mutate_node_type(g, min_dose=10, max_dose=75, dose_interval=5):
    # select random part from list
    old_node = np.random.choice(g.part_list)
    if old_node[0] == 'Z':
        # when replacing node, I with
        # same number will need to have
        # same out edges if in circuit 
        # (accounted for later)
        same_type = 'I'
        # node_avail contains all tfs not
        # already in circuit
        node_avail = []
        for item in tf_list:
            if item not in g.part_list:
                node_avail.append(item)
    else:
        # when replacing node, Z with
        # same number will need to have
        # same out edges if in circuit 
        # (accounted for later)
        same_type = 'Z'
        # node_avail contains all inhibitors
        # not already in circuit
        node_avail = []
        for item in inhibitor_list:
            if item not in g.part_list:
                node_avail.append(item)

    # randomly choose new node from node_avail
    new_node = np.random.choice(node_avail)
    # Z or I with same number stored to make sure 
    # it has same out edges as new_node if in
    # circuit
    same = same_type + new_node[1:]
    # add new node with same dose as old node and 
    # remove old node/dose from doses
    g.dose.update({new_node: g.dose[old_node]})
    g.dose.pop(old_node)
    # in each edge, replace old node with new node
    # if in edge and maintain other edges
    new_edges = switch_node(g, old_node, new_node)
    # remove old node and associated edges from
    # graph
    g.graph.remove_node(old_node)
    # update topo and associated attributes
    g.update(new_edges)
    # if Z or I with same number as new node is
    # in circuit, get out edges from that node
    # and new node
    if same in g.part_list:
        same_out = list(g.graph.out_edges(same))
        new_node_out = list(g.graph.out_edges(new_node))
        if same_type == 'Z':
            # if same_type is tf, use its out edges for
            # new node and remove others
            new_node_out_new = [(new_node, k[1]) for k in same_out]
            g.graph.remove_edges_from(new_node_out)
            g.graph.add_edges_from(new_node_out_new)
        else:
            # if same_type is inhibitor, use out edges 
            # from new node and remove ones from
            # same_type
            same_out_new = [(same, k[1]) for k in new_node_out]
            g.graph.remove_edges_from(same_out)
            g.graph.add_edges_from(same_out_new)
        # update topo and associated edges
        g.update(list(g.graph.edges))

# add a node to the circuit 
def add_node(g, circuit_tf_list, min_dose=10, 
             max_dose=75, dose_interval=5, inhibitor=False):
    
    # node_avail is tfs not in circuit if not using inhibitors
    if not inhibitor:
        node_avail = [k for k in tf_list if k not in circuit_tf_list]
    # if using inhibitors, node_avail is all parts not in circuit
    else:
        node_avail = [k for k in parts.keys() if k not in g.part_list]
    #choose random node from list
    new_node = np.random.choice(node_avail)
    new_node_type = new_node[0]
    # add random dose within interval for new node
    g.dose.update({new_node: get_dose(min_dose,
                                      max_dose, dose_interval, 1)[0]})

    # Create list of new edges starting with
    # edges currently in circuit
    new_edges = []
    for k in g.edge_list:
        new_edges.append(k)

    # Add new edges from new_node_in to new_edges 
    # by iterating through the list (uses tfs currently)
    # in topology
    new_node_in = get_in_path(new_node, g.promo_node, circuit_tf_list)
    for edge in new_node_in:
        if edge not in new_edges:
            new_edges.append(edge)

    if new_node_type == 'I':
        # if new node is I, check if Z with same number
        # is in circuit
        same_type = 'Z'
        same = same_type + new_node[1:]
        if same in g.part_list:
            # if same Z is in part list, use it's out edges
            # and tf regulation of new_node if in new_node_in
            same_out = list(g.graph.out_edges(same))
            all_out_edges_same = [edge for edge in new_node_in if edge[0] == same] + same_out

            # Update same out edges to use new node by 
            # iterating through the list (must be the 
            # same for out edges for I and Z of same 
            # number)
            for k in all_out_edges_same:
                possible_edge = (new_node, k[1])
                if possible_edge not in new_edges:
                    new_edges.append(possible_edge)
        else:
            # if same type Z is not in part list,
            # Add new edges from new_node_out to new_edges by 
            # iterating through the list
            new_node_out = get_out_path(new_node, g.part_list)
            for edge in new_node_out:
                if edge not in new_edges:
                    new_edges.append(edge)
    else:
        # if new node is Z, check if I with same number
        # is in circuit
        same_type = 'I'
        same = same_type + new_node[1:]
        if same in g.part_list:
            # if same I is in part list, get its out
            # edges
            same_out = list(g.graph.out_edges(same))

            # Remove shared edges between new_edges and 
            # same_out from the new_edges array
            # (getting new edges for the I of same type)
            for i in range(len(new_edges)):
                if new_edges[i] in same_out:
                    new_edges.pop(i)

            # get our edges for new node
            new_node_out = get_out_path(new_node, g.part_list)

            # Add new edges from new_node_out to 
            # new_edges by iterating through the list
            for edge in new_node_out:
                if edge not in new_edges:
                    new_edges.append(edge)

            # get all out edges for new node, including self
            # regulation
            all_out_edges_new = [edge for edge in new_node_in if edge[0] == new_node] + new_node_out

            # Update new node out edges to use same node 
            # by iterating through the list (must be the 
            # same for out edges for I and Z of same 
            # number)
            for k in all_out_edges_new:
                possible_edge = (same, k[1])
                if possible_edge not in new_edges:
                    new_edges.append(possible_edge)
        else:
            # if same type I is not in part list,
            # Add new edges from new_node_out to new_edges by 
            # iterating through the list
            new_node_out = get_out_path(new_node, g.part_list)
            for edge in new_node_out:
                if edge not in new_edges:
                    new_edges.append(edge)

    # update topo and associated edges
    g.update(new_edges)

# remove a random node from the circuit
def remove_node(g, circuit_tf_list):
    # get list of inhibitors in circuit
    circuit_in_list = [k for k in g.part_list if k[0] == 'I']

    # if more than 1 tf or more than 1 inhibitor
    # if only 1 tf and 1 inhibitor, does nothing
    if len(circuit_tf_list) > 1 | len(circuit_in_list) > 1:
        # if 1 or fewer inhibitor, choose from tf list
        if len(circuit_in_list) <= 1:
            old_node = np.random.choice(circuit_tf_list)
        # if 1 or fewer tf, choose from inhibitor list
        elif len(circuit_tf_list) <= 1:
            old_node = np.random.choice(circuit_in_list)
        # otherwise, choose from entire part list
        else:
            old_node = np.random.choice(g.part_list)

        # remove old_node from part list, doses, and
        # graph
        g.part_list.remove(old_node)
        g.dose.pop(old_node)
        g.graph.remove_node(old_node)

        # ensure child is valid circuit by updating if
        # needed
        validate(g)

# mutate number of nodes in the circuit
def mutate_node_num(g, max_part, min_dose=10, max_dose=75, 
                    dose_interval=5, inhibitor=False):
    # get list of tfs in circuit
    circuit_tf_list = [k for k in g.part_list if k[0] == 'Z']
    # only add or remove node if max_part > 1
    if max_part > 1:
        # if one part in circuit, add a node
        if len(g.part_list) == 1:
            add_node(g, circuit_tf_list, min_dose, max_dose, 
                     dose_interval, inhibitor)
        # if parts in circuit less than max part,
        # randomly choose to add or remove node
        elif len(g.part_list) < max_part:
            if np.random.uniform() < 0.5:
                add_node(g, circuit_tf_list, min_dose, max_dose, 
                         dose_interval, inhibitor)
            else:
                remove_node(g, circuit_tf_list)
        # otherwise, remove node
        else:
            remove_node(g, circuit_tf_list)

# get edges for part_list corresponding to all possible
# edges for circuit
def get_full_connected(part_list, promo_node):

    # all nodes have P1->node->rep 
    edge_list = [(promo_node, k) for k in part_list]
    edge_list.extend([(k, 'Rep') for k in part_list])
    # all nodes have node->node (self regulation)
    edge_list.extend([(k, k) for k in part_list])
    # all nodes have edges from all permutations from 
    # part_list
    edge_list.extend(list(permutations(part_list, 2)))
    return edge_list

# add or remove edges from circuit
def mutate_edge(g, inhibitor=False):
    # if only 1 part, has max of 3 edges
    if len(g.part_list) == 1:
        if g.graph.size() == 3:
            # if all possible edges, remove
            # self regulation
            g.edge_list.remove((g.part_list[0], g.part_list[0]))
        else:
            # if less than 3 edges, add self
            # regulation (other edges have to
            # be present)
            g.edge_list.append((g.part_list[0], g.part_list[0]))
        g.update(g.edge_list)
    else:
        # if more than 1 part, get full set of possible edges
        edge_full = get_full_connected(g.part_list, g.promo_node)
        # if circuit has all possible edges, remove edge
        if g.graph.size() == len(edge_full):
            if inhibitor:
                num_parts = len(g.part_list)

                # edge_avail = [k for k in g.edge_list]
                # choose random edge to remove from edge_full, 
                # if it makes circuit invalid, remove from 
                # edge_full and try again until run out of 
                # edges or circuit is valid
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
                # if no inhibitors, choose any edge to remove
                # (will still be valid)
                ind = np.random.choice(g.graph.size())
                g.edge_list.remove(edge_full[ind])
                g.update(g.edge_list)

        else:
            # if circuit does not have all possible edges
            # add or remove edges with 50/50 probability
            if np.random.uniform() < 0.5:
                # available dges to add are not in circuit but
                # are in full list of edges
                edge_avail = [k for k in edge_full if k not in g.graph.edges]
                # choose random edge to add
                ind = np.random.choice(len(edge_avail))
                # add to edge list and update topo and
                # associated edges
                g.edge_list.append(edge_avail[ind])
                g.update(g.edge_list)
            else:
                num_parts = len(g.part_list)
                # available edges to remove are in circuit
                edge_avail = [k for k in g.edge_list]
                valid = 0
                # choose random edge to remove from circuit, 
                # if it makes circuit invalid, remove from 
                # edge_avail and try again until run out of 
                # edges or circuit is valid
                while (valid == 0) and (len(edge_avail) > 0):
                    g_graph = deepcopy(g.graph)
                    ind = np.random.choice(len(edge_avail))
                    edge_removed = edge_avail[ind]
                    g_graph.remove_edge(edge_removed[0], edge_removed[1])
                    valid = check_valid(g_graph, g.promo_node, g.part_list)
                    edge_avail.remove(edge_removed)
                if valid == 1:
                    g.update(list(g_graph.edges))

# not used
def is_equal(x1, x2):
    return compare_circuit(x1[0], x2[0])

# not used
def tournament(X, obj, pressure=3):
    throuple = np.random.choice(range(len(X)), size=pressure, replace=False)
    return throuple[obj[throuple].argsort()][1:]

# perform crossover on all circuits in pop
def crossover(X, obj, rank_dict=None, **kwargs):
    # 2 parents per mating
    n_matings = int(len(X) / 2)
    Y = np.full_like(X, None, dtype=object)
    for k in range(n_matings):
        #choose 3 random indices in X
        throuple = np.random.choice(range(len(X)), 3,
                                    replace=False)
        # if single objective optimization
        if rank_dict is None:
            # parents = indices of circuits with top
            # 2 obj func in those indices
            parents = throuple[obj[throuple].argsort()][:-1]
        # if multi-objective optimization
        else:
            # sort by rank_dict for the three 
            # circuits and choose top 2
            throuple_rank = np.asarray([rank_dict[i]['rank'] for i in throuple])
            parents = throuple[throuple_rank.argsort()][:-1]
        # if strategy == "tournament":
        #     parents = tournament(X, obj)
        # elif strategy == "random":
        #     parents = np.random.choice(len(X), size=2, replace=False)
        # for mating 0, add to results to row 0, row 1 of Y; for
        # mating 1 add to row 2, row 3 of Y, etc
        # get topo from X with parents indices and column 0 (X 
        # # has shape (num_circuit, 1))
        Y[2 * k, 0], Y[2 * k + 1, 0] = crossover_structure(X[parents[0], 0], X[parents[1], 0])
    return Y

# perform mutation on all circuits in pop
def mutate(problem, X, prob, dose=False, **kwargs):
    for i in range(len(X)):
        # only mutation if < probability
        # for mutation
        if np.random.uniform(0, 1) < prob:
            # if dose, selecting randomly
            # from 4 possible mutations
            # (0, 1, 2, 3)
            if dose:
                r = np.random.choice(4)
            # if not dose, selecting
            # randomly from 3 possible
            # mutations (0, 1, 2)
            else:
                r = np.random.choice(3)

            # perform specified mutation
            # depending on random number
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