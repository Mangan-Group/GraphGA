def is_valid(g, n):
    circuit_tf_list = []
    same_list = []
    for k in g.part_list:
        if k[0] == 'Z':
            circuit_tf_list.append(k)
            if ('I'+k[1:]) in g.part_list:
                same_list.append([k, ('I'+k[1:])])

    if len(circuit_tf_list) == 0:
        raise Exception("Something's wrong. No TFs in the circuit.")

    if ('Rep' not in g.graph.nodes) | (len([k for k in g.graph.predecessors('Rep') if k[0] == 'Z']) == 0):
        return 0



    viable_type = [k[-2][0] for k in nx.all_simple_paths(g.graph, g.promo_node, n)]
            if len(viable_type) == 0:
                return 0
            else:
                if ('I' in viable_type) & ('Z' not in viable_type):
                    return 0

            if len(list(nx.all_simple_paths(g.graph, n, 'Rep'))) == 0:
                return 0

    #if set(g.graph.edges) != set(g.edge_list):
     #   g.update(list(g.graph.edges))

    for z, i in same_list:
        z_succ = set(g.graph.successors(z))
        i_succ = set(g.graph.successors(i))
        if z_succ != i_succ:
            return 0

    #if set(g.graph.edges) != set(g.edge_list):
    #    g.update(list(g.graph.edges))


def mutate_edge(g):
    if len(g.part_list) == 1:
        if (g.part_list[0], g.part_list[0]) in g.edge_list:
            g.graph.remove_edge(g.part_list[0], g.part_list[0])
        else:
            g.graph.add_edge(g.part_list[0], g.part_list[0])
        g.update(list(g.graph.edges))
    # elif (len(g.part_list) == 2) and ('I' in [k[0] for k in g.part_list]):

    else:
        edge_full = get_full_connected(g.part_list, g.promo_node)
        if g.graph.size() == len(edge_full):
            remove_edge(g)
        else:
            if np.random.uniform() < 0.5:
                remove_edge(g)
            else:
                add_edge(g)

def remove_edge(g):
    node = np.random.choice(g.part_list)
    # circuit_tf_list = [k for k in g.part_list if (k[0] == 'Z') and (k != node)]
    node_type = node[0]
    # new_edges = set([k for k in g.edge_list])
    if np.random.uniform() < 0.5:
        neighbor = np.random.choice(list(g.graph.predecessors(node)))
        g.graph.remove_edge(neighbor, node)
        if not list(nx.all_simple_paths(g.graph, 'P1', node)):
            potential_in_list = [k for k in g.part_list if (k[0] == 'Z') and (k != node) and (k != neighbor)]
            potential_in_list += ['P1']
            g.graph.add_edge.append((np.random.choice(potential_in_list), node))
        if not list(nx.all_simple_paths(g.graph, node, 'Rep')):
            potential_out_list = [k for k in g.part_list if (k != node) and ]

    else:
        neighbor = np.random.choice(list(g.graph.successors(node)))
        new_edges.remove_edge(node, neighbor)



    g.update(list(g.graph.edges))

def add_edge(g, edge_full):
    edge_avail = [k for k in edge_full if k not in g.graph.edges]
    ind = np.random.choice(len(edge_avail))
    g.graph.add_edge(edge_avail[ind])
    g.update(list(g.graph.edges))