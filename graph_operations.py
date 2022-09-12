from define_circuit import *
from itertools import combinations_with_replacement

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
    in_node = np.random.choice(circuit_tf_list + [promo_node])
    edges = [(in_node, n)]
    if in_node == n:
        new_tf_list = [k for k in circuit_tf_list if k != n]
        in_node = np.random.choice(new_tf_list + [promo_node])
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

def sample_circuit(promo_node, num_circuit, max_part=2, max_dose=200, min_dose=20, inhibitor=False):
    circuits = []
    if not inhibitor:
        for i in range(num_circuit):
            num_part = np.random.randint(1, max_part+1)
            part_list = np.random.choice(tf_list, num_part, replace=False)
            dose_list = dict(zip(part_list, np.random.randint(min_dose, max_dose, size=num_part)))
            edge_list = get_edges(promo_node, part_list)
            circuits.append(Topo(edge_list, dose_list, promo_node))
    else:
        for i in range(num_circuit):
            num_tf = np.random.randint(1, max_part)
            num_in = np.random.randint(1, max_part-num_tf+1)
            part_list = np.append(np.random.choice(tf_list, num_tf), np.random.choice(in_list, num_in))
            dose_list = dict(zip(part_list, np.random.randint(min_dose, max_dose, size=(num_tf+num_in))))
            edge_list = get_edges(promo_node, part_list)
            circuits.append(Topo(edge_list, dose_list, promo_node))

    return circuits

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

def validate(g):
    new_edges = set([k for k in g.graph.edges])
    circuit_tf_list = [k for k in g.part_list if k[0] == 'Z']
    for n in g.part_list:
        if g.graph.in_degree(n) <= len(g.in_dict[n]['I']):
            new_edges.update(get_in_path(n, g.promo_node, circuit_tf_list))

        if len(list(nx.all_simple_paths(g.graph, n, 'Rep'))) == 0:
            new_edges.update(get_out_path(n, g.part_list))

    return list(new_edges)

def compare_circuit(g1, g2):
        ind = (g1.dose == g2.dose) & (set(g1.edge_list) == set(g2.edge_list))
        return ind

def crossover_naive(g1, g2):
    pt1, pt2 = get_crosspt(g1.part_list, g2.part_list,)

    child1_edge = switch_node(g1, pt1, pt2)
    child2_edge = switch_node(g2, pt2, pt1)
    child1_dose = {k: g1.dose[k] for k in g1.part_list if k != pt1}
    child1_dose.update({pt2: g2.dose[pt2]})
    child2_dose = {k: g2.dose[k] for k in g2.part_list if k != pt2}
    child2_dose.update({pt1: g1.dose[pt1]})

    child1 = Topo(child1_edge, child1_dose, g1.promo_node)
    child2 = Topo(child2_edge, child2_dose, g2.promo_node)
    # validate(child1)
    # child1.update(child1.edge_list)
    # validate(child2)
    # child2.update(child2.edge_list)
    # child1.check_valid()
    # child2.check_valid()
    # if child1.valid == 0 | child2.valid == 0:
    #     print("error")
    return child1, child2

def mutate_dose(g, min_dose=10, max_dose=75):
    n = np.random.choice(g.part_list)
    g.dose.update({n: np.random.randint(min_dose, max_dose)})

def mutate_node_type(g, min_dose=10, max_dose=75):
    old_node = np.random.choice(g.part_list)
    if old_node[0] == 'Z':
        node_avail = list(set(tf_list).difference(set(g.part_list)))
    else:
        node_avail = list(set(in_list).difference(set(g.part_list)))
    new_node = np.random.choice(node_avail)
    g.dose.pop(old_node)
    g.dose.update({new_node: np.random.randint(min_dose, max_dose)})
    new_edges = switch_node(g, old_node, new_node)
    g.graph.remove_node(old_node)
    g.update(new_edges)

def mutate_node_num(g, max_part, min_dose=10, max_dose=75, inhibitor=False):
    circuit_tf_list = [k for k in g.part_list if k[0] == 'Z']
    if len(g.part_list) < max_part:
        # print(0)
        if not inhibitor:
            node_avail = [k for k in tf_list if k not in circuit_tf_list]
        else:
            node_avail = [k for k in parts.keys() if k not in g.part_list]
        new_node = np.random.choice(node_avail)
        g.dose.update({new_node: np.random.randint(min_dose, max_dose)})
        new_edges = set([k for k in g.edge_list])
        new_edges.update(get_in_path(new_node, g.promo_node, circuit_tf_list))
        new_edges.update(get_out_path(new_node, g.part_list))
        g.update(list(new_edges))

    else:
        circuit_in_list = [k for k in g.part_list if k[0] == 'I']
        if len(circuit_tf_list) > 1 & len(circuit_in_list) <= 1:
            old_node = np.random.choice(circuit_tf_list)
        elif len(circuit_tf_list) == 1 & len(circuit_in_list) > 1:
            old_node = np.random.choice(circuit_in_list)
        elif len(circuit_tf_list) > 1 & len(circuit_in_list) > 1:
            old_node = np.random.choice(g.part_list)

        g.part_list.remove(old_node)
        g.dose.pop(old_node)
        g.graph.remove_node(old_node)
        new_edges = validate(g)
        g.update(new_edges)

def get_full_connected(part_list, promo_node):
    edge_list = [(promo_node, k) for k in part_list]
    edge_list.extend([(k, 'Rep') for k in part_list])
    edge_list.extend(list(combinations_with_replacement(X[0].part_list, 2)))
    return edge_list

def mutate_edge(g):
    pass
