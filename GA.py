from define_circuit import *
import copy

def get_circuit(promo_node, part_list):
    tf_list = [k for k in part_list if k[0] == 'Z']
    edge_list = set([(promo_node, np.random.choice(tf_list)), (np.random.choice(tf_list), 'Rep')])
    for n in part_list:
        edge_list.update([(np.random.choice(tf_list + [promo_node]), n)])
        out_list = [k for k in part_list if k != n]
        edge_list.update([(n, np.random.choice(out_list + ['Rep']))])

    return list(edge_list)


def sample_circuit(promo_node, tf_list, num_circuit, max_part=2, max_dose=200, min_dose=20, inhibitor=False, in_list=None):
    circuits = []

    if not inhibitor:
        for i in range(num_circuit):
            num_part = np.random.randint(1, max_part+1)
            part_list = np.random.choice(tf_list, num_part)
            dose_list = np.random.randint(min_dose, max_dose, size=num_part)
            edge_list = get_circuit(promo_node, part_list)
            circuits.append(Topo(edge_list, part_list, dose_list, promo_node))

    elif in_list != None:
        for i in range(num_circuit):
            num_tf = np.random.randint(1, max_part)
            num_in = np.random.randint(1, max_part-num_tf+1)
            part_list = np.append(np.random.choice(tf_list, num_tf), np.random.choice(in_list, num_in))
            dose_list = np.random.randint(min_dose, max_dose, size=(num_tf+num_in))
            edge_list = get_circuit(promo_node, part_list)
            circuits.append(Topo(edge_list, part_list, dose_list, promo_node))
    else:
        raise Exception("Need a list of inihibitors")

    return circuits

def node_switch(g, old_node, new_node):
    child_edge = []
    for edge in g.edge_list:
        if old_node in edge:
            edge = list(edge)
            ind = edge.index(old_node)
            edge[ind] = new_node
            child_edge.append(tuple(edge))
        else:
            child_edge.append(edge)
    return child_edge

def crossover(g1, node1, g2, node2):
    child1_edge = node_switch(g1, node1, node2)
    child2_edge = node_switch(g2, node2, node1)

    return child1_edge, child2_edge

# edge_list = [("P1", "Z1"), ("Z1", "Rep")]
# edge_list2 = [("P1", "Z2"), ("Z2", "Rep")]
#
# g1 = Topo(edge_list, {'Z1': 200}, 'P1')
# g2 = Topo(edge_list2, {'Z2': 200}, 'P1')
# t = np.arange(0, 48+1, 1)

def get_crosspt(list_1, list_2):
    same = set(list_1).intersection(set(list_2))
    if len(same) > 0:
        crosspt_1 = np.random.choice(same)
        crosspt_2 = crosspt_1
    else:
        crosspt_1 = np.random.choice(list_1)
        crosspt_2 = np.random.choice(list_2)
    return crosspt_1, crosspt_2




a = 2