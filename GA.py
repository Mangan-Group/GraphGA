from define_circuit import *
import copy

def get_circuit(promo_node, part_list):
    tf_list = [k for k in part_list if k[0] == 'Z']
    edge_list = set([(promo_node, np.random.choice(tf_list)), (np.random.choice(tf_list), 'Rep')])
    for n in part_list:
        edge_list.update([(np.random.choice(tf_list + [promo_node]), n)])
        out_list = [k for k in part_list if k != n]
        num_connect = np.random.randint(len(out_list))
        out_path = [n]
        out_path.extend(np.random.choice(out_list, num_connect, replace=False))
        out_path.append('Rep')
        edges = [(i, j) for i, j in zip(out_path[:-1], out_path[1:])]
        edge_list.update(edges)

    return list(edge_list)

def sample_circuit(promo_node, tf_list, num_circuit, max_part=2, max_dose=200, min_dose=20, inhibitor=False, in_list=None):
    circuits = []
    if not inhibitor:
        for i in range(num_circuit):
            num_part = np.random.randint(1, max_part+1)
            part_list = np.random.choice(tf_list, num_part)
            dose_list = dict(zip(part_list, np.random.randint(min_dose, max_dose, size=num_part)))
            edge_list = get_circuit(promo_node, part_list)
            circuits.append(Topo(edge_list, dose_list, promo_node))

    elif in_list == None:
        raise Exception("Need a list of inihibitors")

    else:
        for i in range(num_circuit):
            num_tf = np.random.randint(1, max_part)
            num_in = np.random.randint(1, max_part-num_tf+1)
            part_list = np.append(np.random.choice(tf_list, num_tf), np.random.choice(in_list, num_in))
            dose_list = dict(zip(part_list, np.random.randint(min_dose, max_dose, size=(num_tf+num_in))))
            edge_list = get_circuit(promo_node, part_list)
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
        pt2 = np.random.choice(list2)

    return pt1, pt2

def switch_naive(g, old_node, new_node):
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

def crossover(g1, g2, naive=False):
    pt1, pt2 = get_crosspt(g1.part_list, g2.part_list)
    if naive:
        child1_edge = switch_naive(g1, pt1, pt2)
        child2_edge = switch_naive(g2, pt2, pt1)
        child1_dose = {k:g1.dose[k] for k in g1.part_list if k != pt1}
        child1_dose.update({pt2: g2.dose[pt2]})
        child2_dose = {k: g2.dose[k] for k in g2.part_list if k != pt2}
        child2_dose.update({pt1: g1.dose[pt1]})

    else:
        pass

    child1 = Topo(child1_edge, child1_dose, g1.promo)
    child2 = Topo(child2_edge, child2_dose, g2.promo)
    return child1, child2









tf_list = [k for k in parts.keys() if k[0]=='Z']
in_list = [k for k in parts.keys() if k[0]=='I']

circuits = sample_circuit('P1', tf_list, num_circuit=20, max_part=3, inhibitor=True, in_list=in_list)


a = 2