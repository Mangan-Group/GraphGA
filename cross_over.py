from define_circuit import *
import copy
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

edge_list = [("P1", "Z1"), ("Z1", "Rep")]
edge_list2 = [("P1", "Z2"), ("Z2", "Rep")]

g1 = Topo(edge_list, {'Z1': 200}, 'P1')
g2 = Topo(edge_list2, {'Z2': 200}, 'P1')
t = np.arange(0, 48+1, 1)

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