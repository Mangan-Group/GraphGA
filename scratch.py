import numpy as np
import pandas as pd
import pickle
import networkx as nx
import sympy as sm

import matplotlib.pyplot as plt
with open("promo.pkl", "rb") as fid:
    promo = pickle.load(fid)
with open("parts.pkl", "rb") as fid:
    parts = pickle.load(fid)

class Topo():
    def __init__(self, edge_list):
        self.graph = nx.DiGraph()
        self.graph.add_edges_from(edge_list)
        self.nodes = self.graph.nodes
        # self.in_edges = self.graph.in_edges
        # self.out_edges = self.graph.out_edges

    def plot_graph(self):
        plt.tight_layout()
        nx.draw_networkx(self.graph, arrows=True)
        plt.show()

    def get_in_nodes(self):
        self.node_dict = dict()
        # state = 0
        for n in self.nodes:
            if n[0] != "P":
                edge = self.graph.in_edges(n)
                self.node_dict.update({n: [i[0] for i in edge]})

    def make_equations(self):
        self.get_in_nodes()
        self.num_states = len(self.node_dict.keys())
        # states_sm = ''
        # for i in range(self.num_states):
        #     states_sm += ('x' + str(i))
        #     if i == self.num_states - 1:
        #         break
        #     states_sm += ','
        # X = sm.Matrix(sm.symbols(states_sm))
        self.eq_dict = dict(zip(self.node_dict.keys(), np.arange(self.num_states)))

edge_list = [("P1", "Z1"), ("Z1", "Rep")]
g1 = Topo(edge_list)

def fun(x, state):
    system = []
    for n in g1.node_dict.keys():
        eq = -2.7 * x[2 * g1.eq_dict[n]]
        b = 0
        num = 0
        denom = 1
        for k in eq_dict['P']:
            eq += 200 * promo[k][state]
        for k in eq_dict['Z']:
            b += parts[k][0]
            num += parts[k][1] * parts[k][2] * x[2 * g1.eq_dict[k] + 1]
            denom += parts[k][2] * x[2 * g1.eq_dict[k] + 1]
        for k in eq_dict['I']:
            denom += parts[k][0] * x[2 * g1.eq_dict[k] + 1]
        if len(eq_dict['Z']) > 0:
            b /= len(eq['Z'])
        eq += 200*(b + num) / denom
        system.extend([eq, -0.35 * x[2 * g1.eq_dict[n] + 1] + x[2 * g1.eq_dict[n]]])
    return system

a = 1