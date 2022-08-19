import numpy as np
# import pandas as pd
import pickle
import networkx as nx
from scipy.integrate import odeint
import matplotlib.pyplot as plt

with open("promo.pkl", "rb") as fid:
    promo = pickle.load(fid)
with open("parts.pkl", "rb") as fid:
    parts = pickle.load(fid)
with open("Ref.pkl", "rb") as fid:
    Ref = pickle.load(fid)

class Topo():
    def __init__(self, edge_list):
        self.graph = nx.DiGraph()
        self.graph.add_edges_from(edge_list)
        self.nodes = self.graph.nodes
        self.dose = {'Z': 200, 'I': 200, 'R': 9}
        self.protein_deg = {'Z': 0.35, 'I': 0.35, 'R': 0.029}
        self.node_dict = None
        self.num_states = None
        self.var_dict = None

    def plot_graph(self):
        plt.tight_layout()
        nx.draw_networkx(self.graph, arrows=True)
        plt.show()

    def get_nodes(self):
        self.node_dict = dict()
        for n in self.nodes:
            if n[0] != "P":
                edge = self.graph.in_edges(n)
                # self.node_dict.update({n: [i[0] for i in edge]})
                self.node_dict.update({n:
                                         {'P': [i[0] for i in edge if i[0][0] == 'P'],
                                          'Z': [i[0] for i in edge if i[0][0] == 'Z'],
                                          'I': [i[0] for i in edge if i[0][0] == 'I']}})
        self.num_states = len(self.node_dict.keys())
        self.var_dict = dict(zip(self.node_dict.keys(), np.arange(self.num_states)))

    def system_equations(self, x, t, state):
        system = []
        for n in self.node_dict.keys():
            eq = -2.7 * x[2 * self.var_dict[n]]
            b = 0
            num = 0
            denom = 1
            for k in self.node_dict[n]['P']:
                eq += self.dose[n[0]] * promo[k][state]
            for k in self.node_dict[n]['Z']:
                b += parts[k][0]
                num += parts[k][1] * parts[k][2] * x[2 * self.var_dict[k] + 1]
                denom += parts[k][2] * x[2 * self.var_dict[k] + 1]
            for k in self.node_dict[n]['I']:
                denom += parts[k][0] * x[2 * self.var_dict[k] + 1]
            if len(self.node_dict[n]['Z']) > 0:
                b /= len(self.node_dict[n]['Z'])
            eq += self.dose[n[0]] * (b + num) / denom
            system.extend([eq, -self.protein_deg[n[0]] * x[2 * self.var_dict[n] + 1] + x[2 * self.var_dict[n]]])
        return system

    def simulate(self):
        t = np.arange(0,48,1)
        sol = odeint(self.system_equations, np.zeros(self.num_states*2), t, args=('off',))
        rep_off = sol[-1,-1]
        sol = odeint(self.system_equations, np.zeros(self.num_states*2), t, args=('on',))
        rep_on = sol[-1,-1]
        return rep_off, rep_on

    def get_fitness(self):




edge_list = [("P1", "Z1"), ("Z1", "Rep")]
g1 = Topo(edge_list)
g1.get_nodes()
# g1.simulate()
# def fun(x, state):
#     system = []
#     for n in g1.node_dict.keys():
#         eq = -2.7 * x[2 * g1.var_dict[n]]
#         b = 0
#         num = 0
#         denom = 1
#         for k in g1.node_dict[n]['P']:
#             eq += 200 * promo[k][state]
#         for k in g1.node_dict[n]['Z']:
#             b += parts[k][0]
#             num += parts[k][1] * parts[k][2] * x[2 * g1.var_dict[k] + 1]
#             denom += parts[k][2] * x[2 * g1.var_dict[k] + 1]
#         for k in g1.node_dict[n]['I']:
#             denom += parts[k][0] * x[2 * g1.var_dict[k] + 1]
#         if len(g1.node_dict[n]['Z']) > 0:
#             b /= len(g1.node_dict[n]['Z'])
#         eq += 200*(b + num) / denom
#         system.extend([eq, -0.35 * x[2 * g1.var_dict[n] + 1] + x[2 * g1.var_dict[n]]])
#     return system

a = 1