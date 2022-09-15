import numpy as np
import networkx as nx
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from load_files import *


class Topo():
    def __init__(self, edge_list, dose_list, promo_node):
        self.edge_list = edge_list
        self.graph = nx.DiGraph()
        self.graph.add_edges_from(self.edge_list) # Create graph object from edge_list

        self.promo_node = promo_node  # Get promoter nodes

        self.dose = dose_list
        self.dose.update({'Rep': 9})
        self.part_list = [k for k in self.dose.keys() if k != 'Rep']

        self.protein_deg = {'Z': 0.35, 'I': 0.35, 'R': 0.029}

        self.in_dict = dict() # Classify nodes
        for n in (self.part_list + ['Rep']):
            # in_edge = self.graph.in_edges(n)
            # self.in_dict.update({n:
            #                             {'P': [i[0] for i in in_edge if i[0][0] == 'P'],
            #                             'Z': [i[0] for i in in_edge if i[0][0] == 'Z'],
            #                             'I': [i[0] for i in in_edge if i[0][0] == 'I']}})
            pre = list(self.graph.predecessors(n))
            self.in_dict.update({n:
                                     {'P': [i for i in pre if i[0] == 'P'],
                                      'Z': [i for i in pre if i[0] == 'Z'],
                                      'I': [i for i in pre if i[0] == 'I']}})

        self.num_states = len(self.in_dict.keys())
        self.var_dict = dict(zip((self.in_dict.keys()), np.arange(self.num_states)))
        self.valid = None

    def check_valid(self):
        self.valid = 1
        for n in self.part_list:
            # if (len(self.in_dict[n]['I']) > 0) & (len(self.in_dict[n]['Z']) == 0):
            #     self.valid = 0
            # elif (n != 'Rep') & ((len(self.graph.in_edges(n)) == 0) | (len(self.graph.out_edges(n)) == 0)):
            #     self.valid = 0
            if self.graph.in_degree(n) <= len(self.in_dict[n]['I']):
                self.valid = 0
            if len(list(nx.all_simple_paths(self.graph, n, 'Rep'))) == 0:
                self.valid = 0

    def update(self, edge_list):
        self.edge_list = edge_list
        self.graph = nx.DiGraph()
        self.graph.add_edges_from(self.edge_list)
        self.part_list = [k for k in self.dose.keys() if k != 'Rep']

        self.in_dict = dict()  # Classify nodes
        for n in (self.part_list + ['Rep']):
            pre = list(self.graph.predecessors(n))
            self.in_dict.update({n:
                                     {'P': [i for i in pre if i[0] == 'P'],
                                      'Z': [i for i in pre if i[0] == 'Z'],
                                      'I': [i for i in pre if i[0] == 'I']}})

        self.num_states = len(self.in_dict.keys())
        self.var_dict = dict(zip((self.in_dict.keys()), np.arange(self.num_states)))

    def system_equations(self, x, t, state):
        system = []
        for n in self.in_dict.keys():
            eq = -2.7 * x[2 * self.var_dict[n]]
            b = 0
            num = 0
            denom = 1
            for k in self.in_dict[n]['P']:
                eq += self.dose[n] * promo[k][state]
            for k in self.in_dict[n]['Z']:
                # print(k)
                b += parts[k][0]
                num += parts[k][1] * parts[k][2] * x[2 * self.var_dict[k] + 1]
                denom += parts[k][2] * x[2 * self.var_dict[k] + 1]
            for k in self.in_dict[n]['I']:
                denom += parts[k][0] * x[2 * self.var_dict[k] + 1]
            if len(self.in_dict[n]['Z']) > 0:
                b /= len(self.in_dict[n]['Z'])
            eq += self.dose[n] * (b + num) / denom
            system.extend([eq, -self.protein_deg[n[0]] * x[2 * self.var_dict[n] + 1] + x[2 * self.var_dict[n]]])
        return system

    def simulate(self, max_time=48):
        t = np.arange(0, max_time+1, 1)
        rep_off = odeint(self.system_equations, np.zeros(self.num_states*2), t, args=('off',))[-1,-1]
        rep_on = odeint(self.system_equations, np.zeros(self.num_states*2), t, args=('on',))[-1,-1]
        return rep_off, rep_on

    def plot_graph(self):
        plt.figure()
        plt.tight_layout()
        nx.draw_networkx(self.graph, arrows=True, arrowsize=15, node_size=600, node_shape='s')
        plt.show()


