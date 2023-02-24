import numpy as np
import networkx as nx
# from scipy.integrate import odeint
import matplotlib.pyplot as plt
from load_files import *


class Topo():
    def __init__(self, edge_list, dose_list, promo_node):
        self.edge_list = edge_list
        self.graph = nx.DiGraph()
        self.graph.add_edges_from(self.edge_list) # Create graph object from edge_list

        self.promo_node = promo_node  # Get promoter nodes

        self.dose = dose_list
        self.dose.update({'Rep': 1})
        self.part_list = [k for k in self.dose.keys() if k != 'Rep']

        self.protein_deg = {'Z': 0.35, 'I': 0.35, 'R': 0.029}

        self.in_dict = dict() # Classify nodes
        self.pool = dict()
        for n in (self.part_list + ['Rep']):
            pre = list(self.graph.predecessors(n))
            self.in_dict.update({n:
                                     {'P': [i for i in pre if i[0] == 'P'],
                                      'Z': [i for i in pre if i[0] == 'Z'],
                                      'I': [i for i in pre if i[0] == 'I']}})
            self.pool.update({n: (self.in_dict[n]['P'] != []) + (self.in_dict[n]['Z'] != [])})
        if 0 in list(self.pool.values()):
            raise Exception("Something's wrong. No activator in the circuit.")

        self.num_states = len(self.in_dict.keys())
        self.var_dict = dict(zip((self.in_dict.keys()), np.arange(self.num_states)))
        self.valid = None

    def check_valid(self):
        self.valid = 1
        for n in self.part_list:
            if self.graph.in_degree(n) <= len(self.in_dict[n]['I']):
                self.valid = 0
            if len(list(nx.all_simple_paths(self.graph, n, 'Rep'))) == 0:
                self.valid = 0

    def update(self, edge_list):
        self.edge_list = edge_list
        self.graph = nx.DiGraph(self.edge_list)
        # self.graph.add_edges_from(self.edge_list)
        self.part_list = [k for k in self.dose.keys() if k != 'Rep']

        self.in_dict = dict()  # Classify nodes
        self.pool = dict()
        for n in (self.part_list + ['Rep']):
            pre = list(self.graph.predecessors(n))
            self.in_dict.update({n:
                                     {'P': [i for i in pre if i[0] == 'P'],
                                      'Z': [i for i in pre if i[0] == 'Z'],
                                      'I': [i for i in pre if i[0] == 'I']}})
            self.pool.update({n: (self.in_dict[n]['P'] != []) + (self.in_dict[n]['Z'] != [])})
        if 0 in list(self.pool.values()):
            raise Exception("Something's wrong. No activator in the circuit.")

        self.num_states = len(self.in_dict.keys())
        self.var_dict = dict(zip((self.in_dict.keys()), np.arange(self.num_states)))


    def plot_graph(self):
        plt.figure()
        plt.tight_layout()
        nx.draw_networkx(self.graph, arrows=True, arrowsize=15, node_size=600, node_shape='s')
        plt.show()


