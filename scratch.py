import numpy as np
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
        self.graph.add_edges_from(edge_list) # Create graph object from edge_list
        self.node_list = [n for n in self.graph.nodes if n[0] != 'P'] # Get nodes that are not promoters
        self.promo = set(self.graph.nodes).difference(set(self.node_list))  # Get promoter nodes

        self.dose = {'Z': 200, 'I': 200, 'R': 9}
        self.protein_deg = {'Z': 0.35, 'I': 0.35, 'R': 0.029}
        self.node_dict = dict() # Classify nodes
        for n in self.node_list:
            edge = self.graph.in_edges(n)
            self.node_dict.update({n:
                                        {'P': [i[0] for i in edge if i[0][0] == 'P'],
                                        'Z': [i[0] for i in edge if i[0][0] == 'Z'],
                                        'I': [i[0] for i in edge if i[0][0] == 'I']}})
        self.num_states = len(self.node_dict.keys())
        self.var_dict = dict(zip(self.node_dict.keys(), np.arange(self.num_states)))
        self.valid = None

    def check_valid(self):
        self.valid = 1
        for n in self.node_list:
            if (len(self.node_dict[n]['I']) > 0) & (len(self.node_dict[n]['P']) > 0) & (len(self.node_dict[n]['Z']) == 0):
                self.valid = 0
            elif (n != 'Rep') & ((len(self.graph.in_edges('Z1')) == 0) | (len(self.graph.out_edges('Z1')) == 0)):
                self.valid = 0
        if len(self.promo) != 1:
            self.valid = 0

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
        t = np.arange(0, 48+1, 1)
        rep_off = odeint(self.system_equations, np.zeros(self.num_states*2), t, args=('off',))[-1,-1]
        rep_on = odeint(self.system_equations, np.zeros(self.num_states*2), t, args=('on',))[-1,-1]
        return rep_off, rep_on

    def get_fitness(self):
        rep_off, rep_on = self.simulate()
        ON_ratio = rep_on/Ref[next(iter(self.promo))]['on']
        return ON_ratio

    def plot_graph(self):
        plt.tight_layout()
        nx.draw_networkx(self.graph, arrows=True)
        plt.show()



edge_list = [("P1", "Z1"), ("Z1", "Rep")]
g1 = Topo(edge_list)
g1.get_fitness()
# g1.simulate()

a = 1