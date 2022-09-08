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
tf_list = [k for k in parts.keys() if k[0]=='Z']
in_list = [k for k in parts.keys() if k[0]=='I']

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
        # self.out_dict = dict()
        for n in (self.dose.keys()):
            in_edge = self.graph.in_edges(n)
            self.in_dict.update({n:
                                        {'P': [i[0] for i in in_edge if i[0][0] == 'P'],
                                        'Z': [i[0] for i in in_edge if i[0][0] == 'Z'],
                                        'I': [i[0] for i in in_edge if i[0][0] == 'I']}})

            # out_edge = self.graph.out_edges(n)
            # self.out_dict.update({n:
            #                          {'R': [i[1] for i in out_edge if i[1][0] == 'R'],
            #                           'Z': [i[1] for i in out_edge if i[1][0] == 'Z'],
            #                           'I': [i[1] for i in out_edge if i[1][0] == 'I']}})

        self.num_states = len(self.dose.keys())
        self.var_dict = dict(zip(self.in_dict.keys(), np.arange(self.num_states)))
        self.valid = None

    def check_valid(self):
        self.valid = 1
        for n in self.dose.keys():
            if (len(self.in_dict[n]['I']) > 0) & (len(self.in_dict[n]['Z']) == 0):
                self.valid = 0
            elif (n != 'Rep') & ((len(self.graph.in_edges(n)) == 0) | (len(self.graph.out_edges(n)) == 0)):
                self.valid = 0
        # if len(self.promo) != 1:
        #     self.valid = 0

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

    def get_fitness(self):
        rep_off, rep_on = self.simulate()
        ON_ratio = rep_on/Ref[self.promo_node]['on']
        return ON_ratio

    def plot_graph(self):
        plt.figure()
        plt.tight_layout()
        nx.draw_networkx(self.graph, arrows=True, arrowsize=15, node_size=600, node_shape='s')
        plt.show()

def amplifier_obj(g):
    rep_off, rep_on = g.simulate()
    ON_ratio = rep_on / Ref[g.promo_node]['on']
    FIn = rep_on / rep_off
    return np.array([ON_ratio, FIn])
