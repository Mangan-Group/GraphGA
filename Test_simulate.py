from scipy.integrate import odeint
from load_files_pop import *
from get_system_equations_pop import system_equations_pop
from define_circuit import Topo

def simulate_cell(topology, Z_list, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_on = odeint(system_equations_pop, np.zeros(topology.num_states * 2), t, args=('on', Z_list, topology,))[-1, -1]
    return rep_on

def simulate_pop(topology, Z, max_time=42):
    pop_rep = []
    for i in range(0, len(Z)):
        Z_list = Z[i, :]
        rep_on = simulate_cell(topology, Z_list, max_time)
        pop_rep.append(rep_on)

    pop_rep_mean = np.mean(pop_rep)

    return pop_rep_mean

topology = Topo([('P1', 'Z6'), ('Z6', 'Rep')], {'Z6': 75}, 'P1')

pop_rep_mean = simulate_pop(topology, Z_200, 42)

print(pop_rep_mean)
