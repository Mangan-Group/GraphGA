from scipy.integrate import odeint
from multiprocessing import Pool
from load_files_pop import *
from get_system_equations_pop import system_equations_pop
from define_circuit import Topo

def simulate_cell(topology, Z_list, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_on = odeint(system_equations_pop, np.zeros(topology.num_states * 2), t, args=('on', Z_list, topology,))[-1, -1]
    return rep_on

def simulate_pop(topology, Z, num_processes, max_time=42):
    nc = len(Z)
    zipped_args = list(zip([topology]*nc, Z, [max_time]*nc))
    pop_rep_on = []
    with Pool(num_processes) as pool:
        pop_rep_on = pool.starmap(
            simulate_cell,
            zipped_args,
        )

    pool.close()
    pool.join()

    rep_on_mean = np.mean(pop_rep_on)
    return rep_on_mean


if __name__ == '__main__':

    topology = Topo([('P1', 'Z6'), ('P1', 'I6'), ('Z6', 'Rep'), ('I6', 'Rep')], {'Z6': 75, 'I6': 5}, 'P1')
    rep_on_mean = simulate_pop(topology, Z_200, 42)

    print(rep_on_mean)
