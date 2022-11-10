from scipy.integrate import odeint
from multiprocessing import Pool
from load_files_pop import *
from get_system_equations_pop import system_equations_pop
from define_circuit import Topo

def simulate_cell(topology, Z_list, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_off_end = odeint(system_equations_pop, np.zeros(topology.num_states * 2), t, args=('off', Z_list, topology,))[-1, -1]
    rep_on_ts = odeint(system_equations_pop, np.zeros(topology.num_states * 2), t, args=('on', Z_list, topology,))[:, -1]
    rep_on_end = rep_on_ts[-1]
    return rep_off_end, rep_on_end, rep_on_ts

def simulate_pop(topology, Z, num_processes, max_time=42):
    nc = len(Z)
    zipped_args = list(zip([topology]*nc, Z, [max_time]*nc))
    with Pool(num_processes) as pool:
        results = pool.starmap(
            simulate_cell,
            zipped_args,
        )

    pool.close()
    pool.join()

    pop_rep_off, pop_rep_on, pop_rep_on_ts = zip(*results)
    rep_off_mean = np.mean(pop_rep_off)
    rep_on_mean = np.mean(pop_rep_on)
    rep_on_ts_mean = [np.mean(k) for k in zip(*pop_rep_on_ts)]

    return rep_off_mean, rep_on_mean, rep_on_ts_mean


if __name__ == '__main__':

    topology = Topo([('P1', 'Z6'), ('P1', 'I6'), ('Z6', 'Rep'), ('I6', 'Rep')], {'Z6': 5, 'I6': 5}, 'P1')
    rep_off_mean, rep_on_mean, rep_on_ts_mean = simulate_pop(topology, Z_20, 8, 42)

    # print(rep_off_mean)
    # print(rep_on_mean)
    print(rep_on_ts_mean)
