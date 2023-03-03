import numpy as np
from scipy.integrate import odeint
from multiprocessing import Pool
from get_system_equations_pop import system_equations_pop
from load_files_pop import (
    Ref,
    Z_20,
    Ref_pop20,
)

class Amplifier:
    def __init__(
            self,
            promo_node: str,
            max_part: int, 
            min_dose: int, 
            max_dose: int, 
            dose_interval: int, 
            inhibitor: bool,
            num_dict: dict, 
            n_gen: int,
            pop: bool=False,
            num_processes: int=None, 
            combinatorial: bool=False,
            ) -> None:
        
        self.promo_node = promo_node
        self.max_part = max_part
        self.min_dose = min_dose
        self.max_dose = max_dose
        self.dose_interval = dose_interval
        self.inhibitor = inhibitor
        self.num_dict = num_dict
        self.n_gen = n_gen
        self.num_processes = num_processes
        self.combinatorial = combinatorial

        if pop:
            # set ref = simulation for 20-cell population
            self.ref = Ref_pop20
            # set Z = 20-cell population matrix np.array(20, 5) one row/cell, 1 columm/plasmid
            self.Z = Z_20
            # set simulate function for population using multiprocessing
            self.simulate = simulate_pop
        else:
            # set ref = simulation for single cell population
            self.ref = Ref
            self.Z = None
            # set simulate function for single cell
            self.simulate = simulate_cell

@staticmethod
def simulate_cell(
    topology: object,
    max_time: int =42,
    Z_row: np.ndarray = np.ones(5)
):

    t = np.arange(0, max_time + 1, 1)
    rep_on = odeint(system_equations_pop, np.zeros(topology.num_states * 2), t, args=('on', Z_row, topology,))[-1, -1]
    return rep_on

def simulate_pop(
    self, 
    topology: object, 
    max_time: int =42
):

    nc = len(self.Z)
    zipped_args = list(zip([topology]*nc, self.Z, [max_time]*nc))
    with Pool(self.num_processes) as pool:
        pop_rep_on = pool.starmap(
            simulate_cell,
            zipped_args,
        )
    rep_on_mean = np.mean(pop_rep_on)
    return rep_on_mean

def func(
    self,
    topology: object
):
    return -self.simulate(topology) / self.ref[topology.promo_node]['on']
