import numpy as np
from scipy.integrate import odeint
from multiprocessing import Pool
from get_system_equations_pop import (
    system_equations_pop,
    system_equations_DsRed_pop
)
from load_files_pop import (
    Ref,
    Z_20,
    Ref_pop20,
)

class SignalConditioner:
    def __init__(
            self,
            promo_node: str,
            # max_part: int, 
            min_dose: int, 
            max_dose: int, 
            dose_interval: int, 
            inhibitor: bool,
            DsRed_inhibitor: bool,
            num_dict: dict, 
            n_gen: int,
            pop: bool=False,
            num_processes: int=None, 
            ) -> None:
        
        self.promo_node = promo_node
        # self.max_part = max_part
        self.min_dose = min_dose
        self.max_dose = max_dose
        self.dose_interval = dose_interval
        self.inhibitor = inhibitor
        self.num_dict = num_dict
        self.n_gen = n_gen
        self.num_processes = num_processes
        self.system_eqs = system_equations_pop
        
        if inhibitor:
            if DsRed_inhibitor:
                self.system_eqs = system_equations_DsRed_pop

        if pop:
            # set ref = simulation for 20-cell population
            self.ref = Ref_pop20
            # set Z = 20-cell population matrix np.array(20, 5) one row/cell, 1 columm/plasmid
            self.Z = Z_20
            # set simulate function for population using multiprocessing
            self.simulate = self.simulate_pop
        else:
            # set ref = simulation for single cell population
            self.ref = Ref
            self.Z = None
            # set simulate function for single cell
            self.simulate = self.simulate_cell

    @staticmethod
    def simulate_cell(
        topology: object,
        max_time: int =42,
        Z_row: np.ndarray = np.ones(5)
    ):

        t = np.arange(0, max_time + 1, 1)
        rep_off = odeint(
            system_equations_pop,
            np.zeros(topology.num_states * 2), 
            t,
            args=('off', Z_row, topology,)
        )[-1, -1]
        rep_on = odeint(
            system_equations_pop,
            np.zeros(topology.num_states * 2),
            t, 
            args=('on', Z_row, topology,)
        )[-1, -1]
        return rep_off, rep_on

    def simulate_pop(
        self, 
        topology: object, 
        max_time: int =42
    ):

        nc = len(self.Z)
        zipped_args = list(zip([topology]*nc, [max_time]*nc, self.Z))
        with Pool(self.num_processes) as pool:
            results = pool.starmap(
                self.simulate_cell,
                zipped_args,
            )
        pop_rep_off, pop_rep_on = zip(*results)
        rep_off_mean = np.mean(pop_rep_off)
        rep_on_mean = np.mean(pop_rep_on)
        return rep_off_mean, rep_on_mean

    def calc_ON_rel(self, topology, rep_on):
        reference_on = self.ref[topology.promo_node]['on']
        ON_rel = rep_on/reference_on
        return ON_rel

    @staticmethod
    def calc_FI(off, on):
        FI = on/off
        return FI

    def calc_FI_rel(self, topology, FI):
        reference_off = self.ref[topology.promo_node]['off']
        reference_on = self.ref[topology.promo_node]['on']
        FI_ref = self.calc_FI(reference_off, reference_on)
        FI_rel = FI/FI_ref
        return FI_rel

    def func(
        self,
        topology: object
    ):
        
        rep_off, rep_on = self.simulate(topology)
        ON_rel = self.calc_ON_rel(topology, rep_on)
        FI_sc = self.calc_FI(rep_off, rep_on)
        FI_rel = self.calc_FI_rel(topology, FI_sc)

        return [-ON_rel, -FI_rel]
