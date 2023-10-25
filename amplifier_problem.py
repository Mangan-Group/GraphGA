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

class Amplifier:
    def __init__(
            self,
            promo_node: str,
            dose_specs: list,
            max_part: int,
            inhibitor: bool,
            DsRed_inhibitor: bool,
            num_dict: dict, 
            n_gen: int,
            probability_crossover: float, 
            probability_mutation: float,
            mutate_dose: bool=False,
            pop: bool=False,
            CI: float=None,
            Z_mat: np.ndarray=Z_20,
            num_processes: int=None, 
            ) -> None:
        
        self.promo_node = promo_node
        self.min_dose = dose_specs[0]
        self.max_dose = dose_specs[1]
        self.dose_interval = dose_specs[2]
        self.max_part = max_part
        self.inhibitor = inhibitor
        self.num_dict = num_dict
        self.n_gen = n_gen
        self.prob_crossover = probability_crossover
        self.prob_mutation = probability_mutation
        self.mutate_dose = mutate_dose
        self.pop = pop
        self.CI = CI
        self.num_processes = num_processes
        self.system_eqs = system_equations_pop
        
        if inhibitor:
            if DsRed_inhibitor:
                self.system_eqs = system_equations_DsRed_pop

        if pop:
            # set ref = simulation for 20-cell population
            self.ref = Ref_pop20
            # set Z = 20-cell population matrix np.array(20, 5) one row/cell, 1 columm/plasmid
            self.Z = Z_mat
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
        Z_row: np.ndarray =np.ones(5)
    ):

        t = np.arange(0, max_time + 1, 1)
        rep_on = odeint(
            system_equations_pop,
            np.zeros(topology.num_states * 2),
            t,
            args=('on', Z_row, topology,)
        )[-1, -1]
        return rep_on

    def simulate_pop(
        self, 
        topology: object, 
        max_time: int =42
    ):
        pop_rep_on = []
        nc = len(self.Z)
        zipped_args = list(zip([topology]*nc, [max_time]*nc, self.Z))
        for cell in range(0, nc):
            rep_on = self.simulate_cell(
                zipped_args[cell][0],
                zipped_args[cell][1],
                zipped_args[cell][2])
            pop_rep_on.append(rep_on)
        # with Pool(self.num_processes) as pool:
        #     pop_rep_on = pool.starmap(
        #         self.simulate_cell,
        #         zipped_args,
        #     )
        rep_on_mean = np.mean(pop_rep_on)
        return rep_on_mean

    def calc_ON_rel(self, topology, rep_on):
        reference_on = self.ref[topology.promo_node]['on']
        ON_rel = rep_on/reference_on
        return ON_rel

    @staticmethod
    def calc_FI(off, on):
        FI = on/off
        return FI

    @staticmethod
    def calc_FI_rel(ref_FI, FI):
        FI_rel = FI/ref_FI
        return FI_rel

    def func(
        self,
        topology: object
    ):
        
        rep_on = self.simulate(topology)
        ON_rel = self.calc_ON_rel(topology, rep_on)

        return -ON_rel
