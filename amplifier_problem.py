import numpy as np
from scipy.integrate import odeint
import pandas as pd
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
            obj_labels: list=["ON_rel"],
            max_time: int=42,
            single_cell_tracking: bool=False
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
        self.obj_labels = obj_labels
        self.max_time = max_time
        self.system_eqs = system_equations_pop
        
        if inhibitor:
            if DsRed_inhibitor:
                self.system_eqs = system_equations_DsRed_pop

        if pop:
            # set ref = simulation for 20-cell population
            self.ref = Ref_pop20
            # set Z = 20-cell population matrix np.array(20, 5) one row/cell, 1 columm/plasmid
            self.Z = Z_mat
            # set simulate function for population based on whether to track single cell
            # outputs
            if single_cell_tracking:
                self.simulate = self.simulate_pop_single_cell_tracking
                self.func = self.func_single_cell_tracking
            else:
                self.simulate = self.simulate_pop
                self.func = self.func_obj
        else:
            # set ref = simulation for single cell population
            self.ref = Ref
            self.Z = None
            # set simulate function for single cell
            self.simulate = self.simulate_cell
            self.func = self.func_obj

    def simulate_cell(
        self,
        topology: object,
        Z_row: np.ndarray =np.ones(5)
    ):

        max_time = self.max_time
        t = np.arange(0, max_time + 1, 1)
        rep_on = odeint(
            self.system_eqs,
            np.zeros(topology.num_states * 2),
            t,
            args=('on', Z_row, topology,)
        )[-1, -1]
        return rep_on

    def simulate_pop_single_cell_tracking(
        self, 
        topology: object, 
    ):
        pop_rep_on = []
        nc = len(self.Z)
        zipped_args = list(zip([topology]*nc, self.Z))
        for cell in range(0, nc):
            rep_on = self.simulate_cell(
                zipped_args[cell][0],
                zipped_args[cell][1]
            )
            pop_rep_on.append(rep_on)

        pop_ON_rel = self.calc_all_cell_metrics(
            topology, pop_rep_on
        )
        all_cells_dict = {"Topology": topology, 
                          "ON_rel for each cell": [pop_ON_rel],
                          "Rep ON state for each cell": [pop_rep_on]}

        rep_on_mean = np.mean(pop_rep_on)

        return rep_on_mean, all_cells_dict
    
    def simulate_pop(
        self, 
        topology: object, 
    ):
        pop_rep_on = []
        nc = len(self.Z)
        zipped_args = list(zip([topology]*nc, self.Z))
        for cell in range(0, nc):
            rep_on = self.simulate_cell(
                zipped_args[cell][0],
                zipped_args[cell][1]
            )
            pop_rep_on.append(rep_on)

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

    def func_obj(
        self,
        topology: object
    ):
        
        rep_on = self.simulate(topology)
        ON_rel = self.calc_ON_rel(topology, rep_on)

        return -ON_rel
    
    def func_single_cell_tracking(
        self,
        topology: object
    ):
        
        rep_on, all_cells_dict  = self.simulate(topology)
        ON_rel = self.calc_ON_rel(topology, rep_on)

        return [-ON_rel, all_cells_dict]

    def calc_all_cell_metrics(
                self, topology,
                pop_rep_on
    ):
        
        pop_ON_rel = []
        for i in range(len(pop_rep_on)):
            ON_rel = self.calc_ON_rel(topology, pop_rep_on[i])
            pop_ON_rel.append(ON_rel)

        return pop_ON_rel
