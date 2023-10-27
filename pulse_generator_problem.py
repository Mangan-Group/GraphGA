import numpy as np
from scipy.integrate import odeint
from scipy.signal import find_peaks, peak_prominences
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

class PulseGenerator:
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
            CI: list=None,
            Z_mat: np.ndarray=Z_20,
            num_processes: int=None,
            obj_labels: list=["peak_rel",
                "prominence_rel"
            ]
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

    def simulate_cell(
        self,
        topology: object,
        max_time: int =42,
        Z_row: np.ndarray = np.ones(5)
    ):

        t = np.arange(0, max_time + 1, 1)
        rep_on_ts = odeint(
            self.system_eqs,
            np.zeros(topology.num_states * 2),
            t, 
            args=('on', Z_row, topology,)
        )[:, -1]
        return rep_on_ts

    def simulate_pop(
        self, 
        topology: object, 
        max_time: int =42
    ):
        rep_on_ts_all = []
        nc = len(self.Z)
        zipped_args = list(zip([topology]*nc, [max_time]*nc, self.Z))
        for cell in range(0, nc):
            rep_on_ts = self.simulate_cell(                
                zipped_args[cell][0],
                zipped_args[cell][1],
                zipped_args[cell][2])
            rep_on_ts_all.append(rep_on_ts)
        rep_on_ts_means = [np.mean(k) for k in zip(*rep_on_ts_all)]
        # with Pool(self.num_processes) as pool:
        #     results = pool.starmap(
        #         self.simulate_cell,
        #         zipped_args,
        #     )

        # rep_on_ts_means = [np.mean(k) for k in zip(*results)]
        # rep_on_ts_all = results
        return rep_on_ts_means #, rep_on_ts_all

    def calc_rep_rel(self, topology, rep_on_ts):
        reference_on = self.ref[topology.promo_node]['on']
        rep_on_ts_rel = [i/reference_on for i in rep_on_ts]

        return rep_on_ts_rel

    @staticmethod
    def calc_peak_rel(rep_on_ts_rel):
        return max(rep_on_ts_rel)

    @staticmethod
    def calc_prominence_rel(rep_on_ts_rel, peak_rel):
        peaks_rep, _ = find_peaks(rep_on_ts_rel, prominence=0.1*peak_rel)
        prominence_rep_list = peak_prominences(rep_on_ts_rel, peaks_rep)[0]
        if len(prominence_rep_list) == 0:
            prominence_rel = 0
        else:
            prominence_rel = prominence_rep_list[0]
        return prominence_rel

    # def calc_num_pulses(self, topology, rep_on_ts_all):
    #     count = 0
    #     pulse_cell_ts = []

    #     for cell_ts in rep_on_ts_all:
    #         rep_on_ts_rel = self.calc_rep_rel(topology, cell_ts)
    #         peak_rel = self.calc_peak_rel(rep_on_ts_rel)
    #         prominence_cell = self.calc_prominence(rep_on_ts_rel, peak_rel)

    #         if prominence_cell != 0:
    #             pulse_cell_ts.append(rep_on_ts_rel)
    #             count += 1
    #     return count, pulse_cell_ts

    def func(
        self,
        topology: object
    ):
        rep_on_ts = self.simulate(topology)
        rep_on_ts_rel = self.calc_rep_rel(topology, rep_on_ts)
        # want peak_rel >= 0.25
        peak_rel = self.calc_peak_rel(rep_on_ts_rel)
        # want prominence_rel >= 0.2
        prominence_rel = self.calc_prominence_rel(rep_on_ts_rel, peak_rel)
        return [-peak_rel, -prominence_rel]
