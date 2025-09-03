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
    Z_200,
    Ref_pop200,
)
# create dictionary of reference simulations
# for different population models for setting
# ref attribute
pop_model_ref = {"20 cell": Ref_pop20,
                 "200 cell": Ref_pop200,
}

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
            # this argument isn't used here but
            # is set in run_test_case for the pulse 
            Ref_pop: dict=None, 
            num_processes: int=None,
            obj_labels: list=["ON_rel"],
            max_time: int=42,
            single_cell_tracking: bool=False
    ) -> None:
        
        # set attributes based on input arguments
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
            # change system equations if using DsRed
            # inhibitors in design space
            if DsRed_inhibitor:
                self.system_eqs = system_equations_DsRed_pop

        if pop:
            # set Z matrix and reference simulation according
            # to Z matrix
            self.Z = Z_mat
            self.ref = pop_model_ref[str(len(self.Z))+" cell"]

            # set simulate function for population based
            # on whether to track single cell outputs
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
        Z_row: np.ndarray=np.ones(5)
    ):
        """Solves the ODEs for a given 
        topology for a single cell."""

        # set simulation time (h)
        max_time = self.max_time
        t = np.arange(0, max_time + 1, 1)
        # solve ODEs in the ON state for
        # pEnd to calculate amplifier objective
        # funcion (save final time point reporter
        # expression)
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
        """Solves the ODEs for a given topology for
        each cell in the population using simulate_cell(),
        and calculates ON_rel metric for each cell in the 
        population."""

        pop_rep_on = []
        # get number of cells in population
        nc = len(self.Z)
        zipped_args = list(zip([topology]*nc, self.Z))
        # simulate each cell in the population with
        # corresponding Z matrix row
        for cell in range(0, nc):
            rep_on = self.simulate_cell(
                zipped_args[cell][0],
                zipped_args[cell][1]
            )
            # save simulation output (reporter
            # expression) for each cell in
            # pop_rep_on list
            pop_rep_on.append(rep_on)

        # calculate ON_rel metric and add to
        # dict along with reporter ON state
        # for each cell in the population
        pop_ON_rel = self.calc_all_cell_metrics(
            topology, pop_rep_on
        )
        all_cells_dict = {"Topology": topology, 
                          "ON_rel for each cell": [pop_ON_rel],
                          "Rep ON state for each cell": [pop_rep_on]}

        # calculate mean reporter expression
        # across cells in the population
        rep_on_mean = np.mean(pop_rep_on)
        return rep_on_mean, all_cells_dict
    

    def simulate_pop(
        self, 
        topology: object, 
    ):
        """Solves the ODEs for a given topology
        for each cell in the population using 
        simulate_cell()."""

        pop_rep_on = []
        # get number of cells in population
        nc = len(self.Z)
        zipped_args = list(zip([topology]*nc, self.Z))
        # simulate each cell in the population with
        # corresponding Z matrix row
        for cell in range(0, nc):
            rep_on = self.simulate_cell(
                zipped_args[cell][0],
                zipped_args[cell][1]
            )
            # save simulation output (reporter
            # expression) for each cell in
            # pop_rep_on list
            pop_rep_on.append(rep_on)

        # calculate mean reporter expression
        # across cells in the population
        rep_on_mean = np.mean(pop_rep_on)
        return rep_on_mean


    def calc_ON_rel(
        self,
        topology: object,
        rep_on: float
    ):
        """Calculates ON_rel metric for the given
        topology."""

        reference_on = self.ref[topology.promo_node]['on']
        ON_rel = rep_on/reference_on
        return ON_rel


    def func_obj(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function and calculates ON_rel metric."""

        rep_on = self.simulate(topology)
        ON_rel = self.calc_ON_rel(topology, rep_on)
        # return negative ON_rel for minimization in GA
        # (actually want to maximize ON_rel)
        return -ON_rel
    

    def func_single_cell_tracking(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function, calculates ON_rel metric,
        and saves dict of ON_rel and reporter ON state
        for each cell in the population (relevant for
        population model only)."""

        rep_on, all_cells_dict  = self.simulate(topology)
        ON_rel = self.calc_ON_rel(topology, rep_on)

        return [-ON_rel, all_cells_dict]


    def calc_all_cell_metrics(
        self,
        topology: object,
        pop_rep_on: list
    ):
        """Calculates ON_rel metric for each cell in the
        population (relevant for population model only)."""

        pop_ON_rel = []
        for i in range(len(pop_rep_on)):
            ON_rel = self.calc_ON_rel(topology, pop_rep_on[i])
            pop_ON_rel.append(ON_rel)

        return pop_ON_rel
