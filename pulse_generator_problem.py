import numpy as np
from scipy.integrate import odeint
from scipy.signal import find_peaks, peak_prominences
from scipy.stats.mstats import gmean
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
            mean: str="arithmetic",
            Z_mat: np.ndarray=Z_20,
            Ref_pop: dict=None, 
            num_processes: int=None,
            obj_labels: list=["t_pulse (hr)",
                            "prominence_rel"
            ],
            max_time: float=42,
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
        self.mean = mean,
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
            # if reference is specified (e.g. corresponding
            # to Z matrix sample, not one of the imported Z 
            # matrices & Ref, use that as Ref & Z matrix
            if Ref_pop is not None:
                self.ref = Ref_pop
            else:
                self.ref = pop_model_ref[str(len(self.Z))+" cell"]

            # set simulate function for population based
            # on whether to track single cell outputs
            if single_cell_tracking:
                self.simulate = self.simulate_pop_single_cell_tracking

                # specify objective function based on objective
                # labels argument
                if len(self.obj_labels) == 3:
                    self.func = self.func_single_cell_tracking_3obj
                elif ("t_pulse" in '\t'.join(self.obj_labels) and
                      "frac_pulse" in '\t'.join(self.obj_labels)):
                    self.func = self.func_single_cell_tracking_frac_t_pulse
                    
                elif ("t_pulse" in '\t'.join(self.obj_labels) and
                      len(self.obj_labels) == 3):
                    self.func = self.func_single_cell_tracking_frac_pulse_3obj

                elif ("t_pulse" in '\t'.join(self.obj_labels) and
                      "prominence_rel" in '\t'.join(self.obj_labels)):
                    self.func = self.func_single_cell_tracking_t_pulse

                elif "peak_rel" in '\t'.join(self.obj_labels):
                    self.func = self.func_single_cell_tracking_peak_rel

                elif ("frac_pulse" in '\t'.join(self.obj_labels) and 
                      "prominence_rel" in '\t'.join(self.obj_labels)):
                    self.func = self.func_single_cell_tracking_frac_pulse
                    
            else:
                self.simulate = self.simulate_pop

                # specify objective function based on objective
                # labels argument
                if (len(self.obj_labels) == 3 and 
                    "peak_rel" in '\t'.join(self.obj_labels)):
                    self.func = self.func_obj_3obj

                elif ("frac_pulse" in '\t'.join(self.obj_labels) and
                      len(self.obj_labels) == 3):
                    self.simulate = self.simulate_pop_single_cell_tracking
                    self.func = self.func_obj_frac_pulse_3obj

                elif ("t_pulse" in '\t'.join(self.obj_labels)  and
                      "frac_pulse" in '\t'.join(self.obj_labels)):
                    self.simulate = self.simulate_pop_single_cell_tracking
                    self.func = self.func_obj_frac_t_pulse

                elif ("t_pulse" in '\t'.join(self.obj_labels) and
                      "prominence_rel" in '\t'.join(self.obj_labels)):
                    self.func = self.func_obj_t_pulse

                elif "peak_rel" in '\t'.join(self.obj_labels):
                    self.func = self.func_obj_peak_rel

                elif ("frac_pulse" in '\t'.join(self.obj_labels) and
                      "prominence_rel" in '\t'.join(self.obj_labels)):
                    self.simulate = self.simulate_pop_single_cell_tracking
                    self.func = self.func_obj_frac_pulse
        else:
            # set ref = simulation for single cell population
            self.ref = Ref
            self.Z = None
            # set simulate function for single cell
            self.simulate = self.simulate_cell

            # specify objective function based on objective
            # labels argument
            if len(self.obj_labels) == 3:
                self.func = self.func_obj_3obj
            elif "t_pulse" in '\t'.join(self.obj_labels):
                self.func = self.func_obj_t_pulse
            elif "peak_rel" in '\t'.join(self.obj_labels):
                self.func = self.func_obj_peak_rel


    def simulate_cell(
        self,
        topology: object,
        Z_row: np.ndarray = np.ones(5)
    ):
        """Solves the ODEs for a given 
        topology for a single cell."""
        
        # set simulation time (h)
        max_time = self.max_time
        t = np.arange(0, max_time + 1, 1)
        # solve ODEs in the ON state for
        # pEnd to calculate pulse objective
        # funcion (save reporter expression
        # time series and t for t_pulse
        # metric calculation)
        rep_on_ts = odeint(
            self.system_eqs,
            np.zeros(topology.num_states * 2),
            t, 
            args=('on', Z_row, topology,)
        )[:, -1]
        return t, rep_on_ts


    def simulate_pop_single_cell_tracking(
        self, 
        topology: object, 
    ):
        """Solves the ODEs for a given topology for
        each cell in the population using simulate_cell(),
        and calculates reporter_rel for each cell in the 
        population, mean reporter_rel at each time point
        in t."""

        rep_on_ts_all = []
        # get number of cells in population
        nc = len(self.Z)
        zipped_args = list(zip([topology]*nc, self.Z))
        # simulate each cell in the population with
        # corresponding Z matrix row
        for cell in range(0, nc):
            t, rep_on_ts = self.simulate_cell(                
                zipped_args[cell][0],
                zipped_args[cell][1]
            )
            # save simulation output (reporter
            # expression) for each cell in
            # rep_on_ts_all list
            rep_on_ts_all.append(rep_on_ts)

        rep_on_ts_rel_all = self.calc_all_cell_rep_rel(
            topology,
            rep_on_ts_all
        )
        # calculate mean reporter expression at each time in t
        rep_on_ts_means = [np.mean(k) for k in zip(*rep_on_ts_all)]
        # calculate reporter_rel from mean reporter at each time
        # in t
        rep_on_ts_rel_mean = self.calc_rep_rel(
            topology,
            rep_on_ts_means)

        all_cells_dict = {"Topology": topology, 
                          "Rep_rel time series for each cell": rep_on_ts_rel_all,
                          "Rep_rel time series mean": rep_on_ts_rel_mean}

        return t, rep_on_ts_means, all_cells_dict, rep_on_ts_all
    

    def simulate_pop(
        self, 
        topology: object, 
    ):
        """Solves the ODEs for a given topology
        for each cell in the population using 
        simulate_cell()."""

        rep_on_ts_all = []
        # get number of cells in population
        nc = len(self.Z)
        zipped_args = list(zip([topology]*nc, self.Z))
        # simulate each cell in the population with
        # corresponding Z matrix row
        for cell in range(0, nc):
            t, rep_on_ts = self.simulate_cell(                
                zipped_args[cell][0],
                zipped_args[cell][1]
            )
            # save simulation output (reporter
            # expression) for each cell in
            # rep_on_ts_all list
            rep_on_ts_all.append(rep_on_ts)
        # calculate mean reporter expression at each time in t
        rep_on_ts_means = [np.mean(k) for k in zip(*rep_on_ts_all)]

        return t, rep_on_ts_means


    def calc_rep_rel(
        self,
        topology: object,
        rep_on_ts: list
    ):
        """Calculates relative reporter expression
        for the given topology."""

        reference_on = self.ref[topology.promo_node]['on']
        rep_on_ts_rel = [i/reference_on for i in rep_on_ts]

        return rep_on_ts_rel


    @staticmethod
    def calc_peak_rel(
        rep_on_ts_rel: list
    ):
        """Calculates peak_rel metric
        from the given time series
        list of reporter_rel."""

        return max(rep_on_ts_rel)


    @staticmethod
    def calc_prominence_rel(
        rep_on_ts_rel: list,
        peak_rel: float
    ):
        """Calculates the prominence_rel metric
        from the given time series list of 
        reporter_rel and peak_rel."""

        peaks_rep, _ = find_peaks(
            rep_on_ts_rel,
            prominence=0.1*peak_rel
        )
        prominence_rep_list = peak_prominences(
            rep_on_ts_rel,
            peaks_rep
        )[0]
        if len(prominence_rep_list) == 0:
            prominence_rel = 0
        else:
            prominence_rel = prominence_rep_list[0]
        return prominence_rel
    

    def calc_t_pulse(
        self,
        t: np.ndarray,
        rep_on_ts_rel: list,
        peak_rel: float,
        prominence_rel: float
    ):
        """Calculates the t_pulse metric
        from the given time series list of 
        reporter_rel, peak_rel, and array
        of t values."""

        if prominence_rel != 0:
            idx_peak_rel = rep_on_ts_rel.index(peak_rel)
            t_pulse = t[idx_peak_rel]
        else:
            t_pulse = self.max_time
        
        return t_pulse


    def calc_frac_pulse(
        self,
        topology: object,
        rep_on_ts_all: list
    ):
        """Calculates the frac_pulse metric
        from the given time series list of 
        reporter_rel."""

        count = 0
        pulse_cell_ts = []
        for cell_ts in rep_on_ts_all:
            rep_on_ts_rel = self.calc_rep_rel(topology, cell_ts)
            peak_rel = self.calc_peak_rel(rep_on_ts_rel)
            prominence_cell = self.calc_prominence_rel(rep_on_ts_rel, peak_rel)
            if prominence_cell != 0:
                pulse_cell_ts.append(cell_ts)
                count += 1
        pulse_cell_ts_means = [np.mean(k) for k in zip(*pulse_cell_ts)]
        frac_pulse = count/len(rep_on_ts_all)
        return frac_pulse, pulse_cell_ts_means


    def func_obj_peak_rel(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function and calculates peak_rel 
        and prominence_rel metrics."""

        _, rep_on_ts = self.simulate(topology)
        rep_on_ts_rel = self.calc_rep_rel(topology, rep_on_ts)
        # want peak_rel >= 0.25
        peak_rel = self.calc_peak_rel(rep_on_ts_rel)
        # want prominence_rel >= 0.2
        prominence_rel = self.calc_prominence_rel(rep_on_ts_rel, peak_rel)
        # return negative peak_rel and prominence_rel for 
        # minimization in GA (actually want to maximize
        # peak_rel and prominence_rel)
        return [-peak_rel, -prominence_rel]


    def func_obj_t_pulse(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function and calculates t_pulse and 
        prominence_rel metrics."""

        t, rep_on_ts = self.simulate(topology)
        rep_on_ts_rel = self.calc_rep_rel(topology, rep_on_ts)
        peak_rel = self.calc_peak_rel(rep_on_ts_rel)
        prominence_rel = self.calc_prominence_rel(rep_on_ts_rel, peak_rel)

        t_pulse = self.calc_t_pulse(
            t, rep_on_ts_rel, peak_rel, prominence_rel)

        # return negative prominence_rel for 
        # minimization in GA (actually want to
        # maximize prominence_rel)
        return [t_pulse, -prominence_rel]
    

    def func_obj_3obj(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function and calculates t_pulse, peak_rel, 
        and prominence_rel metrics."""

        t, rep_on_ts = self.simulate(topology)
        rep_on_ts_rel = self.calc_rep_rel(topology, rep_on_ts)
        peak_rel = self.calc_peak_rel(rep_on_ts_rel)
        prominence_rel = self.calc_prominence_rel(
            rep_on_ts_rel,
            peak_rel
        )
        t_pulse = self.calc_t_pulse(
            t, rep_on_ts_rel, peak_rel, prominence_rel
        )
        # return negative peak_rel and prominence_rel for 
        # minimization in GA (actually want to maximize
        # peak_rel and prominence_rel)
        return [t_pulse, -peak_rel, -prominence_rel]
    

    def func_obj_frac_pulse(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function and calculates frac_pulse and 
        prominence_rel (for pulses) metrics."""

        _, _, _, rep_on_ts_all = self.simulate(topology)
        frac_pulse, pulse_cell_ts_means = self.calc_frac_pulse(
            topology,
            rep_on_ts_all
        )
        if len(pulse_cell_ts_means) > 0:
            pulse_rep_ts_rel = self.calc_rep_rel(topology, pulse_cell_ts_means)
            pulse_peak_rel = self.calc_peak_rel(pulse_rep_ts_rel)
            pulse_prom_rel = self.calc_prominence_rel(
                pulse_rep_ts_rel,
                pulse_peak_rel
            )
        else:
            pulse_prom_rel = 0
        # return negative frac_pulse and prominence_rel for 
        # minimization in GA (actually want to maximize
        # frac_pulse and prominence_rel)
        return [-frac_pulse, -pulse_prom_rel]
   

    def func_obj_frac_t_pulse(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function and calculates frac_pulse and 
        t_pulse metrics."""

        t, rep_on_ts_means, _, rep_on_ts_all = self.simulate(topology)
        frac_pulse, pulse_cell_ts_means = self.calc_frac_pulse(
            topology,
            rep_on_ts_all
        )
        if len(pulse_cell_ts_means) > 0:
            pulse_rep_ts_rel = self.calc_rep_rel(
                topology,
                pulse_cell_ts_means
            )
            pulse_peak_rel = self.calc_peak_rel(pulse_rep_ts_rel)
            pulse_prom_rel = self.calc_prominence_rel(
                pulse_rep_ts_rel,
                pulse_peak_rel
            )
            # t_pulse for PULSE cells
            t_pulse = self.calc_t_pulse(
                t, 
                pulse_rep_ts_rel, 
                pulse_peak_rel, 
                pulse_prom_rel
            )
        else:
            t_pulse = 126

        # return negative frac_pulse for minimization
        # in GA (actually want to maximize frac_pulse)
        return [-frac_pulse, t_pulse]


    def func_obj_frac_pulse_3obj(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function and calculates frac_pulse t_pulse,
        peak_rel, and prominence_rel (for pulses) metrics."""

        t, rep_on_ts_means, _, rep_on_ts_all = self.simulate(topology)
        frac_pulse, pulse_cell_ts_means = self.calc_frac_pulse(topology, rep_on_ts_all)
        if len(pulse_cell_ts_means) > 0:
            pulse_rep_ts_rel = self.calc_rep_rel(topology, pulse_cell_ts_means)
            pulse_peak_rel = self.calc_peak_rel(pulse_rep_ts_rel)
            pulse_prom_rel = self.calc_prominence_rel(pulse_rep_ts_rel, pulse_peak_rel)
            # t_pulse for PULSE cells
            t_pulse = self.calc_t_pulse(
                t, pulse_rep_ts_rel, pulse_peak_rel, pulse_prom_rel)
        else:
            pulse_prom_rel = 0
            t_pulse = 126
            
        # return negative frac_pulse and prominence_rel for 
        # minimization in GA (actually want to maximize
        # frac_pulse and prominence_rel)
        return [-frac_pulse, t_pulse, -pulse_prom_rel]


    def func_single_cell_tracking_peak_rel(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function, calculates peak_rel and 
        prominence_rel metrics, and saves dict of metrics
        and reporter time series for each cell in the
        population (relevant for population model only)."""

        _, rep_on_ts, all_cells_dict, _ = self.simulate(topology)
        rep_on_ts_rel = self.calc_rep_rel(topology, rep_on_ts)
        peak_rel = self.calc_peak_rel(rep_on_ts_rel)
        prominence_rel = self.calc_prominence_rel(rep_on_ts_rel, peak_rel)

        return [[-peak_rel, -prominence_rel], all_cells_dict]


    def func_single_cell_tracking_t_pulse(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function, calculates t_pulse and 
        prominence_rel metrics, and saves dict of metrics
        and reporter time series for each cell in the
        population (relevant for population model only)."""

        t, rep_on_ts, all_cells_dict, _ = self.simulate(topology)
        rep_on_ts_rel = self.calc_rep_rel(topology, rep_on_ts)
        peak_rel = self.calc_peak_rel(rep_on_ts_rel)
        prominence_rel = self.calc_prominence_rel(rep_on_ts_rel, peak_rel)

        t_pulse = self.calc_t_pulse(
            t, rep_on_ts_rel, peak_rel, prominence_rel)

        return [[t_pulse, -prominence_rel], all_cells_dict]
    

    def func_single_cell_tracking_3obj(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function, calculates t_pulse, peak_rel, 
        and prominence_rel metrics, and saves dict of metrics
        and reporter time series for each cell in the
        population (relevant for population model only)."""

        t, rep_on_ts, all_cells_dict, _ = self.simulate(topology)
        rep_on_ts_rel = self.calc_rep_rel(topology, rep_on_ts)
        peak_rel = self.calc_peak_rel(rep_on_ts_rel)
        prominence_rel = self.calc_prominence_rel(rep_on_ts_rel, peak_rel)

        t_pulse = self.calc_t_pulse(
            t, rep_on_ts_rel, peak_rel, prominence_rel
        )

        return [[t_pulse, -peak_rel, -prominence_rel], all_cells_dict]


    def func_single_cell_tracking_frac_pulse(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function, calculates frac_pulse and 
        prominence_rel (for pulses) metrics, and saves
        dict of metrics and reporter time series for 
        each cell in the population (relevant for 
        population model only)."""

        _, _, all_cells_dict, rep_on_ts_all = self.simulate(topology)
        frac_pulse, pulse_cell_ts_means = self.calc_frac_pulse(topology, rep_on_ts_all)
        if len(pulse_cell_ts_means) > 0:
            pulse_rep_ts_rel = self.calc_rep_rel(topology, pulse_cell_ts_means)
            pulse_peak_rel = self.calc_peak_rel(pulse_rep_ts_rel)
            pulse_prom_rel = self.calc_prominence_rel(pulse_rep_ts_rel, pulse_peak_rel)
        else:
            pulse_prom_rel = 0

        return [[-frac_pulse, -pulse_prom_rel], all_cells_dict]


    def func_single_cell_tracking_frac_t_pulse(
        self,
        topology: object
    ):
        """Simulates the given topology using specified
        simulate function, calculates frac_pulse and 
        t_pulse metrics, and saves dict of metrics
        and reporter time series for each cell in the
        population (relevant for population model only)."""

        t, rep_on_ts_means, all_cells_dict, rep_on_ts_all  = self.simulate(topology)
        frac_pulse, pulse_cell_ts_means = self.calc_frac_pulse(topology, rep_on_ts_all)
        if len(pulse_cell_ts_means) > 0:
            pulse_rep_ts_rel = self.calc_rep_rel(topology, pulse_cell_ts_means)
            pulse_peak_rel = self.calc_peak_rel(pulse_rep_ts_rel)
            pulse_prom_rel = self.calc_prominence_rel(pulse_rep_ts_rel, pulse_peak_rel)
            rep_on_ts_rel = self.calc_rep_rel(topology, rep_on_ts_means)
            peak_rel = self.calc_peak_rel(rep_on_ts_rel)
            prominence_rel = self.calc_prominence_rel(rep_on_ts_means, peak_rel)
            t_pulse = self.calc_t_pulse(
                t, rep_on_ts_rel, peak_rel, prominence_rel)
        else:
            pulse_prom_rel = 0
            t_pulse = 126

        return [[-frac_pulse, t_pulse], all_cells_dict]


    def calc_all_cell_rep_rel(
        self,
        topology: object,
        rep_on_ts_all: list
    ):
        """Calculates relative reporter expression
        for each cell in the population for the 
        given topology."""

        rep_on_ts_rel_all = []
        for i in range(len(rep_on_ts_all)):
            rep_on_ts_rel = self.calc_rep_rel(
                topology, rep_on_ts_all[i]
            )
            rep_on_ts_rel_all.append(
                rep_on_ts_rel
            )
        return rep_on_ts_rel_all
