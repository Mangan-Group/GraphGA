import numpy as np
from scipy.integrate import odeint
from scipy.signal import find_peaks, peak_prominences
import pandas as pd
from define_circuit import Topo
from get_system_equations_pop import (
    system_equations_pop,
    system_equations_DsRed_pop
)
from load_files_pop import (
    Ref
)


######################################################
######### function to simulate each topology #########
######################################################

## this variable will need to be set to False for a
## combinatorial search with the weaker inhibitors
DsRed_inhibitor=True,

if DsRed_inhibitor:
        system_eqs = system_equations_DsRed_pop
else:
        system_eqs = system_equations_pop

def simulate_ciruit(
    topology: object,
    system_eqs,
    max_time_off: int=42,
    max_time_on: int=126,
):
    """
    Simulates one circuit topology for the combinatorial
    search

    Parameters
    ----------

    topology
        an instance of the Topo class defining the topology
        to be simulated
    
    system_eqs
        a function that compiles the system of ODEs for the
        given topology
    
    max_time_off
        an int specifying the max time point for simulating
        the circuit in the OFF state
        the default value should be used here - the signal
        conditioner requires the OFF state 42h reporter
        expression for calculating FI_rel
    
    max_time_on
        an int specifying the max time point for simulating
        the circuit in the ON state.
        the default value should be used here - the pulse
        generator requires the full time series for 126h,
        since this time point was used in the GA

    Returns
    -------

    circuit_rep_vals
        a dict containing the neccessary simulation values
        for calculating the objectives for each topology.
        the dict includes the 42h rep_off value, th 42h
        rep_on value, and the full 126h rep_on time series
    """
    t_off = np.arange(0, max_time_off + 1, 1)
    rep_off = odeint(
        system_eqs,
        np.zeros(topology.num_states * 2), 
        t_off,
        args=('off', np.ones(5), topology)
    )[-1, -1]

    t_on = np.arange(0, max_time_on + 1, 1)
    rep_on_ts = odeint(
            system_eqs,
            np.zeros(topology.num_states * 2),
            t_on, 
            args=('on', np.ones(5), topology)
        )[:, -1]
    rep_on = rep_on_ts[42]

    circuit_rep_vals = {"rep_off": rep_off, "rep_on": rep_on, "rep_on_ts": rep_on_ts}

    return circuit_rep_vals


######################################################
### function to calculate metrics for each topology ###
######################################################

def calculate_metrics(topology:object,
                      circuit_rep_vals:dict,
                      max_time_on:int=126
):
    """
    Calculates ON_rel for the amplifier and signal
    conditioner, FI_rel for the signal conditoner,
    and peak_rel, prom_rel, and t_pulse for the
    pulse

    Parameters
    ----------

    topology
        an instance of the Topo class defining the topology
        to be simulated

    Returns
    -------

    circuit_metrics
        a dict containing ON_rel, FI_rel, peak_rel,
        prom_rel, and t_pulse for the given topology
    """
    
    # circuit_rep_vals = {"rep_off": rep_off, "rep_on": rep_on, rep_on_ts: rep_on_ts}
    t_on = np.arange(0, max_time_on + 1, 1)
    rep_off = circuit_rep_vals["rep_off"]
    rep_on = circuit_rep_vals["rep_on"]
    rep_on_ts = circuit_rep_vals["rep_on_ts"]

    ON_rel = calc_ON_rel(topology, rep_on)
    FI = calc_FI(rep_off, rep_on)
    FI_rel = calc_FI_rel(topology, FI)
    rep_on_ts_rel = calc_rep_rel(topology, rep_on_ts)
    peak_rel = calc_peak_rel(rep_on_ts_rel)
    prominence_rel = calc_prominence_rel(rep_on_ts_rel,
                                         peak_rel
    )
    t_pulse = calc_t_pulse(t_on, rep_on_ts_rel,
                           peak_rel,
                           prominence_rel
    )
    circuit_metrics = {"ON_rel": ON_rel,
                       "FI_rel": FI_rel,
                       "peak_rel": peak_rel, 
                       "prominence_rel": prominence_rel,
                       "t_pulse": t_pulse
    }
    return circuit_metrics



######################################################
######### functions to calculate each metric #########
######################################################
def calc_ON_rel(topology, rep_on):
    reference_on = Ref[topology.promo_node]['on']
    ON_rel = rep_on/reference_on
    return ON_rel


def calc_FI(rep_off, rep_on):
    FI = rep_on/rep_off
    return FI


def calc_FI_rel(topology, FI):
    FI_ref = Ref[topology.promo_node]['fi']
    FI_rel = FI/FI_ref
    return FI_rel


def calc_rep_rel(topology, rep_on_ts):
    reference_on = Ref[topology.promo_node]['on']
    rep_on_ts_rel = [i/reference_on for i in rep_on_ts]
    return rep_on_ts_rel


def calc_peak_rel(rep_on_ts_rel):
    return max(rep_on_ts_rel)


def calc_prominence_rel(rep_on_ts_rel, peak_rel):
    peaks_rep, _ = find_peaks(rep_on_ts_rel, prominence=0.1*peak_rel)
    prominence_rep_list = peak_prominences(rep_on_ts_rel, peaks_rep)[0]
    if len(prominence_rep_list) == 0:
        prominence_rel = 0
    else:
        prominence_rel = prominence_rep_list[0]
    return prominence_rel


def calc_t_pulse(
    t, rep_on_ts_rel, peak_rel, prominence_rel
):
    if prominence_rel != 0:
        idx_peak_rel = rep_on_ts_rel.index(peak_rel)
        t_pulse = t[idx_peak_rel]
    else:
        t_pulse = 0
    return t_pulse