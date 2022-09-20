# -*- coding: utf-8 -*-
"""
Created on Mon Sep 12 14:20:18 2022

@author: Katie_Dreyer
"""

import numpy as np
from load_files import Ref

class amplifier:
    
    def __init__(self):
        
        self.promo_node = 'P1'
        self.max_part = 3
        self.min_dose = 10
        self.max_dose = 75
        self.dose_interval = 5
        self.inhibitor = False
        self.n_obj = 1
        self.n_ieq_constr = 0
        self.num_circuit = 20
        self.n_gen = 15
        
    def objective(self, topology):
        
        rep_off, rep_on = topology.simulate()
        ON_rel = rep_on / Ref[topology.promo_node]['on']
        # FI_amp = rep_on / rep_off
        # objectives = np.array([-ON_rel, topology.num_states-1]) # second one is the number of regulators
        # objectives = np.array([-ON_rel, topology.graph.number_of_edges()]) # we can also try number of regulations
        objectives = -ON_rel
                
        return objectives
    
class signal_conditioner:
    
    def __init__(self):
        
        self.promo_node = 'P1'
        self.max_part = 3
        self.min_dose = 10
        self.max_dose = 75
        self.dose_interval = 5
        self.inhibitor = True
        self.n_obj = 2
        self.n_ieq_constr = 1
        self.num_circuit = 50
        self.n_gen = 5
        
    def objective(self, topology):
        
        rep_off, rep_on = topology.simulate()
        ON_rel = rep_on / Ref[topology.promo_node]['on']
        OFF_rel = rep_off / Ref[topology.promo_node]['off']
        
        objectives = np.array([-ON_rel, OFF_rel])
                
        return objectives

    def constr(self, topology):
        # rep_off, rep_on = topology.simulate()
        # ON_rel = rep_on / Ref[topology.promo_node]['on']
        return topology.graph.size()-7

