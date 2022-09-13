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
        self.max_part = 4
        self.min_dose = 10
        self.max_dose = 75
        self.inhibitor = False
        self.n_obj = 2
        self.num_circuit = 50
        self.n_gen = 5
        
    def objective(self, topology):
        
        rep_off, rep_on = topology.simulate()
        ON_rel = rep_on / Ref[topology.promo_node]['on']
        FI_amp = rep_on / rep_off
        
        objectives = -np.array([ON_rel, FI_amp])
                
        return objectives
    
class signal_conditioner:
    
    def __init__(self):
        
        self.promo_node = 'P1'
        self.max_part = 4
        self.min_dose = 10
        self.max_dose = 75
        self.inhibitor = True
        self.n_obj = 2
        self.num_circuit = 50
        self.n_gen = 5
        
    def objective(self, topology):
        
        rep_off, rep_on = topology.simulate()
        ON_rel = rep_on / Ref[topology.promo_node]['on']
        OFF_rel = rep_off / Ref[topology.promo_node]['off']
        
        objectives = np.array([-ON_rel, OFF_rel])
                
        return objectives   


