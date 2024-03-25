import numpy as np
import pickle

# path_amp_const = "/Users/kdreyer/Documents/Github/GraphGA/hyperparameter_optimization/Results/Amplifier_single_cell_const_dose/"
# with open(path_amp_const + "pareto_front.pkl", "rb") as fid:
#     pareto_front = pickle.load(fid)
# print(pareto_front)

# path_amp_vary = "/Users/kdreyer/Documents/Github/GraphGA/hyperparameter_optimization/Results/Amplifier_single_cell_vary_dose/run0/"
# with open(path_amp_vary + "pareto_front.pkl", "rb") as fid:
#     pareto_front = pickle.load(fid)
# print(pareto_front)

# path_sig_cond = "/Users/kdreyer/Documents/Github/GraphGA/hyperparameter_optimization/Results/Signal_conditioner_single_cell/run0_ngen50/"
# with open(path_sig_cond + "pareto_front.pkl", "rb") as fid:
#     pareto_front = pickle.load(fid)
# print(pareto_front)

path_pulse = "/Users/kdreyer/Documents/Github/GraphGA/hyperparameter_optimization/Results/Pulse_generator_single_cell/"
with open(path_pulse + "pareto_front.pkl", "rb") as fid:
    pareto_front = pickle.load(fid)
print(pareto_front)