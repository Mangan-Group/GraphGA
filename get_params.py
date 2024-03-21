import pickle
import numpy as np

promo = {}

promo['P0'] = {}
promo['P1'] = {}
promo['P2'] = {}
promo["P_exp_amp"] = {}
promo["P_exp_sc"] = {}

promo['P0']['off'] = 0.00001
promo['P0']['on'] = 0.42
promo['P1']['off'] = 0.42
promo['P1']['on'] = 0.93
promo['P2']['off'] = 0.93
promo['P2']['on'] = 3.5
promo["P_exp_amp"]["off"] = 0.66
promo["P_exp_amp"]["on"] = 3.2
promo["P_exp_sc"]["off"] = 0.45
promo["P_exp_sc"]["on"] = 2.0
promo['k_txn'] = 8.

with open('promo.pkl', 'wb') as fid:
    pickle.dump(promo, fid)

parts = {}
parts['Z1'] = np.array([0.08, 33., 0.036])
parts['Z2'] = np.array([2.5e-01, 5.4e+01, 1.8e-02])
parts['Z6'] = np.array([2.0e-02, 5.8e+01, 4.3e-02])
parts['Z7'] = np.array([1.1e-01, 4.6e+01, 2.5e-02])
parts['Z8'] = np.array([7.0e-02, 4.3e+01, 4.1e-02])
parts['Z9'] = np.array([0.46, 33., 0.096])
parts['Z10'] = np.array([1.0e-02, 3.1e+01, 3.7e-02])
parts['Z11'] = np.array([8.0e-02, 3.2e+01, 2.5e-02])
parts['Z12'] = np.array([0.15, 33., 0.065])
parts['Z13'] = np.array([4.0e-02, 4.1e+01, 1.2e-02])
parts['Z14'] = np.array([0.2, 30., 0.069])
parts['Z15'] = np.array([1.8e-01, 3.3e+01, 7.0e-03])
parts['I1'] = np.array([0.036])
parts['I2'] = np.array([0.018])
parts['I6'] = np.array([0.043])
parts['I7'] = np.array([0.025])
parts['I8'] = np.array([0.041])
parts['I9'] = np.array([0.096])
parts['I10'] = np.array([0.037])
parts['I11'] = np.array([0.025])
parts['I12'] = np.array([0.065])
parts['I13'] = np.array([0.012])
parts['I14'] = np.array([0.069])
parts['I15'] = np.array([0.007])

# with open('parts.pkl', 'wb') as fid:
#     pickle.dump(parts, fid)