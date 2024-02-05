import numpy as np
import pickle
from scipy.integrate import odeint

with open("promo.pkl", "rb") as fid:
    promo = pickle.load(fid)
# with open("Z_200.npy", 'rb') as fid:
#     Z_200 = np.load(fid)
Z_path = "/Users/kdreyer/Documents/Github/GraphGA/Z_matrices/"
with open(Z_path + "Z_mat_20_cell0.npy", 'rb') as fid:
    Z_20 = np.load(fid)

def reference(y,t,k_end, Z):
    return np.array([k_end*Z - 2.7*y[0],
                    y[0] - 0.029*y[1]])


def simulate_reference(Z, filename):
    Ref = dict()
    for k in list(promo.keys())[:3]:
        ref_off = []
        ref_on = []
        for i in range(0, len(Z)):
            off = odeint(reference, np.zeros(2), np.arange(0, 42 + 1), args=(promo[k]['off']*promo['k_txn'], Z[i, 0]))[-1,-1]
            on = odeint(reference, np.zeros(2), np.arange(0, 42 + 1), args=(promo[k]['on']*promo['k_txn'], Z[i, 0]))[-1, -1]
            ref_off.append(off)
            ref_on.append(on)
        ref_off_mean = np.mean(ref_off)
        ref_on_mean = np.mean(ref_on)
        Ref.update({k: {'off': ref_off_mean, 'on': ref_on_mean, 'fi': ref_on_mean/ref_off_mean}})

    with open(filename, "wb") as fid:
        pickle.dump(Ref, fid)
    
    return Ref
    
# filename_200 = "Ref_pop200.pkl"
filename_20 = "Ref_pop20.pkl"

# Ref_200 = simulate_reference(Z_200, filename_200)
Ref_20 = simulate_reference(Z_20, filename_20)

print(Ref_20)
