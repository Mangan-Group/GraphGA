import numpy as np
import pickle
from scipy.integrate import odeint
from load_Z_mat_samples import Z_mat_list_20, Z_mat_list_200

repo_path = "/Users/kdreyer/Documents/Github/GraphGA/"

# load promoter parameters
with open(repo_path+"promo.pkl", "rb") as fid:
    promo = pickle.load(fid)

# load Z
Z_path = "/Users/kdreyer/Documents/Github/GraphGA/Z_matrix_samples/"
with open(Z_path + "Z_mat_20_cell0.npy", 'rb') as fid:
    Z_20 = np.load(fid)

with open(Z_path + "Z_mat_200_cell.npy", 'rb') as fid:
    Z_200 = np.load(fid)


def reference(y,t,k_end, Z):
    """Defines the ODEs for the
    reference case, for the
    population model."""

    return np.array([k_end*Z - 2.7*y[0],
                    y[0] - 0.029*y[1]])


def simulate_reference(Z, filename):
    """Simulates the reference for the
    specified promoter parameters and saves
    the final time point, mean reporter 
    expression across the population of cells
    for the OFF and ON states."""
    Ref = dict()
    for k in list(promo.keys())[:5]:
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


def simulate_reference_time_series(promo_list, Z):
    Ref_all_cells = dict()
    for k in promo_list:
        ref_off_time_series = []
        ref_on_time_series = []
        for i in range(0, len(Z)):
            off_time_series = odeint(reference, np.zeros(2), np.arange(0, 46 + 1), args=(promo[k]['off']*promo['k_txn'], Z[i, 0]))[:,-1]
            on_time_series = odeint(reference, np.zeros(2), np.arange(0, 46 + 1), args=(promo[k]['on']*promo['k_txn'], Z[i, 0]))[:, -1]
            ref_off_time_series.append(off_time_series)
            ref_on_time_series.append(on_time_series)

        Ref_all_cells.update({k: {'off all cells': ref_off_time_series, 'on all cells': ref_on_time_series}})
    
    return Ref_all_cells


# simulate the reference with the 20-cell model and
# the 200-cell model 
# filename_20 = "Ref_pop20.pkl"
# Ref_20 = simulate_reference(Z_20, filename_20)

# filename_200 = "Ref_pop200.pkl"
# Ref_200 = simulate_reference(Z_200, filename_200)

# simulate the reference with the 10 manifestations
# of the 20-cell model
for i, z_mat in enumerate(Z_mat_list_20):
    z_mat = Z_mat_list_20[i]
    file_name = "Z_matrix_samples/Ref_pop20_"+str(i)+".pkl"
    Ref_20 = simulate_reference(z_mat, file_name)

# simulate the reference with the 10 manifestations
# of the 200-cell model
for i, z_mat in enumerate(Z_mat_list_200):
    z_mat = Z_mat_list_200[i]
    file_name = "Z_matrix_samples/Ref_pop200_"+str(i)+".pkl"
    Ref_200 = simulate_reference(z_mat, file_name)