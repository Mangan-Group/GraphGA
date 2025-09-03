import pickle
import numpy as np

"""Loads the 10 manifestations of the Z matrices (for the 20-cell and
200-cell models), and the corresponding reference case simulations."""

repo_path = "/Users/kdreyer/Documents/Github/GraphGA/"
Z_path = repo_path+"Z_matrix_samples/"

Z_mat_list_20 = []
Ref_list_20 = []
for i in range(10):
    Z_fname = "Z_mat_20_cell"+str(i)+".npy"
    with open(Z_path+Z_fname, "rb") as fid:
        Z_20 = np.load(fid)
    Z_mat_list_20.append(Z_20)
    Ref_fname = "Ref_pop20_"+str(i)+".pkl"
    with open(Z_path+Ref_fname, "rb") as fid:
        Ref_20 = pickle.load(fid)
    Ref_list_20.append(Ref_20)

Z_mat_list_200 = []
Ref_list_200 = []
for i in range(10):
    Z_fname = "Z_mat_200_cell_"+str(i)+".npy"
    with open(Z_path+Z_fname, "rb") as fid:
        Z_200 = np.load(fid)
    Z_mat_list_200.append(Z_200)
    Ref_fname = "Ref_pop200_"+str(i)+".pkl"
    with open(Z_path+Ref_fname, "rb") as fid:
        Ref_200 = pickle.load(fid)
    Ref_list_200.append(Ref_200)
