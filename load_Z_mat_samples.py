import pickle
import numpy as np

quest = False
if quest:
    Z_path_20 = "/home/ksd844/GCAD/GA_hyperparameter_opt/Z_matrix_samples/"
else:
    Z_path_20 = "/Users/kdreyer/Documents/Github/GraphGA/Z_matrix_samples/"

Z_path_200 = "/Users/kdreyer/Documents/Github/GraphGA/Z_matrix_samples/"

Z_mat_list_20 = []
Ref_list_20 = []
for i in range(10):
    Z_fname = "Z_mat_20_cell"+str(i)+".npy"
    with open(Z_path_20+Z_fname, "rb") as fid:
        Z_20 = np.load(fid)
    Z_mat_list_20.append(Z_20)
    Ref_fname = "Ref_pop20_"+str(i)+".pkl"
    with open(Z_path_20+Ref_fname, "rb") as fid:
        Ref_20 = pickle.load(fid)
    Ref_list_20.append(Ref_20)

Z_mat_list_200 = []
Ref_list_200 = []
for i in range(10):
    Z_fname = "Z_mat_200_cell_"+str(i)+".npy"
    with open(Z_path_200+Z_fname, "rb") as fid:
        Z_200 = np.load(fid)
    Z_mat_list_200.append(Z_200)
    Ref_fname = "Ref_pop200_"+str(i)+".pkl"
    with open(Z_path_200+Ref_fname, "rb") as fid:
        Ref_200 = pickle.load(fid)
    Ref_list_200.append(Ref_200)
