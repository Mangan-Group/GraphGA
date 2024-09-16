import pickle
import numpy as np

quest = False
if quest:
    Z_path = "/home/ksd844/GCAD/GA_hyperparameter_opt/Z_matrix_samples/"
else:
    Z_path = "/Users/kdreyer/Documents/Github/GraphGA/Z_matrix_samples/"

Z_mat_list = []
Ref_list = []
for i in range(10):
    Z_fname = "Z_mat_20_cell"+str(i)+".npy"
    with open(Z_path+Z_fname, "rb") as fid:
        Z_20 = np.load(fid)
    Z_mat_list.append(Z_20)
    # print(Z_20[0,0])
    Ref_fname = "Ref_pop20_"+str(i)+".pkl"
    with open(Z_path+Ref_fname, "rb") as fid:
        Ref_20 = pickle.load(fid)
    Ref_list.append(Ref_20)
    # print(Ref_20["P1"]["off"])
