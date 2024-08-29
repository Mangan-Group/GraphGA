import pickle
import numpy as np

Z_path = "/Users/kdreyer/Documents/Github/GraphGA/Z_matrices/"

Z_mat_list = []

for i in range(10):
    fname = "Z_mat_20_cell"+str(i)+".npy"
    with open(Z_path+fname, "rb") as fid:
        Z_20 = np.load(fid)
    Z_mat_list.append(Z_20)

# print(Z_mat_list)

