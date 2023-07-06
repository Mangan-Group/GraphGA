import pickle
import numpy as np

Z_path = "/Users/kdreyer/Desktop/Github/GCAD_Test_Cases/20_cell/"

Z_mat_list = []

for i in range(10):
    fname = "Z_mat20_new"+str(i)+".npy"
    with open(Z_path+fname, "rb") as fid:
        Z_20 = np.load(fid)
    Z_mat_list.append(Z_20)

# with open(Z_path+"Z_mat_nc20_p5_1.npy", 'rb') as fid:
#     Z_20_1 = np.load(fid)

# with open(Z_path+"Z_mat_nc20_p5_1.npy", 'rb') as fid:
#     Z_20_1 = np.load(fid)

# Z_mat_list = [
#     Z_20_1, 
    # Z_20_2, 
    # Z_20_3, 
    # Z_20_4, 
    # Z_20_5,
    # Z_20_6,
    # Z_20_7,
    # Z_20_8,
    # Z_20_9,
    # Z_20_10
# ]

for i in range(len(Z_mat_list)):
    print(Z_mat_list[i])
# for i in range(len(Z_row1)):
#     print(Z_row1[i])

# print(Z_mat_list[1][0,0]- Z_mat_list[0][0,0])