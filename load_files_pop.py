import pickle
import numpy as np

quest = False
if quest:
    repo_path = "/home/ksd844/GCAD/GA_hyperparameter_opt/"
else:
    repo_path = "/Users/kdreyer/Documents/Github/GraphGA/"
with open(repo_path + "promo.pkl", "rb") as fid:
    promo = pickle.load(fid)
with open(repo_path + "parts.pkl", "rb") as fid:
    parts = pickle.load(fid)
with open(repo_path + "Ref.pkl", "rb") as fid:
    Ref = pickle.load(fid)
with open(repo_path + "Ref_pop20.pkl", "rb") as fid:
    Ref_pop20 = pickle.load(fid)

with open(repo_path + "Ref_pop200.pkl", "rb") as fid:
    Ref_pop200 = pickle.load(fid)

with open(repo_path + "Ref_pop2000.pkl", "rb") as fid:
    Ref_pop2000 = pickle.load(fid)

Z_path = repo_path + "Z_matrix_samples/"
with open(Z_path + "Z_mat_20_cell0.npy", 'rb') as fid:
    Z_20 = np.load(fid)

with open(Z_path + "Z_mat_200_cell.npy", 'rb') as fid:
    Z_200 = np.load(fid)

with open(Z_path + "Z_mat_2000_cell.npy", 'rb') as fid:
    Z_2000 = np.load(fid)

tf_list = [k for k in parts.keys() if k[0] == 'Z']
inhibitor_list = [k for k in parts.keys() if k[0] == 'I']