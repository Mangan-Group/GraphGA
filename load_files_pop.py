import pickle
import numpy as np

with open("promo.pkl", "rb") as fid:
    promo = pickle.load(fid)
with open("parts.pkl", "rb") as fid:
    parts = pickle.load(fid)
with open("Ref.pkl", "rb") as fid:
    Ref = pickle.load(fid)
with open("Z_200.npy", 'rb') as fid:
    Z_200 = np.load(fid)
with open("Z_20.npy", 'rb') as fid:
    Z_20 = np.load(fid)
with open("Ref_pop200.pkl", "rb") as fid:
    Ref_pop200 = pickle.load(fid)
with open("Ref_pop20.pkl", "rb") as fid:
    Ref_pop20 = pickle.load(fid)

tf_list = [k for k in parts.keys() if k[0] == 'Z']
inhibitor_list = [k for k in parts.keys() if k[0] == 'I']