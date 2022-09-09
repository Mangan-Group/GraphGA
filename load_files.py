import pickle

with open("promo.pkl", "rb") as fid:
    promo = pickle.load(fid)
with open("parts.pkl", "rb") as fid:
    parts = pickle.load(fid)
with open("Ref.pkl", "rb") as fid:
    Ref = pickle.load(fid)
tf_list = [k for k in parts.keys() if k[0]=='Z']
in_list = [k for k in parts.keys() if k[0]=='I']