import pickle

repo_path = "/Users/kdreyer/Documents/Github/GraphGA/"
with open(repo_path + "promo.pkl", "rb") as fid:
    promo = pickle.load(fid)
with open(repo_path + "parts.pkl", "rb") as fid:
    parts = pickle.load(fid)
with open(repo_path + "Ref.pkl", "rb") as fid:
    Ref = pickle.load(fid)

tf_list = [k for k in parts.keys() if k[0] == 'Z']
inhibitor_list = [k for k in parts.keys() if k[0] == 'I']