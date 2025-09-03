import numpy as np
import pickle
from scipy.integrate import odeint

# load promoter parameters
with open("promo.pkl", "rb") as fid:
    promo = pickle.load(fid)


def reference(y,t,k_end):
    """Defines the ODEs for the
    reference case, for the single
    cell model."""

    return np.array([k_end - 2.7*y[0],
                    y[0] - 0.029*y[1]])



# simulate the reference for the desired promoter
# parameters and save the final time point reporter
# expression for the OFF and ON states.
Ref = dict()
for k in list(promo.keys())[:2]:
    off = odeint(reference, np.zeros(2), np.arange(0, 42 + 1), args=(promo[k]['off']*promo['k_txn'],))[-1,-1]
    on = odeint(reference, np.zeros(2), np.arange(0, 42 + 1), args=(promo[k]['on']*promo['k_txn'],))[-1, -1]
    Ref.update({k: {'off': off, 'on': on, 'fi': on/off}})

with open("Ref.pkl", "wb") as fid:
    pickle.dump(Ref, fid)
