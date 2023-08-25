import pandas as pd

with open('Results/SetDose.pkl', 'rb') as f:
    circuit = pd.read_pickle(f)

print(circuit)