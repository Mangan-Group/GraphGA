import pandas as pd
from main_function_sigcond import*
import matplotlib.pyplot as plt





def simulate(topology, max_time=42):
    t = np.arange(0, max_time + 1, 1)
    rep_off = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('off', topology,))[-1, -1]
    rep_on = odeint(system_equations, np.zeros(topology.num_states * 2), t, args=('on', topology,))[-1, -1]
    return rep_off, rep_on


def func(toplogy):
    rep_off, rep_on = simulate(toplogy)
    ON_rel = rep_on / Ref[toplogy.promo_node]['on']
    FI_rel = (rep_on / rep_off) / (Ref[toplogy.promo_node]['fi'])
    return [-ON_rel, -FI_rel]



with open('SigCond/SigCond_combo_2.pkl', 'rb') as f:
    sigcond = pd.read_pickle(f)

print("SigCond:", len(sigcond))

with open('Amplifier/Amplifier_combo_1.pkl', 'rb') as f:
    amp1 = pd.read_pickle(f)

print("Amp1:", len(amp1))

with open('Amplifier/Amplifier_combo_2.pkl', 'rb') as f:
    amp2 = pd.read_pickle(f)

print("Amp2:", len(amp2))


total_runs = len(sigcond) + len(amp1) + len(amp2)
print("Total Runs:", total_runs)
current_run = 0

num_circuits = len(sigcond)

full = []

for i in sigcond:
    current_run += 1
    if current_run % 1000 == 0:
        print(current_run,"-------------------------------------------------------------------")
    result = func(i)
    print(result)
    full.append(result)

# for i in amp1:
#     current_run += 1
#     if current_run % 1000 == 0:
#         print(current_run,"-------------------------------------------------------------------")
#     result = func(i)
#     print(result)
#     full.append(result)
#
# for i in amp2:
#     current_run += 1
#     if current_run % 1000 == 0:
#         print(current_run,"-------------------------------------------------------------------")
#     result = func(i)
#     print(result)
#     full.append(result)

full = np.asarray(full)

fronts = NonDominatedSorting().do(full)

first = full[fronts[0],:]

plt.scatter(first[:,0], first[:,1], color = "blue")
plt.show()

# for i in amp1:
#     current_run += 1
#     if current_run % 1000 == 0:
#         print(current_run,"-------------------------------------------------------------------")
#     result = func(i)
#     print(result)
#     full.append(result)
#
# for i in amp2:
#     current_run += 1
#     if current_run % 1000 == 0:
#         print(current_run,"-------------------------------------------------------------------")
#     result = func(i)
#     print(result)
#     full.append(result)
