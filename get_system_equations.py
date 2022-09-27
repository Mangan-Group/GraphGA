import numpy as np
from load_files import *

def system_equations(x, t, state, topology):
    system = []
    for n in topology.in_dict.keys():
        eq = -2.7 * x[2 * topology.var_dict[n]]
        b = []
        num = 0
        denom = 1
        for k in topology.in_dict[n]['P']:
            eq += float(topology.dose[n]) / topology.pool[n] * promo[k][state]
        for k in topology.in_dict[n]['Z']:
            b.append(parts[k][0])
            num += parts[k][1] * parts[k][2] * x[2 * topology.var_dict[k] + 1]
            denom += parts[k][2] * x[2 * topology.var_dict[k] + 1]
        for k in topology.in_dict[n]['I']:
            if ('Z' + k[1:]) not in topology.in_dict[n]['Z']:
                b.append(parts['Z' + k[1:]][0])
            denom += parts[k][0] * x[2 * topology.var_dict[k] + 1]
        if len(b) == 0:
            b = 0
        else:
            b = np.mean(b)
        eq += float(topology.dose[n]) / topology.pool[n] * (b + num) / denom
        system.extend([eq, -topology.protein_deg[n[0]] * x[2 * topology.var_dict[n] + 1] + x[2 * topology.var_dict[n]]])
    return system

def system_equations_DsRed(x, t, state, topology):
    system = []
    for n in topology.in_dict.keys():
        eq = -2.7 * x[2 * topology.var_dict[n]]
        b = []
        num = 0
        denom = 1
        for k in topology.in_dict[n]['P']:
            eq += float(topology.dose[n]) / topology.pool[n] * promo[k][state]

        if len(topology.in_dict[n]['I']) == 0:
            for k in topology.in_dict[n]['Z']:
                b.append(parts[k][0])
                num += parts[k][1] * parts[k][2] * x[2 * topology.var_dict[k] + 1]
                denom += parts[k][2] * x[2 * topology.var_dict[k] + 1]
        else:
            zfa = topology.in_dict[n]['Z'][0]
            zfi = topology.in_dict[n]['I'][0]
            b.extend([parts[zfa][0], parts[zfi][0]])
            u = 2.
            l = 0.25
            m = (parts[zfa][1] + parts['Z' + zfi[1:]][1]) / 2
            A = ((4 * parts[zfi][0] * x[2 * topology.var_dict[zfi] + 1]) /
                 (1e-10 + parts[zfa][2] * x[2 * topology.var_dict[zfa] + 1]) - 4 * l) * (1 - m) / (4 * u - 4 * l) + m
            num += np.piecewise(A, [A < 1, A >= m], [1., m, A])
            denom += parts[zfa][2] * x[2 * topology.var_dict[zfa] + 1] + 4 * parts[zfi][0] * x[2 * topology.var_dict[zfi] + 1]
        if len(b) == 0:
            b = 0
        else:
            b = np.mean(b)
        eq += float(topology.dose[n]) / topology.pool[n] * (b + num) / denom
        system.extend([eq, -topology.protein_deg[n[0]] * x[2 * topology.var_dict[n] + 1] + x[2 * topology.var_dict[n]]])
    return system

