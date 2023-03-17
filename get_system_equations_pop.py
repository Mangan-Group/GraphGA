import numpy as np
from load_files_pop import (
    promo, parts, parts_order
)

def system_equations_pop(x, t, state, Z_list, topology):
    system = []
    index = 0
    for n in topology.in_dict.keys():
        eq = -2.7 * x[2 * topology.var_dict[n]]
        b = []
        num = 0
        denom = 1
        Z_ZF = 0
        for k in topology.in_dict[n]['P']:
            eq += float(topology.dose[n]) / topology.pool[n] * promo[k][state] * promo['k_txn'] * Z_list[index]
            index += 1
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
            Z_ZF = Z_list[index]
            index += 1
        eq += float(topology.dose[n]) / topology.pool[n] * (b + num) / denom * 9. * Z_ZF
        system.extend([eq, -topology.protein_deg[n[0]] * x[2 * topology.var_dict[n] + 1] + x[2 * topology.var_dict[n]]])
    return system

def system_equations_DsRed_pop(x, t, state, Z_list, topology):
    system = []
    index = 0
    for n in topology.in_dict.keys():
        eq = -2.7 * x[2 * topology.var_dict[n]]
        b = []
        num = 0
        denom = 1
        Z_ZF = 0
        for k in topology.in_dict[n]['P']:
            eq += float(topology.dose[n]) / topology.pool[n] * promo[k][state] * promo['k_txn'] * Z_list[index]
            index += 1
        if len(topology.in_dict[n]['I']) == 0:
            for k in topology.in_dict[n]['Z']:
                b.append(parts[k][0])
                num += parts[k][1] * parts[k][2] * x[2 * topology.var_dict[k] + 1]
                denom += parts[k][2] * x[2 * topology.var_dict[k] + 1]
        else:
            zfa = topology.in_dict[n]['Z'][0]
            zfi = topology.in_dict[n]['I'][0]
            b.extend([parts[zfa][0], parts['Z' + zfi[1:]][0]])
            A = parts[zfa][1] + (1 - parts[zfa][1])/6 * (4 * parts[zfi][0] * x[2 * topology.var_dict[zfi] + 1]/(1e-10 + parts[zfa][2] * x[2 * topology.var_dict[zfa] + 1]) - 2)
            m = np.piecewise(A, [A < 1, A >= parts[zfa][1]], [1., parts[zfa][1], A])
            num += m * parts[zfa][2] * x[2 * topology.var_dict[zfa] + 1]
            denom += parts[zfa][2] * x[2 * topology.var_dict[zfa] + 1] + 4 * parts[zfi][0] * x[2 * topology.var_dict[zfi] + 1]
        if len(b) == 0:
            b = 0
        else:
            b = np.mean(b)
            Z_ZF = Z_list[index]
            index += 1
        eq += float(topology.dose[n]) / topology.pool[n] * (b + num) / denom * 9. * Z_ZF
        system.extend([eq, -topology.protein_deg[n[0]] * x[2 * topology.var_dict[n] + 1] + x[2 * topology.var_dict[n]]])
    return system

