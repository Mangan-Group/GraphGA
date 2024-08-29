import numpy as np
from load_files import *

def system_equations(x, t, state, topology):
    system = []
    # for each part and reporter
    for n in topology.in_dict.keys():
        # add mRNA degradation term (index
        # is 2*state number assigned in
        # define_circuit because 2 eqs per
        # state)
        eq = -2.7 * x[2 * topology.var_dict[n]]
        b = []
        num = 0
        denom = 1
        # add promoter regulation terms if
        # promoter regulation (pool is dose
        # pool if part is regulated by tf and
        # promoter (2 plasmids))
        for k in topology.in_dict[n]['P']:
            eq += ((float(topology.dose[n]) / topology.pool[n])/200)**0.5 * promo[k][state] * promo['k_txn']
        # add b to list, and add numerator and
        # denominator for each tf (index
        # is 2*state number assigned in
        # define_circuit +1 for protein)
        for k in topology.in_dict[n]['Z']:
            b.append(parts[k][0])
            num += parts[k][1] * parts[k][2] * x[2 * topology.var_dict[k] + 1]
            denom += parts[k][2] * x[2 * topology.var_dict[k] + 1]
        # add b if different inhibitor than 
        # tf and add inhibitor denominator
        # for each ingibitor
        for k in topology.in_dict[n]['I']:
            if ('Z' + k[1:]) not in topology.in_dict[n]['Z']:
                b.append(parts['Z' + k[1:]][0])
            denom += parts[k][0] * x[2 * topology.var_dict[k] + 1]
        if len(b) == 0:
            b = 0
        # take mean of b values (needed
        # if more than 1)
        else:
            b = np.mean(b)
        # full production term
        eq += ((float(topology.dose[n]) / topology.pool[n])/200)**0.5 * (b + num) / denom * 9.
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
            eq += ((float(topology.dose[n]) / topology.pool[n])/200)**0.5 * promo[k][state] * promo['k_txn']

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
        eq += ((float(topology.dose[n]) / topology.pool[n])/200)**0.5 * (b + num) / denom * 9.
        # add mRNA equation and protein equation to
        # system
        system.extend([eq, -topology.protein_deg[n[0]] * x[2 * topology.var_dict[n] + 1] + x[2 * topology.var_dict[n]]])
    return system
