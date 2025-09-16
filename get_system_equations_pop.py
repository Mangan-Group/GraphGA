import numpy as np
from load_files_pop import (
    promo, parts
)

def system_equations_pop(
        x,
        t,
        state: str,
        Z_list: list,
        topology: object
):
    """Builds the system of ODEs for
    the given topology, for the synTF-R.
    inhibitors. When using the
    single cell model, Z_list is a list
    of ones."""

    system = []
    # need to initialize index=0 for indexing
    # z_list
    index = 0
    # loop through in_dict keys (contains
    # each part)-- need ODEs for each part
    for n in topology.in_dict.keys():
        # add RNA degradation term
        eq = -2.7 * x[2 * topology.var_dict[n]]
        b = []
        num = 0
        denom = 1
        Z_ZF = 0
        # loop through promoter regulation in in_dict 
        # to get RNA ODE term for promoter regulation  
        for k in topology.in_dict[n]['P']:
            eq += ((float(topology.dose[n]) / topology.pool[n])/200)**0.5 * promo[k][state] * promo['k_txn'] * Z_list[index]
            # increase index by 1 for indexing z_list
            index += 1
        # loop through Z (synTF-A) regulation in in_dict
        # to build RNA ODE term for synTF-A regulation
        for k in topology.in_dict[n]['Z']:
            b.append(parts[k][0])
            num += parts[k][1] * parts[k][2] * x[2 * topology.var_dict[k] + 1]
            denom += parts[k][2] * x[2 * topology.var_dict[k] + 1]
        # loop through I (synTF-R) regulation in in_dict
        # to build RNA ODE term for synTF-R regulation
        for k in topology.in_dict[n]['I']:
            if ('Z' + k[1:]) not in topology.in_dict[n]['Z']:
                b.append(parts['Z' + k[1:]][0])
            denom += parts[k][0] * x[2 * topology.var_dict[k] + 1]
        # calculate b term (0 if no synTF regulation, otherwise
        # mean of b values for synTFs)
        if len(b) == 0:
            b = 0
        else:
            b = np.mean(b)
            Z_ZF = Z_list[index]
            index += 1
        # define term for all synTF regulation from 
        # synTF-A/R regulation
        eq += ((float(topology.dose[n]) / topology.pool[n])/200)**0.5 * (b + num) / denom * 9. * Z_ZF
        system.extend([eq, -topology.protein_deg[n[0]] * x[2 * topology.var_dict[n] + 1] + x[2 * topology.var_dict[n]]])
    return system

def system_equations_DsRed_pop(
        x,
        t,
        state,
        Z_list,
        topology
):
    """Builds the system of ODEs for
    the given topology, for the synTF-R.
    inhibitors. When using the
    single cell model, Z_list is a list
    of ones.
    Here, the process is the same as in
    system_equations_pop(), but the synTF-R-DsRED
    term is used to define synTF-R regulation.
    See system_equations_pop() for comments."""

    system = []
    index = 0
    for n in topology.in_dict.keys():
        eq = -2.7 * x[2 * topology.var_dict[n]]
        b = []
        num = 0
        denom = 1
        Z_ZF = 0
        for k in topology.in_dict[n]['P']:
            eq += ((float(topology.dose[n]) / topology.pool[n])/200)**0.5 * promo[k][state] * promo['k_txn'] * Z_list[index]
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
        eq += ((float(topology.dose[n]) / topology.pool[n])/200)**0.5 * (b + num) / denom * 9. * Z_ZF
        # eq += float(topology.dose[n]) / topology.pool[n] * (b + num) / denom * 9. * Z_ZF

        system.extend([eq, -topology.protein_deg[n[0]] * x[2 * topology.var_dict[n] + 1] + x[2 * topology.var_dict[n]]])
    return system

