import numpy as np
import networkx as nx
import pickle
import matplotlib.pyplot as plt
from load_files_pop import (
    Z_20,
    Z_200
)

def plot_obj_distribution(
        obj_list: np.ndarray, 
        x_label: str, 
        n_bins: int =30,
        log_scale =False
):
    fig, axs = plt.subplots(1, 2, figsize= (6, 4))
    axs.ravel()
    
    if log_scale:
        hist, bins, _ = axs[1].hist(obj_list, bins=n_bins)
        logbins = np.logspace(np.log10(bins[0]),np.log10(bins[-1]),len(bins))
        axs[0].hist(obj_list, bins=logbins)
        axs[0].set_xscale('log')
    else:
        axs[0].hist(obj_list, bins=n_bins)

    axs[0].set_xlabel(x_label)
    axs[0].set_ylabel("Count")
    plt.show()
    # plt.savefig("230316_amp_1part_ON_rel_dist.svg")
    
def plot_metric(metric_list: np.ndarray, label: str):
    generations = np.arange(len(metric_list))

    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    ax.plot(generations, metric_list)
    ax.set_xlabel("Generation")
    ax.set_ylabel(label)
    plt.savefig("Amplifier_"+label+"_"+str(len(generations)-1)+".svg")

with open("Results/GA_results/230316_Amplifier_n_gen20_pop.pkl", "rb") as fid:
    amp_GA = pickle.load(fid)

# print(type(amp_GA["genotype"]))
# print(amp_GA["phenotype"])
# print(np.arange(len(amp_GA["phenotype"])))
# print(amp_GA["obj_min"][-1])
# print(amp_GA["circuit_min"][-1][0].edge_list)

# plot_metric(amp_GA["phenotype"], "Phenotype")



plot_obj_distribution(
    Z_200[:, 0],
    "mean-centered fluor. (log10 a.u.)",
    n_bins=20,
    log_scale=True
)

# print(Z_20[:, 1])