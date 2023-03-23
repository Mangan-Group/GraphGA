import numpy as np
import networkx as nx
import pickle
import matplotlib.pyplot as plt

with open("Results/GA_results/230316_Amplifier_n_gen20_pop.pkl", "rb") as fid:
    amp_GA = pickle.load(fid)

# print(type(amp_GA["genotype"]))
# print(amp_GA["phenotype"])
# print(np.arange(len(amp_GA["phenotype"])))
# print(amp_GA["obj_min"][-1])
# print(amp_GA["circuit_min"][-1][0].edge_list)

def plot_obj_distribution(obj_list: np.ndarray, label: str, n_bins: int =30):
    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    ax.hist(obj_list, bins=n_bins)
    ax.set_xlabel(label)
    ax.set_ylabel("Count")
    # plt.show()
    # plt.savefig("230316_amp_1part_ON_rel_dist.svg")
    
def plot_metric(metric_list: np.ndarray, label: str):
    generations = np.arange(len(metric_list))

    fig, ax = plt.subplots(1, 1, figsize= (4, 4))
    ax.plot(generations, metric_list)
    ax.set_xlabel("Generation")
    ax.set_ylabel(label)
    plt.savefig("Amplifier_"+label+"_"+str(len(generations)-1)+".svg")

plot_metric(amp_GA["phenotype"], "Phenotype")

