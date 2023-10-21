import numpy as np 
import pickle
from copy import deepcopy
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from pymoo.indicators.hv import HV
from rankcrowding import RankAndCrowding
from GA import (
    crossover,
    mutate 
)
from diversity_metrics import (
    geno_diversity,
    pheno_diversity
)

# set up GA for single objective
def single_obj_GA(
        folder_path: str,
        problem: object,
        n_gen: int,
        population: np.ndarray,
        num_circuits: int, 
        obj: np.ndarray,
        probability_crossover: float,
        probability_mutation: float,
        mutate_dose: bool,
        metrics: bool =False
):
    # create list to store min obj function, 
    # circuits with min obj function for each
    # generation
    obj_min = np.zeros(n_gen + 1)

    # create list to store all obj functions 
    # and circuits for initial population 
    # and all generations 
    all_obj = []
    all_obj.append(obj)
    all_circuits = []
    all_circuits.append(population)

    # index of list that contains min obj function
    ind_min = np.argmin(obj)
    obj_min[0] = obj[ind_min]
    circuit_min = []
    circuit_min.append(population[ind_min])

    geno = None
    pheno = None
    # if storing metrics, create lists for those
    # metrics and store initial population value
    if metrics:
        geno = np.zeros(n_gen+1)
        geno[0] = geno_diversity(population)

        pheno = np.zeros_like(geno)
        pheno[0] = pheno_diversity(obj)

    for gen in range(n_gen):
        # perform crossover to generate new
        # population (children) from parent 
        # circuits if randomly generated float  
        # is less than probability_crossover
        if np.random.uniform() < probability_crossover:
            children = crossover(population, obj)
        else:
            children = deepcopy(population)

        # perform mutation on children if 
        # randomly generated float is less 
        # than probability_mutation (used 
        # in mutate function)
        mutate(problem, children, 
               probability_mutation, 
               dose=mutate_dose
        )

        # simulate topology and calculate obj
        # function for each circuit in children
        # and append to obj array
        obj_children = np.asarray([problem.func(g[0]) for g in children])
        obj = np.append(obj, obj_children)

        # append obj of children and children
        # to respective lists
        all_obj.append(obj_children)
        all_circuits.append(children)

        # add children to population array
        population = np.vstack((population, children))
        # return array of indices that would sort
        # obj
        S = np.lexsort([obj])
        # select top num_circuits obj from obj array 
        # and top num_circuits from population
        # (initial population + children of each gen)
        obj = obj[S[:num_circuits]]
        population = population[S[:num_circuits], :]

        # return index of minimum obj function
        # from obj
        ind_min = np.argmin(obj)

        # add min obj to obj_min and circuit with
        # min obj to circuit_min
        obj_min[gen + 1] = obj[ind_min]
        circuit_min.append(population[ind_min])

        # calculate metrics for population
        if metrics:
             geno[gen+1] = geno_diversity(population)
             pheno[gen+1] = pheno_diversity(obj)

        print("generation "+ str(gen) + "  complete")
    
    # print in which gen the min obj first appeared
    print(first_seen(obj_min))
    print(circuit_min[-1][0].dose)

    # reshape all_obj and all_circuits to be 
    # 1 column arrays
    all_obj = np.asarray(all_obj).reshape(
        num_circuits*(1 + n_gen), 1)
    all_circuits = np.asarray(all_circuits).reshape(
        num_circuits*(1 + n_gen), 1)
    
    # print minimum objectives from each
    # generation
    # print("minimum objectives: ", obj_min)
    # unique objectives and circuits for all
    # objectives and all_circuits (all gens)
    unique_obj, unique_indices = np.unique(all_obj,
                                            return_index=True)
    unique_circuits = all_circuits[unique_indices]

    print(len(unique_circuits))

    # save results as .pkl files
    file_name = "minimum_obj_all_gens.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(obj_min, fid)

    file_name = "min_obj_circuit_all_gens.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(circuit_min, fid)

    file_name = "all_objectives.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(all_obj, fid)

    file_name = "all_circuits.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(all_circuits, fid)

    file_name = "all_unique_obj.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(unique_obj, fid)

    file_name = "all_unique_circuits.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(unique_circuits, fid)

    # plot graph of circuit with final min 
    # objective
    graph_file_name = "circuit_with_min_obj"
    plot_graph(circuit_min[-1][0], 
                folder_path + "/" + graph_file_name)


def multi_obj_GA(
        folder_path: str,
        problem: object, 
        n_gen: int,
        population: np.ndarray,
        num_circuits: int,
        obj: np.ndarray,
        probability_crossover: float,
        probability_mutation: float,
        mutate_dose: bool,
        plot: str=False
):
    
    # define reference point and class
    # instance of hypervolume calculator
    ref_point = np.array([0, 0])
    hv = HV(ref_point=ref_point)

    # store the progression of hypervolumes
    hypervolumes = []
    
    # create list to store all obj functions 
    # anbd circuits for initial population 
    # and all generations 
    all_obj = []
    all_obj.append(obj)
    all_circuits = []
    all_circuits.append(population)

    # create class instance of non-dominated
    # sorting class (to sort multi-objective
    # and determine pareto front)
    nds = RankAndCrowding()
    for gen in range(n_gen):
        # sort objectives using non-dominated
        # sorting algorithm and return ranks 
        # for each circuit index in population
        _, rank_dict = nds.do(obj, num_circuits, return_rank=True)

        # perform crossover to generate new
        # population (children) from parent 
        # circuits if randomly generated float  
        # is less than probability_crossover
        if np.random.uniform() < probability_crossover:
            children = crossover(population, obj, rank_dict)
        else:
            children = deepcopy(population)

        # perform mutation on children if 
        # randomly generated float is less 
        # than probability_mutation (used 
        # in mutate function)
        mutate(problem, children, 
                probability_mutation, 
                dose=mutate_dose
        )

        # simulate topology and calculate obj
        # function for each circuit in children
        # and append to all_obj list and obj
        # array
        obj_children = [problem.func(g[0]) for g in children]
        all_obj.append(obj_children)
        obj_children = np.asarray(obj_children)

        # append children to all circuits
        all_circuits.append(children)

        obj = np.vstack((obj, obj_children))
        # add children to population array
        population = np.vstack((population, children))

        # sort objectives using non-dominated
        # sorting algorithm and return indices
        # of num_circuits highest rank
        S = nds.do(obj, num_circuits)

        # select top num_circuits obj from obj array 
        # and top num_circuits from population
        # (initial population + children of each gen)
        obj = obj[S]
        population = population[S, :]

        # append hypervolume to list
        hypervolumes.append(hv(obj))

        print("generation "+ str(gen) + "  complete")

    fronts = NonDominatedSorting().do(obj)

    # create a list based on whether the circuit 
    # has an inhibitor
    types = []
    for topo in population:
        inhib = "Activators"
        for part in topo[0].part_list:
            if part[0] == "I":
                inhib = "Inhibitors"
        types.append(inhib)

    obj_df = pd.DataFrame(obj, columns=["ON_rel", "FI_rel"])
    obj_df["type"] = types

    # reshape all_obj to be 2 column array and
    # all_circuits to be to be 1 column array
    all_obj = np.asarray(
        all_obj).reshape(num_circuits*(1 + n_gen), 2)
    all_circuits = np.asarray(all_circuits).reshape(
        num_circuits*(1 + n_gen), 1)
    
    # print final objectives
    print("final objectivces: ", obj)
    # unique objectives and circuits for all
    # objectives and all_circuits (all gens)
    unique_obj, unique_indices = np.unique(all_obj,
                                           axis=0,
                                            return_index=True)
    unique_circuits = all_circuits[unique_indices]

    # save results as .pkl files
    file_name = "final_objectives_df_with_type.pkl"
    obj_df.to_pickle(folder_path + "/" + file_name)
 
    file_name = "final_population.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(population)

    file_name = "all_objectives.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(all_obj, fid)

    file_name = "all_circuits.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(all_circuits, fid)

    file_name = "all_unique_obj.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(unique_obj, fid)

    file_name = "all_unique_circuits.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(unique_circuits, fid)

    file_name = "hypervolumes.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(hypervolumes, fid)

    if plot:
        for i, circuit in enumerate(population):
            graph_file_name = ("pareto_front_circuit_" 
                               + str(i))
            plot_graph(circuit[0],
                    folder_path + "/" + graph_file_name)
                   

def plot_graph(topology, file_name):
    plt.figure()
    plt.tight_layout()
    nx.draw_networkx(topology.graph, arrows=True, arrowsize=15, 
                     node_size=600, node_shape='s')
    plt.savefig(file_name+".svg")


def first_seen(progression):

    looking = True
    gen_num = 0

    # for each value in reversed list
    # of minimum objective function for
    # each generation
    for gen in reversed(progression):
        # if the value does not equal the
        # final value in reversed list
        # return length of progression
        # (number of generations) - gen_num
        # (number of generations to end)
        if progression[-1] != gen:
            return len(progression) - gen_num
        gen_num += 1

    return 0