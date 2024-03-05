import numpy as np 
import pickle
from copy import deepcopy
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from multiprocessing import Pool
from alive_progress import alive_bar
from pymoo.util.nds.non_dominated_sorting import NonDominatedSorting
from pymoo.indicators.hv import HV
from rankcrowding import RankAndCrowding
from pulse_generator_problem import PulseGenerator
from GA import (
    crossover,
    mutate 
)
from diversity_metrics import (
    geno_diversity,
    pheno_diversity,
    first_seen
)
from plot_search_results import(
    plot_graph,
    plot_1D_obj_scatter,
    plot_pareto_front,
    plot_pareto_front3D,
    plot_hypervolume
)

# set up GA for single objective
def single_obj_GA(
        folder_path: str,
        problem: object,
        population: np.ndarray,
        num_circuits: int, 
        obj: np.ndarray,
        get_unique: bool=False,
        metrics: bool =False
):
    
    # create list to store min obj function, 
    # circuits with min obj function for each
    # generation
    obj_min = np.zeros(problem.n_gen + 1)

    # create list to store all obj functions 
    # and circuits for initial population 
    # and all generations 
    all_obj = []
    all_obj.append(obj)
    all_circuits = []
    all_circuits.append(population)

    # index of obj list that contains min 
    # obj function
    ind_min = np.argmin(obj)
    obj_min[0] = obj[ind_min]
    circuit_min = []
    circuit_min.append(population[ind_min])

    geno = None
    pheno = None
    # if storing metrics, create lists for those
    # metrics and store initial population value
    if metrics:
        geno = np.zeros(problem.n_gen+1)
        geno[0] = geno_diversity(population)

        pheno = np.zeros_like(geno)
        pheno[0] = pheno_diversity(obj)
    with alive_bar(problem.n_gen) as bar:
        for gen in range(problem.n_gen):
            # perform crossover to generate new
            # population (children) from parent 
            # circuits if randomly generated float  
            # is less than probability_crossover
            if np.random.uniform() < problem.prob_crossover:
                children = crossover(population, obj)
            else:
                children = deepcopy(population)

            # perform mutation on children if 
            # randomly generated float is less 
            # than probability_mutation (used 
            # in mutate function)
            mutate(
                problem, children, 
                problem.prob_mutation, 
                dose=problem.mutate_dose
            )

            # simulate topology and calculate obj
            # function for each circuit in children
            # and append to obj array
            if problem.pop:
                child_topologies = [g[0] for g in children]
                with Pool(problem.num_processes) as pool:
                    obj_list = pool.imap(problem.func, child_topologies)

                    pool.close()
                    pool.join()
                obj_list = list(obj_list)
                obj_children = np.asarray(obj_list)
                
            else:
                obj_children = np.asarray(
                    [problem.func(g[0]) for g in children])
                
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

            # print("generation "+ str(gen) + " complete")
            bar()
    
    # print in which gen the min obj first appeared
    # print(first_seen(obj_min))
    # print doses and edge list for opt circuit
    # print(circuit_min[-1][0].dose)
    # print(circuit_min[-1][0].edge_list)

    # reshape all_obj and all_circuits to be 
    # 1 column arrays
    all_obj = np.asarray(all_obj).reshape(
        num_circuits*(1 + problem.n_gen), 1)
    all_circuits = np.asarray(all_circuits).reshape(
        num_circuits*(1 + problem.n_gen), 1)
    
    # save results for either single cell or population
    file_name = "all_objectives.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(all_obj, fid)

    file_name = "all_circuits.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(all_circuits, fid)
    
    if problem.pop:
        # save all_cell results df (see problem class)
        file_name = "All_circuits_all_cell_results.pkl"
        problem.all_cells.to_pickle(folder_path + "/" + file_name)

        # get edge lists for each circuit sampled in GA
        # and add dose as an "edge"
        circuit_edge_lists = []
        for circuit in all_circuits:
            circuit_edges = circuit[0].edge_list
            for key, val in circuit[0].dose.items():
                circuit_edges.append((key, str(val)))
            circuit_edge_lists.append(circuit_edges)


        combo_edges_lists = []
        for edges in circuit_edge_lists:
            edge_combos = [edge[0] + edge[1] for
                edge in edges]
            combo_edges_lists.append(edge_combos)

        unique_edge_combo = []
        index_list = []
        seen = set()
        for i, combo_list in enumerate(
            combo_edges_lists):
            combo_set = frozenset(combo_list)
            if combo_set not in seen:
                seen.add(combo_set)
                index_list.append(i)
                unique_edge_combo.append(combo_list)

        unique_obj = all_obj[index_list]
        unique_circuits = all_circuits[index_list]

        file_name = "unique_objectives.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(unique_obj, fid)

        file_name = "unique_circuits.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(unique_circuits, fid)

        graph_file_name = "objs_scatter_plot.svg"
        plot_1D_obj_scatter(
            folder_path + "/" + graph_file_name,
            unique_obj, 
            problem.obj_labels,
        )

    else:
        # save single cell specific results as .pkl files
        file_name = "minimum_obj_all_gens.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(obj_min, fid)

        file_name = "min_obj_circuit_all_gens.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(circuit_min, fid)
    
        # plot graph of circuit with final min 
        # objective
        graph_file_name = "circuit_with_min_obj.svg"
        plot_graph(
            folder_path + "/" + graph_file_name,
            circuit_min[-1][0]
        )

        if get_unique:
            # unique objectives and circuits for all
            # objectives and all_circuits (all gens)
            unique_obj, unique_indices = np.unique(all_obj,
                                            return_index=True)
            unique_circuits = all_circuits[unique_indices]

            print(len(unique_circuits))

            file_name = "all_unique_obj.pkl"
            with open(folder_path + "/" + file_name, "wb") as fid:
                pickle.dump(unique_obj, fid)

            file_name = "all_unique_circuits.pkl"
            with open(folder_path + "/" + file_name, "wb") as fid:
                pickle.dump(unique_circuits, fid)


def multi_obj_GA(
        folder_path: str,
        problem: object, 
        population: np.ndarray,
        num_circuits: int,
        obj: np.ndarray,
        get_unique: bool=False,
        plot: str=False
):
    
    # set plotting function for pareto
    # front and
    # define reference point and class
    # instance of hypervolume calculator
    if "t_pulse" in '\t'.join(problem.obj_labels):
        if len(problem.obj_labels) == 3:
            problem.pareto_plot = plot_pareto_front3D
            ref_point = np.array([problem.max_time, 0, 0])
        else:
            problem.pareto_plot = plot_pareto_front
            ref_point = np.array([problem.max_time, 0])
    else:
        problem.pareto_plot = plot_pareto_front
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
    with alive_bar(problem.n_gen) as bar:
        for gen in range(problem.n_gen):
            # sort objectives using non-dominated
            # sorting algorithm and return ranks 
            # for each circuit index in population
            _, rank_dict = nds.do(obj, num_circuits, return_rank=True)

            # perform crossover to generate new
            # population (children) from parent 
            # circuits if randomly generated float  
            # is less than probability_crossover
            if np.random.uniform() < problem.prob_crossover:
                children = crossover(population, obj, rank_dict)
            else:
                children = deepcopy(population)

            # perform mutation on children if 
            # randomly generated float is less 
            # than probability_mutation (used 
            # in mutate function)
            mutate(problem, children, 
                    problem.prob_mutation, 
                    dose=problem.mutate_dose
            )

            # simulate topology and calculate obj
            # function for each circuit in children
            # and append to all_obj list and obj
            # array
            if problem.pop:
                child_topologies = [g[0] for g in children]
                with Pool(problem.num_processes) as pool:
                    obj_list = pool.imap(problem.func, child_topologies)

                    pool.close()
                    pool.join()
                obj_list = list(obj_list)
                obj_children = np.asarray(obj_list)
                
            else:
                obj_children = np.asarray(
                    [problem.func(g[0]) for g in children])
                    
            all_obj.append(obj_children)
            
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

            # print("generation "+ str(gen) + " complete")
            bar()

    fronts = NonDominatedSorting().do(obj)

    obj_df = pd.DataFrame(obj, columns=problem.obj_labels)

    # don't get types for pulse (all use inhibitors)
    if not isinstance(problem, PulseGenerator):
        # create a list based on whether the circuit 
        # has an inhibitor
        types = []
        for topo in population:
            inhib = "Activators"
            for part in topo[0].part_list:
                if part[0] == "I":
                    inhib = "Inhibitors"
            types.append(inhib)

        obj_df["type"] = types
        types_ = True
    else:
        types_ = False

    # reshape all_obj to be array with num columns = # objs and
    # all_circuits to be to be 1 column array
    all_obj = np.asarray(
        all_obj).reshape(num_circuits*(1 + problem.n_gen), len(problem.obj_labels))
    all_circuits = np.asarray(all_circuits).reshape(
        num_circuits*(1 + problem.n_gen), 1)

    # save results as .pkl files
    file_name = "final_objectives_df.pkl"
    obj_df.to_pickle(folder_path + "/" + file_name)
 
    file_name = "final_population.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(population, fid)

    file_name = "all_objectives.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(all_obj, fid)

    file_name = "all_circuits.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(all_circuits, fid)

    file_name = "hypervolumes.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(hypervolumes, fid)

    graph_file_name = "final_population_pareto_front.svg"
    problem.pareto_plot(
        folder_path + "/" + graph_file_name,
        obj_df,
        problem.obj_labels,
        types=types_
    )
    hv_progression_file_name = "hv_progression.svg"
    plot_hypervolume(
        folder_path + "/" + hv_progression_file_name,
        problem.n_gen,
        hypervolumes
    )

    if problem.pop:
        # save all_cell results df (see problem class)
        file_name = "All_circuits_all_cell_results.pkl"
        problem.all_cells.to_pickle(folder_path + "/" + file_name)

        # S_all = nds.do(all_obj, len(all_obj))
        # sorted_all_obj = all_obj[S_all, :]
        # sorted_all_circuits = all_circuits[S_all]

        circuit_edge_lists = []
        for circuit in all_circuits:
            circuit_edges = circuit[0].edge_list
            for key, val in circuit[0].dose.items():
                circuit_edges.append((key, str(val)))
            circuit_edge_lists.append(circuit_edges)

        combo_edges_lists = []
        for edges in circuit_edge_lists:
            edge_combos = [edge[0] + edge[1] for
                edge in edges]
            combo_edges_lists.append(edge_combos)

        unique_edge_combo = []
        index_list = []
        seen = set()
        for i, combo_list in enumerate(
            combo_edges_lists):
            combo_set = frozenset(combo_list)
            if combo_set not in seen:
                seen.add(combo_set)
                index_list.append(i)
                unique_edge_combo.append(combo_list)

        unique_obj = all_obj[index_list]
        unique_circuits = all_circuits[index_list]

        # unique_types = []
        # for topo in unique_circuits:
        #     inhib = "Activators"
        #     for part in topo[0].part_list:
        #         if part[0] == "I":
        #             inhib = "Inhibitors"
        #     unique_types.append(inhib)

        unique_obj_df = pd.DataFrame(unique_obj, columns=problem.obj_labels)
        # unique_obj_df["type"] = unique_types

        file_name = "unique_objectives_df.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(unique_obj_df, fid)

        file_name = "unique_circuits.pkl"
        with open(folder_path + "/" + file_name, "wb") as fid:
            pickle.dump(unique_circuits, fid)

        # scatter plot of all obj values for all unique
        # circuits in GA run
        graph_file_name = "unique_obj_scatter_plot.svg"
        plot_pareto_front(
            folder_path + "/" + graph_file_name,
            unique_obj_df,
            problem.obj_labels,
            types=False
        )

    else:
        # get unique circuits and obj for all 
        # generations of GA
        if get_unique:
            # unique objectives and circuits for all
            # objectives and all_circuits (all gens),
            # for determining population model CI
            unique_obj, unique_indices = np.unique(all_obj,
                                            axis=0,
                                            return_index=True)
            unique_circuits = all_circuits[unique_indices]

            print(len(unique_circuits))

            file_name = "all_unique_obj.pkl"
            with open(folder_path + "/" + file_name, "wb") as fid:
                pickle.dump(unique_obj, fid)

            file_name = "all_unique_circuits.pkl"
            with open(folder_path + "/" + file_name, "wb") as fid:
                pickle.dump(unique_circuits, fid)

    # plot all circuits along pareto front
    if plot:
        for i, circuit in enumerate(population):
            graph_file_name = ("pareto_front_circuit_" 
                               + str(i))
            plot_graph(
                folder_path + "/" + graph_file_name,
                circuit[0]
            )
                   