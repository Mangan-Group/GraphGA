import sys
import json
import numpy as np
import pickle
from multiprocessing import Pool
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from pulse_generator_problem import PulseGenerator
from load_files_pop import Z_200
from load_Z_mat_samples import Z_mat_list, Ref_list
from saving import make_main_directory
from GA import sampling
from GA_setup import (
    single_obj_GA,
    multi_obj_GA
)

def run(
        testcase: object,
        settings: dict,
    ):
    '''Run the genetic algorithm for a given test case'''

    folder_path = make_main_directory(settings)
    # if a reference needs to be specified for the
    # corresponding Z matrix (one of the Z_20 samples
    # that isn't usually used), both ref and Z will
    # need to be specified. otherwise, can just specify 
    # a different size Z matrix (Z_200 or Z_2000) & the
    # corresponding reference will be used
    # else, Z_mat = Z_20, the population model default
    # (first Z matrix in Z_mat_list)
    if "reference" in settings:
        Ref_pop = settings["reference"]
    else:
        Ref_pop = None
    if "Z_matrix" in settings:
        Z_mat = settings["Z_matrix"]
    else:
        Z_mat = Z_mat_list[0]

    # if "mean" in settings:
    #     problem = testcase(
    #         promo_node=settings["promo_node"],
    #         dose_specs=settings["dose_specs"],
    #         max_part=settings["max_part"],
    #         inhibitor=settings["inhibitor"],
    #         DsRed_inhibitor=settings["DsRed_inhibitor"],
    #         num_dict=settings["num_dict"],
    #         n_gen=settings["n_gen"],
    #         probability_crossover=settings["probability_crossover"],
    #         probability_mutation=settings["probability_mutation"],
    #         mutate_dose=settings["mutate_dose"],
    #         pop=settings["pop"],
    #         mean=settings["mean"],
    #         Z_mat=Z_mat,
    #         Ref_pop=Ref_pop,
    #         num_processes=settings["num_processes"],
    #         obj_labels=settings["obj_labels"],
    #         max_time=settings["max_time"]
    # )

    # else:
    problem = testcase(
        promo_node=settings["promo_node"],
        dose_specs=settings["dose_specs"],
        max_part=settings["max_part"],
        inhibitor=settings["inhibitor"],
        DsRed_inhibitor=settings["DsRed_inhibitor"],
        num_dict=settings["num_dict"],
        n_gen=settings["n_gen"],
        probability_crossover=settings["probability_crossover"],
        probability_mutation=settings["probability_mutation"],
        mutate_dose=settings["mutate_dose"],
        pop=settings["pop"],
        Z_mat=Z_mat,
        Ref_pop=Ref_pop,
        num_processes=settings["num_processes"],
        obj_labels=settings["obj_labels"],
        max_time=settings["max_time"]
)

    population = sampling(
        problem.promo_node,
        problem.num_dict,
        problem.min_dose,
        problem.max_dose,
        problem.dose_interval,
        inhibitor=problem.inhibitor
    )

    # save initial population
    file_name = "initial_population.pkl"
    with open(folder_path + "/" + file_name, "wb") as fid:
        pickle.dump(population, fid)

    num_circuits = len(population)

    # raise an exception if num_circuits is odd
    if num_circuits % 2 != 0:
        raise Exception("Population size must be an even number")
    
    # calculate objective for each circuit 
    # in initial population
    if settings["pop"]:
        topologies = [g[0] for g in population]
        with Pool(settings["num_processes"]) as pool:
            obj_list = pool.imap(problem.func, topologies)
            # obj_all_cells_dict_list = pool.imap(problem.func, topologies)

            pool.close()
            pool.join()
        obj_list = list(obj_list)
        obj = np.asarray(obj_list)
    else:
        obj = np.asarray(
            [problem.func(g[0]) for g in population])

    # run multi-objective GA if multiple objectives
    # (determined by type of first first objective 
    # in array- will be array if more than one)
    if isinstance (obj[0], np.ndarray):
        print('objective is an array- using '+
              'multi-objective optimization')
        multi_obj_GA(
            folder_path,
            problem,
            population,
            num_circuits,
            obj,
            settings["plot"]
        )
        
    else:
        print('objective is not an array- using '+
              'single objective optimization')
        single_obj_GA(
            folder_path,
            problem,
            population,
            num_circuits,
            obj,
        )


if __name__ == "__main__":
    settings_path = sys.argv[1]
    settings_file = open(settings_path, encoding="utf-8")
    settings = json.load(settings_file)
    if "pop_size" in settings:
        total_population = settings["pop_size"]*2
        one_part_circuits = int(total_population*settings["pop_ratio"])
        two_part_circuits = int(total_population - one_part_circuits)
    else:
        one_part_circuits = settings["num_dict"]["1"]
        two_part_circuits = settings["num_dict"]["2"]
    settings["num_dict"] = {1: one_part_circuits, 2: two_part_circuits}

    if settings["test_case"] == "Amplifier":
        test_case = Amplifier
    elif settings["test_case"] == "SignalConditioner":
        test_case = SignalConditioner
    elif settings["test_case"] == "PulseGenerator":
        test_case = PulseGenerator
    else:
        raise Exception("Error: test case not defined")

    for seed in range(0, 1):
        settings["folder_name"] = settings["folder_name"].removesuffix("_seed_" + str(seed-1))
        np.random.seed(seed)
        settings["seed"] = seed
        settings["folder_name"] = settings["folder_name"] + "_seed_" + str(seed)

        run(test_case, settings)
        print("seed "+str(seed)+" complete")

    # seed = 0
    # np.random.seed(seed)
    # settings["seed"] = seed
    # for i, zmat in enumerate(Z_mat_list[1:]):
    #     settings["folder_name"] = settings["folder_name"].removesuffix("_Z20_" + str(i))
    #     settings["folder_name"] = settings["folder_name"] + "_Z20_" + str(i+1)
    #     Ref_pop = Ref_list[i+1]
    #     settings["reference"] = Ref_pop
    #     settings["Z_matrix"] = zmat
    #     run(test_case, settings)
        # print("Z_20 "+str(i+1)+" run complete")

    # settings["Z_matrix"] = Z_200
    # run(test_case, settings)