from GA import *
from itertools import combinations, permutations, product
import json
from amplifier_problem import Amplifier
from signal_conditioner_problem import SignalConditioner
from load_Z_mat_samples import Z_mat_list, Ref_list
from copy import deepcopy
import sys

index = int(sys.argv[1])

promo_node = 'P1'
min_dose = 5
max_dose = 75
dose_step = 5
dose_range = [[k] for k in range(min_dose, max_dose+1, dose_step)]
dose_range_2 = [tuple([i, j]) for i in range(min_dose, max_dose+1, dose_step) for j in range(min_dose, max_dose+1, dose_step)]

with open('./combo_search/main_combo_2_inhibitor.pkl', 'rb') as fid:
    combo = pickle.load(fid)

###################################
test_case = Amplifier

with open("./settings_amp_initial.json", "r") as fid:
    settings = json.load(fid)

if "Z_matrix" in settings:
    Z_mat = settings["Z_matrix"]
else:
    Z_mat = Z_mat_list[0]

if "reference" in settings:
    Ref_pop = settings["reference"]
else:
    Ref_pop = None


problem = test_case(
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
        mutate_scheme_prob=settings["mutate_scheme_prob"],
        pop=settings["pop"],
        Z_mat=Z_mat,
        Ref_pop=Ref_pop,
        num_processes=settings["num_processes"],
        obj_labels=settings["obj_labels"],
        max_time=settings["max_time"]
)


obj = []
combo_expanded = []
start = index * 150
cnt = 0
for topo in combo[start:]:
    if cnt == 150:
        break
    for dose in dose_range_2:
        for i, part in enumerate(topo.part_list):
            topo.dose[part] = dose[i]
        combo_expanded.append(deepcopy(topo))
        obj.append(problem.simulate_cell(topo))
    cnt += 1
obj = np.asarray(obj)

with open('./combo_search/amplifier_topo_2_inhibitor_p%d.pkl'%index, 'wb') as fid:
    pickle.dump(combo_expanded, fid)
with open('./combo_search/amplifier_obj_2_inhibitor_p%d.pkl'%index, 'wb') as fid:
    pickle.dump(obj, fid)
    
    