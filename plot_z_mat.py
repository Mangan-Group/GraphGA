import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics as stats
from load_files_pop import Z_20

plt.style.use('/Users/kdreyer/Documents/Github/GraphGA/paper.mplstyle.py')
# print(Z_20)

repo_path = "/Users/kdreyer/Library/CloudStorage/OneDrive-NorthwesternUniversity/KatieD_LL/GCAD_Collab/GA_results/"
#3 obj [:5]
results_path = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/Pulse_pop_DsRED_inhibitor_3obj_126h_ZF1_ZF2_seed_0/Results_analysis/"
file_name = "all_cell_metrics_low_t_pulse.csv"

#t pulse [1:2]
# results_path = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/Pulse_pop_DsRED_inhibitor_t_pulse_126h_ZF1_ZF2_seed_0/Results_analysis_sub_opt/"
# file_name = "all_cell_metrics_sub_opt.csv"

all_cell_results = pd.read_csv(repo_path+results_path+file_name)
all_cell_time_series_opt = all_cell_results.copy().iloc[:5]
# print(all_cell_time_series_opt["single_cell_prominence"].to_list)
all_cell_prom_list = eval(all_cell_time_series_opt["single_cell_prominence"].tolist()[0])
# print(all_cell_prom_list)


Z_20_df = pd.DataFrame(data = Z_20, columns = ["plasmid_" + str(i) for i in range(5)])
Z_20_df_log = Z_20_df.copy()
for col in ["plasmid_" + str(i) for i in range(5)]:
#     max_val = max(Z_20_df[col])
    Z_20_df_log[col] = np.log10(Z_20_df_log[col])
# print(Z_20_df_log)

cells = np.arange(20)

####plot all plasmids overlapping
# path_save = "/Users/kdreyer/Documents/Github/GraphGA/GA_results/Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/"
# color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']
# fig1, axs1 = plt.subplots(1, 1, figsize=(6, 3))
# fig2, axs2 = plt.subplots(1, 1, figsize=(6, 4))
# for plasmid in range(5):
#     axs1.plot(cells, Z_20[:, plasmid])
#     # axs[plasmid, 1].hist(np.log10(Z_20[:, plasmid]), bins=20)
#     axs1.set_xticks(cells)
#     axs1.set_xlabel("Cell number")
#     axs1.set_ylabel("Mean norm. amount \n of plasmid uptaken")

#     Z_20_df_log["plasmid_" + str(plasmid)].plot.density(ind=20, color=color_list[plasmid], ax=axs2)
#     axs2.set_xlim([-3, 1.5])
#     axs2.set_xlabel("Mean norm. amount of \n plasmid uptaken (log10 a.u.)")
#     axs2.set_ylabel("Probability density")
#     axs2.set_box_aspect(1)
# pulse_cell_list = [1, 7, 10, 13, 16, 17, 18, 19]
# label_list = ["cell with pulse"] + [""]*(len(pulse_cell_list)-1)
# for i, pulse_cell in enumerate(pulse_cell_list):
#     axs1.axvline(x=pulse_cell, color="k", linewidth=1, label=label_list[i])
# axs1.legend()
# # plt.show()
# fig1.savefig(path_save + "plasmid_uptake_per_cell.svg")
# fig2.savefig(path_save + "plasmid_uptake_distribution.svg")



#### stdev plasmid uptake vs prominence
plasmids = ["plasmid_0", "plasmid_1", "plasmid_2", "plasmid_3", "plasmid_4"]
# plasmid_uptake_maxes = [max(Z_20_df[plasmid]) for plasmid in plasmids]
# # print(plasmid_uptake_maxes)
# stdev_list = []
# stdev_list_norm = []
# for cell in range(20):
#     uptake_list = Z_20_df.iloc[cell].tolist()
#     uptake_list_norm = [uptake_list[i]/plasmid_uptake_maxes[i] for i in range(len(uptake_list))]
#     # uptake_list_norm = 
#     stdev_uptake = stats.stdev(uptake_list)
#     stdev_uptake_norm = stats.stdev(uptake_list_norm)
#     stdev_list.append(stdev_uptake)
#     stdev_list_norm.append(stdev_uptake_norm)
#     print(stdev_uptake, stdev_uptake_norm, all_cell_prom_list[cell])
# plt.plot(stdev_list, all_cell_prom_list, linestyle="none", marker="o")
# plt.show()
# plt.plot(stdev_list_norm, all_cell_prom_list, linestyle="none", marker="o")
# plt.show()


#### plasmid uptake vs prominence
# fig = plt.figure(figsize=(4, 3.5))
# for plasmid in plasmids:
#     idx = Z_20_df.sort_values(plasmid)[plasmid].index.tolist()
#     all_cell_prom_list = eval(all_cell_time_series_opt["single_cell_prominence"].tolist()[0])
#     all_cell_prom_list_sorted = [all_cell_prom_list[i] for i in idx]
#     plt.plot(Z_20_df.sort_values(plasmid)[plasmid], all_cell_prom_list_sorted, linestyle="none", marker="o")
# plt.axvline(x=1.8, color="k", linewidth=2, linestyle="dotted", label="selected threshold")
# plt.xlabel("plasmid uptake")
# plt.ylabel("prominence_rel")
# plt.ylim(bottom=0)
# plt.legend()
# plt.show()

#pulses = 1, 7, 10, 13, 16-19
# best pulses = 10, 13, 16, 18
#plasmid 0: top 8 = 8 pulses 1, 7, 10, 13, 16-19
#plasmid 1: top 8 has 1 miss (cell 15 = #7, 17 = #8, and 7= #9)
#plasmid 2: top 8 has 2 misses (cell 12 = #7, 15 = #8, and 7= #9, 1 = #12)
#plasmid 3: top 8 has 1 miss (cell 12 = # 6, 19 = #7, 1 = #8, 7 = #9)
#plasmid 4: top 8 has has 2 misses (cell 12 = #5, 17 = #6, and 15= #7, 10 = #8, 1 = #11, 7 = #15)



#### percentile plasmid uptake vs. % cells that pulse
# fig = plt.figure(figsize=(4, 3.5))
all_cell_prom_list = eval(all_cell_time_series_opt["single_cell_prominence"].tolist()[0])
plasmids = ["plasmid_0", "plasmid_1", "plasmid_2", "plasmid_3", "plasmid_4"]
percentile = np.arange(0, 100, 10).tolist()
percentile_fractions = [i/100 for i in percentile]
fraction_included = [round(1.0-i, 2) for i in percentile_fractions]
# print(fraction_included)
for plasmid in plasmids:
    plasmid_uptake_list = Z_20_df[plasmid].tolist()
    plasmid_uptake_arr = np.array(plasmid_uptake_list)
    # print("full list: ", plasmid_uptake_arr)
    max_uptake = max(Z_20_df[plasmid])
    threshold = np.mean(plasmid_uptake_arr)*1.8
    threshold_percent = (1 - threshold/max_uptake)*100
    print("threshold percent : ", threshold_percent, threshold)
    percent_pulse_cells_list = []
    for percent in fraction_included:
        print(percent, ": ")
        uptake_min = max_uptake*percent

        # print("threshold: ", uptake_min)
        idx_uptake_included = np.where(plasmid_uptake_arr >= uptake_min)[0]
        print("cells included: ", idx_uptake_included)
        uptake_included = plasmid_uptake_arr[idx_uptake_included]
        prom_included = [all_cell_prom_list[i] for i in idx_uptake_included]
        # print("prom_rel cells included: ", prom_included)
        pulse_cell_prom = [prom for prom in prom_included if prom > 0]
        percent_pulse_cells = (len(pulse_cell_prom)/len(prom_included))*100
        print("percent pulse cells", percent_pulse_cells, "%")
        percent_pulse_cells_list.append(percent_pulse_cells)
    pct_plot = [90, 80, 70, 60, 50, 40, 30, 20, 10]
    pct_plot_str = [str(i) for i in pct_plot]#reversed(percentile[1:]).tolist()
    plt.plot(pct_plot_str, percent_pulse_cells_list[1:], label=str(round(threshold_percent, 1)))

plt.legend(title="threshold %")
plt.xlabel("percentile plasmid uptake included")
plt.ylabel("% cells with a pulse")
plt.show()