import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics as stats
from load_files_pop import Z_20, Z_200, Z_2000

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
all_cell_prom_list = eval(all_cell_time_series_opt["single_cell_prominence"].tolist()[4])
# print(all_cell_prom_list)


Z_20_df = pd.DataFrame(data = Z_20, columns = ["plasmid_" + str(i) for i in range(5)])
Z_20_df_log = Z_20_df.copy()
Z_200_df = pd.DataFrame(data = Z_200, columns = ["plasmid_" + str(i) for i in range(9)])
Z_200_df_log = Z_200_df.copy()
Z_2000_df = pd.DataFrame(data = Z_2000, columns = ["plasmid_" + str(i) for i in range(5)])
Z_2000_df_log = Z_2000_df.copy()
for col in ["plasmid_" + str(i) for i in range(5)]:
#     max_val = max(Z_20_df[col])
    Z_20_df_log[col] = np.log10(Z_20_df_log[col])
    Z_200_df_log[col] = np.log10(Z_200_df_log[col])
    Z_2000_df_log[col] = np.log10(Z_2000_df_log[col])
# print(Z_20_df_log)

cells = np.arange(2000)

####plot all plasmids overlapping
path_save = repo_path + "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/"
color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']
fig1, axs1 = plt.subplots(1, 1, figsize=(6, 3))
# fig2, axs2 = plt.subplots(1, 1, figsize=(6, 4))
for plasmid in range(4):
    axs1.plot(cells, Z_2000[:, plasmid], label="plasmid "+str(plasmid))
#     # axs[plasmid, 1].hist(np.log10(Z_20[:, plasmid]), bins=20)
    axs1.set_xticks(cells)
    # axs1.set_xlabel("Cell number")
    axs1.set_xticks([])
    axs1.set_ylabel("Mean norm. amount \n of plasmid uptaken")
    axs1.set_ylim(top=20)
    # Z_20_df_log["plasmid_" + str(plasmid)].plot.density(ind=20, color=color_list[plasmid], ax=axs1)

    # axs1.set_xlim([-4, 2])
    # axs1.set_ylim([0, 0.8])
    # axs1.set_xlabel("Mean norm. amount of \n plasmid uptaken (log10 a.u.)")
    # axs1.set_ylabel("Probability density")
    # axs1.set_box_aspect(1)
# pulse_cell_list = [1, 7, 10, 13, 16, 17, 18, 19]
# label_list = ["cell with pulse"] + [""]*(len(pulse_cell_list)-1)
# for i, pulse_cell in enumerate(pulse_cell_list):
#     axs1.axvline(x=pulse_cell, color="k", linewidth=1, label=label_list[i])
axs1.legend()
plt.show()
# fig1.savefig(path_save + "plasmid_uptake_per_cell_2000_cell_zoomed.svg")
# fig1.savefig(path_save + "plasmid_uptake_distribution_20_cell.svg")



#### avg plasmid uptake vs prominence
# plasmids = ["plasmid_0", "plasmid_1", "plasmid_2", "plasmid_3", "plasmid_4"]
# plasmid_uptake_maxes = [max(Z_20_df[plasmid]) for plasmid in plasmids]
# # print(plasmid_uptake_maxes)
# fig, ax = plt.subplots(1, 1, figsize= (3, 3))
# avg_list = []
# for cell in range(20):
#     uptake_list = Z_20_df.iloc[cell].tolist()
#     avg_uptake = np.mean(uptake_list)
#     avg_list.append(avg_uptake)
# max_avg_uptake = max(avg_list)
# threshold = max_avg_uptake*0.6
# # print(max_avg_uptake*0.6)
# ax.plot(avg_list, all_cell_prom_list, linestyle="none", marker="o")
# ax.axvline(x=threshold, color="k", linewidth=1, linestyle="dotted", label="selected threshold")
# ax.set_xlabel("average plasmid uptake")
# ax.set_ylabel("prominence_rel")
# ax.set_ylim(bottom=0)
# plt.legend()
# plt.show()

#### avg plasmid uptake percentile vs prominence
# fig, ax = plt.subplots(1, 1, figsize= (3.5, 3))
# Z_20_df["average"] = np.mean(Z_20_df, axis=1)
# plasmid_uptake_list = Z_20_df["average"].tolist()
# plasmid_uptake_arr = np.array(plasmid_uptake_list)
# max_uptake = max(Z_20_df["average"])
# threshold = 60
# percentile = np.arange(10, 100, 10).tolist()
# percentile_fractions = [round(i/100, 2) for i in percentile]
# # fraction_included = [round(1.0-i, 2) for i in percentile_fractions]
# for i, percent in enumerate(percentile_fractions):
#     print(percent, ": ")
#     uptake_min = max_uptake*percent
#     print("threshold: ", uptake_min)
#     idx_uptake_included = np.where(plasmid_uptake_arr >= uptake_min)[0]
#     print("cells included: ", idx_uptake_included)
#     uptake_included = plasmid_uptake_arr[idx_uptake_included]
#     prom_included = [all_cell_prom_list[i] for i in idx_uptake_included]
#     print("prom_rel cells included: ", prom_included)
#     # pct_plot = [90, 80, 70, 60, 50, 40, 30, 20, 10]
#     ax.plot([percentile[i]]*len(prom_included), prom_included, linestyle="none", marker="o", color="grey")


# # print(max_avg_uptake*0.6)
# ax.axvline(x=threshold, color="k", linewidth=1, linestyle="dotted", label="selected threshold")
# ax.set_xlabel("percentile average plasmid uptake")
# ax.set_xticks(percentile)
# ax.set_ylabel("prominence_rel")
# ax.set_ylim(bottom=0)
# plt.legend(loc="center right")
# plt.show()

#pulses = 1, 7, 10, 13, 16-19
# best pulses = 10, 13, 16, 18
#plasmid 0: top 8 = 8 pulses 1, 7, 10, 13, 16-19
#plasmid 1: top 8 has 1 miss (cell 15 = #7, 17 = #8, and 7= #9)
#plasmid 2: top 8 has 2 misses (cell 12 = #7, 15 = #8, and 7= #9, 1 = #12)
#plasmid 3: top 8 has 1 miss (cell 12 = # 6, 19 = #7, 1 = #8, 7 = #9)
#plasmid 4: top 8 has has 2 misses (cell 12 = #5, 17 = #6, and 15= #7, 10 = #8, 1 = #11, 7 = #15)


#### percentile average plasmid uptake vs. % cells that pulse
# fig = plt.figure(figsize=(4, 3.5))
# Z_20_df["average"] = np.mean(Z_20_df, axis=1)
# print(Z_20_df)
# all_cell_prom_list = eval(all_cell_time_series_opt["single_cell_prominence"].tolist()[0])
# percentile = np.arange(10, 100, 10).tolist()
# percentile_fractions = [round(i/100, 2) for i in percentile]
# # fraction_included = [round(1.0-i, 2) for i in percentile_fractions]
# # print(fraction_included)

# plasmid_uptake_list = Z_20_df["average"].tolist()
# plasmid_uptake_arr = np.array(plasmid_uptake_list)
# print("full list: ", plasmid_uptake_arr)
# max_uptake = max(Z_20_df["average"])
# threshold_percent = 60
# print("threshold percent : ", threshold_percent)
# percent_pulse_cells_list = []
# for percent in percentile_fractions:
#     print(percent, ": ")
#     uptake_min = max_uptake*percent
#     print("threshold: ", uptake_min)
#     idx_uptake_included = np.where(plasmid_uptake_arr >= uptake_min)[0]
#     print("cells included: ", idx_uptake_included)
#     uptake_included = plasmid_uptake_arr[idx_uptake_included]
#     prom_included = [all_cell_prom_list[i] for i in idx_uptake_included]
#     # print("prom_rel cells included: ", prom_included)
#     pulse_cell_prom = [prom for prom in prom_included if prom > 0]
#     percent_pulse_cells = (len(pulse_cell_prom)/len(prom_included))*100
#     print("percent pulse cells", percent_pulse_cells, "%")
#     percent_pulse_cells_list.append(percent_pulse_cells)
# plt.plot(percentile, percent_pulse_cells_list)
# plt.axvline(x=threshold_percent, color="k", linewidth=1,linestyle="dotted", label="selected threshold")
# plt.legend()
# plt.xlabel("percentile average plasmid uptake")
# plt.xticks(percentile)
# plt.ylabel("% cells with a pulse")
# plt.show()