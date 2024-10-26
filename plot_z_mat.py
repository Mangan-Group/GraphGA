import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statistics as stats
import scipy.stats as sstats
from copy import deepcopy
import matplotlib.mlab as mlab
from load_files_pop import Z_20, Z_200, Z_2000

plt.style.use('/Users/kdreyer/Documents/Github/GraphGA/paper.mplstyle.py')
# print(Z_20)

repo_path = "/Users/kdreyer/Library/CloudStorage/OneDrive-NorthwesternUniversity/KatieD_LL/GCAD_Collab/GA_results/"
#3 obj [:5]
results_path = "Pulse_seed_pop_DsRED_inhibitor/ZF1_ZF2_only/Pulse_pop_DsRED_inhibitor_3obj_126h_ZF1_ZF2_seed_0/Results_analysis/"
file_name = "all_cell_metrics_low_t_pulse.csv"
path = "/Users/kdreyer/Library/CloudStorage/OneDrive-NorthwesternUniversity/KatieD_LL/GCAD_Collab/Selected_GA_results_paper/Pulse_pop/Experiment_synTF1&2/200_cell_all_cell_results_df.pkl"
pulse_exp_200 = pd.read_pickle(path)
print(pulse_exp_200)
Z_20_df = pd.DataFrame(data = Z_20, columns = ["plasmid_" + str(i) for i in range(5)])
Z_200_df = pd.DataFrame(data = Z_200, columns = ["plasmid_" + str(i) for i in range(9)])
Z_200_df = Z_200_df.drop(labels=["plasmid_" + str(i) for i in range(5,9)], axis=1)
Z_2000_df = pd.DataFrame(data = Z_2000, columns = ["plasmid_" + str(i) for i in range(5)])

####plot all plasmids overlapping
def plot_all_plasmid_uptake(Z_mat1, Z_mat2, path_save):
    num_cells1 = len(Z_mat1)
    cells1 = np.arange(0, num_cells1)
    num_plasmids1 = np.shape(Z_mat1)[-1]
    color_list = plt.rcParams['axes.prop_cycle'].by_key()['color']
    num_cells2 = len(Z_mat2)
    cells2 = np.arange(0, num_cells2)
    num_plasmids2 = np.shape(Z_mat2)[-1]
    # fig1, axs1 = plt.subplots(1, 1, figsize=(6, 3))
    fig2, axs2 = plt.subplots(1, 1, figsize=(6, 4))
    # for plasmid in range(num_plasmids):
    #     axs1.plot(cells, Z_mat[:, plasmid], label="plasmid "+str(plasmid))
    # axs1.set_xticks(cells)
    # axs1.set_xlabel("Cell number")
    # axs1.set_xticks([])
    # axs1.set_ylabel("Mean norm. amount \n of plasmid uptaken")
    # axs1.legend()
    Z_df1 = pd.DataFrame(data = Z_mat1, columns = ["plasmid_" + str(i) for i in range(num_plasmids1)])
    Z_df1_log = Z_df1.copy()
    Z_df2 = pd.DataFrame(data = Z_mat2, columns = ["plasmid_" + str(i) for i in range(num_plasmids2)])
    Z_df2_log = Z_df2.copy()
    for i, column in enumerate(Z_df1.columns):
        Z_df1_log[column] = np.log10(Z_df1_log[column])
    #     mu1, sigma1 = sstats.norm.fit(Z_df1_log[column])
    #     points1 = np.linspace(sstats.norm.ppf(0.01,loc=mu1,scale=sigma1),
    #                     sstats.norm.ppf(0.9999,loc=mu1,scale=sigma1),1000)
    #     pdf1 = sstats.norm.pdf(points1,loc=mu1,scale=sigma1)
    #     axs2.plot(points1, pdf1, color='k')
        
    #     # Z_df1_log[column].plot.density(ind=1000, color="k", ax=axs2)
    #     axs2 = Z_df1_log[column].plot.hist(bins=100, density=True, color="k")

        axs2.hist(Z_df1_log[column], bins=20, color="k", density=True)
    for i, column in enumerate(Z_df2.columns):
        Z_df2_log[column] = np.log10(Z_df2_log[column])
        n, bins, patches = axs2.hist(Z_df2_log[column], bins=10, color="red", density=True)

        mu2, sigma2 = sstats.norm.fit(Z_df2_log[column])
        # points2 = np.linspace(sstats.norm.ppf(0.01,loc=mu2,scale=sigma2),
        #                 sstats.norm.ppf(0.9999,loc=mu2,scale=sigma2),1000)
        # pdf2 = sstats.norm.pdf(points2,loc=mu2,scale=sigma2) * sum(n * np.diff(bins))
        # y = sstats.norm.pdf(bins, mu2, sigma2) * sum(n * np.diff(bins))
        # axs2.plot(bins, y, '-o', linewidth=2)
        # axs2.plot(points2, pdf2, color='k')

    #     Z_df2_log[column].plot.density(ind=1000, color="red", ax=axs2) #color=color_list[i+num_plasmids1])
        # axs2 = Z_df2_log[column].plot.hist(bins=20, density=True, color="red")
    

    # Z_df1_log["average"] = Z_df1_log.mean(axis=1)
    # axs2 = Z_df1_log["average"].plot.hist(bins=50, density=True, color="k")
    # mu1, sigma1 = sstats.norm.fit(Z_df1_log["average"])
    # points1 = np.linspace(sstats.norm.ppf(0.01,loc=mu1,scale=sigma1),
    #                 sstats.norm.ppf(0.9999,loc=mu1,scale=sigma1),100)
    # pdf1 = sstats.norm.pdf(points1,loc=mu1,scale=sigma1)
    # axs2.plot(points1, pdf1, color='k')

    # Z_df2_log["average"] = Z_df2_log.mean(axis=1)
    # axs2 = Z_df2_log["average"].plot.hist(bins=5, density=True, color="red")
    # mu2, sigma2 = sstats.norm.fit(Z_df2_log["average"])
    # points2 = np.linspace(sstats.norm.ppf(0.01,loc=mu2,scale=sigma2),
    #                 sstats.norm.ppf(0.9999,loc=mu2,scale=sigma2),100)
    # pdf2 = sstats.norm.pdf(points2,loc=mu2,scale=sigma2)
    # axs2.plot(points2, pdf2, color='r')


    axs2.set_xlim([-4, 2])
    # axs2.set_ylim([0, 0.8])
    axs2.set_xlabel("Mean norm. amount of \n plasmid uptaken (log10 a.u.)")
    axs2.set_ylabel("Probability density")
    axs2.set_box_aspect(1)
    plt.show()
    # fig1.savefig(path_save + "plasmid_uptake_per_cell_" +str(num_cells)+ "_cell.svg")
    # fig2.savefig(path_save + "plasmid_uptake_distribution_" +str(num_cells)+ "_cell.svg")

# plot_all_plasmid_uptake(Z_200, Z_20, "")

# pulse_cell_list = [1, 7, 10, 13, 16, 17, 18, 19]
# label_list = ["cell with pulse"] + [""]*(len(pulse_cell_list)-1)
# for i, pulse_cell in enumerate(pulse_cell_list):
#     axs1.axvline(x=pulse_cell, color="k", linewidth=1, label=label_list[i])



#### avg plasmid uptake vs prominence
# def plot_avg_uptake_prom_rel():
# plasmids = ["plasmid_0", "plasmid_1", "plasmid_2", "plasmid_3", "plasmid_4"]
# plasmid_uptake_maxes = [max(Z_20_df[plasmid]) for plasmid in plasmids]
# # print(plasmid_uptake_maxes)
# fig, ax = plt.subplots(1, 1, figsize= (3, 3))
# avg_list = []
# for cell in range(20):
#     uptake_list = Z_20_df.iloc[cell].tolist()a
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
# def plot_avg_uptake_pct_prom_rel():
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
# Z_200_df["average"] = np.mean(Z_200_df, axis=1)
# # print(Z_200_df)
# all_cell_prom_list = pulse_exp_200["single_cell_prominence"].tolist()[4]
# percentile = np.arange(10, 100, 10).tolist()
# percentile_fractions = [round(i/100, 2) for i in percentile]
# # fraction_included = [round(1.0-i, 2) for i in percentile_fractions]
# # # print(fraction_included)

# plasmid_uptake_list = Z_200_df["average"].tolist()
# plasmid_uptake_arr = np.array(plasmid_uptake_list)
# # print("full list: ", plasmid_uptake_arr)
# max_uptake = max(Z_200_df["average"])
# threshold_percent = 70
# # print("threshold percent : ", threshold_percent)
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
# fig, ax = plt.subplots(1, 1, figsize=(1.4, 1.4))

# ax.plot(percentile, percent_pulse_cells_list, color="k")
# ax.axvline(x=threshold_percent, color="k",linestyle="dashed", label="selected threshold")
# # plt.legend()
# ax.set_xlabel("pct avg plasmid uptake")
# ax.set_xticks([0, 20, 40, 60, 80])
# ax.set_ylabel("% cells with a pulse")
# ax.set_ylim(bottom=0)
# ax.set_yticks([0, 20, 40, 60, 80, 100])
# ax.set_box_aspect(1)
# # plt.show()
# plt.savefig("/Users/kdreyer/Desktop/200_cell_plasmid_pct_paper.svg")
