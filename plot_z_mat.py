import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from load_files_pop import Z_20

plt.style.use('/Users/kdreyer/Documents/Github/GraphGA/paper.mplstyle.py')
# print(Z_20)

Z_20_df = pd.DataFrame(data = Z_20, columns = ["plasmid_" + str(i) for i in range(5)])
# print(Z_20_df.sort_values("plasmid_4"))
Z_20_df_log = Z_20_df.copy()
for col in ["plasmid_" + str(i) for i in range(5)]:
    max_val = max(Z_20_df[col])
    # print(max_val*0.6)
    # print(Z_20_df[Z_20_df[col] >= 0.6*max_val])
    print(Z_20_df[Z_20_df[col] >= 1.8])
    # Z_20_df_log[col] = np.log10(Z_20_df_log[col])
# print(Z_20_df_log)

cells = np.arange(20)

# plot each plasmid separately
# fig1, axs1 = plt.subplots(5, 1, figsize=(6, 12))
# fig2, axs2 = plt.subplots(5, 1, figsize=(6, 12))
# for plasmid in range(5):
#     axs1[plasmid].plot(cells, Z_20[:, plasmid])
#     # axs[plasmid, 1].hist(np.log10(Z_20[:, plasmid]), bins=20)
#     axs1[plasmid].set_xticks(cells)
#     axs1[plasmid].set_xlabel("Cell number")
#     axs1[plasmid].set_ylabel("Mean norm. amount \n of plasmid uptaken")

#     Z_20_df_log["plasmid_" + str(plasmid)].plot.density(ind=20, color='k', alpha=0.5, ax=axs2[plasmid])
#     axs2[plasmid].set_xlim([-3, 1.5])
#     axs2[plasmid].set_xlabel("Mean norm. amount of \n plasmid uptaken (log10 a.u.)")
#     axs2[plasmid].set_ylabel("Probability density")
#     axs2[plasmid].set_box_aspect(1)
# plt.show()

#plot all plasmids overlapping
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