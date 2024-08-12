import numpy as np
import matplotlib.pyplot as plt
from load_files_pop import Z_20


cells = np.arange(20)
# print(cells)
fig, axs = plt.subplots(2, 1, figsize=(5, 7))
axs.ravel()
for plasmid in range(5):
    axs[0].plot(cells, Z_20[:, plasmid])
    axs[1].hist(np.log10(Z_20[:, plasmid]), bins=10)
axs[0].set_xticks(cells)
# plt.xticks(cells)
plt.show()