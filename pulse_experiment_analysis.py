import numpy as np
from flow_cytometry_calculations import (
    run_flow_cytometry_calculations
)

path_data = "/Users/kdreyer/Library/CloudStorage/OneDrive-NorthwesternUniversity/KatieD_LL/GCAD_Collab/Experimental_data_&_planning/240806_Pulse_flow_cytometry/formatted_single_cell_data/FITC_PacBlue/"
path_save = "/Users/kdreyer/Library/CloudStorage/OneDrive-NorthwesternUniversity/KatieD_LL/GCAD_Collab/Experimental_data_&_planning/240806_Pulse_flow_cytometry/results/all_cells_deg3_test/"#filtered_tfx70_paper_deg3/"
times = ["14", "18", "22", "26", "38", "42", "46"]
conditions = ["pulse1", "pulse3", "pulse4", "pulse5", "reference"]
background_subtractions = [
    np.mean([90.8, 103, 93.1]), 
    np.mean([65, 70.5, 67]), 
    np.mean([22, 21.3, 23.7]),  
    np.mean([23.1, 22.4, 25]),
    np.mean([30.1, 28.3, 29.8]),
    np.mean([38, 33.3, 34.7]),
    np.mean([49.9, 50.4, 51.7])
]
background_subtractions_stderr = [
    np.std([90.8, 103, 93.1], ddof=1)/np.sqrt(3),  
    np.std([65, 70.5, 67], ddof=1)/np.sqrt(3), 
    np.std([22, 21.3, 23.7], ddof=1)/np.sqrt(3),
    np.std([23.1, 22.4, 25], ddof=1)/np.sqrt(3),
    np.std([30.1, 28.3, 29.8], ddof=1)/np.sqrt(3),
    np.std([38, 33.3, 34.7], ddof=1)/np.sqrt(3),
    np.std([49.9, 50.4, 51.7], ddof=1)/np.sqrt(3)
]
MEFL_conversions = [
    323.345243705024, 
    631.120965591862, 
    2543.43441686963,
    3121.36550948891,
    6829.80404249601,
    6812.01906915132,
    6880.88864816803
]
MEFL_conversions_stderr = [
    2.41597846744345, 
    5.29854152370918, 
    29.5262709531106,
    38.1100625157131,
    92.0650011641531,
    99.1214182171212,
    101.84815344387
]

run_flow_cytometry_calculations(
    path_data,
    path_save,
    conditions,
    times,
    background_subtractions,
    background_subtractions_stderr,
    MEFL_conversions,
    MEFL_conversions_stderr,
    percentile=None,#70 
    save_data=False,#False,
    log_scale=False,
    split_plot=False
)