import pickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import spearmanr, ks_2samp

plt.style.use('/Users/kdreyer/Documents/Github/GCAD_Test_Cases/paper.mplstyle.py')
orange_ = [i/255 for i in [230, 159, 0]]
sky_blue = [i/255 for i in [86, 180, 233]]
pink_ = [i/255 for i in [204, 121, 167]]
bluish_green = [i/255 for i in [0, 158, 115]]
vermillion = [i/255 for i in [213, 94, 0]]

def load_files(
        base_path: str,
        conditions: list,
        times: list):
    """Loads the experimental data files
    and saves as a dict of dicts of dfs 
    (one dict per condition with one df per
    time point)."""

    experiment_dfs_dict = {}
    # loop through conditions (pulses and reference)
    for condition in conditions:
        experiment_dfs_dict[condition] = {}
        # loop through each time point
        for time in times:
            # load data for conditino and time point
            data_path = base_path + condition + "_time" + time + ".csv"
            df = pd.read_csv(data_path)
            experiment_dfs_dict[condition][time] = df

    return experiment_dfs_dict


def filter_tfx_single(
        df_condition: pd.DataFrame,
        percentile: float
):
    """Filter dataset (one condition at one time 
    point) to only keep cells with Pacific Blue 
    (transfection control) expression >= threshold."""

    df_to_filter = df_condition.copy()

    # calculate percentile as a fraction (decimal)
    percentile_fraction = percentile/100
    replicates = ["rep1", "rep2", "rep3"]
    df_filtered_list = []
    # loop through replicates
    for rep in replicates:
        # subset data for specific replicate
        df_rep = df_to_filter.filter(regex=rep)
        # calculate maximum Pacific Blue expression
        max_pac_blue = df_rep[rep+"_Pac_Blue"].max()
        # subset replicate data to only keep cells with
        # Pacific Blue expression >= maximum*percentile
        df_rep_filtered = df_rep[
            df_rep[rep+"_Pac_Blue"] >=
            max_pac_blue*percentile_fraction
        ].reset_index(drop=True)
        df_filtered_list.append(df_rep_filtered)

    # concatenate replicate data
    df_all_filtered = pd.concat(df_filtered_list, axis=1)

    return df_all_filtered


def filter_tfx_condition(
        df_condition_dict: dict,
        time_points: list,
        percentile: float
):
    """Filter dataset (one condition all time 
    point) to only keep cells with Pacific Blue 
    (transfection control) expression >= threshold."""
    
    df_condition_dict_filtered = {}
    # loop through time points/ data for each
    # condition at each time point
    for i, time in enumerate(time_points):
        df_time_point = df_condition_dict[time].copy()
        # filter data for condition at time point
        df_time_point = filter_tfx_single(
            df_time_point,
            percentile
        )
        # add filtered data to dict
        df_condition_dict_filtered[time] = (
            df_time_point
        )

    return df_condition_dict_filtered


def filter_tfx_experiment(
        experiment_dfs_dict: dict,
        conditions: list,
        times: list,
        percentile: float
):
    """Filter dataset (all conditions all time 
    point) to only keep cells with Pacific Blue 
    (transfection control) expression >= threshold."""
    
    experiment_dfs_filtered_dict = {}
    # loop through conditions 
    for condition in conditions:
        df_condition_dict = experiment_dfs_dict[condition]
        # filter data for condition for all time points
        df_condition_dict_filtered = filter_tfx_condition(
            df_condition_dict,
            times, 
            percentile
        )
        # add filtered data to dict
        experiment_dfs_filtered_dict[condition] = (
            df_condition_dict_filtered
        )

    return experiment_dfs_filtered_dict


def background_subtract_single(
        df_condition: pd.DataFrame, 
        background_subtraction: float,
        background_subtraction_stderr: float
):
    """Subtract background expression from FITC (MFI) (one
    condition one time point) and add standard error to data
    (since there is no error associated with the individual 
    cell data, the standard error is the same as that of the 
    background expression)."""
    
    df_condition_background_subtracted = df_condition.copy()

    replicates = ["rep1", "rep2", "rep3"]
    # loop through replicates 
    for rep in replicates:
        # subtract the background expression (averaged)
        df_condition_background_subtracted[rep + "_FITC_background_subtracted"] = (
            df_condition_background_subtracted[rep+"_FITC"] - background_subtraction
        )
        # add standard error value
        df_condition_background_subtracted[rep + "_FITC_background_subtracted_stderr"] = (
        background_subtraction_stderr)
    
    return df_condition_background_subtracted


def background_subtract_condition(
        df_condition_dict: dict,
        time_points: list,
        background_subtractions: list,
        background_subtractions_stderr: list
):
    """Subtract background expression from FITC (MFI) (one
    condition all time points) and add standard error to data."""
    
    df_condition_dict_background_subtracted = {}
    # loop through time points and dfs for each time point
    for i, time in enumerate(time_points):
        df_time_point = df_condition_dict[time].copy()
        # subtract background expression
        df_time_point = background_subtract_single(
            df_time_point,
            background_subtractions[i],
            background_subtractions_stderr[i]
        )
        # add dataframe to dict
        df_condition_dict_background_subtracted[time] = (
          df_time_point  
        )

    return df_condition_dict_background_subtracted


def convert_to_MEFLs_single(
        df_condition: pd.DataFrame,
        MEFL_conversion: float,
        MEFL_conversion_stderr: float
):
    """Convert FITC (MFI) values to MEFLs (one condition 
    one time point) and calculate standard error."""
    
    df_condition_MEFLs = df_condition.copy()

    replicates = ["rep1", "rep2", "rep3"] 
    # loop through replicates 
    for rep in replicates:
        # convert to MEFLs
        df_condition_MEFLs[rep + "_MEFLs"] = (
            df_condition_MEFLs[rep + "_FITC_background_subtracted"]*MEFL_conversion
        )
        # calculate standard error
        df_condition_MEFLs[rep + "_MEFLs_stderr"] = (
            df_condition_MEFLs[rep + "_MEFLs"]*np.sqrt(
                (df_condition_MEFLs[rep + "_FITC_background_subtracted_stderr"]/
                 df_condition_MEFLs[rep + "_FITC_background_subtracted"])**2 + 
                 (MEFL_conversion_stderr/MEFL_conversion)**2)
        )

    return df_condition_MEFLs


def convert_to_MEFLs_condition(
        df_condition_dict,
        time_points,
        MEFL_conversions,
        MEFL_conversions_stderr
):
    """Convert FITC (MFI) values to MEFLs (one condition 
    all time points) and calculate standard error."""
    
    df_condition_dict_MEFLs = {}
    # loop through time points and dfs for each time point
    for i, time in enumerate(time_points):
        df_time_point = df_condition_dict[time].copy()
        # convert to MEFLs
        df_time_point = convert_to_MEFLs_single(
            df_time_point,
            MEFL_conversions[i],
            MEFL_conversions_stderr[i]
        )
        # add dataframe to dict
        df_condition_dict_MEFLs[time] = (
            df_time_point
        )

    return df_condition_dict_MEFLs


def calculate_MEFLs_experiment(
        experiment_dfs_dict: dict,
        conditions: list,
        times: list,
        background_subtractions: list,
        background_subtractions_stderr: list,
        MEFL_conversions: list,
        MEFL_conversions_stderr: list
):
    """Convert FITC (MFI) values to MEFLs (all conditions 
    all time points) and calculate standard error."""

    experiment_dfs_calculations_dict = {}
    # loop through conditions
    for condition in conditions:
        df_condition_dict = experiment_dfs_dict[condition]
        # subtract background expression
        df_condition_dict_background_subtracted = background_subtract_condition(
            df_condition_dict,
            times,
            background_subtractions,
            background_subtractions_stderr
        )
        # convert to MEFLs
        df_condition_dict_MEFLs = convert_to_MEFLs_condition(
            df_condition_dict_background_subtracted,
            times,
            MEFL_conversions,
            MEFL_conversions_stderr
        )
        # add to dict
        experiment_dfs_calculations_dict[condition] = df_condition_dict_MEFLs

    return experiment_dfs_calculations_dict


def format_and_save_data(
        path: str, 
        experiment_dfs_calculations_dict: dict,
        conditions: list,
        time_points: list,
        save_data: bool
):
    """Format dataframes and save data: all calculations for each
    time point/condition and and just MEFLs for each condition."""

    replicates = ["rep1", "rep2", "rep3"]
    experiment_dfs_MEFLs_dict = {}
    experiment_dfs_MEFLs_combined_dict = {}
    # loop through conditions 
    for condition in conditions:
        df_condition_MEFLs_series_list = []
        df_condition_MEFLs_combined_series_list = []
        # loop through time points
        for time in time_points:
            # get df for condition and time point
            df_single = experiment_dfs_calculations_dict[condition][time]
            # save full df
            if save_data:
                df_single.to_csv(path + condition + "_calculations_" + time + "h.csv")
            # subset data to only keep MEFLs/MEFLs standard error columns for 
            # each replicate at the time point and add to list
            df_single_MEFLs = df_single.filter(regex="MEFLs")
            df_single_MEFLs = df_single_MEFLs.add_suffix("_" + time + "h")
            df_condition_MEFLs_series_list.append(df_single_MEFLs)

            # create list of dfs with MEFLs only for each replicate for the 
            # time point (get rid of error)
            df_single_MEFLs_reps_list = [df_single_MEFLs[rep + "_MEFLs_" + time + "h"]
                                    for rep in replicates
            ]
            df_single_MEFLs_reps_list = [df.dropna() for df in df_single_MEFLs_reps_list]
            # concatenate dfs for each replicate (MEFLs only) and 
            # add to list for condition
            df_single_MEFLs_combined = pd.concat(df_single_MEFLs_reps_list, ignore_index=True, axis=0)
            df_condition_MEFLs_combined_series_list.append(df_single_MEFLs_combined)

        # concatenate MEFLs/MEFLs standard error dfs for each 
        # time point and add to list
        df_condition_MEFLs_series = pd.concat(
            df_condition_MEFLs_series_list, axis=1
        )
        if save_data:
            # save df with MEFLs and MEFLs standard error
            df_condition_MEFLs_series.to_csv(path + condition + "_MEFLs.csv")
            
        # save to dict
        experiment_dfs_MEFLs_dict[condition] = df_condition_MEFLs_series

        # concatenate MEFLs only dfs for each time point
        df_condition_MEFLs_combined_series = pd.concat(
            df_condition_MEFLs_combined_series_list, axis=1,
            keys=[t+"h" for t in time_points]
        )
        # save to dict
        experiment_dfs_MEFLs_combined_dict[condition] = df_condition_MEFLs_combined_series

    return experiment_dfs_MEFLs_dict, experiment_dfs_MEFLs_combined_dict


def plot_condition_combined(
        path_figure: str,
        condition: str,
        time_points: list,
        df_condition_MEFLs_series: pd.DataFrame,
):
    """Plots individual cell MEFLs at each time point for
    a given condition (all 3 replicates on the same plot/same
    color)."""

    fig, ax = plt.subplots(1, 1, figsize=(1.75, 1.75))# (3, 2.75)
    # loop through time points 
    for time in time_points:
        # combine data from all 3 replicates for time point
        condition_time_MEFLs = (df_condition_MEFLs_series["rep1_MEFLs_"+time+"h"].tolist() + 
                                df_condition_MEFLs_series["rep2_MEFLs_"+time+"h"].tolist() +
                                df_condition_MEFLs_series["rep3_MEFLs_"+time+"h"].tolist()
        )
        time_list = [time]*len(condition_time_MEFLs)
        time_list = [float(i) for i in time_list]
        grey_ = [(i/255) for i in [150, 150, 150]]
        # plot MEFLs for each cell vs time point
        ax.plot(
            time_list, condition_time_MEFLs, linestyle="none", 
            marker="o", markersize=1.5, color=grey_
        )

    ax.set_xlabel("Time (h)")
    ax.set_xticks([0, 20, 40])
    ax.set_yticks([0, 0.5E9, 1E9]) #scaling set for plotting all cells; for plotting threshold cells only, use [0, 0.5E8, 1E8, 1.5E8]
    ax.set_ylabel("MEFLs")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0, top=1.04E9)
    ax.set_box_aspect(1)
    # plt.show()
    plt.savefig(path_figure+"time_series_"+condition+"_combined.svg")


def plot_reference(
        path_figure: str,
        time_points: list,
        reference_df: pd.DataFrame
):
    """Plots individual cell MEFLs at each time point for
    the reference condition (all 3 replicates on the same 
    plot/same color). Separate from plotting other conditions
    to enable different scaling."""

    fig, ax = plt.subplots(1, 1, figsize=(1.75, 1.75))#(3, 2.75)
    # loop through time points
    for j, time in enumerate(time_points):
        # combine data from all 3 replicates for time point
        reference_time_MEFLs = (reference_df["rep1_MEFLs_"+time+"h"].tolist() +
                                reference_df["rep2_MEFLs_"+time+"h"].tolist() + 
                                reference_df["rep3_MEFLs_"+time+"h"].tolist()
    )
        reference_time_list = [time]*len(reference_time_MEFLs)
        reference_time_list = [float(i) for i in reference_time_list]
        grey_ = [(i/255) for i in [150, 150, 150]]
        # plot MEFLs for each cell vs time point
        ax.plot(
            reference_time_list, reference_time_MEFLs, linestyle="none",
            marker="o", markersize=1.5, color=grey_
        )

        ax.set_xlabel("Time (h)")
        ax.set_xticks([0, 20, 40])
        ax.set_yticks([0, 0.5E9, 1E9]) #scaling set for plotting all cells; for plotting threshold cells only, use [0, 0.5E8, 1E8, 1.5E8]
        ax.set_ylabel("MEFLs")
        ax.set_xlim(left=0, right=48)
        ax.set_ylim(bottom=0, top=1.04E9)
        ax.set_box_aspect(1)
    # plt.show()
    plt.savefig(path_figure+"time_series_reference.svg")


def plot_experiment(
        path_figure: str,
        conditions_plot: list,
        time_points: list,
        experiment_dfs_MEFLs_dict: dict,
):
    """Plots the individual cell data for each
    condition."""

    # loop through conditions to plot
    for condition in conditions_plot:
        # get df for condition from the df dict
        condition_df = experiment_dfs_MEFLs_dict[condition]
        # plot condition
        plot_condition_combined(
            path_figure,
            condition,
            time_points,
            condition_df
        )


def fit_spline_condition(
        df_condition_MEFLs_series: pd.DataFrame
):
    """Fits a spline polynomial to the experimental
    data for the given condition and calculates the 
    derivative of the spline."""

    # reformat the condition df for spline fitting
    df_condition_MEFLs_series_T = df_condition_MEFLs_series.transpose().copy()
    df_condition_MEFLs_series_T.index = (
        df_condition_MEFLs_series_T.index.str.strip('h')
    )
    df_condition_MEFLs_series_T["Time (h)"] = df_condition_MEFLs_series_T.index
    df_condition_MEFLs_series_T_plot = pd.melt(
        frame=df_condition_MEFLs_series_T,
        id_vars="Time (h)",
        var_name="column_name",
        value_name="MEFLs"
    )
    df_condition_MEFLs_series_T_plot["Time (h)"] = (
        df_condition_MEFLs_series_T_plot["Time (h)"].astype(float)
    )
    df_condition_MEFLs_series_T_sorted = (
        df_condition_MEFLs_series_T_plot.sort_values("Time (h)")
    )
    # only keep the cells with MEFLs >= 0 (mostly needed when fitting full dataset)
    df_condition_MEFLs_series_T_no_neg = (
        df_condition_MEFLs_series_T_sorted[
            df_condition_MEFLs_series_T_sorted["MEFLs"] >=0
        ]
    )
    # drop nans (replicates have different number of cells, so some
    # rows are filled with nans in original df)
    df_condition_sorted_no_nan = df_condition_MEFLs_series_T_no_neg.dropna()
    # fit spline with specified degree to data
    spline_condition = np.polyfit(
        df_condition_sorted_no_nan["Time (h)"].tolist(),
        df_condition_sorted_no_nan["MEFLs"].tolist(), deg=3
    )
    poly_spline = np.poly1d(spline_condition)
    # calculate the derivative of the spline fit
    spline_derivative = poly_spline.deriv()
    times = np.arange(14, 47, 0.1)
    # evaluate the spline and derivative at each time point (use
    # smaller dt value than experiment for a smooth curve)
    spline_eval = poly_spline(times)
    derivative_eval = spline_derivative(times)
    
    return spline_eval, derivative_eval, df_condition_MEFLs_series_T_plot


def fit_spline_experiment(
        path_save: str,
        conditions: list,
        experiment_dfs_MEFLs_combined_dict:dict,
):
    """Fits a spline polynomial to the experimental
    data for all conditions, calculates the 
    derivative of the spline, and saves calculations."""

    spline_eval_list = []
    derivative_eval_list = []
    df_condition_MEFLs_series_list = []
    # loop through conditions
    for condition in conditions:
        condition_df = experiment_dfs_MEFLs_combined_dict[condition]
        # fit spline and calculate derivative
        spline_eval, derivative_eval, df_condition_MEFLs_series_T_plot = (
            fit_spline_condition(condition_df, condition)
        )
        spline_eval_list.append(spline_eval)
        derivative_eval_list.append(derivative_eval)
        df_condition_MEFLs_series_list.append(df_condition_MEFLs_series_T_plot)

    # save calculations
    fname_spline = "spline_eval_list.pkl"
    with open(path_save+fname_spline, "wb") as fid:
        pickle.dump(spline_eval_list, fid)
    fname_deriv = "derivative_eval_list.pkl"
    with open(path_save+fname_deriv, "wb") as fid:
        pickle.dump(derivative_eval_list, fid)

    return spline_eval_list, derivative_eval_list, df_condition_MEFLs_series_list


def plot_spline_and_deriv_condition(
        path_figure: str,
        spline_eval: np.ndarray,
        derivative_eval: np.ndarray, 
        spline_eval_all_cells: np.ndarray,
        deriv_eval_all_cells: np.ndarray,
        condition: str
):
    """Plots the spline and spline derivative
    for a given condition."""
    
    # define times used in spline and derivative
    # calculations
    times = np.arange(14, 47, 0.1)

    # plot the spline fit
    fig1, ax1 = plt.subplots(1, 1, figsize=(1.75, 1.75))
    # this code can be uncommented to set specific
    # y-axis limits and ticks for the reference
    # vs. the other conditions
    # if condition == "reference":
    #     ylim = 2.2E7
    #     yticks = [0, 1E7, 2E7]
    # else:
    #     ylim = 6E7
    #     yticks = [0, 2E7, 4E7, 6E7]
    ax1.plot(times, spline_eval, color="k", label="subpop")
    ax1.plot(times, spline_eval_all_cells, color="grey", label="fullpop")
    ax1.set_xticks(np.arange(0, 47, 10))
    ax1.set_xlabel("time (h)")
    ax1.set_ylabel("polynomial fit")  
    ax1.set_ylim(bottom=0)#,top=ylim)  
    # ax1.set_yticks(yticks)
    ax1.set_xticks([0, 20, 40])
    ax1.legend()
    ax1.set_box_aspect(1)
    plt.savefig(path_figure+"exp_spline_"+condition+"_paper.svg")

    # plot the spline derivative fit
    fig2, ax2 = plt.subplots(1, 1, figsize=(1.75, 1.75))
    # this code can be uncommented to set specific
    # y-axis limits and ticks for the reference
    # if condition == "reference":
    #     ylim_ = [0, 2.2E6]
    #     yticks_ = [0, 1E6, 2E6]
    # else:
    #     ylim_ = [-1E6, 4E6]
    #     yticks_ = [-1E6, 0, 1E6, 2E6, 3E6, 4E6]
    ax2.plot(times, derivative_eval, color="k", label="subpop")
    ax2.plot(times, deriv_eval_all_cells, color="grey", label="fullpop")
    ax2.set_xticks(np.arange(0, 47, 10))
    ax2.set_xlabel("time (h)")
    ax2.set_ylabel("fit deriv")    
    # ax2.set_ylim(ylim_)  
    # ax2.set_yticks(yticks_)
    ax2.set_xticks([0, 20, 40])
    ax2.set_yticks([-2E6, 0, 2E6, 4E6])
    ax2.set_title(condition + " deriv")
    ax2.legend()
    ax2.set_box_aspect(1)
    plt.savefig(path_figure+"exp_deriv_"+condition+"_paper.svg")


def plot_spline_and_deriv_experiment(
        path_figure: str,
        conditions: list,
        spline_eval_list: list,
        derivative_eval_list: list
):
    
    #load in data for all cells manually with path
    path_all_cells = ("/Users/kdreyer/Library/CloudStorage/"
    "OneDrive-NorthwesternUniversity/KatieD_LL/GCAD_Collab/"
    "Experimental_data_&_planning/240806_Pulse_flow_cytometry/"
    "results/all_cells_paper_deg3/"
    )
    with open(path_all_cells+"spline_eval_list.pkl", "rb") as fid:
        spline_eval_all_cells = pickle.load(fid)

    with open(path_all_cells+"derivative_eval_list.pkl", "rb") as fid:
        deriv_eval_all_cells = pickle.load(fid)
    
    # loop through conditions
    for i, condition in enumerate(conditions):
        # plot spline and derivative for condition
        plot_spline_and_deriv_condition(path_figure, spline_eval_list[i],
                                        derivative_eval_list[i], 
                                        spline_eval_all_cells[i], 
                                        deriv_eval_all_cells[i], condition)


def run_flow_cytometry_calculations(
        path_data: str,
        path_save: str,
        conditions: list,
        time_points: list,
        background_subtractions: list,
        background_subtractions_stderr: list,
        MEFL_conversions: list,
        MEFL_conversions_stderr: list,
        percentile: int=None,
        save_data: bool=True,
):
    """Runs the above functions to do the time
    series analysis for the flow cytometry experiment."""

    # load the data
    experiment_dfs_dict = load_files(path_data, conditions, time_points)

    # if a percentile threshold is set, filter the data
    # before doing the MEFLs calculations
    if percentile:
        experiment_dfs_filtered_dict = filter_tfx_experiment(
            experiment_dfs_dict, 
            conditions, 
            time_points,
            percentile)

        experiment_dfs_calculations_dict = calculate_MEFLs_experiment(
            experiment_dfs_filtered_dict,
            conditions,
            time_points, 
            background_subtractions, 
            background_subtractions_stderr,
            MEFL_conversions,
            MEFL_conversions_stderr
    )
        
    else:
        experiment_dfs_calculations_dict = calculate_MEFLs_experiment(
            experiment_dfs_dict,
            conditions,
            time_points, 
            background_subtractions, 
            background_subtractions_stderr,
            MEFL_conversions,
            MEFL_conversions_stderr,
    )
    
    # format the data for subsequent calculations and plotting
    # and save the data if save_data=True
    (experiment_dfs_MEFLs_dict, 
     experiment_dfs_MEFLs_combined_dict) = format_and_save_data(
        path_save,
        experiment_dfs_calculations_dict,
        conditions,
        time_points,
        save_data
    )

    # plot the individual cell data for each time point for
    # the reference and the other conditions
    conditions_plot = [i for i in conditions if i != "reference"]
    plot_experiment(
        path_save,
        conditions_plot,
        time_points,
        experiment_dfs_MEFLs_dict,
    )

    reference_df = experiment_dfs_MEFLs_dict["reference"]
    plot_reference(
        path_save,
        time_points,
        reference_df,
    )

    # fit splines to each condition and calculate
    # spline derivative
    (spline_eval_list, derivative_eval_list,
    df_condition_MEFLs_series_list) = fit_spline_experiment(
        path_save,
        conditions,
        experiment_dfs_MEFLs_combined_dict
    )
    # plot the spline and derivative for each condition
    plot_spline_and_deriv_experiment(
        path_save,
        conditions,
        spline_eval_list,
        derivative_eval_list,
        df_condition_MEFLs_series_list
    )