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

def load_files(base_path, conditions, times):

    experiment_dfs_dict = {}
    for condition in conditions:
        experiment_dfs_dict[condition] = {}
        for time in times:
            data_path = base_path + condition + "_time" + time + ".csv"
            df = pd.read_csv(data_path)
            experiment_dfs_dict[condition][time] = df

    return experiment_dfs_dict


def filter_tfx_single(df_condition, percentile):

    df_to_filter = df_condition.copy()

    if percentile != "average":
        percentile_fraction = percentile/100
        replicates = ["rep1", "rep2", "rep3"]
        df_filtered_list = []
        for rep in replicates:
            df_rep = df_to_filter.filter(regex=rep)
            max_pac_blue = df_rep[rep+"_Pac_Blue"].max()
            df_rep_filtered = df_rep[
                df_rep[rep+"_Pac_Blue"] >=
                max_pac_blue*percentile_fraction
            ].reset_index(drop=True)
            df_filtered_list.append(df_rep_filtered)
    else:
        replicates = ["rep1", "rep2", "rep3"]
        df_filtered_list = []
        for rep in replicates:
            df_rep = df_to_filter.filter(regex=rep)
            mean_pac_blue = df_rep[rep+"_Pac_Blue"].mean()
            df_rep_filtered = df_rep[
                df_rep[rep+"_Pac_Blue"] >=
                mean_pac_blue*1.8
            ].reset_index(drop=True)
            df_filtered_list.append(df_rep_filtered)

    df_all_filtered = pd.concat(df_filtered_list, axis=1)

    return df_all_filtered


def filter_tfx_condition(
        df_condition_dict,
        time_points,
        percentile
):
    
    df_condition_dict_filtered = {}
    for i, time in enumerate(time_points):
        df_time_point = df_condition_dict[time].copy()
        df_time_point = filter_tfx_single(
            df_time_point,
            percentile
        )
        df_condition_dict_filtered[time] = (
            df_time_point
        )

    return df_condition_dict_filtered


def filter_tfx_experiment(
        experiment_dfs_dict,
        conditions,
        times,
        percentile
):
    
    experiment_dfs_filtered_dict = {}
    for condition in conditions:
        df_condition_dict = experiment_dfs_dict[condition]
        df_condition_dict_filtered = filter_tfx_condition(
            df_condition_dict,
            times, 
            percentile
        )
        experiment_dfs_filtered_dict[condition] = (
            df_condition_dict_filtered
        )

    return experiment_dfs_filtered_dict


def background_subtract_single(
        df_condition, 
        background_subtraction,
        background_subtraction_stderr
):
    
    df_condition_background_subtracted = df_condition.copy()

    replicates = ["rep1", "rep2", "rep3"] 
    for rep in replicates:
        df_condition_background_subtracted[rep + "_FITC_background_subtracted"] = (
            df_condition_background_subtracted[rep+"_FITC"] - background_subtraction
        )
        df_condition_background_subtracted[rep + "_FITC_background_subtracted_stderr"] = (
        background_subtraction_stderr)
    
    return df_condition_background_subtracted


def background_subtract_condition(
        df_condition_dict,
        time_points,
        background_subtractions,
        background_subtractions_stderr
):
    
    df_condition_dict_background_subtracted = {}
    for i, time in enumerate(time_points):
        df_time_point = df_condition_dict[time].copy()
        df_time_point = background_subtract_single(
            df_time_point,
            background_subtractions[i],
            background_subtractions_stderr[i]
        )
        df_condition_dict_background_subtracted[time] = (
          df_time_point  
        )

    return df_condition_dict_background_subtracted


def convert_to_MEFLs_single(
        df_condition,
        MEFL_conversion,
        MEFL_conversion_stderr
):
    
    df_condition_MEFLs = df_condition.copy()

    replicates = ["rep1", "rep2", "rep3"] 
    for rep in replicates:
        df_condition_MEFLs[rep + "_MEFLs"] = (
            df_condition_MEFLs[rep + "_FITC_background_subtracted"]*MEFL_conversion
        )
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
    
    df_condition_dict_MEFLs = {}
    for i, time in enumerate(time_points):
        df_time_point = df_condition_dict[time].copy()
        df_time_point = convert_to_MEFLs_single(
            df_time_point,
            MEFL_conversions[i],
            MEFL_conversions_stderr[i]
        )
        df_condition_dict_MEFLs[time] = (
            df_time_point
        )

    return df_condition_dict_MEFLs


def calculate_MEFLs_experiment(
        experiment_dfs_dict,
        conditions,
        times,
        background_subtractions,
        background_subtractions_stderr,
        MEFL_conversions,
        MEFL_conversions_stderr
):
    
    experiment_dfs_calculations_dict = {}
    for condition in conditions:
        df_condition_dict = experiment_dfs_dict[condition]
        df_condition_dict_background_subtracted = background_subtract_condition(
            df_condition_dict,
            times,
            background_subtractions,
            background_subtractions_stderr
        )
        df_condition_dict_MEFLs = convert_to_MEFLs_condition(
            df_condition_dict_background_subtracted,
            times,
            MEFL_conversions,
            MEFL_conversions_stderr
        )
        experiment_dfs_calculations_dict[condition] = df_condition_dict_MEFLs

    return experiment_dfs_calculations_dict


def format_and_save_data(
        path, 
        experiment_dfs_calculations_dict,
        conditions,
        time_points,
        save_data
):
    replicates = ["rep1", "rep2", "rep3"]
    experiment_dfs_MEFLs_dict = {}
    experiment_dfs_MEFLs_combined_dict = {}
    for condition in conditions:
        df_condition_MEFLs_series_list = []
        df_condition_MEFLs_combined_series_list = []
        for time in time_points:
            df_single = experiment_dfs_calculations_dict[condition][time]
            if save_data:
                df_single.to_csv(path + condition + "_calculations_" + time + "h.csv")
            df_single_MEFLs = df_single.filter(regex="MEFLs")
            df_single_MEFLs = df_single_MEFLs.add_suffix("_" + time + "h")
            df_condition_MEFLs_series_list.append(df_single_MEFLs)

            df_single_MEFLs_reps_list = [df_single_MEFLs[rep + "_MEFLs_" + time + "h"]
                                    for rep in replicates
            ]
            df_single_MEFLs_reps_list = [df.dropna() for df in df_single_MEFLs_reps_list]
            df_single_MEFLs_combined = pd.concat(df_single_MEFLs_reps_list, ignore_index=True, axis=0)
            df_condition_MEFLs_combined_series_list.append(df_single_MEFLs_combined)

        df_condition_MEFLs_series = pd.concat(
            df_condition_MEFLs_series_list, axis=1
        )
        if save_data:
            df_condition_MEFLs_series.to_csv(path + condition + "_MEFLs.csv")
        experiment_dfs_MEFLs_dict[condition] = df_condition_MEFLs_series

        df_condition_MEFLs_combined_series = pd.concat(
            df_condition_MEFLs_combined_series_list, axis=1,
            keys=[t+"h" for t in time_points]
        )
        experiment_dfs_MEFLs_combined_dict[condition] = df_condition_MEFLs_combined_series

    return experiment_dfs_MEFLs_dict, experiment_dfs_MEFLs_combined_dict


def plot_condition(
        path_figure,
        condition,
        time_points,
        df_condition_MEFLs_series,
):
    fig, axs = plt.subplots(1, 3, figsize=(7, 2.75))#, sharey=True)
    axs = axs.ravel()
    labels_condition = [condition] + [""]*(len(time_points)-1)
    for i, rep in enumerate(["rep1", "rep2", "rep3"]):
        for j, time in enumerate(time_points):
            condition_time_MEFLs = df_condition_MEFLs_series[rep+"_MEFLs_"+time+"h"].tolist()
            time_list = [time]*len(condition_time_MEFLs)
            time_list = [float(i) for i in time_list]
            axs[i].plot(
                time_list, condition_time_MEFLs, linestyle="none", 
                marker="o", markersize=4, color=orange_, label=labels_condition[j] #, fillstyle="none")
            )

        axs[i].set_xlabel("Time (h)")
        axs[i].set_xticks([float(i) for i in time_points])
        axs[i].set_ylabel("MEFLs")
        axs[i].set_title(rep+" time series")
        axs[i].set_box_aspect(1)
    axs[0].legend()
    # plt.show()
    plt.savefig(path_figure+"time_series_"+condition+".svg")


def plot_condition_combined(
        path_figure,
        condition,
        time_points,
        df_condition_MEFLs_series,
):
    fig, ax = plt.subplots(1, 1, figsize=(1.75, 1.75))# (3, 2.75)
    for time in time_points:
        condition_time_MEFLs = (df_condition_MEFLs_series["rep1_MEFLs_"+time+"h"].tolist() + 
                                df_condition_MEFLs_series["rep2_MEFLs_"+time+"h"].tolist() +
                                df_condition_MEFLs_series["rep3_MEFLs_"+time+"h"].tolist()
        )
        time_list = [time]*len(condition_time_MEFLs)
        time_list = [float(i) for i in time_list]
        grey_ = [(i/255) for i in [150, 150, 150]]
        ax.plot(
            time_list, condition_time_MEFLs, linestyle="none", 
            marker="o", markersize=1.5, color=grey_
        )

    ax.set_xlabel("Time (h)")
    ax.set_xticks([0, 20, 40])
    ax.set_yticks([0, 0.5E9, 1E9]) #[0, 0.5E8, 1E8, 1.5E8]
    ax.set_ylabel("MEFLs")
    ax.set_xlim(left=0)
    ax.set_ylim(bottom=0, top=1.04E9)
    ax.set_box_aspect(1)
    # plt.show()
    plt.savefig(path_figure+"time_series_"+condition+"_combined.svg")


def plot_condition_combined_max_values(
        path_figure,
        condition,
        time_points,
        df_condition_MEFLs_series,
        threshold_value
):
    if condition == "reference":
        plot_color = bluish_green
    else:
        plot_color = orange_
    fig, ax = plt.subplots(1, 1, figsize=(3, 2.75))#, sharey=True)
    for time in time_points:
        condition_time_MEFLs = (df_condition_MEFLs_series["rep1_MEFLs_"+time+"h"].tolist() + 
                                df_condition_MEFLs_series["rep2_MEFLs_"+time+"h"].tolist() +
                                df_condition_MEFLs_series["rep3_MEFLs_"+time+"h"].tolist()
        )
        condition_time_MEFLs_max_vals = [j for j in condition_time_MEFLs if j >= threshold_value]
        print(condition, " ", time, " fraction of cells above threshold: ", round(len(condition_time_MEFLs_max_vals)/len(condition_time_MEFLs), 5))
        time_list = [time]*len(condition_time_MEFLs_max_vals)
        time_list = [float(i) for i in time_list]

        ax.plot(
            time_list, condition_time_MEFLs_max_vals, linestyle="none", 
            marker="o", markersize=4, color=plot_color, fillstyle="none"
        )

    ax.set_xlabel("Time (h)")
    ax.set_xticks([float(i) for i in time_points])
    ax.set_ylabel("MEFLs >= threshold value")
    ax.set_box_aspect(1)
    # plt.show()
    plt.savefig(path_figure+"time_series_"+condition+"_zoomed.svg")


def plot_reference(
        path_figure,
        time_points,
        reference_df
):
    fig, ax = plt.subplots(1, 1, figsize=(1.75, 1.75))#(3, 2.75)
    labels_reference = ["reference"] + [""]*(len(time_points)-1)
    for j, time in enumerate(time_points):
        reference_time_MEFLs = (reference_df["rep1_MEFLs_"+time+"h"].tolist() +
                                reference_df["rep2_MEFLs_"+time+"h"].tolist() + 
                                reference_df["rep3_MEFLs_"+time+"h"].tolist()
    )
        reference_time_list = [time]*len(reference_time_MEFLs)
        reference_time_list = [float(i) for i in reference_time_list]
        grey_ = [(i/255) for i in [150, 150, 150]]
        ax.plot(
            reference_time_list, reference_time_MEFLs, linestyle="none",
            marker="o", markersize=1.5, color=grey_
        )

        ax.set_xlabel("Time (h)")
        ax.set_xticks([0, 20, 40])
        ax.set_yticks([0, 0.5E9, 1E9]) #[0, 0.5E8, 1E8, 1.5E8]
        ax.set_ylabel("MEFLs")
        ax.set_xlim(left=0, right=48)
        ax.set_ylim(bottom=0, top=1.04E9)
        ax.set_box_aspect(1)
    # plt.show()
    plt.savefig(path_figure+"time_series_reference.svg")


def plot_experiment(
        path_figure,
        conditions_plot,
        time_points,
        experiment_dfs_MEFLs_dict,
):
    
    for condition in conditions_plot:
        condition_df = experiment_dfs_MEFLs_dict[condition]
        plot_condition_combined(
            path_figure,
            condition,
            time_points,
            condition_df
        )


def fit_spline_condition(df_condition_MEFLs_series, condition):
    df_condition_MEFLs_series_T = df_condition_MEFLs_series.transpose().copy()
    df_condition_MEFLs_series_T.index = df_condition_MEFLs_series_T.index.str.strip('h')
    df_condition_MEFLs_series_T["Time (h)"] = df_condition_MEFLs_series_T.index
    df_condition_MEFLs_series_T_plot = pd.melt(
        frame=df_condition_MEFLs_series_T,
        id_vars="Time (h)",
        var_name="column_name",
        value_name="MEFLs")
    df_condition_MEFLs_series_T_plot["Time (h)"] = df_condition_MEFLs_series_T_plot["Time (h)"].astype(float)
    df_condition_MEFLs_series_T_sorted = df_condition_MEFLs_series_T_plot.sort_values("Time (h)")
    df_condition_MEFLs_series_T_no_neg = df_condition_MEFLs_series_T_sorted[df_condition_MEFLs_series_T_sorted["MEFLs"] >=0]
    df_condition_sorted_no_nan = df_condition_MEFLs_series_T_no_neg.dropna()
    # print(df_condition_sorted_no_nan)
    spline_condition = np.polyfit(df_condition_sorted_no_nan["Time (h)"].tolist(), df_condition_sorted_no_nan["MEFLs"].tolist(), deg=3)
    poly_spline = np.poly1d(spline_condition)
    spline_derivative = poly_spline.deriv()
    # print(spline_condition)
    times = np.arange(14, 47, 0.1)
    # times_exp = [14, 18, 22, 26, 38, 42, 46]
    spline_eval = poly_spline(times)
    # spline_eval_exp = poly_spline(times_exp)
    derivative_eval = spline_derivative(times)
    # derivative_eval_exp = spline_derivative(times_exp)
    
    return spline_eval, derivative_eval, df_condition_MEFLs_series_T_plot


def fit_spline_experiment(
        path_save,
        conditions,
        experiment_dfs_MEFLs_combined_dict,
):

    spline_eval_list = []
    derivative_eval_list = []
    df_condition_MEFLs_series_list = []
    for condition in conditions:
        condition_df = experiment_dfs_MEFLs_combined_dict[condition]
        spline_eval, derivative_eval, df_condition_MEFLs_series_T_plot = fit_spline_condition(condition_df, condition)
        spline_eval_list.append(spline_eval)
        derivative_eval_list.append(derivative_eval)
        df_condition_MEFLs_series_list.append(df_condition_MEFLs_series_T_plot)

    fname_spline = "spline_eval_list.pkl"
    with open(path_save+fname_spline, "wb") as fid:
        pickle.dump(spline_eval_list, fid)

    fname_deriv = "derivative_eval_list.pkl"
    with open(path_save+fname_deriv, "wb") as fid:
        pickle.dump(derivative_eval_list, fid)

    return spline_eval_list, derivative_eval_list, df_condition_MEFLs_series_list

def plot_spline_and_deriv_condition(
        path_figure, spline_eval, derivative_eval, 
        spline_eval_all_cells, deriv_eval_all_cells,
        condition):
    
    times = np.arange(14, 47, 0.1)
    fig1, ax1 = plt.subplots(1, 1, figsize=(1.75, 1.75))

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

    fig2, ax2 = plt.subplots(1, 1, figsize=(1.75, 1.75))

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
    # plt.show()
    plt.savefig(path_figure+"exp_deriv_"+condition+"_paper.svg")


def plot_spline_and_deriv_experiment(path_figure, conditions,
        spline_eval_list, derivative_eval_list, df_condition_MEFLs_series_list):
    
    #load in manually with path
    path_all_cells = "/Users/kdreyer/Library/CloudStorage/OneDrive-NorthwesternUniversity/KatieD_LL/GCAD_Collab/Experimental_data_&_planning/240806_Pulse_flow_cytometry/results/all_cells_paper_deg3/"
    with open(path_all_cells+"spline_eval_list.pkl", "rb") as fid:
        spline_eval_all_cells = pickle.load(fid)

    with open(path_all_cells+"derivative_eval_list.pkl", "rb") as fid:
        deriv_eval_all_cells = pickle.load(fid)
    
    for i, condition in enumerate(conditions):
        plot_spline_and_deriv_condition(path_figure, spline_eval_list[i],
                                        derivative_eval_list[i], 
                                        spline_eval_all_cells[i], 
                                        deriv_eval_all_cells[i], condition)


def run_flow_cytometry_calculations(
        path_data,
        path_save,
        conditions,
        time_points,
        background_subtractions,
        background_subtractions_stderr,
        MEFL_conversions,
        MEFL_conversions_stderr,
        percentile=None,
        save_data=True,
        log_scale=False,
        split_plot=False
):
    
    experiment_dfs_dict = load_files(path_data, conditions, time_points)

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
    
    (experiment_dfs_MEFLs_dict, 
     experiment_dfs_MEFLs_combined_dict) = format_and_save_data(
        path_save,
        experiment_dfs_calculations_dict,
        conditions,
        time_points,
        save_data
    )


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

    (spline_eval_list, derivative_eval_list,
    df_condition_MEFLs_series_list) = fit_spline_experiment(path_save, conditions, experiment_dfs_MEFLs_combined_dict)
    plot_spline_and_deriv_experiment(path_save, conditions, spline_eval_list, derivative_eval_list, df_condition_MEFLs_series_list)