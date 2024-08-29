### this code is based off of saving.py from the GAMES workflow and 
### is mostly the same: https://github.com/leonardlab/GAMES

import os
import json
from datetime import date

def make_main_directory(settings: dict, custom_folder_path: str=None) -> str:
    """Makes main results folder

    Parameters
    ----------
    settings
        a dictionary of settings

    Returns
    -------
    folder_path
        path leading to main results folder
    """
    # make results folder and change directories
    if custom_folder_path:
        results_folder_path = custom_folder_path
    else:
        results_folder_path = settings["repository_path"] + "GA_results/"
    date_today = date.today()
    folder_path = results_folder_path + str(date_today) + "_" + settings["folder_name"]
    try:
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
    except FileExistsError:
        print("Directory already exists")
    os.chdir(folder_path)

    # save settings
    with open("./settings.json", "w", encoding="utf-8") as fid:
        json.dump(settings, fid)
    return folder_path

def make_results_analysis_directory(full_results_path: str,
                                    settings
) -> str:
    """Makes sub-folder for results analysis within main results
        folder

    Parameters
    ----------
    full_results_path
        an absolute path to the main results folder

    Returns
    -------
    folder_path
        path leading to the results analysis sub-folder
    """
    # make results folder and change directories
    date_today = date.today()
    folder_path = full_results_path + str(date_today) + "_results_analysis/"
    try:
        if not os.path.exists(folder_path):
            os.makedirs(folder_path)
    except FileExistsError:
        print("Directory already exists")
    os.chdir(folder_path)

    # save settings
    with open("./settings.json", "w", encoding="utf-8") as fid:
        json.dump(settings, fid)

    return folder_path
