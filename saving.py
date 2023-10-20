### this code is based off of saving.py from the GAMES workflow and 
### is mostly the same: https://github.com/leonardlab/GAMES

import os
import json
from datetime import date

def make_main_directory(settings: dict) -> str:
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
