import os
import pandas as pd


def get_git_directory(config):
    if config.get("git_dir") is not None:
        git_dir = config["git_dir"]
    else:
        git_dir = os.getenv("GIT_DIR", "/users/jsun1/SpaceHack2023")

    if not git_dir.endswith("/"):
        git_dir += "/"
    return git_dir


def get_sample_dirs(data_dir):
    return [f.path for f in os.scandir(data_dir) if f.is_dir()]


def check_files_in_folder(folder_path, file_list):
    # Get a list of all files in the folder
    files_in_folder = os.listdir(folder_path)
    # Check each file in the file_list
    for file in file_list:
        if file not in files_in_folder:
            return False
    return True


def get_ncluster(file_path, sample, default_value=7):
    if not os.path.exists(file_path):
        return default_value
    try:
        df = pd.read_csv(file_path, sep="\t", index_col=0)
        df_filtered = df[df["directory"] == sample]
        return int(df_filtered["n_clusters"].mean())
    except:
        return default_value
