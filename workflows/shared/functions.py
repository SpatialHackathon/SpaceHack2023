import os
import pandas as pd


def get_sample_dirs(data_dir):
    return [ f.path for f in os.scandir(data_dir) if f.is_dir() ]


def check_files_in_folder(folder_path, file_list):
    # Get a list of all files in the folder
    files_in_folder = os.listdir(folder_path)
    # Check each file in the file_list
    for file in file_list:
        if file not in files_in_folder:
            return False
    return True


def generate_metrics_results(data_dir, metrics_name, methods, file_ext, configfiles):
    result_files = []
    for sample_dir in get_sample_dirs(data_dir):
        for method in methods:
            method_dir = sample_dir + "/" + method
            if check_files_in_folder(method_dir, ["domains.tsv"]):
                if configfiles:
                    for config_file_name in configfiles.keys():
                        result_files.append(method_dir + "/" + metrics_name + "/" + config_file_name + "/results." + file_ext)
                else:
                    result_files.append(method_dir + "/" + metrics_name + "/results." + file_ext)
            else:
                for config_dir in os.listdir(method_dir):
                    if check_files_in_folder(method_dir + "/" + config_dir, ["domains.tsv"]):
                        if "config_files" in config[methods].keys():
                            for config_file_name in config[methods]["config_files"].keys():
                                result_files.append(method_dir + "/" + config_dir + "/" + metrics_name + "/" + config_file_name + "/results." + file_ext)
                        else: result_files.append(method_dir + "/" + config_dir + "/" + metrics_name + "/results." + file_ext)
    return(result_files)




#configfiles = {"name1":"file1"}
#methods = ["spaGCN"] 
#data_dir = "/home/ubuntu/tmp_data/libd_dlpfc"
#metrics_name = "test_metrics"
#file_ext = "json"
#generate_metrics_results(data_dir, metrics_name, methods, file_ext,configfiles)
#generate_metrics_results(data_dir, metrics_name, methods, file_ext, None)

def get_ncluster(file_path, sample, default_value=7):
    if not os.path.exists(file_path):
        return default_value
    try:
        df = pd.read_csv(file_path, sep='\t', index_col=0)
        df_filtered = df[df["directory"] == sample]
        return int(df_filtered["n_clusters"].mean())
    except:
        return default_value

#sample = "Br5595_151669"
#file_path = "/home/ubuntu/tmp_data/xenium_mouse_brain_SergioSalas/samples.tsv"
#get_ncluster(file_path, sample)
