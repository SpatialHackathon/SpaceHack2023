import os


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


def generate_metrics_results(data_dir, metrics_name, methods, file_ext, configfiles=None):
    result_files = []
    for sample_dir in get_sample_dirs(data_dir):
        for method in methods:
            method_dir = os.path.join(sample_dir, method)
            dirs_to_check = [method_dir] if check_files_in_folder(method_dir, ["domains.tsv"]) else os.listdir(method_dir)
            for dir_to_check in dirs_to_check:
                if check_files_in_folder(os.path.join(method_dir, dir_to_check), ["domains.tsv"]):
                    config_files = configfiles.keys() if configfiles else [""]
                    for config_file_name in config_files:
                        result_files.append(os.path.join(method_dir, dir_to_check, metrics_name, config_file_name, "results." + file_ext))
    return result_files

#configfiles = {"name1":"file1"}
#methods = ["spaGCN"] 
#data_dir = "/home/ubuntu/tmp_data/libd_dlpfc"
#metrics_name = "test_metrics"
#file_ext = "json"
#generate_metrics_results(data_dir, metrics_name, methods, file_ext,configfiles)
#generate_metrics_results(data_dir, metrics_name, methods, file_ext, None)
