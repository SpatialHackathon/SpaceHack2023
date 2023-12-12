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

