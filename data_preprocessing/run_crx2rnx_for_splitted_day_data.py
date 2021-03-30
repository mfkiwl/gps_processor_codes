import os
from datetime import datetime

from tqdm import tqdm


def get_subfolders_with_source(root_path):
    return [os.path.join(root_path, f) for f in os.listdir(root_path) if os.path.isdir(os.path.join(root_path, f))]


def run_CRX2RNX(day_dir):
    day_index = os.path.split(day_dir)[-1]
    # zipped_files = str(os.path.normpath(os.path.join(day_dir, "*.gz")))
    # unzipp = "gunzip " + zipped_files
    # os.system(unzipp)
    input_files = [os.path.join(day_dir, f) for f in os.listdir(day_dir) if
                   (os.path.isfile(os.path.join(day_dir, f)) and str(f)[-2:] == "0d")]
    # print(day_dir)
    print(len(input_files))
    for file in input_files:  # [:1]:
        file = os.path.normpath(file)
        tool = os.path.normpath(r"C:\SzabolcsKelemen\PhD\GPS\RTKLIB\RTKLIB_bin-rtklib_2.4.3\RTKLIB_bin-rtklib_2.4.3\bin\CRX2RNX.exe")
        command1 = r"{} {} -f".format(tool, file)
        # print(file)
        os.system(command1)

    # command2 = "rm *gz"
    # os.system(command2)


def run_converting(root_path):
    binexes = [os.path.join(root_path, f) for f in os.listdir(root_path) if os.path.isdir(os.path.join(root_path, f))]
    for folder in tqdm(binexes):  # [:1]:

        print(folder)
        # day_data_dir = os.path.join(folder, "hkks/1s")
        day_data_dir = folder
        run_CRX2RNX(day_data_dir)


# input_ = r"C:\SzabolcsKelemen\PhD\GPS\HKKS\downloaded"
# input_ = r"D:\GPS\raw_data\TIDV\raw_data\extracted_collected\2020"
input_ = r"D:\GPS\raw_data\IISC\extracted_collected\2020"

run_converting(input_)
