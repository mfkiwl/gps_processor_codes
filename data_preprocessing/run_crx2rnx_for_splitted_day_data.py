import os
from datetime import datetime 



def get_subfolders_with_source(root_path):
    return [os.path.join(root_path, f) for f in os.listdir(root_path) if os.path.isdir(os.path.join(root_path, f))]



def run_CRX2RNX(day_dir):
    day_index = os.path.split(day_dir)[-1]

    input_files = [os.path.join(day_dir, f) for f in os.listdir(day_dir) if (os.path.isfile(os.path.join(day_dir, f)) and  str(f)[-2:] != "gz" )]

    
    for file in input_files:
        command1 = "./CRX2RNX {} -f -d".format(file)
        
        # print(file)
        os.system(command1)
        
    command2 = "rm *gz"
    os.system(command2)


    
def run_converting(root_path):
    binexes = [os.path.join(root_path, f) for f in os.listdir(root_path) if os.path.isdir(os.path.join(root_path, f))]
    for folder in binexes:
        day_data_dir = os.path.join(folder, "hkks/1s")
        run_CRX2RNX(day_data_dir)





input_ = r"/Users/kelemensz/Documents/Research/GPS/keep_from_qsync/HKG/HKKS_rinex/januar"
# output = r"C:\SzabolcsKelemen\PhD\GPS\perth\december"

run_converting(input_)





