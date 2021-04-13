import os
from numpy import *
import pandas as pd


def create_dir(root_path, dir_name):
    results_dir = os.path.join(root_path, dir_name)
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    return results_dir


def create_generalinfo_path(directory):
    generalinfo_path = os.path.join(directory, 'allsatellites')
    if not os.path.exists(generalinfo_path):
        os.makedirs(generalinfo_path)
    return generalinfo_path


def get_obs_files_simple(directory):
    obs_ext = ".obs"
    obs_files = []
    files = next(os.walk(directory))[2]
    for file in files:
        if os.path.splitext(file)[1] == obs_ext:
            obs_files.append(os.path.abspath(os.path.join(directory, file)))
    return obs_files


def get_pos_files_simple(directory):
    obs_ext = ".pos"
    pos_files = []
    files = next(os.walk(directory))[2]
    for file in files:
        if os.path.splitext(file)[1] == obs_ext:
            pos_files.append(os.path.abspath(os.path.join(directory, file)))
    return pos_files


def posparser(file):
    a = []
    with open(file) as in_file:
        lineList = [line.rstrip('\n').split(" ") for line in in_file]
        for line in lineList[15:]:
            # print(len(line))
            # print(line)
            # print(line[3], line[6], line[8])
            # out_file.write("{},{},{}\n".format(line[3], line[6], line[8]))
            a.append([line[3], line[6], line[8]])
    return array(a)


def find_corresponding_dirs_in_different_roots(root_A, root_B, day=False):
    dir_pairs = []
    A_subfolders_with_paths = [f.path for f in os.scandir(root_A) if f.is_dir()]
    B_subfolders_with_paths = [f.path for f in os.scandir(root_B) if f.is_dir()]
    if not day:
        for A_dir in A_subfolders_with_paths:
            A_dirname = os.path.split(A_dir)[-1]
            sign = 0
            for B_dir in B_subfolders_with_paths:
                B_dirname = os.path.split(B_dir)[-1]
                if sign == 1:
                    sign = 0
                    break
                cond_for_days = (len(B_dirname) == 12 and len(A_dirname) == 12 and A_dirname[-6:] == B_dirname[-6:])
                if (A_dirname == B_dirname) or cond_for_days:
                    dir_pairs.append([A_dir, B_dir])
                    sign = 1

        return dir_pairs


def too_bad_positions_std(path, filename, std_limit=30.0):
    path = os.path.join(path, "allsatellites")
    file = os.path.join(path, filename)
    if os.path.isfile(file):
        user_ = pd.read_csv(file, skiprows=1).values  # .transpose()
        std_pos = mean(std(user_, axis=0))
        # print("                             STD: ", std_pos)
        if std_pos > std_limit:
            return True
        return False
    return False


def is_reliable(A_day, B_day, needed_file):
    stdUA = too_bad_positions_std(A_day, needed_file)
    stdUB = too_bad_positions_std(B_day, needed_file)
    if stdUA or stdUB:
        return False
    return True


def are_reliable(A_day, B_day, needed_file):
    if is_all_data(A_day, needed_file, add_allsatellites=True) and is_all_data(B_day, needed_file, add_allsatellites=True):
        stdUA = too_bad_positions_std(A_day, needed_file[0])
        stdUB = too_bad_positions_std(B_day, needed_file[0])
        if not stdUA and not stdUB:
            return True
    return False


def is_all_data(path, needed_files, add_allsatellites=False):
    if add_allsatellites:
        path = os.path.join(path, "allsatellites")
    try:

        list_files_with_paths = [f.path for f in os.scandir(path) if f.is_file()]

        count = 0
        for needed in needed_files:
            for file in list_files_with_paths:
                if os.path.basename(file) == needed:
                    count += 1
                    if count == len(needed_files):
                        return True
    except:
        pass
    return False
