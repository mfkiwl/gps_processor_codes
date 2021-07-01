import time

from numpy import *
import pandas as pd
from itertools import chain
import os
import math

import sys

sys.path.insert(1, '/Users/kelemensz/Documents/Research/GPS/gps_processor_codes/data_postprocessing')
from data_locations_handler import AllGPSDataLocations
from triangle_method_rework.rework_lib import v_length, parse_sats_data, Defaults, \
    get_mean_user_position, get_common, get_proportion_v3, unit, normvec, transform_matrix, ecef_to_gcs, \
    get_user_and_sat_positions, cart2sph, get_ij_on_map, get_global_stars, raw_results_to_GCS, \
    get_mean_direction_in_system_over_time, consecutively_group_results, save_matrices, plot_mollweid, create_dir, \
    plot_save_imshow_3_maps, rotate2darray, input_fileset_of_the_day, get_star_directions_in_GCS

sys.path.insert(1, '/Users/kelemensz/Documents/Research/GPS/gps_processor_codes/utility')
from frecvently_used_functions import is_all_data, is_reliable, are_reliable


def extract_common_sats(satA, satB, based_on_id=True):
    common_parts = {}
    for kA, vA in satA.items():
        vB = satB.get(kA, None)
        if vB is not None:
            common = get_common(vA, vB, is_id=based_on_id)
            common_parts[kA] = common
    return common_parts


def calc_for_one_satellite_beta(posA, posB, sat_pr_both):
    S_A = sat_pr_both[0][:3]
    S_B = sat_pr_both[1][:3]
    t1 = sat_pr_both[0][3]
    t2 = sat_pr_both[1][3]
    AS = v_length(S_A, posA)
    BS = v_length(S_B, posB)

    A_SA_i = unit(array(posA) - array(S_A))
    B_SB_j = unit(array(posB) - array(S_B))
    n = A_SA_i - B_SB_j

    mod_n = v_length(array([0, 0, 0]), n)
    return unit(n), get_proportion_v3(AS, t1, BS, t2), mod_n


# ===============================================================================================================================================================================


def process_one_epoch(posA, posB, data_from_common_sats):
    direction_value = []
    for sats in data_from_common_sats:
        direction_value.append(calc_for_one_satellite_beta(posA, posB, sats))
    return direction_value


# ===============================================================================================================================================================================

def process_all(posA, posB, common_data):
    raw_results = {}
    for k, v in common_data.items():  # k = epoch index
        raw_results[k] = process_one_epoch(posA, posB, v)
    return raw_results


# =================================================================================================
# =================================================================================================


def get_raw_results(pathA, pathB, based_on_id=False):
    posA, dataA = get_user_and_sat_positions(pathA)
    posB, dataB = get_user_and_sat_positions(pathB)
    comonAB = extract_common_sats(dataA, dataB, based_on_id=based_on_id)
    raw_results = process_all(posA, posB, comonAB)
    return raw_results


def get_raw_results_using_current_days_mean_positions(pathA, pathB, mean_positions, based_on_id=False):
    posA = mean_positions[0]
    posB = mean_positions[1]
    dataA = get_user_and_sat_positions(pathA, False)
    dataB = get_user_and_sat_positions(pathB, False)
    commonAB = extract_common_sats(dataA, dataB, based_on_id=based_on_id)
    raw_results = process_all(posA, posB, commonAB)  # dictionary, keys are the index of the epochs
    return raw_results


def process_raw_GCS_data(raw_results_GCS, resolution):
    theta_max = math.pi
    phi_max = math.pi / 2.0
    rot_theta = arange(-theta_max, theta_max, resolution)
    rot_phi = arange(-phi_max, phi_max, resolution)
    l_theta = int(len(rot_theta) / 2.0)
    l_phi = int(len(rot_phi) / 2.0)
    cmap_v = zeros((len(rot_theta), len(rot_phi)))
    cmap_count = zeros((len(rot_theta), len(rot_phi)))
    cmap_mod_n = zeros((len(rot_theta), len(rot_phi)))
    for direction, value, mod_n in raw_results_GCS:
        # i, j = get_ij_on_map_(direction, resolution)
        i, j = get_ij_on_map(direction, l_theta, l_phi, resolution)
        cmap_v[i][j] += value
        cmap_mod_n[i][j] += mod_n
        cmap_count[i][j] += 1
    cmap_count[cmap_count < 1] = 0
    return cmap_v, cmap_count, cmap_mod_n


def process_one_day_symmetrized(pathA, pathB, star_dir, resolution, mean_positions=None, root=None, fill_out=0.0,
                                id_based=False):  # fill_out=0.0 in case when the "r * (r-1)^2" values are in the data matrix
    if root is None:
        print("No day_result directory is given...process stops here!")
        sys.exit()

    Nx, Ny, Nz, S, D = get_global_stars(star_dir, os.path.dirname(pathB))
    nr_of_time_updates = len(Nx)
    star_directions_in_GCS = get_star_directions_in_GCS(Nx, Ny, Nz, S, D)
    # print(star_directions_in_GCS)
    # =================================================================================================================
    GCS_all = [array([Nx[i], Ny[i], Nz[i]]) for i in range(nr_of_time_updates)]

    # -------------------------------Symmetrize the calculated measure-------------------------------------------------
    raw_results_ECEF_AB = get_raw_results_using_current_days_mean_positions(pathA, pathB, mean_positions,
                                                                            based_on_id=id_based)
    # raw_results_ECEF_AB = get_raw_results(pathA, pathB, based_on_id=id_based)
    groupped_raw_results_ECEF_AB = consecutively_group_results(raw_results_ECEF_AB, nr_of_time_updates)
    groupped_raw_results_GCS_AB = raw_results_to_GCS(groupped_raw_results_ECEF_AB, GCS_all)
    raw_results_GCS_AB = list(chain(*groupped_raw_results_GCS_AB))

    raw_results_ECEF_BA = get_raw_results_using_current_days_mean_positions(pathB, pathA, mean_positions[::-1],
                                                                            based_on_id=id_based)
    # raw_results_ECEF_BA = get_raw_results(pathB, pathA, based_on_id=id_based)
    groupped_raw_results_ECEF_BA = consecutively_group_results(raw_results_ECEF_BA, nr_of_time_updates)
    groupped_raw_results_GCS_BA = raw_results_to_GCS(groupped_raw_results_ECEF_BA, GCS_all)
    raw_results_GCS_BA = list(chain(*groupped_raw_results_GCS_BA))
    # -----------------------------------------------------------------------------------------------------------------
    raw_results_GCS = raw_results_GCS_AB + raw_results_GCS_BA

    print("Raw results in GCS are in! (data size)", len(raw_results_GCS))
    raw_n_v_nmod = pd.DataFrame(raw_results_GCS)
    raw_n_v_nmod.to_csv(os.path.join(root, "GCS_Ndir_measure_Nmod_AB_and_BA.csv"), index=True)
    del raw_n_v_nmod

    day_data, day_count, day_cmap_n_mod = process_raw_GCS_data(raw_results_GCS, resolution)

    day_data = nan_to_num(day_data, nan=fill_out)
    day_cmap_n_mod = nan_to_num(day_cmap_n_mod, nan=fill_out)
    day_count = nan_to_num(day_count, nan=fill_out)

    save_matrices([day_data, day_count, day_cmap_n_mod], ["measure", "histogram", "n_mod"], root, resolution)
    plot_mollweid(day_count, star_directions_in_GCS, root, "histogram", str(int(degrees(resolution))), anot=True)
    plot_mollweid(divide(day_data, day_count), star_directions_in_GCS, root, "measure", str(int(degrees(resolution))),
                  anot=True)
    plot_mollweid(divide(day_cmap_n_mod, day_count), star_directions_in_GCS, root, "n_mod",
                  str(int(degrees(resolution))), anot=True)
    # plot_save_imshow(divide(day_cmap_n_mod, day_count), root, "n_mod")
    return day_data, day_count, day_cmap_n_mod


def find_corresponding_dirs_in_different_roots(root_A, root_B):
    dir_pairs = []
    A_subfolders_with_paths = [f.path for f in os.scandir(root_A) if f.is_dir()]
    B_subfolders_with_paths = [f.path for f in os.scandir(root_B) if f.is_dir()]
    for A_dir in A_subfolders_with_paths:
        A_dirname = os.path.split(A_dir)[-1]
        # sign = 0
        for B_dir in B_subfolders_with_paths:
            B_dirname = os.path.split(B_dir)[-1]
            # if sign == 1:
            #     break
            cond_for_days = (len(B_dirname) == 12 and len(A_dirname) == 12 and A_dirname[-6:] == B_dirname[-6:])
            if (A_dirname == B_dirname) or cond_for_days:
                dir_pairs.append([A_dir, B_dir])
                # sign = 1
                break
    return dir_pairs


def find_same_days_and_process(path_A, path_B, result_path, star_dir, resolution, month_names):
    d = 0
    all_hist = []
    all_value = []
    all_n_mod = []
    # print(path_B, '\n', path_A)
    if os.path.isdir(path_A) and os.path.isdir(path_B) and os.path.isdir(result_path):
        month_pairs = find_corresponding_dirs_in_different_roots(path_A, path_B)
        # print('Common months:  ', len(month_pairs))
        for A_month, B_month in month_pairs[:]:
            month_name = os.path.split(A_month)[-1]
            condition = month_name in month_names
            if condition:
                print(month_name)
                day_pairs = find_corresponding_dirs_in_different_roots(A_month, B_month)
                print("Number of days: ", len(day_pairs))
                for A_day, B_day in day_pairs[:]:
                    start = time.time()
                    date = str(os.path.split(B_day)[-1])[-8:]
                    # print('Date: ', date)
                    fileset_A_day = input_fileset_of_the_day(A_day, True)
                    fileset_B_day = input_fileset_of_the_day(B_day, True)
                    # print(fileset_A_day, fileset_B_day)
                    if fileset_A_day * fileset_B_day != 0:
                        if are_reliable(A_day, B_day, Defaults.USER_POSITIONS_FILENAME.get('user'),
                                        std_limits=[10.0, 10.0]):
                            result_month = create_dir(result_path, month_name)
                            result_day = create_dir(result_month, date)

                            if len(os.listdir(result_day)) == 0:
                                print(" Data will be processed from: ", os.path.split(A_day)[-1], "    ",
                                      os.path.split(B_day)[-1],
                                      "\n", "Index of the process: ", d, "\n")
                                A_day = os.path.join(A_day, "allsatellites")
                                B_day = os.path.join(B_day, "allsatellites")

                                posUA = get_mean_user_position(A_day)
                                posUB = get_mean_user_position(B_day)
                                based_on_id = True if fileset_A_day == 2 and fileset_B_day == 2 else False
                                value, hist, n_mod = process_one_day_symmetrized(A_day, B_day, star_dir, resolution,
                                                                                 mean_positions=[posUA, posUB],
                                                                                 root=result_day, fill_out=0.0,
                                                                                 id_based=based_on_id)
                                all_hist.append(hist)
                                all_value.append(value)
                                all_n_mod.append(n_mod)
                                d += 1
                                print('Elapsed time of the current day: ', time.time() - start, date)
                            else:
                                print("{} already processed!".format(date))
                        else:
                            print('Positions are not reliable: {}'.format(date))

                    else:
                        print("\n Data not found for: ", date, "\n")

        all_hist = sum(array(all_hist), axis=0)
        all_value = sum(array(all_value), axis=0)
        all_n_mod = sum(array(all_n_mod), axis=0)
        # save_matrices([all_value, all_hist, all_n_mod], ["measure", "histogram", "n_mod"], result_path, resolution)

        all_value_av = divide(all_value, all_hist)
        all_n_mod_av = divide(all_n_mod, all_hist)
        # print(all_value, all_hist, all_n_mod)

        all_hist = rotate2darray(nan_to_num(all_hist, nan=0.0))
        all_value_av = rotate2darray(nan_to_num(all_value_av, nan=0.0))
        all_n_mod_av = rotate2darray(nan_to_num(all_n_mod_av, nan=0.0))

        plot_save_imshow_3_maps([all_hist, all_value_av, all_n_mod_av], ["Histogram", "(|1-r|/|n|)^(1/4)", "<n_mod>"],
                                result_path, resolution="5", logplot=False)

    print("Total nr of days: ", d)


# =================================================================================================
# =================================================================================================


star_dir = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020"
resolution = radians(5.0)
# # needed_files = ["user_pos_allsatellites.csv", satellite_positions]
# needed_files = [Defaults.USER_POSITIONS_FILENAME, Defaults.SAT_POS_FILENAMES.get('no_ID'),
#                 Defaults.SAT_POS_FILENAMES.get('ID')]
#
#
# place_B = AllGPSDataLocations.user_and_satellites.get('CUTA')
# place_A = AllGPSDataLocations.user_and_satellites.get('CUTB')
# # results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_test/CUTA_CUTB"
# # results_root = create_dir(results_root, r'r_inv_r_symmetrized')
# results_root = r'/Users/kelemensz/Documents/Research/GPS/process/triangle_test/CUTA_CUTB/r_inv_r_symmetrized'
#
#
# find_same_days_and_process(place_A, place_B, result_path=results_root, star_dir=star_dir, resolution=resolution, month_names=AllGPSDataLocations.all_months)


def process_for_all(main_root, result_roots, star_dir, resolution):
    for result_root in result_roots[5:6]:
        root = os.path.join(main_root, result_root)
        print(root)
        if not os.path.isdir(root):
            os.makedirs(root)

        pair = result_root.split("/")[-2]
        # print(pair, root)
        name_A = pair.split('_')[0]
        name_B = pair.split('_')[1]
        place_A = AllGPSDataLocations.user_and_satellites.get(name_A)
        place_B = AllGPSDataLocations.user_and_satellites.get(name_B)
        # print(place_A, '\n', place_B, '\n', root, '\n\n', )
        find_same_days_and_process(place_A, place_B, root, star_dir, resolution,
                                   month_names=AllGPSDataLocations.all_months)


process_for_all(main_root=AllGPSDataLocations.main_root, result_roots=AllGPSDataLocations.result_roots,
                  star_dir=star_dir, resolution=resolution)
