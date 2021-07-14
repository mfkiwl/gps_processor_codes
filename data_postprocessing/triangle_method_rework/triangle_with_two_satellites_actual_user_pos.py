import copy
import time
import sys
from numpy import *
from itertools import chain
import os
import math
from scipy import signal
from tqdm import tqdm

sys.path.insert(1, '/Users/kelemensz/Documents/Research/GPS/gps_processor_codes/data_postprocessing')
from data_locations_handler import AllGPSDataLocations
from triangle_method_rework.rework_lib import v_length, parse_sats_data, Defaults, \
    get_mean_user_position, get_common, get_proportion_v3, unit, normvec, transform_matrix, ecef_to_gcs, \
    get_user_and_sat_positions, cart2sph, get_ij_on_map, get_global_stars, raw_results_to_GCS, \
    get_mean_direction_in_system_over_time, consecutively_group_results, save_matrices, plot_mollweid, create_dir, \
    plot_save_imshow_3_maps, rotate2darray, input_fileset_of_the_day, is_string_in_filenames, filter_sats_within_solid_angle_above_the_sky
sys.path.insert(1, '/Users/kelemensz/Documents/Research/GPS/gps_processor_codes/utility')
from frecvently_used_functions import is_all_data, is_reliable, too_bad_positions_std


def calc_for_two_satellites(user_position, satA, satB):
    S_A = satA[:3]
    S_B = satB[:3]
    t1 = satA[3]
    t2 = satB[3]
    AS = v_length(S_A, user_position)
    BS = v_length(S_B, user_position)

    # A_SA_i = unit(array(user_position) - array(S_A))
    # B_SB_j = unit(array(user_position) - array(S_B))
    A_SA_i = unit(array(S_A) - array(user_position))
    B_SB_j = unit(array(S_B) - array(user_position))
    n = A_SA_i - B_SB_j

    mod_n = v_length(array([0, 0, 0]), n)
    return unit(n), get_proportion_v3(AS, t1, BS, t2), mod_n


# ===============================================================================================================================================================================


def process_one_epoch(user_position, sats_data):
    direction_value = []
    nr_of_sats = len(sats_data)
    # nr_of_sats = int(len(sats_data)/4)
    for i in range(nr_of_sats-1):
        for j in range(i+1, nr_of_sats):
            direction_value.append(calc_for_two_satellites(user_position, sats_data[i], sats_data[j]))
    return direction_value


def process_one_epoch_symmetrized(user_position, sats_data):
    direction_value = []
    nr_of_sats = len(sats_data)
    for i in range(nr_of_sats-1):
        for j in range(i+1, nr_of_sats):
            direction_value.append(calc_for_two_satellites(user_position, sats_data[i], sats_data[j]))
            direction_value.append(calc_for_two_satellites(user_position, sats_data[j], sats_data[i]))
    return direction_value
# ===============================================================================================================================================================================


def process_all_(user_position, sat_data, symmetrized):
    raw_results = {}
    if symmetrized:
        for epoch, sats_dat in tqdm(sat_data.items()):  # k = epoch index
            raw_results[epoch] = process_one_epoch_symmetrized(user_position, sats_dat)
    else:
        for epoch, sats_dat in tqdm(sat_data.items()):  # k = epoch index
            raw_results[epoch] = process_one_epoch(user_position, sats_dat)
    return raw_results


def smooth_user_positions(user_positions):
    user_positionsT = user_positions.T
    x_smooth = signal.savgol_filter(user_positionsT[0], 53, 3)
    y_smooth = signal.savgol_filter(user_positionsT[1], 53, 3)
    z_smooth = signal.savgol_filter(user_positionsT[2], 53, 3)
    return asarray(list(zip(x_smooth, y_smooth, z_smooth)))


def process_all(user_position, sat_data, symmetrized):
    raw_results = {}

    smooth_positions = smooth_user_positions(user_position)
    if symmetrized:
        for epoch, sats_dat in tqdm(sat_data.items()):  # k = epoch index
            raw_results[epoch] = process_one_epoch_symmetrized(smooth_positions[int(epoch)], sats_dat)
    else:
        for epoch, sats_dat in tqdm(sat_data.items()):  # k = epoch index
            raw_results[epoch] = process_one_epoch(smooth_positions[int(epoch)], sats_dat)
    return raw_results


# =================================================================================================


def get_raw_results_of_day(sat_data, mean_position, symmetrized):
    raw_results = process_all(mean_position, sat_data, symmetrized=symmetrized)  # dictionary, keys are the index of the epochs
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


def process_one_day(day_root, day_result_location, star_dir, resolution, filenames, fill_out=0.0, n_filter=False, symmetrized=True):
    # fill_out=0.0 in case when the "r * (r-1)^2" values are in the data matrix

    if day_result_location is None:
        print("No day_result directory is given...process stops here!")
        sys.exit()
    Nx, Ny, Nz, S, D = get_global_stars(star_dir, os.path.dirname(day_root))
    if isinstance(Nx, int):
        sys.exit()
    nr_of_galactic_frames = len(Nx)
    star_directions_in_GCS = Defaults.get_star_directions_in_GCS(Nx, Ny, Nz, D, S)

    GCS_all = [array([Nx[i], Ny[i], Nz[i]]) for i in range(nr_of_galactic_frames)]
    # =================================================================================================================
    # =================================================================================================================
    sat_data = get_user_and_sat_positions(day_root, get_u=False)
    mean_position, user_positions = get_mean_user_position(day_root, both=True)

    sat_data = filter_sats_within_solid_angle_above_the_sky(sat_data, mean_position, angular_distance=70.0)


    raw_results_ECEF_AB = get_raw_results_of_day(sat_data, user_positions, symmetrized=symmetrized)
    groupped_raw_results_ECEF_AB = consecutively_group_results(raw_results_ECEF_AB, nr_of_galactic_frames)
    groupped_raw_results_GCS_AB = raw_results_to_GCS(groupped_raw_results_ECEF_AB, GCS_all)
    # groupped_raw_results_GCS_AB = groupped_raw_results_ECEF_AB

    raw_results_GCS = list(chain(*groupped_raw_results_GCS_AB))
    if isinstance(n_filter, float):
        raw_results_GCS = asarray(raw_results_GCS)
        print(len(raw_results_GCS))
        raw_results_GCS = raw_results_GCS[raw_results_GCS[:, -1] < n_filter]
        print(len(raw_results_GCS))
    print("Raw results in GCS are in! (data size)", len(raw_results_GCS))

    day_data, day_count, day_cmap_n_mod = process_raw_GCS_data(raw_results_GCS, resolution)

    day_data = nan_to_num(day_data, nan=fill_out)
    day_cmap_n_mod = nan_to_num(day_cmap_n_mod, nan=fill_out)
    day_count = nan_to_num(day_count, nan=fill_out)

    # =================================================================================================================
    # =================================================================================================================

    save_matrices([day_data, day_count, day_cmap_n_mod], filenames, day_result_location, resolution)
    plot_mollweid(day_count, star_directions_in_GCS, day_result_location, filenames[1],
                  str(int(degrees(resolution))), anot=True)
    plot_mollweid(divide(day_data, day_count), star_directions_in_GCS, day_result_location, filenames[0],
                  str(int(degrees(resolution))), anot=True)
    plot_mollweid(divide(day_cmap_n_mod, day_count), star_directions_in_GCS, day_result_location, filenames[2],
                  str(int(degrees(resolution))), anot=True)
    # plot_save_imshow(divide(day_cmap_n_mod, day_count), root, "n_mod")

    return day_data, day_count, day_cmap_n_mod
    # return 0,0,0


def find_days_and_process(user_sat_data_root, result_root, star_dir, resolution, month_names=AllGPSDataLocations.all_months, symmetrized=True):
    filenames = Defaults.triangle_result_names_n_filter_notsym if not symmetrized else Defaults.triangle_result_names_n_filter_sym
    print(filenames)
    d = 0
    all_hist = []
    all_value = []
    all_n_mod = []

    if os.path.isdir(user_sat_data_root) and os.path.isdir(result_root):
        months = [f.path for f in os.scandir(user_sat_data_root) if f.is_dir()]
        for month in months[:]:
            month_name = os.path.split(month)[-1]
            condition = month_name in month_names
            # condition = d < 1
            # condition = month_name in ['aprilis']
            if condition:
                print(month_name)
                day_folders = [f.path for f in os.scandir(month) if f.is_dir()]  # and len(str(f)) == 12 and str(f)[-8:].isnumeric()]
                print("Number of days: ", len(day_folders))
                for day_folder in day_folders[:]:
                    try:
                        start = time.time()
                        date = str(os.path.split(day_folder)[-1])[-8:]
                        fileset_day = input_fileset_of_the_day(day_folder, True)

                        if fileset_day != 0 and not too_bad_positions_std(day_folder, Defaults.USER_POSITIONS_FILENAME.get('user')):
                            result_month = create_dir(result_root, month_name)
                            result_day = create_dir(result_month, date)
                            dayfiles = [f for f in os.listdir(result_day) if os.path.isfile(os.path.join(result_day, f))]
                            # if len(os.listdir(result_day)) == 0:
                            if not is_string_in_filenames(filenames[0], dayfiles):
                                print(" Data will be processed from: ", os.path.split(day_folder)[-1], "    ",
                                      "\n", "Index of the process: ", d, "\n")
                                day_folder = os.path.join(day_folder, "allsatellites")

                                value, hist, n_mod = process_one_day(day_folder, result_day, star_dir, resolution, filenames, n_filter=1.0, symmetrized=symmetrized)
                                all_hist.append(hist)
                                all_value.append(value)
                                all_n_mod.append(n_mod)
                                d += 1
                                print('Elapsed time of the current day: ', time.time() - start, date)
                            else:
                                print("{} already processed!".format(date))
                    except:
                        pass
                    else:
                        print("\n Data not found for: ", date, "\n")

        print("Total nr of days: ", d)
        all_hist = sum(array(all_hist), axis=0)
        all_value = sum(array(all_value), axis=0)
        all_n_mod = sum(array(all_n_mod), axis=0)

        save_matrices([all_value, all_hist, all_n_mod], filenames, result_root, resolution)

        # all_value_av = divide(all_value, all_hist)
        # all_n_mod_av = divide(all_n_mod, all_hist)

        # all_hist = rotate2darray(nan_to_num(all_hist, nan=0.0))
        # all_value_av = rotate2darray(nan_to_num(all_value_av, nan=0.0))
        # all_n_mod_av = rotate2darray(nan_to_num(all_n_mod_av, nan=0.0))

        # plot_save_imshow_3_maps([all_hist, all_value_av, all_n_mod_av], ["Histogram", "(|1-r|/|n|)^(1/4)", "<n_mod>"],
        #                         result_root, resolution="5", logplot=False)
    # plot_save_imshow(all_value_av, result_path, "(|1-r|/|n|)^(1/4)")

    print("Nr of days: ", d)

# =================================================================================================
# =================================================================================================


star_dir = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020"
resolution = radians(5.0)
needed_files = [Defaults.USER_POSITIONS_FILENAME, Defaults.SAT_POS_FILENAMES.get('no_ID'),
                Defaults.SAT_POS_FILENAMES.get('ID')]


user_sat_data_root = AllGPSDataLocations.user_and_satellites.get('CUTB')

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_method_two_sats/CUTB_vertical_cone70"
results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_method_two_sats_smart/CUTB_vertical_cone70"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_method_two_sats/HKKS"
results_root = create_dir(results_root, r'r_inv_r_twoSats')



find_days_and_process(user_sat_data_root, result_root=results_root, star_dir=star_dir, resolution=resolution, symmetrized=False)
