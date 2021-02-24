import time
from data_postprocessing.postprocess_utility.general_functions import *
from numpy import *
import pandas as pd
from itertools import chain, islice
import os
import math
import copy as cp

s_o_l = 1.0  # 3.0*10**8
# satellite_positions = "all_sats_pos_time.csv"
satellite_positions = "sats_pos_time_id.csv"


def group_results(results_dict, l):
    groupped = []
    group = []
    n = len(list(results_dict.keys()))
    l_groups = int(n / l)
    i = 0
    if l_groups > 1:
        for k, v in results_dict.items():
            if i >= l_groups:
                groupped.append(list(chain(*group)))
                group = [v]
                i = 0
                continue
            group.append(v)
            i += 1
        groupped.append(list(chain(*group)))
    return groupped


def parse_sats_data(sat_data_file):
    data = {}
    epoch = []
    with open(sat_data_file) as in_file:
        lineList = [line.rstrip('\n').split(",") for line in in_file]
        for line in lineList[1:]:
            if line[0][0] == "E":
                data[line[-2]] = epoch
                epoch = []
                continue
            epoch.append([float(line[0]), float(line[1]), float(line[2]), float(line[3]), line[4]])
    return data


def get_comon(vA, vB):
    data = []
    for a in vA:
        # print(a)
        for b in vB:
            if a[4] == b[4]:
                # print(dr_lenght(a[:3], b[:3]))
                data.append([a, b])
                break
    return data


def prepare_u_and_s_positions(path, get_u=True):
    sat_data = parse_sats_data(os.path.join(path, satellite_positions))
    if get_u:
        user_ = pd.read_csv(path + '/user_pos_allsatellites.csv', skiprows=1).values  # .transpose()
        user_mean = mean(user_, axis=0).astype('float64')
        return user_mean, sat_data
    return sat_data


def get_mean_position(path, positions_file):
    user_ = pd.read_csv(os.path.join(path, positions_file), skiprows=1).values  # .transpose()
    return mean(user_, axis=0).astype('float64')


def extract_common_sats(satA, satB):
    comon_parts = {}
    for kA, vA in satA.items():
        vB = satB.get(kA, None)
        if vB:
            comon = get_comon(vA, vB)
            comon_parts[kA] = comon
    return comon_parts


def calc_for_one_sattelite_satID(posA, posB, sat_pr_both):
    S_A = sat_pr_both[0][:3]
    S_B = sat_pr_both[1][:3]
    t1 = sat_pr_both[0][3]
    t2 = sat_pr_both[1][3]
    AS = dr_lenght(S_A, posA)
    BS = dr_lenght(S_B, posB)

    A_SA_i = unit(array(posA) - array(S_A))
    B_SB_j = unit(array(posB) - array(S_B))
    n = A_SA_i - B_SB_j

    mod_n = dr_lenght(array([0, 0, 0]), n)
    return {sat_pr_both[0][4]: [unit(n), get_proportion_v3(AS, t1, BS, t2), mod_n]}


def get_proportion_v3(AS, t1, BS, t2):
    r = float(AS) * float(t2) / float(BS) / float(t1)
    return r - 1.0 / float(r)


def process_one_epoch(posA, posB, data_from_common_sats):
    dirrection_value = []
    for sats in data_from_common_sats:
        dirrection_value.append(calc_for_one_sattelite_satID(posA, posB, sats))
    return dirrection_value


def process_all(posA, posB, common_data):
    raw_results = {}
    # for k, v in islice(common_data.items(), 100):  # k = epoch index
    for k, v in common_data.items():  # k = epoch index
        raw_results[k] = process_one_epoch(posA, posB, v)
    return raw_results


def reorganize_dictlist_dictlist_structure():
    """Initial structure: {e1:[{id1:[n, m, |n|]}, {id2:[n, m, |n|]},..], e2:[{id1:[n, m, |n|]}, {id2:[n, m, |n|]},..],
                            ...}
    Final structure:    [[e1, id1, n, m, |n|],
                        [e1, id2, n, m, |n|],
                        ...
                        [ek, id1, n, m, |n|]]
    """




def raw_results_to_GCS_satID(results_ECEF, GCS):
    results_GCS = []
    nr_groups = len(results_ECEF)
    nr_axis = len(GCS)
    print("Systems and groups: ", nr_axis, nr_groups)
    # if nr_axis == nr_groups:
    for i in range(nr_groups):
        results_group = results_ECEF[i]
        GCS_group = GCS[i]
        results_group_GCS = []
        for line in results_group:
            satID = list(line.keys())[0]
            direction, value, mod_n = list(line.values())
            # results_group_GCS.append([ecef_to_gcs(GCS_group, direction), value, mod_n])
            results_GCS.append({satID: [ecef_to_gcs(GCS_group, direction), value, mod_n]})
    return results_GCS


# =================================================================================================
# =================================================================================================

def get_raw_results_using_mean_positions(pathA, pathB, mean_positions):
    posA = mean_positions[0]
    posB = mean_positions[1]
    dataA = prepare_u_and_s_positions(pathA, False)
    dataB = prepare_u_and_s_positions(pathB, False)
    comonAB = extract_common_sats(dataA, dataB)
    raw_results = process_all(posA, posB, comonAB)  # dictionary, keys are the index of the epochs
    return raw_results


def process_raw_GCS_data(raw_results_GCS, resolution):
    theta_max = math.pi
    phi_max = math.pi / 2.0
    rot_theta = arange(-theta_max, theta_max, resolution)
    rot_phi = arange(-phi_max, phi_max, resolution)
    cmap_v = zeros((len(rot_theta), len(rot_phi)))
    cmap_count = zeros((len(rot_theta), len(rot_phi)))
    cmap_mod_n = zeros((len(rot_theta), len(rot_phi)))
    for direction, value, mod_n in raw_results_GCS:
        i, j = get_ij_on_map(direction, resolution)
        cmap_v[i][j] += value
        cmap_mod_n[i][j] += mod_n
        cmap_count[i][j] += 1
    cmap_count[cmap_count < 1] = 0
    return cmap_v, cmap_count, cmap_mod_n


def test_one_day_n_dirrections(pathA, pathB, star_dir, resolution, mean_positions=None, root=None, fill_out=0.0,
                               day=None):

    Nx, Ny, Nz, S, D = get_global_stars(star_dir, os.path.dirname(pathB))
    l = len(Nx)
    star_directions_in_GCS = {'GNP': get_mean_direction_over_time(array([Nx, Ny, Nz]), Nz),
                              'GC': get_mean_direction_over_time(array([Nx, Ny, Nz]), Nx),
                              'OY of GCS': get_mean_direction_over_time(array([Nx, Ny, Nz]), Ny),
                              'Sadalmelic': get_mean_direction_over_time(array([Nx, Ny, Nz]), S),
                              'Denebola': get_mean_direction_over_time(array([Nx, Ny, Nz]), D)}
    # =================================================================================================================
    GCS_all = [array([Nx[i], Ny[i], Nz[i]]) for i in range(l)]

    raw_results_ECEF_AB = get_raw_results_using_mean_positions(pathA, pathB, mean_positions)
    # groupped_raw_results_ECEF_AB = group_results(raw_results_ECEF_AB, l)
    # groupped_raw_results_GCS_AB = raw_results_to_GCS_satID(groupped_raw_results_ECEF_AB, GCS_all)
    sat_id = "R21"
    # sat_ids = ["R21", "R22", "R23", "R09", "R08", "R07", "R01", "G02", "G05", "G13"]
    # filter_collected_triangles_many(raw_results_ECEF_AB, sat_identifiers=sat_ids)

    sat_data = filter_collected_triangles(raw_results_ECEF_AB, sat_identifier=sat_id)
    spherical_data = sat_data_to_spherical(sat_data)
    plot_phi_theta_n(spherical_data, sat_id, day)
    # -----------------------------------------------------------------------------------------------------------------

    # plot_mollweid(day_count, star_directions_in_GCS, root, "histogram", str(int(degrees(resolution))), anot=True)

    return 0, 0, 0


def find_same_days_and_process(path_A, path_B, result_path, needed_files, star_dir, resolution):
    d = 0
    if os.path.isdir(path_A) and os.path.isdir(path_B) and os.path.isdir(result_path):
        month_pairs = find_corresponding_dirs_in_different_roots(path_A, path_B)
        # mean_pos_A = get_mean_pos_from_root(path_A, needed_files[0], max_deviations=5)
        # mean_pos_B = get_mean_pos_from_root(path_B, needed_files[0], max_deviations=0.5)  # NZLD eseten 0.2
        mean_pos_A = [0.0, 0.0, 0.0]
        mean_pos_B = [0.0, 0.0, 0.0]
        for A_month, B_month in month_pairs:
            month_name = os.path.split(A_month)[-1]
            condition = month_name in ["januar"]  # , "marcius", "aprilis", "majus", "junius", "november"]
            if condition:
                print(month_name)
                day_pairs = find_corresponding_dirs_in_different_roots(A_month, B_month)
                print("Number of days: ", len(day_pairs))
                for A_day, B_day in day_pairs[:1]:
                    start = time.time()
                    date = str(os.path.split(B_day)[-1])[-8:]
                    cond2 = is_all_data(A_day, needed_files[1:], True) and is_all_data(B_day, needed_files[1:], True)
                    if cond2:
                        result_month = create_dir(result_path, month_name)
                        result_day = create_dir(result_month, date)
                        print(" Data will be processed from: ", os.path.split(A_day)[-1], "    ",
                              os.path.split(B_day)[-1], "\n")

                        A_day = os.path.join(A_day, "allsatellites")
                        B_day = os.path.join(B_day, "allsatellites")

                        posUA = cp.deepcopy(mean_pos_A)
                        posUB = cp.deepcopy(mean_pos_B)
                        if is_all_data(A_day, needed_files[:1]):
                            print("Actual position considered for: {}".format(str(A_day).split("/")[-2]))
                            posUA = get_mean_position(A_day, needed_files[0])
                        if is_all_data(B_day, needed_files[:1]):
                            print("Actual position considered for: {}".format(str(B_day).split("/")[-2]))
                            posUB = get_mean_position(B_day, needed_files[0])

                        value, hist, n_mod = test_one_day_n_dirrections(A_day, B_day, star_dir, resolution,
                                                                        mean_positions=[posUA, posUB],
                                                                        root=result_day, fill_out=0.0, day=date)
                    else:
                        print("\n Data not found for: ", date, "\n")
                    print('Elapsed time of the current day: ', time.time() - start, date)

    print("Nr of days: ", d)


# =================================================================================================
# =================================================================================================


star_dir = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020"
resolution = radians(5.0)
needed_files = ["user_pos_allsatellites.csv", satellite_positions]

# --------------------------------------------PERTH-Hong-Kong--------------------------------------------
place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/PERTH_daily_measurements"
place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_HKKS"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_PERTH/r_inv_r_symmetrized"
results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/pairs_by_identifier/HKKS_PERTH/r_inv_r_symmetrized"


# --------------------------------------------Hong-Kong-India--------------------------------------------
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_HKKS"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_IIGC"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/pairs_by_identifier/HKKS_IIGC/r_inv_r_symmetrized"


find_same_days_and_process(place_A, place_B, results_root, needed_files, star_dir, resolution)