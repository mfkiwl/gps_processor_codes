from numpy import *
import pandas as pd
import matplotlib.pyplot as plt

from itertools import chain
import pymap3d as pm
import os
from vpython import rotate
import pylab as pl
import math
from scipy.ndimage import rotate

s_o_l = 1.0  # 3.0*10**8


def rotateAntiClockwise(array):
    return rotate(array, 90)


def getSTARS_for_galactic_system(path):
    # D = normvec(pd.read_csv(path+'/star_d_positions.csv',skiprows=0).values)
    # S = normvec(pd.read_csv(path+'/star_s_positions.csv',skiprows=0).values)
    D = normvec(pd.read_csv(path + '/Denebola_positions.csv', skiprows=0).values)
    S = normvec(pd.read_csv(path + '/Sadalmelik_positions.csv', skiprows=0).values)

    sun = pd.read_csv(path + '/SUN_positions.csv', skiprows=0).values
    galactic_n_p = pd.read_csv(path + '/GNP_positions.csv', skiprows=0).values
    galactic_center = pd.read_csv(path + '/GC_positions.csv', skiprows=0).values

    # star_1 = pd.read_csv(path+'/star_Mirphak_positions.csv',skiprows=0).values
    # star_2 = pd.read_csv(path+'/star_Vega_positions.csv',skiprows=0).values

    # nunki = pd.read_csv(path+'/star_Nunki_positions.csv',skiprows=0).values
    # capella = pd.read_csv(path+'/star_Capella_positions.csv',skiprows=0).values
    nunki = pd.read_csv(path + '/Nunki_positions.csv', skiprows=0).values
    capella = pd.read_csv(path + '/Capella_positions.csv', skiprows=0).values
    return sun, galactic_n_p, galactic_center, nunki, capella, S, D


def get_mean_direction_over_time(systems, directions):
    l = min(len(directions), len(systems[0]))
    phi = 0.
    theta = 0.
    for i in range(l):
        out = cartesian_to_galactic(array([systems[0][i], systems[1][i], systems[2][i]]), directions[i])
        theta += out[0]
        phi += out[1]
    return theta / float(l), phi / float(l)


def normvec(vec):
    a = empty(shape(vec))
    for i in range(len(vec)):
        a[i] = unit(vec[i])
    return a


def scal(v1, v2):
    return vdot(array(v1), array(v2))


def unit(v):
    return 1 / sqrt(scal(v, v)) * array(v)


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
                data[line[-1]] = epoch
                epoch = []
                continue
            epoch.append([float(line[0]), float(line[1]), float(line[2]), float(line[3])])
    return data


def dr_lenght(v1, v2):
    return sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def get_comon(vA, vB):
    data = []
    for a in vA:
        # print(a)
        for b in vB:
            if dr_lenght(a[:3], b[:3]) < 500:
                # print(dr_lenght(a[:3], b[:3]))
                data.append([a, b])
    return data


def prepare_data(path):
    user_ = pd.read_csv(path + '/user_pos_allsatellites.csv', skiprows=1).values  # .transpose()
    user_mean = mean(user_, axis=0).astype('float64')
    sat_data = parse_sats_data(os.path.join(path, "all_sats_pos_time.csv"))
    return user_mean, sat_data


def extract_common_sats(satA, satB):
    comon_parts = {}
    for kA, vA in satA.items():
        vB = satB.get(kA, None)
        if vB:
            comon = get_comon(vA, vB)
            comon_parts[kA] = comon
    return comon_parts


# def calc_for_one_sattelite(posA, posB, sat_pr_both):
#     S_A = sat_pr_both[0][:3]
#     S_B = sat_pr_both[1][:3]
#     t1 = sat_pr_both[0][3]  # / s_o_l
#     t2 = sat_pr_both[1][3]  # / s_o_l
#     AS = dr_lenght(S_A, posA)
#     BS = dr_lenght(S_B, posB)
#     return unit((array(posA) + array(posB)) / 2.0 - array(S_A)), get_proportion_v0(AS, t1, BS, t2)


def calc_for_one_sattelite_beta(posA, posB, sat_pr_both):
    S_A = sat_pr_both[0][:3]
    S_B = sat_pr_both[1][:3]
    t1 = sat_pr_both[0][3]  # / s_o_l
    t2 = sat_pr_both[1][3]  # / s_o_l
    AS = dr_lenght(S_A, posA)
    BS = dr_lenght(S_B, posB)

    A_SA_i = unit(array(posA) - array(S_A))
    B_SB_j = unit(array(posB) - array(S_B))
    n = A_SA_i - B_SB_j

    mod_n = dr_lenght(array([0, 0, 0]), n)
    return unit(n), get_proportion_v3(AS, t1, BS, t2), mod_n


def get_proportion_v0(AS, t1, BS, t2):
    return (float(AS) * float(t2)) / (float(BS) * float(t1))


def get_proportion_v1(AS, t1, BS, t2):
    r = float(AS) * float(t2) / float(BS) / float(t1)
    return r * (r - 1) ** 2


def get_proportion_v2(AS, t1, BS, t2):
    r = float(AS) * float(t2) / float(BS) / float(t1)
    return abs(r - 1)


def get_proportion_v3(AS, t1, BS, t2):
    r = float(AS) * float(t2) / float(BS) / float(t1)
    return r - 1.0 / float(r)


# ===============================================================================================================================================================================

def process_one_epoch(posA, posB, data_from_common_sats):
    dirrection_value = []
    for sats in data_from_common_sats:
        # dirrection_value.append(calc_for_one_sattelite(posA, posB, sats))
        dirrection_value.append(calc_for_one_sattelite_beta(posA, posB, sats))

    return dirrection_value


# ===============================================================================================================================================================================

def process_all(posA, posB, common_data):
    raw_results = {}
    for k, v in common_data.items():  # k = epoch index
        raw_results[k] = process_one_epoch(posA, posB, v)
    return raw_results


def raw_results_to_GCS(results_ECEF, GCS):
    results_GCS = []
    nr_groups = len(results_ECEF)
    nr_axis = len(GCS)
    print("Systems and groups: ", nr_axis, nr_groups)
    # if nr_axis == nr_groups:
    for i in range(nr_groups):
        results_group = results_ECEF[i]
        GCS_group = GCS[i]
        results_group_GCS = []
        for direction, value, mod_n in results_group:
            results_group_GCS.append([ecef_to_gcs(GCS_group, direction), value, mod_n])
        results_GCS.append(results_group_GCS)
    return results_GCS


def transform_matrix(f1, f2):  # transforms from f1 to f2
    R = array([
        [dot(f2[0], f1[0]), dot(f2[0], f1[1]), dot(f2[0], f1[2])],
        [dot(f2[1], f1[0]), dot(f2[1], f1[1]), dot(f2[1], f1[2])],
        [dot(f2[2], f1[0]), dot(f2[2], f1[1]), dot(f2[2], f1[2])]
    ])
    return R


def get_theta_phi(v):
    if 0.99 < v[0] and 1.01 > v[0] and -0.01 < v[1] and 0.01 > v[1] and -0.01 < v[2] and 0.01 > v[2]:
        return 0.0, 0.0
    if 0.99 < v[1] and 1.01 > v[1] and -0.01 < v[0] and 0.01 > v[0] and -0.01 < v[2] and 0.01 > v[2]:
        return 90.0, 0.0
    if 0.99 < v[2] and 1.01 > v[2] and -0.01 < v[1] and 0.01 > v[1] and -0.01 < v[0] and 0.01 > v[0]:
        return 0.0, 90.0
    return None, None


def ecef_to_gcs(system, vector):
    vector = unit(vector)
    system = normvec(system)
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R = transform_matrix(ECEF, system)
    return around(R.dot(vector), decimals=3)


def cartesian_to_galactic(system, vector):
    vector = unit(vector)
    system = normvec(system)
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R = transform_matrix(ECEF, system)
    vector_in_system = around(R.dot(vector), decimals=3)
    theta, phi = get_theta_phi(vector_in_system)
    if not theta and not phi:
        phi, theta, _ = pm.ecef2geodetic(vector_in_system[0], vector_in_system[1], vector_in_system[2])
        phi = 90 - degrees(arccos(scal(ECEF[2], vector_in_system)))
    return theta, phi


def cart2sph(x, y, z):
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    el = arctan2(z, hxy)
    az = arctan2(y, x)
    return az, el


def cartesian_to_spherical(vector):
    theta, phi = get_theta_phi(vector)
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    if not theta and not phi:
        phi, theta, _ = pm.ecef2geodetic(vector[0], vector[1], vector[2])
        phi = 90 - degrees(arccos(scal(ECEF[2], vector)))
    return radians(theta), radians(phi)


# =================================================================================================
# =================================================================================================


def get_global_stars(star_dir, directory):
    star_dirs = [f.path for f in os.scandir(star_dir) if f.is_dir()]
    for stardir in star_dirs:
        if directory[-8:] == stardir[-8:]:
            return get_stars_for_CMBR_cmap(stardir)


def get_stars_for_CMBR_cmap(path):
    sun, galactic_n_p, galactic_center, nunki, galactic_anti_center, S, D = getSTARS_for_galactic_system(path)

    sun = normvec(sun)
    galactic_n_p = normvec(galactic_n_p)
    galactic_center = normvec(galactic_center)
    galactic_anti_center = normvec(galactic_anti_center)
    nunki = normvec(nunki)
    Nx = normvec(galactic_center).astype('float64')
    Nz = normvec(galactic_n_p).astype('float64')
    Ny = -get_third_dir_by_cross_product(Nz, Nx)

    return Nx, Ny, Nz, S, D


def get_third_dir_by_cross_product(A, B):
    l = min(len(A), len(B))
    C = empty((l, 3))
    for i in range(l):
        C[i] = cross(A[i], B[i])
    return normvec(C)


def get_raw_results(pathA, pathB):
    posA, dataA = prepare_data(pathA)
    posB, dataB = prepare_data(pathB)
    comonAB = extract_common_sats(dataA, dataB)
    raw_results = process_all(posA, posB, comonAB)
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


def get_ij_on_map(direction, resolution):
    theta_max = math.pi
    phi_max = math.pi / 2.0
    rot_theta = arange(-theta_max, theta_max, resolution)
    rot_phi = arange(-phi_max, phi_max, resolution)
    # theta, phi = cartesian_to_spherical(direction)
    theta, phi = cart2sph(direction[0], direction[1], direction[2])
    I_f = 0
    J_f = 0
    l_theta = int(len(rot_theta) / 2.0)
    l_phi = int(len(rot_phi) / 2.0)
    for i in range(-l_theta, l_theta):
        if i * resolution > theta:
            I_f = i + l_theta
            break
    for i in range(-l_phi, l_phi):
        if i * resolution > phi:
            J_f = i + l_phi
            break
    return I_f, J_f


def plot_mollweid(matrix, star_directions, root_directory, name, resolution, anot=True):
    matrix = nan_to_num(matrix, nan=0.0)
    child_dirname = os.path.split(root_directory)[-1] + '_24h'
    plt.clf()
    try:
        cmap_save = pd.DataFrame(matrix)
        cmap_save.to_csv(os.path.join(root_directory, 'GCS_' + child_dirname + "_" + name + "_" + resolution + '.csv'),
                         index=True)
        print('Matrix.csv saved!')
    except:
        pass
    ra = linspace(-math.pi, math.pi, len(matrix))
    dec = linspace(-math.pi / 2, math.pi / 2, len(matrix[0]))

    X, Y = meshgrid(ra, dec)
    Z = matrix.T
    plt.figure()
    ax = pl.subplot(111)  # , projection = 'mollweide')
    fig = ax.contourf(X, Y, Z, 100)
    # fig = ax.imshow(rot90(fliplr(matrix), -1))
    if anot:
        for k, v in star_directions.items():
            add_star_annotated(v[0], v[1], k, ax)

    # ax.set_title('---$(24h)$', fontsize=15)  # , fontweight='bold')
    plt.xlabel(r'$\theta$', fontsize=15)  # Italic font method
    plt.ylabel(r'$\phi$', fontsize=15)  # Bold font method without fontweight parameters
    pl.colorbar(fig)
    ax.grid()
    # ax.contour(X,Y,Z,10,colors='k')
    # pl.show()
    fig_name1 = os.path.join(root_directory, 'GCS_' + child_dirname + "_" + name + "_" + resolution + '.png')
    pl.savefig(fig_name1, bbox_inches='tight')


def plot_save_imshow(matrix, root_directory, name, resolution="5", logplot=False):
    child_dirname = os.path.split(root_directory)[-1] + "_" + name + '_24h' + "_" + resolution
    plt.clf()
    try:
        cmap_save = pd.DataFrame(matrix)
        cmap_save.to_csv(os.path.join(root_directory, child_dirname + '.csv'), index=True)
        print('Matrix.csv saved: ', name, "   ", os.path.split(root_directory)[-1])
    except:
        pass
    if logplot:
        matrix = around(log(matrix + 1.0), decimals=0)

    # matrix = sqrt(sqrt(matrix))
    # a = 6
    # b = a #* 2
    # matrix = rebin(matrix, (int(len(matrix)/b), int(len(matrix[0])/a)))

    plt.imshow(matrix)
    plt.colorbar()

    # plt.title("(r*(r-1)^2)^0.25")
    # plt.show()
    fig_name = os.path.join(root_directory, "imshow_" + child_dirname + '.png')
    plt.savefig(fig_name, bbox_inches='tight')
    plt.clf()


def add_star_annotated(theta, phi, name, ax):
    theta = radians(theta)
    phi = radians(phi)
    ax.text(theta, phi, name, fontsize=12)
    ax.scatter(theta, phi, marker='x', c='k', s=15)


# ax.annotate(name,
#            xy=array(theta, phi),
#            xycoords='data')

#            arrowprops=
#                dict(facecolor='black', shrink=0.05),
#                horizontalalignment='left',
#                verticalalignment='top')


def plot_save_imshow_3_maps(matrices, names, root_directory, resolution="5", logplot=False, show=False, fill_out=0.0):
    matrices = array(matrices)
    # plt.clf()
    if logplot:
        matrices[0] = around(log(matrices[0] + 1.0), decimals=0)

    # matrices[1] = sqrt(matrices[1])
    a = 4
    b = a  # * 2
    # matrices[1] = rebin(matrices[1], (int(len(matrices[1])/b), int(len(matrices[1][0])/a)))
    fig, axis = plt.subplots(len(matrices), 1)
    fig.subplots_adjust(left=0.02, bottom=0.1, right=0.95, top=0.94, wspace=0.8, hspace=0.5)

    for matrix, name, ax in zip(matrices, names, axis):
        nan_to_num(matrix, nan=fill_out)
        sp = ax.imshow(matrix)
        ax.set_title(name)
        fig.colorbar(sp, ax=ax)

    if not show:
        print('Matrix.csv saved: ', os.path.split(root_directory)[-1])
        figname = os.path.split(root_directory)[-1] + "_" + resolution
        fig_name_with_path = os.path.join(root_directory, figname + '.png')
        fig.savefig(fig_name_with_path, bbox_inches='tight')
    else:
        plt.show()
    # plt.clf()


def save_matrices(matrices, names, directory, resolution):
    for matrix, name in zip(matrices, names):
        cmap_save = pd.DataFrame(matrix)
        cmap_save.to_csv(os.path.join(directory, name + "_not_averaged_" + str(int(degrees(resolution))) + '.csv'),
                         index=True)
        print('Matrix.csv saved!: ', name)


def process_one_day(pathA, pathB, star_dir, resolution, root=None,
                    fill_out=0.0):  # fill_out=0.0 in case when the "r * (r-1)^2" values are in the data matrix

    Nx, Ny, Nz, S, D = get_global_stars(star_dir, os.path.dirname(pathB))
    l = len(Nx)
    star_directions_in_GCS = {'GNP': get_mean_direction_over_time(array([Nx, Ny, Nz]), Nz),
                              'GC': get_mean_direction_over_time(array([Nx, Ny, Nz]), Nx),
                              'OY of GCS': get_mean_direction_over_time(array([Nx, Ny, Nz]), Ny),
                              'Sadalmelic': get_mean_direction_over_time(array([Nx, Ny, Nz]), S),
                              'Denebola': get_mean_direction_over_time(array([Nx, Ny, Nz]), D)}
    # =================================================================================================================
    GCS_all = [array([Nx[i], Ny[i], Nz[i]]) for i in range(l)]
    raw_results_ECEF = get_raw_results(pathA, pathB)

    groupped_raw_results_ECEF = group_results(raw_results_ECEF, l)

    groupped_raw_results_GCS = raw_results_to_GCS(groupped_raw_results_ECEF, GCS_all)

    raw_results_GCS = list(chain(*groupped_raw_results_GCS))  # [::10000]
    print("Raw results in GCS are in! (data size)", len(raw_results_GCS))

    day_data, day_count, day_cmap_n_mod = process_raw_GCS_data(raw_results_GCS, resolution)

    day_data = nan_to_num(day_data, nan=fill_out)
    day_cmap_n_mod = nan_to_num(day_cmap_n_mod, nan=fill_out)
    day_count = nan_to_num(day_count, nan=fill_out)

    if root is None:
        root = os.path.dirname(os.path.dirname(pathA))
    save_matrices([day_data, day_count, day_cmap_n_mod], ["measure", "histogram", "n_mod"], root, resolution)
    plot_mollweid(day_count, star_directions_in_GCS, root, "histogram", str(int(degrees(resolution))), anot=True)
    plot_mollweid(divide(day_data, day_count), star_directions_in_GCS, root, "measure", str(int(degrees(resolution))),
                  anot=True)
    plot_mollweid(divide(day_cmap_n_mod, day_count), star_directions_in_GCS, root, "n_mod",
                  str(int(degrees(resolution))), anot=True)
    # plot_save_imshow(divide(day_cmap_n_mod, day_count), root, "n_mod")
    return day_data, day_count, day_cmap_n_mod


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


def get_same_days_folders(root, needed_files):
    list_subfolders_with_paths = [f.path for f in os.scandir(root) if f.is_dir()]
    for folderA in list_subfolders_with_paths:
        folderA = os.path.join(folderA, "allsatellites")
        if not is_all_data(folderA, needed_files):
            continue
        for folderB in list_subfolders_with_paths:
            folderB = os.path.join(folderB, "allsatellites")
            if not is_all_data(folderB, needed_files):
                continue
            print(os.path.basename(os.path.dirname(folderA))[:4], os.path.basename(os.path.dirname(folderB))[:4])
            if os.path.basename(os.path.dirname(folderA))[-8:] == os.path.basename(os.path.dirname(folderB))[
                                                                  -8:] and os.path.basename(os.path.dirname(folderA))[
                                                                           :4] != os.path.basename(
                os.path.dirname(folderB))[:4]:
                return folderA, folderB
    return None, None


# def process_one_day_from_root(path, needed_files, star_dir, resolution):
#     folderA, folderB = get_same_days_folders(path, needed_files)
#
#     if folderA:
#         print("\n\n Data will be processed in: ", path, "\n\n")
#         process_one_day(folderA, folderB, star_dir, resolution)
#         return True
#     print("\n\n Data not found in: ", path, "\n\n")
#     return False


# def process_many_days_from_root(root_directory, needed_files, star_dir, resolution):
#     list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]
#     for day_path in list_subfolders_with_paths:
#         process_one_day_from_root(day_path, needed_files, star_dir, resolution)


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


def create_dir(root_path, dir_name):
    results_dir = os.path.join(root_path, dir_name)
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    return results_dir


def find_same_days_and_process(path_A, path_B, result_path, needed_files, star_dir, resolution):
    d = 0
    all_hist = []
    all_value = []
    all_n_mod = []
    if os.path.isdir(path_A) and os.path.isdir(path_B) and os.path.isdir(result_path):
        month_pairs = find_corresponding_dirs_in_different_roots(path_A, path_B)
        for A_month, B_month in month_pairs:
            month_name = os.path.split(A_month)[-1]
            print(month_name)
            # condition = True  # month_name in ["augusztus", "november"]
            # if condition:
            # print(month_name)

            day_pairs = find_corresponding_dirs_in_different_roots(A_month, B_month)
            for A_day, B_day in day_pairs:

                date = str(os.path.split(B_day)[-1])[-8:]

                if is_all_data(A_day, needed_files, True) and is_all_data(B_day, needed_files, True):
                    result_month = create_dir(result_path, month_name)
                    result_day = create_dir(result_month, date)
                    print(" Data will be processed from: ", os.path.split(A_day)[-1], "    ", os.path.split(B_day)[-1],
                          "\n", "Index of the process: ", d, "\n")
                    A_day = os.path.join(A_day, "allsatellites")
                    B_day = os.path.join(B_day, "allsatellites")
                    value, hist, n_mod = process_one_day(A_day, B_day, star_dir, resolution, root=result_day,
                                                         fill_out=0.0)
                    all_hist.append(hist)
                    all_value.append(value)
                    all_n_mod.append(n_mod)
                    d += 1
                else:
                    print("\n Data not found for: ", date, "\n")
        all_hist = sum(array(all_hist), axis=0)
        all_value = sum(array(all_value), axis=0)
        all_n_mod = sum(array(all_n_mod), axis=0)
        save_matrices([all_value, all_hist, all_n_mod], ["measure", "histogram", "n_mod"], result_path, resolution)

        all_value_av = divide(all_value, all_hist)
        all_n_mod_av = divide(all_n_mod, all_hist)

        all_hist = rotateAntiClockwise(nan_to_num(all_hist, nan=0.0))
        all_value_av = rotateAntiClockwise(nan_to_num(all_value_av, nan=0.0))
        all_n_mod_av = rotateAntiClockwise(nan_to_num(all_n_mod_av, nan=0.0))

        plot_save_imshow_3_maps([all_hist, all_value_av, all_n_mod_av], ["Histogram", "(|1-r|/|n|)^(1/4)", "<n_mod>"],
                                result_path, resolution="5", logplot=False)
    # plot_save_imshow(all_value_av, result_path, "(|1-r|/|n|)^(1/4)")

    print("Nr of days: ", d)


# =================================================================================================
# =================================================================================================


star_dir = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020"
resolution = radians(5.0)
needed_files = ["user_pos_allsatellites.csv", "all_sats_pos_time.csv"]
# --------------------------------------------NZLD-PERTH--------------------------------------------
place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/PERTH_daily_measurements"
place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_HKKS"

results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_PERTH/r_inv_r"


# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/automatic_processing_no_weight_AminusB"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/automatic_processing_no_weight_AplusB"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/automatic_processing_no_weight_NplusDxyz"

# --------------------------------------------KOREA-Hong-Kong--------------------------------------------
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NASA"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_HKKS"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_NASA/R_Rinv"


# find_same_days_and_process(place_A, place_B, results_root, needed_files, star_dir, resolution)


# ==================================================================================================================================================================================================
# ==================================================================================================================================================================================================
# ==================================================================================================================================================================================================


def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


def get_matrix(path):
    return pd.read_csv(path, skiprows=0).values[:, 1:]  # .transpose()[0]#.astype('float64')


def get_csv_file(directory):
    csv_ext = ".csv"
    # files = next(os.walk(directory))[2]
    files = [f.path for f in os.scandir(directory) if f.is_file()]

    files_with_path = []
    for file in files:
        if os.path.splitext(file)[1] == csv_ext:
            files_with_path.append(os.path.abspath(os.path.join(directory, file)))
    # print(files_with_path)
    return files_with_path


def select_cmap_hist_n_mod(file_list):
    hist = None
    cmap = None
    n_mod = None
    for file in file_list:
        if "histogram_not_averaged" in str(os.path.split(file)[-1]):
            hist = file
        if "measure_not_averaged" in str(os.path.split(file)[-1]) or "measure" in str(os.path.split(file)[-1]):
            cmap = file
        if "n_mod_not_averaged" in str(os.path.split(file)[-1]):
            n_mod = file
    # print(hist)
    # print(n_mod)
    return cmap, hist, n_mod


def plot_mollweid_simple(matrix):
    ra = linspace(-math.pi, math.pi, len(matrix))
    dec = linspace(-math.pi / 2, math.pi / 2, len(matrix[0]))

    X, Y = meshgrid(ra, dec)
    Z = matrix.T
    plt.figure()
    ax = pl.subplot(111, projection='mollweide')
    fig = ax.contourf(X, Y, Z, 100)
    # fig = ax.imshow(rot90(fliplr(matrix), -1))

    plt.xlabel(r'$\theta$', fontsize=15)  # Italic font method
    plt.ylabel(r'$\phi$', fontsize=15)  # Bold font method without fontweight parameters
    pl.colorbar(fig)
    ax.grid()
    # ax.contour(X,Y,Z,10,colors='k')
    pl.show()


def plot_more(path):
    csv_s = get_csv_file(path)
    M = []
    for f in csv_s:
        a = get_matrix(f)
        # print(shape(a))
        M.append(a)

    mm = zeros(shape(M[0]))
    for m in M:
        mm += m
    M = mm / float(len(M))
    # M = mean(M, axis=0)

    plot_mollweid_simple(M)
    M = rotateAntiClockwise(M)
    a = 6
    b = a  # * 2
    M = rebin(M, (int(len(M) / b), int(len(M[0]) / a)))
    plt.imshow(M)
    plt.colorbar()
    plt.show()


# path = r"/Users/kelemensz/Documents/Research/GPS/process/24h/triangular_method/data_GS__24h_cmap.csv"
path = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/atlagolt"
path = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/test/20200803/data_GS_20200803_24h_histogram_5.csv"
path = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/test/20200803/data_GS_20200803_24h_cmap_5.csv"


# plot_more(path)

def read_matrix(path):
    M = rotateAntiClockwise(pd.read_csv(path).values)
    # M = M[M<0.2]
    M = M[:-1, :-1]
    nn = zeros(shape(M))
    # M = nn + M[M<0.2]
    M[M > 0.1] = 0
    # M = log(M+1.0)
    # M = around(M, decimals=0)

    plt.imshow(M)
    plt.colorbar()
    plt.show()


# read_matrix(path)

f_h = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/test/20200803/data_GS_20200803_24h_histogram_5.csv"
f_d = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/test/20200803/data_GS_20200803_24h_cmap_5.csv"


# read_matrix(f_d, f_h)


def read_2d_matrix(file):
    l = loadtxt(file, dtype="str", delimiter=',')
    # with open(file, 'r') as f:
    #     l = [[num for num in line.split(',')] for line in f]
    # l = array(l)
    l = l[1:]
    l = l[:, 1:]
    return rotateAntiClockwise(l.astype("float64"))
    # return l.astype("float64")

def get_matrices_from_paths(paths):
    matrices = []
    for matrix_path in paths:
        # matrix = rotateAntiClockwise(pd.read_csv(matrix_path, index_col=0, dtype='str').astype('float64').values)  # [:, 1:])
        matrix = read_2d_matrix(matrix_path)
        #  # matrix = rotateAntiClockwise(load(matrix_path, allow_pickle=True))
        #  # matrix[matrix < 1] = 0
        matrices.append(matrix)
    return matrices


def plot_save_imshow_2_maps(matrix1, matrix2, root_directory, name, resolution="5", logplot=True):
    print('Matrix.csv saved: ', name, "   ", os.path.split(root_directory)[-1])
    child_dirname = os.path.split(root_directory)[-1] + "_" + name + '_24h' + "_" + resolution
    plt.clf()
    if logplot:
        matrix1 = around(log(matrix1 + 1.0), decimals=0)

    fig, (ax1, ax2) = plt.subplots(2, 1)
    fig.subplots_adjust(left=0.02, bottom=0.1, right=0.95, top=0.94, wspace=0.5, hspace=0.3)
    sp1 = ax1.imshow(matrix1)
    # print(root_directory, "\n", shape(matrix1))
    ax1.set_title("Histogram")
    sp2 = ax2.imshow(matrix2)
    ax2.set_title("r*(r-1)^2")
    fig.colorbar(sp1, ax=ax1)
    fig.colorbar(sp2, ax=ax2)
    # plt.show()
    fig_name = os.path.join(root_directory, child_dirname + '.png')
    fig.savefig(fig_name, bbox_inches='tight')
    plt.clf()


def clean_from_nans(matrix):
    m = []
    for i in matrix:
        mm = []
        for j in i:
            mm.append(j)
        m.append(mm)
        mm = []
    nan_to_num(m, nan=0.0)
    m = array(m)
    m[m > 0.035] = 0
    return m


def calc_correct_average(H, M, new_shape):
    shape = (new_shape[0], H.shape[0] // new_shape[0], new_shape[1], H.shape[1] // new_shape[1])
    H = nan_to_num(H, nan=0.0)
    M = nan_to_num(M, nan=0.0)
    MH = multiply(M, H)
    MH = nan_to_num(MH, nan=0.0)
    H_reshaped_summed = H.reshape(shape).sum(-1).sum(1)
    # plt.imshow(MH)
    # plt.colorbar()
    # # plt.title("<|1-r|/|n|>")
    # plt.show()
    MH_reshaped_summed = MH.reshape(shape).sum(-1).sum(1)
    reshaped = divide(MH_reshaped_summed, H_reshaped_summed)
    # plt.clf()
    # plt.imshow(reshaped)
    # plt.colorbar()
    # plt.title("<|1-r|/|n|>")
    # plt.show()
    return reshaped


# =============================================================================================================================
# =============================================================================================================================
# =============================================================================================================================


def create_averaged_plots_from_root(root_0, months=None):
    sum_all_cmap = []
    sum_all_hist = []
    sum_all_n_mod = []
    subfolders_with_paths_months = [f.path for f in os.scandir(root_0) if f.is_dir()]
    for month_root in subfolders_with_paths_months:
        month_name = str(month_root).split("/")[-1]
        if months and month_name in months:

            days_with_paths = [f.path for f in os.scandir(month_root) if f.is_dir()]
            print("Month name: ", month_name, "  nr days: ", len(days_with_paths))
            for day_root in days_with_paths:
                csv_files = get_csv_file(day_root)
                cmap, hist, n_mod = select_cmap_hist_n_mod(csv_files)
                if (cmap and hist and n_mod) and (
                        os.path.isfile(cmap) and os.path.isfile(hist) and os.path.isfile(n_mod)):
                    M, H, N = get_matrices_from_paths([cmap, hist, n_mod])

                    sum_all_cmap.append(M)
                    sum_all_hist.append(H)
                    sum_all_n_mod.append(N)
    print(hist, "\n", list(H[-1]), "\n", shape(array(H)))
    print("Total number of days:  ", len(sum_all_cmap))
    sum_all_cmap = sum(array(sum_all_cmap), axis=0)
    sum_all_hist = sum(array(sum_all_hist), axis=0)
    sum_all_n_mod = sum(array(sum_all_n_mod), axis=0)
    # plot_save_imshow_3_maps([sum_all_cmap, sum_all_hist, sum_all_n_mod], ["|1-1/r|", "Histogram", "<n_mod>"],
    #                         root_0, resolution="5")
    # plot_save_imshow(sum_all_cmap, root_0, "|1-1/r|")
    # print(sum_all_hist)
    return sum_all_cmap, sum_all_hist, sum_all_n_mod


def handle_raw_not_averaged_matrices(M, H, N):
    ind_no_data = array(H < 1)

    M[ind_no_data] = 0
    N[ind_no_data] = 0.0

    H[H < 1] = 0

    # N[N<0]=0

    # nan_to_num(H, nan=0.0)
    # nan_to_num(M, nan=0.0)
    # nan_to_num(N, nan=0.0)

    M = divide(M, H)
    # N = divide(N, H)
    # N = nan_to_num(N, nan=0.0)
    # M = nan_to_num(M, nan=0)
    # M = absolute(M)
    # M[absolute(M)>0.008] = nan

    # M = clean_from_nans(M)
    # M[M>0.005]=0
    a = 1
    b = a  # * 2
    # M = rebin(M, (int(len(M)/b), int(len(M[0])/a)))
    M = calc_correct_average(H, M, (int(len(M) / b), int(len(M[0]) / a)))

    # M = M * -1
    # M[M < 0] = -1
    # M[M>0] = 1

    # M = nan_to_num(M, nan=0)
    # H = log(H)
    # plot_save_imshow_3_maps([H, M, N], ["Histogram", "(|1-r|/|n|)", "<n_mod>"], root_directory=None, resolution="5", logplot=False, show=True)

    plt.imshow(M)
    plt.colorbar()

    # plot_mollweid_simple(M[::-1].T)
    # plt.title("<|1-r|/|n|>")
    plt.show()


results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_NASA/r_inv_r"

all_months = ["julius", "szeptember", "februar", "marcius", "augusztus", "januar", "december2019", "oktober",
               "november", "majus", "aprilis", "junius", "december2020"]
months1 = ["julius", "szeptember", "augusztus", "november", "junius", "december2020"]
months2 = ["majus", "februar", "marcius", "aprilis", "januar"]


m, h, n = create_averaged_plots_from_root(results_root, all_months)
handle_raw_not_averaged_matrices(m, h, n)


def plot_the_three_raw_matrix_from_path(path):
    csv_files = get_csv_file(path)
    cmap, hist, n_mod = select_cmap_hist_n_mod(csv_files)
    M, H = get_matrices_from_paths([cmap, hist])
    N = rotateAntiClockwise(nan_to_num(pd.read_csv(n_mod, na_filter=True).values[:, 1:], nan=0.0))
    ind_no_data = array(H < 1)

    M[ind_no_data] = nan
    N[ind_no_data] = 0.0

    H[H < 1] = 0

    # # N[N<0]=0

    # nan_to_num(H, nan=0.0)
    # nan_to_num(M, nan=0.0)
    # nan_to_num(N, nan=0.0)

    M = divide(M, H)
    # N = divide(N, H)
    # N = nan_to_num(N, nan=0.0)
    # M = nan_to_num(M, nan=0)
    # M = absolute(M)
    # M[absolute(M)>0.008] = nan

    # M = clean_from_nans(M)
    # M[M>0.005]=0
    a = 1
    b = a  # * 2
    # M = rebin(M, (int(len(M)/b), int(len(M[0])/a)))
    M = calc_correct_average(H, M, (int(len(M) / b), int(len(M[0]) / a)))

    M[M == 0] = nan
    M[M > 0] = 1
    # M[M>0]=2
    # H = log(H)
    # plot_save_imshow_3_maps([H, M, N], ["Histogram", "(|1-r|/|n|)", "<n_mod>"], root_directory=None, resolution="5", logplot=False, show=True)

    plt.imshow(M)
    plt.colorbar()
    # plt.title("<|1-r|/|n|>")
    plt.show()


# M = rotate(M, -90)
# plot_mollweid_simple(M)


# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_PERTH/r_inv_r"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_NZLD/r_inv_r"
# plot_the_three_raw_matrix_from_path(results_root)

# =================================================================================================================================================================
# =================================================================================================================================================================
# =================================================================================================================================================================


# =====================================================================================
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from matplotlib import cm


def random_point(r=1):
    ct = 2 * np.random.rand() - 1
    st = np.sqrt(1 - ct ** 2)
    phi = 2 * np.pi * np.random.rand()
    x = r * st * np.cos(phi)
    y = r * st * np.sin(phi)
    z = r * ct
    return np.array([x, y, z])


def near(p, pntList, d0):
    cnt = 0
    for pj in pntList:
        dist = np.linalg.norm(p - pj)
        if dist < d0:
            cnt += 1 - dist / d0
    return cnt


"""
https://stackoverflow.com/questions/22128909/plotting-the-temperature-distribution-on-a-sphere-with-python
"""


def plot_on_sphere(WW):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    # WW = get_mean_matrix(root_dir)

    u = np.linspace(0, 2 * np.pi, len(WW))
    v = np.linspace(0, np.pi, len(WW[0]))

    # create the sphere surface
    XX = 10 * np.outer(np.cos(u), np.sin(v))
    YY = 10 * np.outer(np.sin(u), np.sin(v))
    ZZ = 10 * np.outer(np.ones(np.size(u)), np.cos(v))

    WW = WW + abs(np.amin(WW))
    myheatmap = WW / np.amax(WW)

    # ~ ax.scatter( *zip( *pointList ), color='#dd00dd' )
    ax.plot_surface(XX, YY, ZZ, cstride=1, rstride=1, facecolors=cm.jet(myheatmap))
    # plt.colorbar(cm.jet( myheatmap ))
    plt.show()


def prepear_for_sphere(m, h):
    ind_no_data = array(h < 1)
    m[ind_no_data] = 0
    h[h < 1] = 0
    nan_to_num(h, nan=0.0)
    nan_to_num(m, nan=0.0)
    M = divide(m, h)
    M = nan_to_num(M, nan=0)
    a = 1
    b = a  # * 2
    M = calc_correct_average(h, M, (int(len(M) / b), int(len(M[0]) / a)))
    M = nan_to_num(M, nan=0)
    M[M < 0] = -1
    M[M > 0] = 1
    return M

M = prepear_for_sphere(m, h)


# plot_on_sphere(M)

# =====================================================================================
