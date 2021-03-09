import time

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
from scipy.linalg import norm as nrm
import copy as cp

s_o_l = 1.0  # 3.0*10**8
satellite_positions0 = "all_sats_pos_time.csv"
satellite_positions = "sats_pos_time_id.csv"


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


def parse_sats_data_beta(sat_data_file):
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


def parse_sats_data_alpha(sat_data_file):
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


def parse_sats_data(sat_data_file):
    if "_id" in str(sat_data_file).split("/")[-1]:
        return parse_sats_data_alpha(sat_data_file)
    return parse_sats_data_beta(sat_data_file)


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
                break
    return data


def get_comon_alpha(vA, vB):
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
    sat_data_file = os.path.join(path, satellite_positions)
    if not os.path.isfile(sat_data_file):
        sat_data_file = os.path.join(path, satellite_positions0)
    sat_data = parse_sats_data(sat_data_file)
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
    posA, dataA = prepare_u_and_s_positions(pathA)
    posB, dataB = prepare_u_and_s_positions(pathB)
    comonAB = extract_common_sats(dataA, dataB)
    raw_results = process_all(posA, posB, comonAB)
    return raw_results


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
    plt.close('all')
    plt.clf()
    pl.clf()


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


def process_one_day_symmetrized(pathA, pathB, star_dir, resolution, mean_positions=None, root=None,
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

    # -------------------------------Symmetrize the calculated measure-------------------------------------------------
    raw_results_ECEF_AB = get_raw_results_using_mean_positions(pathA, pathB, mean_positions)
    groupped_raw_results_ECEF_AB = group_results(raw_results_ECEF_AB, l)
    groupped_raw_results_GCS_AB = raw_results_to_GCS(groupped_raw_results_ECEF_AB, GCS_all)
    raw_results_GCS_AB = list(chain(*groupped_raw_results_GCS_AB))  # [::10000]

    raw_results_ECEF_BA = get_raw_results_using_mean_positions(pathB, pathA, mean_positions[::-1])
    groupped_raw_results_ECEF_BA = group_results(raw_results_ECEF_BA, l)
    groupped_raw_results_GCS_BA = raw_results_to_GCS(groupped_raw_results_ECEF_BA, GCS_all)
    raw_results_GCS_BA = list(chain(*groupped_raw_results_GCS_BA))  # [::10000]
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

    if root is None:
        print("No day_result directory is given...process stops here!")
        sys.exit()

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


def get_mean_pos_from_root(root_path, positions_file, max_deviations=0.5, nofilter = False):
    positions = []
    month_folders = [f.path for f in os.scandir(root_path) if f.is_dir()]
    for month_root in month_folders:
        day_folders = [f.path for f in os.scandir(month_root) if f.is_dir()]
        for day_root in day_folders:
            day_root = os.path.join(day_root, "allsatellites")
            if is_all_data(day_root, [positions_file], add_allsatellites=False):
                mean_pos = get_mean_position(day_root, positions_file)
                if nofilter:
                    print("Nofilter true: ", mean_pos)
                    return mean_pos
                # if str(os.path.split(day_root)[-2]).split("/")[-1][:4] in ['NZLD', 'NZDL']:
                if str(os.path.split(day_root)[-2]).split("/")[-1][:4] not in ['BLUF', 'HOKI', 'MAVL', 'LKTA', 'MTJO']:
                    # print(str(os.path.split(day_root)[-2]).split("/")[-1], mean_pos)
                    positions.append(mean_pos)
    positions = array(positions)
    std_pos = std(positions, axis=0)
    print("Number of days with positions determined and std of the positions before filter: ", len(positions), std_pos)

    std_pos_norm = sqrt(std_pos.dot(std_pos))
    mean_ = mean(positions, axis=0)
    distance_from_mean = nrm(positions - mean_, axis=1)
    not_outlier = distance_from_mean < max_deviations * std_pos_norm
    # print(not_outlier)
    no_outliers = positions[not_outlier]
    print("After filter: ", len(no_outliers), std(array(no_outliers), axis=0), "\n ")
    print(mean(array(positions), axis=0), mean(array(no_outliers), axis=0))
    return mean(array(no_outliers), axis=0)


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
        mean_pos_A = get_mean_pos_from_root(path_A, needed_files[0], max_deviations=5)
        mean_pos_B = get_mean_pos_from_root(path_B, needed_files[0], max_deviations=0.2, nofilter=True)  # NZLD eseten 0.2
        for A_month, B_month in month_pairs:
            month_name = os.path.split(A_month)[-1]
            condition = month_name in ["marcius", "aprilis", "majus", "junius", "februar"]
            if condition:
                print(month_name)
                day_pairs = find_corresponding_dirs_in_different_roots(A_month, B_month)
                print("Number of days: ", len(day_pairs))
                for A_day, B_day in day_pairs[:5]:
                    start = time.time()
                    date = str(os.path.split(B_day)[-1])[-8:]
                    cond1 = is_all_data(A_day, needed_files[1:], True) or is_all_data(A_day, [satellite_positions0], True)
                    cond2 = is_all_data(B_day, needed_files[1:], True) or is_all_data(B_day, [satellite_positions0], True)
                    if cond1 and cond2:
                        result_month = create_dir(result_path, month_name)
                        result_day = create_dir(result_month, date)
                        print(" Data will be processed from: ", os.path.split(A_day)[-1], "    ",
                              os.path.split(B_day)[-1],
                              "\n", "Index of the process: ", d, "\n")
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

                        # if is_all_data(A_day, needed_files[:1]) and is_all_data(B_day, needed_files[:1]):
                        #     value, hist, n_mod = process_one_day_symmetrized(A_day, B_day, star_dir, resolution,
                        #                                                      root=result_day, fill_out=0.0)
                        # else:
                        #     # print("Mean positoin will be considered! ", str(B_day).split("/")[-1])
                        value, hist, n_mod = process_one_day_symmetrized(A_day, B_day, star_dir, resolution,
                                                                         mean_positions=[posUA, posUB],
                                                                         root=result_day, fill_out=0.0)
                        all_hist.append(hist)
                        all_value.append(value)
                        all_n_mod.append(n_mod)
                        d += 1
                        print('Elapsed time of the current day: ', time.time() - start, date)
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
needed_files = ["user_pos_allsatellites.csv", satellite_positions]

# ======================================================================================================================

# --------------------------------------------PERTH-Hong-Kong--------------------------------------------
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/PERTH_daily_measurements"
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_HKKS"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_PERTH/r_inv_r_symmetrized"

# --------------------------------------------NZLD-Hong-Kong-------------------------------------------- [*******]
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NZLD"
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_HKKS"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/NZLD_HKKS/r_inv_r_AB"


# --------------------------------------------PERTH-Del-korea-symmetrized--------------------------------------------
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/PERTH_daily_measurements"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NASA"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_NASA/r_inv_r_symmetrized"

# --------------------------------------------PERTH-India-symmetrized--------------------------------------------
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/PERTH_daily_measurements"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_IIGC"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_IIGC/r_inv_r_symmetrized"

# --------------------------------------------NZLD-India-symmetrized--------------------------------------------
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NZLD"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_IIGC"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/NZLD_IIGC/r_inv_r_symmetrized"


# --------------------------------------------NZLD-Del-korea-symmetrized--------------------------------------------
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NZLD"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NASA"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/NZLD_NASA/r_inv_r_symmetrized"

# --------------------------------------------NZLD-PERTH-symmetrized-------------------------------------------
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/PERTH_daily_measurements"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NZLD"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_NZLD/r_inv_r_symmetrized"

# --------------------------------------------KOREA-Hong-Kong--------------------------------------------
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NASA"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_HKKS"
# # results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_NASA/r_inv_r"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_NASA/r_inv_r_over_Nmod_symmetrized"


# # --------------------------------------------KOREA-India--------------------------------------------
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NASA"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_IIGC"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/NASA_IIGC/r_inv_r_symmetrized"


# --------------------------------------------Hong-Kong-India--------------------------------------------
# place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_HKKS"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_IIGC"
# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/pairs_by_identifier/HKKS_IIGC/r_inv_r_symmetrized"


# --------------------------------------------PERTH-TIDV-symmetrized-------------------------------------------
# place_A = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_TIDV"
# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/PERTH_TIDV/r_inv_r_symmetrized"


# --------------------------------------------NZLD-TIDV-symmetrized-------------------------------------------
# place_A = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_NZLD"
# place_B = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_TIDV"
# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/NZLD_TIDV/r_inv_r_symmetrized"

# # # --------------------------------------------PERTH CUTA-NZLD--------------------------------------------
# place_A = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements/CUTA"
# place_B = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_NZLD"
# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/CUTA_NZLD/r_inv_r_symmetrized"



# --------------------------------------------NASA-India------test--------------------------------------
place_A = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_NASA"
place_B = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_IIGC"
results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_test/r_inv_r_symmetrized"


find_same_days_and_process(place_A, place_B, results_root, needed_files, star_dir, resolution)
