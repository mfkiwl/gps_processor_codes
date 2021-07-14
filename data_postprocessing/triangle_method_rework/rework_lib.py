import time

from numpy import *
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from itertools import chain
import pymap3d as pm
import os
import pylab as pl
import math
from scipy.ndimage import rotate

import sys

from tqdm import tqdm

sys.path.insert(1, '/Users/kelemensz/Documents/Research/GPS/gps_processor_codes/utility')
from frecvently_used_functions import is_all_data
# from utility.frecvently_used_functions import is_all_data


class Defaults:
    SAT_POS_FILENAMES = {'no_ID': "all_sats_pos_time.csv", 'ID': "sats_pos_time_id.csv"}
    USER_POSITIONS_FILENAME = {'user': "user_pos_allsatellites.csv"}

    STAR_POS_FILES = {
        'denebola_old': 'star_d_positions.csv',
        'denebola': 'Denebola_positions.csv',
        'sadalmelik_old': 'star_s_positions.csv',
        'sadalmelik': 'Sadalmelik_positions.csv',
        'sun': 'SUN_positions.csv',
        'galactic_nord_pole': 'GNP_positions.csv',
        'galactic_center': 'GC_positions.csv',
        'mirphak': 'star_Mirphak_positions.csv',
        'vega': 'star_Vega_positions.csv',
        'nunki_old': 'star_Nunki_positions.csv',
        'nunki': 'Nunki_positions.csv',
        'capella_old': 'star_Capella_positions.csv',
        'capella': 'Capella_positions.csv',
    }


    triangle_result_names = ["measure", "histogram", "n_mod"]
    triangle_result_names_notsym = ["not_symmetrical_measure", "not_symmetrical_histogram", "not_symmetrical_n_mod"]
    triangle_result_names_n_filter_notsym = ["not_symmetrical_angle_filtered_measure", "not_symmetrical_angle_filtered_histogram", "not_symmetrical_angle_filtered_n_mod"]
    triangle_result_names_n_filter_sym = ["angle_filtered_measure", "angle_filtered_histogram", "angle_filtered_n_mod"]

    @classmethod
    def get_star_directions_in_GCS(cls, Nx, Ny, Nz, D, S):
        return {
            'GNP': get_mean_direction_in_system_over_time(array([Nx, Ny, Nz]), Nz),
            'GC': get_mean_direction_in_system_over_time(array([Nx, Ny, Nz]), Nx),
            'OY of GCS': get_mean_direction_in_system_over_time(array([Nx, Ny, Nz]), Ny),
            'Sadalmelic': get_mean_direction_in_system_over_time(array([Nx, Ny, Nz]), S),
            'Denebola': get_mean_direction_in_system_over_time(array([Nx, Ny, Nz]), D)}

    @classmethod
    def get_star_directions_in_GCS_simple(cls, Nx, Ny, Nz, D, S):
        return {
            'GNP': Nz[0],
            'GC': Nx[0],
            'OY of GCS': Ny[0],
            'Sadalmelic': S[0],
            'Denebola': D[0]}


def rotate2darray(array2d, rot_angle=90):
    """90 = AntiClockwise, -90 = Clockwise"""
    return rotate(array2d, rot_angle)


def scal(v1, v2):
    return vdot(array(v1), array(v2))


def unit(v):
    return 1 / sqrt(scal(v, v)) * array(v)


def normvec(vec):
    a = empty(shape(vec))
    for i in range(len(vec)):
        a[i] = unit(vec[i])
    return a


def angle_between(v1, v2):
    unit_vector_1 = unit(v1)
    unit_vector_2 = unit(v2)
    dot_product = np.dot(unit_vector_1, unit_vector_2)
    return degrees(np.arccos(dot_product))


def verify_angle(v1_list, v2_list):
    if len(v1_list) == len(v2_list):
        for i in range(len(v2_list)):
           print(angle_between(v1_list[i], v2_list[i]))


def getSTARS_for_galactic_system(path):
    D = normvec(pd.read_csv(os.path.join(path, Defaults.STAR_POS_FILES.get('denebola')), skiprows=0).values)
    S = normvec(pd.read_csv(os.path.join(path, Defaults.STAR_POS_FILES.get('sadalmelik')), skiprows=0).values)
    sun = pd.read_csv(os.path.join(path, Defaults.STAR_POS_FILES.get('sun')), skiprows=0).values
    galactic_n_p = pd.read_csv(os.path.join(path, Defaults.STAR_POS_FILES.get('galactic_nord_pole')), skiprows=0).values
    galactic_center = pd.read_csv(os.path.join(path, Defaults.STAR_POS_FILES.get('galactic_center')), skiprows=0).values
    nunki = pd.read_csv(os.path.join(path, Defaults.STAR_POS_FILES.get('nunki')), skiprows=0).values
    capella = pd.read_csv(os.path.join(path, Defaults.STAR_POS_FILES.get('capella')), skiprows=0).values
    return sun, galactic_n_p, galactic_center, nunki, capella, S, D


def get_mean_direction_in_system_over_time(systems, directions):
    # l = min(len(directions), len(systems[0]))
    ld = len(directions)
    ls = len(systems[0])
    if ld == ls:
        phi = 0.
        theta = 0.
        for i in range(ld):
            out = ECEFcartesian_to_GCSspherical(array([systems[0][i], systems[1][i], systems[2][i]]), directions[i])
            theta += out[0]
            phi += out[1]
        return theta / float(ld), phi / float(ld)
    return None


def get_star_directions_in_GCS(Nx, Ny, Nz, S, D):
    return {'GNP': get_mean_direction_in_system_over_time(array([Nx, Ny, Nz]), Nz),
            'GC': get_mean_direction_in_system_over_time(array([Nx, Ny, Nz]), Nx),
            'OY of GCS': get_mean_direction_in_system_over_time(array([Nx, Ny, Nz]), Ny),
            'Sadalmelic': get_mean_direction_in_system_over_time(array([Nx, Ny, Nz]), S),
            'Denebola': get_mean_direction_in_system_over_time(array([Nx, Ny, Nz]), D)}


def consecutively_group_results(results_dict, l):
    """Creates a 2D data container of length l, out of the result_dict dictionary by grouping consecutive entries."""
    groupped = []
    group = []
    n = len(list(results_dict.keys()))
    l_groups = int(n / l)
    i = 0
    if l_groups > 1:
        for v in results_dict.values():
            if i >= l_groups:
                groupped.append(list(chain(*group)))
                group = [v]
                i = 0
                continue
            group.append(v)
            i += 1
        groupped.append(list(chain(*group)))
    return groupped


def v_length(v1, v2):
    return sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def parse_sats_data_no_ID(sat_data_file):
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


def parse_sats_data_with_ID(sat_data_file):
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
        return parse_sats_data_with_ID(sat_data_file)
    return parse_sats_data_no_ID(sat_data_file)


def get_common_no_ID(vA, vB):
    data = []
    for a in vA:
        for b in vB:
            if v_length(a[:3], b[:3]) < 500:
                data.append([a, b])
                break
    return data


def get_common_with_ID(vA, vB):
    data = []
    for a in vA:
        for b in vB:
            if a[4] == b[4]:
                data.append([a, b])
                break
    return data


def get_common(vA, vB, is_id=True):
    if is_id:
        return get_common_with_ID(vA, vB)
    return get_common_no_ID(vA, vB)


def get_mean_user_position(root_path, both=False):
    user_ = pd.read_csv(os.path.join(root_path, Defaults.USER_POSITIONS_FILENAME.get('user')), skiprows=1).values
    if both:
        return mean(user_, axis=0).astype('float64'), user_
    return mean(user_, axis=0).astype('float64')


def transform_matrix(f1, f2):  # transforms from f1 to f2
    R = array([
        [dot(f2[0], f1[0]), dot(f2[0], f1[1]), dot(f2[0], f1[2])],
        [dot(f2[1], f1[0]), dot(f2[1], f1[1]), dot(f2[1], f1[2])],
        [dot(f2[2], f1[0]), dot(f2[2], f1[1]), dot(f2[2], f1[2])]
    ])
    return R


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


def ecef_to_gcs_(system, vector):
    vector = unit(vector)
    system = normvec(system)
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R = transform_matrix(ECEF, system)
    # return around(R.dot(vector), decimals=3)
    return R.dot(vector)


def ecef_to_gcs(system, vector):
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R = transform_matrix(ECEF, system)
    # return around(R.dot(vector), decimals=3)
    return R.dot(vector)


def cart2sph(x, y, z):
    hxy = hypot(x, y)
    # r = hypot(hxy, z)
    el = arctan2(z, hxy)
    az = arctan2(y, x)
    return az, el


def raw_results_to_GCS(results_ECEF, GCS):
    results_GCS = []
    nr_groups = len(results_ECEF)
    nr_axis = len(GCS)
    # print("Systems and groups: ", nr_axis, nr_groups)
    if nr_axis == nr_groups:
        for i in tqdm(range(nr_groups)):
            results_group = results_ECEF[i]
            GCS_group = normvec(GCS[i])
            results_group_GCS = []
            for direction, value, mod_n in results_group:
                results_group_GCS.append([ecef_to_gcs(GCS_group, direction), value, mod_n])
            results_GCS.append(results_group_GCS)
        return results_GCS
    return None


def ECEFcartesian_to_GCSspherical(system, vector):
    vector = unit(vector)
    system = normvec(system)
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R = transform_matrix(ECEF, system)
    vector_in_system = around(R.dot(vector), decimals=3)
    # vector_in_system = R.dot(vector)
    # theta, phi = get_theta_phi(vector_in_system)
    theta, phi = cart2sph(vector_in_system[0], vector_in_system[1], vector_in_system[2])
    # print(vector, vector_in_system, theta, phi)
    if theta is None and phi is None:
        # print('Problem with GCS direction calculation...')
        phi, theta, _ = pm.ecef2geodetic(vector_in_system[0], vector_in_system[1], vector_in_system[2])
        phi = 90 - degrees(arccos(scal(ECEF[2], vector_in_system)))
    return theta, phi


def get_third_dir_by_cross_product(A, B):
    # l = min(len(A), len(B))
    lA = min(len(A), len(B))
    lB = min(len(A), len(B))
    if lA == lB:
        C = empty((lA, 3))
        for i in range(lA):
            C[i] = cross(A[i], B[i])
        return normvec(C)
    return None


def get_global_stars(star_dir, directory):
    star_dirs = [f.path for f in os.scandir(star_dir) if f.is_dir()]
    for stardir in star_dirs:
        if directory[-8:] == stardir[-8:]:
            return get_stars_for_CMBR_cmap(stardir)
    return 0, 0, 0, 0, 0


def get_stars_for_CMBR_cmap(path):
    sun, galactic_n_p, galactic_center, nunki, galactic_anti_center, S, D = getSTARS_for_galactic_system(path)

    # sun = normvec(sun)
    S = normvec(S)
    D = normvec(D)
    galactic_n_p = normvec(galactic_n_p)
    galactic_center = normvec(galactic_center)
    # galactic_anti_center = normvec(galactic_anti_center)
    # nunki = normvec(nunki)
    Nx = normvec(galactic_center).astype('float64')
    Nz = normvec(galactic_n_p).astype('float64')
    Ny = -get_third_dir_by_cross_product(Nz, Nx)
    # print('X-Y:')
    # verify_angle(Nx, galactic_n_p)
    # print('X-Z:')
    # verify_angle(galactic_n_p, Nz)
    # print('Y-Z:')
    # verify_angle(Ny, galactic_n_p)

    return Nx, Ny, Nz, S, D


def create_dir(root_path, dir_name):
    results_dir = os.path.join(root_path, dir_name)
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    return results_dir


def get_user_and_sat_positions(path, get_u=True):
    sat_data_file = os.path.join(path, Defaults.SAT_POS_FILENAMES.get('ID'))
    if not os.path.isfile(sat_data_file):
        sat_data_file = os.path.join(path, Defaults.SAT_POS_FILENAMES.get('no_ID'))
    sat_data = parse_sats_data(sat_data_file)
    if get_u:
        return get_mean_user_position(path), sat_data
    return sat_data


def get_ij_on_map_(direction, resolution):
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


def get_ij_on_map(direction, l_theta, l_phi, resolution):
    theta, phi = cart2sph(direction[0], direction[1], direction[2])
    I_f = 0
    J_f = 0
    for i in range(-l_theta, l_theta):
        if i * resolution > theta:
            I_f = i + l_theta
            break
    for i in range(-l_phi, l_phi):
        if i * resolution > phi:
            J_f = i + l_phi
            break
    return I_f, J_f


#  =====================================================================================================================
#  =====================================================================================================================
#  =====================================================================================================================

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
    # theta = radians(theta)
    # phi = radians(phi)
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
        cmap_save = pd.DataFrame(np.array(matrix))
        cmap_save.to_csv(os.path.join(directory, name + "_not_averaged_" + str(int(degrees(resolution))) + '.csv'),
                         index=True)
    print('Matrices saved! :', names)


def input_fileset_of_the_day(day_path, add_allsatellites):
    user_positions = is_all_data(day_path, [Defaults.USER_POSITIONS_FILENAME.get('user')], add_allsatellites)
    if user_positions:
        if is_all_data(day_path, [Defaults.SAT_POS_FILENAMES.get('ID')], add_allsatellites):
            return 2
        if is_all_data(day_path, [Defaults.SAT_POS_FILENAMES.get('no_ID')], add_allsatellites):
            return 1
    return 0


def is_string_in_filenames(string, filenames):
    for filename in filenames:
        if string in str(filename):
            return True
    return False


def filter_sats_within_solid_angle(sats_dat, surface_normal, angular_distance):
    close_sats = []
    for sat in sats_dat:
        if angle_between(sat[:3], surface_normal) < angular_distance/2.0:
            close_sats.append(sat)
    # print(len(sats_dat), '  -  ', len(close_sats))
    return close_sats


def filter_sats_within_solid_angle_above_the_sky(sat_data, surface_normal, angular_distance):
    filtered = {}
    for epoch, sats_dat in tqdm(sat_data.items()):  # k = epoch index
        filtered[epoch] = filter_sats_within_solid_angle(sats_dat, surface_normal, angular_distance)
    return filtered
