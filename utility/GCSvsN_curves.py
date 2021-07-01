from __future__ import division
import scipy as sci
import scipy.special as sp
from numpy import *
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
from sympy import Point3D, Line3D, Segment3D, Point, Line, Plane, Poly
import pymap3d as pm
import os
from vpython import vector, rotate
import pylab as pl
import math


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


def normvec(vec):
    a = empty(shape(vec))
    for i in range(len(vec)):
        a[i] = unit(vec[i])
    return a


def groupvec(vec, l):
    groups = []
    groups_std = []
    ll = int(len(vec) / l)
    for i in range(l):
        groups.append(vec[i * ll:(i + 1) * ll])
    # groups_std.append(std(vec[i*ll:(i+1)*ll]))
    groups = array(groups)
    # groups_std = array(groups_std)
    return groups


def earth_DS(user_mean, SD):
    proj = []
    for x in SD:
        proj.append(dot(unit(user_mean), x))
    return array(proj)


def projUMEAN_users_mean(groupusers, user_mean, l):
    proj2 = empty(l)
    std_shift_groups = empty(l)
    for i in range(l):
        proj2[i] = dot(mean(groupusers[i], axis=0), unit(user_mean))
        std_shift_groups[i] = mean(std(groupusers[i], axis=0))
    return proj2, std_shift_groups


def scal(v1, v2):
    return vdot(v1, v2)


def unit(v):
    return 1 / sqrt(scal(v, v)) * v


def transform_matrix(f1, f2):  # transforms from f1 to f2
    R = array([
        [dot(f2[0], f1[0]), dot(f2[0], f1[1]), dot(f2[0], f1[2])],
        [dot(f2[1], f1[0]), dot(f2[1], f1[1]), dot(f2[1], f1[2])],
        [dot(f2[2], f1[0]), dot(f2[2], f1[1]), dot(f2[2], f1[2])]
    ])
    return R


def vector_3D_galactic_system(system, theta, phi):
    v_new_in_system = [math.cos(phi) * math.cos(theta),
                       math.cos(phi) * math.sin(theta),
                       math.sin(phi)]
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R2 = transform_matrix(system, ECEF)
    rotated_vector = R2.dot(v_new_in_system)
    return rotated_vector


def get_theta_phi(v):
    if 0.99 < v[0] and 1.01 > v[0] and -0.01 < v[1] and 0.01 > v[1] and -0.01 < v[2] and 0.01 > v[2]:
        return 0.0, 0.0
    if 0.99 < v[1] and 1.01 > v[1] and -0.01 < v[0] and 0.01 > v[0] and -0.01 < v[2] and 0.01 > v[2]:
        return 90.0, 0.0
    if 0.99 < v[2] and 1.01 > v[2] and -0.01 < v[1] and 0.01 > v[1] and -0.01 < v[0] and 0.01 > v[0]:
        return 0.0, 90.0
    return None, None


def cartesian_to_galactic(system, vector):
    vector = unit(vector)
    system = normvec(system)
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R = transform_matrix(ECEF, system)
    vector_in_system = around(R.dot(vector), decimals=3)
    # print(vector_in_system)
    theta, phi = get_theta_phi(vector_in_system)
    if not theta and not phi:
        phi, theta, _ = pm.ecef2geodetic(vector_in_system[0], vector_in_system[1], vector_in_system[2])
        # print(theta, phi)
        phi = 90 - degrees(arccos(scal(ECEF[2], vector_in_system)))
    return theta, phi


# return degrees(theta), degrees(phi)


def get_mean_direction_over_time(systems, directions):
    # systems = systems.T
    # print(systems)
    l = min(len(directions), len(systems[0]))
    phi = 0.
    theta = 0.
    # print('lenght', l)
    for i in range(l):
        out = cartesian_to_galactic(array([systems[0][i], systems[1][i], systems[2][i]]), directions[i])
        # print('out:', out)
        theta += out[0]
        phi += out[1]
    return theta / float(l), phi / float(l)


def direction_3D_galactic_system(Nx, Ny, Nz, theta, phi):
    galactic_direction = empty(shape(Nx))
    l = min(len(Nx), len(Ny), len(Nz))
    for i in range(l):
        system = array([Nx[i], Ny[i], Nz[i]])
        galactic_direction[i] = vector_3D_galactic_system(system, theta, phi)
    return galactic_direction


def plot_save_curves(user_mean, SD, theta, phi, dest_dir):
    l = len(SD)
    d1 = earth_DS(user_mean, SD)  # SD-n angles
    # minange = max(abs(d1))
    # d2, std_d2 = projUMEAN_users_mean(usergrupped, user_mean, l)
    # # print('Mean projections: ', mean(absolute(d2)))
    # # print('Mean projections: ', mean(d2))
    d1 = array(d1).astype('float64')
    # d2 = array(d2).astype('float64')
    # diff = p_m_byv_limits(d1,d2)
    plt.clf()
    plt.plot(d1)
    plt.plot(zeros(shape(d1)))
    plt.title("theta: {} - phi: {}".format(around(degrees(theta), 3), around(degrees(phi), 3)))

    fig_name1 = os.path.join(dest_dir, "{}_{}.png".format(around(degrees(theta), 3), around(degrees(phi), 3)))

    plt.savefig(fig_name1, bbox_inches='tight')


# plt.show()

# =================================================================================================
# =================================================================================================

def get_third_dir_by_cross_product(A, B):
    l = min(len(A), len(B))
    C = empty((l, 3))
    for i in range(l):
        C[i] = cross(A[i], B[i])
    return normvec(C)


def prepare_data_for_CMBR_cmap_global_stars(star_dir, directory):
    star_dirs = [f.path for f in os.scandir(star_dir) if f.is_dir()]
    for stardir in star_dirs:
        if directory[-8:] == stardir[-8:]:
            print('Star data: ', stardir)
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


def prepare_data_for_CMBR_cmap(path):
    user_ = pd.read_csv(path + '/user_pos_allsatellites.csv', skiprows=1).values  # .transpose()
    user_mean = mean(user_, axis=0).astype('float64')
    shifts = user_ - user_mean

    try:
        sun, galactic_n_p, galactic_center, nunki, galactic_anti_center, S, D = getSTARS_for_galactic_system(path)
        sun = normvec(sun)
        galactic_n_p = normvec(galactic_n_p)
        galactic_center = normvec(galactic_center)
        galactic_anti_center = normvec(galactic_anti_center)
        nunki = normvec(nunki)
        Nx = normvec(galactic_center).astype('float64')
        Nz = normvec(galactic_n_p).astype('float64')
        Ny = -get_third_dir_by_cross_product(Nz, Nx)
        return user_mean, shifts, Nx, Ny, Nz, S, D
    except:
        print("No local star data!")
        return user_mean, shifts, None, None, None, None, None


def calculate_using_directions_in_galactic_coor(user_mean, usergrupped, Nx, Ny, Nz, theta, phi, dest_dir):
    rDIRECTION = direction_3D_galactic_system(Nx, Ny, Nz, theta, phi)
    plot_save_curves(user_mean, rDIRECTION, theta, phi, dest_dir)


def vectors_dot_prods(v1, v2):
    v1 = normvec(v1)
    v2 = normvec(v2)
    dot_v1_v2 = []
    for i in range(len(v1)):
        dot_v1_v2.append(dot(v1[i], v2[i]))
    print('Local vs. Global directions: ', degrees(arccos(mean(dot_v1_v2))))
    return dot_v1_v2


def create_dir(root, dirname):
    new_dir = os.path.join(root, dirname)
    if not os.path.isdir(new_dir):
        os.makedirs(new_dir)
    return new_dir


def process_one_day(directory, star_dir, resolution, dest_dir):
    dest_dir = create_dir(dest_dir, str(directory).split("/")[-1])
    print("Dest dir: ", dest_dir)
    theta_max = math.pi  # * 2.0 #+ resolution * 0.1
    phi_max = math.pi / 2.0
    rot_theta = arange(-theta_max, theta_max, resolution)
    rot_phi = arange(-phi_max, phi_max, resolution)

    total_dir = os.path.join(directory, "allsatellites")
    print("Positional data: ", total_dir)
    user_mean, shifts, Nx_local, Ny_local, Nz_local, S_local, D_local = prepare_data_for_CMBR_cmap(total_dir)

    Nx_global, Ny_global, Nz_global, S_global, D_global = prepare_data_for_CMBR_cmap_global_stars(star_dir, directory)
    # =================================================================================================================
    # vectors_dot_prods(Nx_local, Nx_global)
    # vectors_dot_prods(Ny_local, Ny_global)
    # vectors_dot_prods(Nz_local, Nz_global)
    # vectors_dot_prods(S_local, S_global)
    # vectors_dot_prods(D_local, D_global)
    # print("-", get_mean_direction_over_time(array([Nx_global, Ny_global, Nz_global]), Nz_global))
    star_directions_in_GCS = {}
    star_directions_in_GCS['GNP'] = get_mean_direction_over_time(array([Nx_global, Ny_global, Nz_global]), Nz_global)
    star_directions_in_GCS['GC'] = get_mean_direction_over_time(array([Nx_global, Ny_global, Nz_global]), Nx_global)
    star_directions_in_GCS['OY of GCS'] = get_mean_direction_over_time(array([Nx_global, Ny_global, Nz_global]),
                                                                       Ny_global)
    star_directions_in_GCS['Sadalmelic'] = get_mean_direction_over_time(array([Nx_global, Ny_global, Nz_global]),
                                                                        S_global)
    star_directions_in_GCS['Denebola'] = get_mean_direction_over_time(array([Nx_global, Ny_global, Nz_global]),
                                                                      D_global)
    # =================================================================================================================
    print('STD of positions: ', std(shifts, axis=0))
    usergrupped = groupvec(shifts, len(Nx_global))
    del shifts
    cmap = empty((len(rot_theta), len(rot_phi)))
    for j in range(len(rot_phi)):
        # print('Phi index: {}/{}'.format(j+1, len(rot_phi)))
        for i in range(len(rot_theta)):
            theta = rot_theta[i]
            phi = rot_phi[j]
            cmap[i][j] = calculate_using_directions_in_galactic_coor(user_mean, usergrupped, Nx_global, Ny_global,
                                                                     Nz_global, theta, phi, dest_dir)


# =================================================================================================
# =================================================================================================

# root_directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/process_NASA/januar"

# root_directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/december2019"

star_dir = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020"
nr = radians(5.0)

directory = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements/januar/CUTB20200101"

dest_dir = r"/Users/kelemensz/Documents/Research/GPS/process/show_user_inGCS/CUTB"
process_one_day(directory, star_dir, nr, dest_dir)
