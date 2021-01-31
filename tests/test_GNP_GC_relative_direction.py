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


def getSTARS_for_galactic_system_alpha(path):
    D = normvec(pd.read_csv(path + '/star_d_positions.csv', skiprows=0).values)
    S = normvec(pd.read_csv(path + '/star_s_positions.csv', skiprows=0).values)

    sun = pd.read_csv(path + '/SUN_positions.csv', skiprows=0).values
    galactic_n_p = pd.read_csv(path + '/GNP_positions.csv', skiprows=0).values
    # galactic_center = pd.read_csv(path+'/GC_positions.csv',skiprows=0).values

    # star_1 = pd.read_csv(path+'/star_Mirphak_positions.csv',skiprows=0).values
    # star_2 = pd.read_csv(path+'/star_Vega_positions.csv',skiprows=0).values

    nunki = pd.read_csv(path + '/star_Nunki_positions.csv', skiprows=0).values
    capella = pd.read_csv(path + '/star_Capella_positions.csv', skiprows=0).values
    # print(galactic_center)
    return sun, galactic_n_p, nunki, capella, S, D


def getSTARS_for_galactic_system(path):
    print('±±±±±±±±±±±±±±±±±±±±±±±')
    S = normvec(pd.read_csv(path + '/Denebola_positions.csv', skiprows=0).values)
    D = normvec(pd.read_csv(path + '/Sadalmelik_positions.csv', skiprows=0).values)

    sun = pd.read_csv(path + '/SUN_positions.csv', skiprows=0).values
    galactic_n_p = pd.read_csv(path + '/GNP_positions.csv', skiprows=0).values
    galactic_center = pd.read_csv(path + '/GC_positions.csv', skiprows=0).values

    # star_1 = pd.read_csv(path+'/star_Mirphak_positions.csv',skiprows=0).values
    # star_2 = pd.read_csv(path+'/star_Vega_positions.csv',skiprows=0).values

    nunki = pd.read_csv(path + '/Nunki_positions.csv', skiprows=0).values
    capella = pd.read_csv(path + '/Capella_positions.csv', skiprows=0).values
    # print(galactic_center)
    return sun, galactic_n_p, galactic_center, nunki, capella, S, D


def normvec(vec):
    a = empty(shape(vec))
    for i in range(len(vec)):
        a[i] = unit(vec[i])
    return a


def scal(v1, v2):
    return vdot(v1, v2)


def unit(v):
    return 1 / sqrt(scal(v, v)) * v


def gram_schmidt(L1):
    L2 = [unit(L1[0])]
    for i in range(1, len(L1)):
        s = L1[i]
        for j in range(i):
            s = s - scal(L1[i], L2[j]) * L2[j]
        L2.append(unit(s))
    return array(L2)


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


def cartesian_to_galactic(system, vector):
    s_phi = round(dot(vector, system[2]), 4)
    s_theta = round(dot(vector, system[0]), 4) / (sqrt(1 - s_phi ** 2))
    return round(s_theta, 4), round(s_phi, 4)


def direction_3D_galactic_system(Nx, Ny, Nz, SD, theta, phi):
    galactic_direction = empty(shape(SD))
    l = min(len(Nx), len(Ny), len(Nz), len(SD))
    SD_lat_lon = empty((l, 2))
    # print('Effectuated star position updates: ', l)
    for i in range(l):
        system = gram_schmidt(array([Nx[i], Nz[i], Ny[i]]))
        system = normvec(array([system[0], system[2], system[1]]))
        galactic_direction[i] = vector_3D_galactic_system(system, theta, phi)
        SD_lat_lon[i] = cartesian_to_galactic(system, SD[i])
    return galactic_direction, SD_lat_lon


# =================================================================================================
# =================================================================================================
def get_third_dir_by_cross_product(X, Z):
    l = min(len(X), len(Z))
    Y = empty((l, 3))
    for i in range(l):
        Y[i] = unit(cross(X[i], Z[i]))
    return Y


def prepare_data_for_CMBR_cmap(path):
    user_ = pd.read_csv(path + '/user_pos_allsatellites.csv', skiprows=1).values  # .transpose()
    user_mean = mean(user_, axis=0).astype('float64')
    shifts = user_ - user_mean
    sun, galactic_n_p, galactic_center, nunki, galactic_anti_center, S, D = getSTARS_for_galactic_system(path)
    S = normvec(S)
    D = normvec(D)
    SD = normvec(S - D).astype('float64')

    sun = normvec(sun)
    galactic_n_p = normvec(galactic_n_p)
    galactic_center = normvec(galactic_center)
    galactic_anti_center = normvec(galactic_anti_center)
    nunki = normvec(nunki)

    # Nx = normvec(galactic_anti_center - galactic_center).astype('float64')
    # Nz = normvec(galactic_n_p - sun).astype('float64')

    Nx = normvec(galactic_center).astype('float64')
    Nz = normvec(galactic_n_p).astype('float64')

    Ny = get_third_dir_by_cross_product(Nx, Nz)

    # Nx = normvec(- sun + galactic_center).astype('float64')
    # Nz = normvec(galactic_n_p - sun).astype('float64')

    for i in range(len(SD) - 1):
        print(degrees(arccos(dot(Nx[i], Nz[i]))))

    return user_mean, shifts, Nx, Ny, Nz, SD


def prepare_data_for_CMBR_cmap_alpha(path):
    user_ = pd.read_csv(path + '/user_pos_allsatellites.csv', skiprows=1).values  # .transpose()
    user_mean = mean(user_, axis=0).astype('float64')
    shifts = user_ - user_mean

    sun, galactic_n_p, galactic_center, galactic_anti_center, S, D = getSTARS_for_galactic_system_alpha(path)
    sun = normvec(sun)
    galactic_n_p = normvec(galactic_n_p)
    galactic_center = normvec(galactic_center)
    galactic_anti_center = normvec(galactic_anti_center)

    Nx = normvec(galactic_anti_center - galactic_center).astype('float64')
    Nz = normvec(array(galactic_n_p)).astype('float64')

    SD = normvec(-array(normvec(D) - normvec(S)).astype('float64'))
    for i in range(len(SD) - 1):
        print(degrees(arccos(dot(Nx[i], Nz[i]))))
    return user_mean, shifts, Nx, Nz, SD


def process_calc_meancumulative_fastest(root_directory, resolution):
    list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]

    theta_max = math.pi  # * 2.0 #+ resolution * 0.1
    phi_max = math.pi / 2.0
    rot_theta = arange(-theta_max, theta_max, resolution)
    rot_phi = arange(-phi_max, phi_max, resolution)
    print(len(rot_theta), len(rot_phi))
    data = []
    SD_lat_lon = []
    for directory in list_subfolders_with_paths[:1]:
        total_dir = os.path.join(directory, "allsatellites")
        print(total_dir)
        # user_mean, shifts, Nx, Nz, SD = prepare_data_for_CMBR_cmap_alpha(total_dir)
        user_mean, shifts, Nx, Ny, Nz, SD = prepare_data_for_CMBR_cmap(total_dir)


# =================================================================================================
# =================================================================================================


nr = radians(15.0)
root_directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/perth_sept_6_12"
# root_directory = r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/mean_projection/perth_sept_6_12"
list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]

process_calc_meancumulative_fastest(root_directory, nr)
