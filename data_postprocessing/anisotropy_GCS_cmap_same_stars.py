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


def p_m_byv_beta(d1, d2, v):
    pn = 0
    pp = 0
    nr_p = 0
    nr_m = 0
    for x in range(len(d1)):
        if float(d1[x]) < v and float(d1[x]) >= 0:
            pp += d2[x]
            nr_p += 1
        if float(d1[x]) > -v and float(d1[x]) <= 0:
            pn += d2[x]
            nr_m += 1
    try:
        pp = pp / nr_p
    except:
        pass
    try:
        pn = pn / nr_m
    except:
        pass
    return pn, pp


def p_m_byv_limits(d1, d2):
    d1 = array(d1)
    d2 = array(d2)
    d1_n_idx = where(d1 < 0.0)
    d1_p_idx = where(d1 > 0.0)
    d2_poz = d2[d1_p_idx]
    d2_neg = d2[d1_n_idx]
    pn = 0.0
    pp = 0.0
    if not d2_poz.size == 0:
        pp = mean(d2_poz)

    if not d2_neg.size == 0:
        pn = mean(d2_neg)
    return pn, pp


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


def plot_mollweid_SD(matrix, root_directory, name):
    plt.clf()
    try:
        cmap_save = pd.DataFrame(matrix.T)
        cmap_save.to_csv(os.path.join(root_directory, name + '_GS_SD_data.csv'), index=False)
        print('Matrix.csv saved!')
    except:
        print()

    ra = linspace(-math.pi, math.pi, len(matrix))
    dec = linspace(-math.pi / 2, math.pi / 2, len(matrix[0]))
    X, Y = meshgrid(ra, dec)
    Z = matrix
    fig = plt.figure()
    ax = pl.subplot(111, projection='mollweide')
    fig = ax.contourf(X, Y, Z, 100)
    ax.set_title('$(24h)$', fontsize=20, fontweight='bold')

    plt.xlabel(r'$\theta$', fontsize=20)  # Italic font method
    plt.ylabel(r'$\phi$', fontsize=20)  # Bold font method without fontweight parameters
    pl.colorbar(fig)
    # ax.grid()
    # ax.contour(X,Y,Z,10,colors='k')
    # pl.show()

    fig_name1 = os.path.join(root_directory, name + '_GS_SD_fig.png')

    try:
        pl.savefig(fig_name1, bbox_inches='tight')
    except:
        print()


def add_star_annotated(theta, phi, name, ax):
    # print(S[0],S[1])
    # print(D[0],D[1])
    # y = radians(array([S[1], D[1]]))
    # x = radians(array([S[0], D[0]]))
    # print(x, y)
    theta = radians(theta)
    phi = radians(phi)
    ax.text(theta, phi, name, fontsize=12)
    ax.scatter(theta, phi, marker='x', c='k', s=15)


# print(theta, phi)
# ax.annotate(name,
#            xy=array(theta, phi),
#            xycoords='data')
# ,

#            arrowprops=
#                dict(facecolor='black', shrink=0.05),
#                horizontalalignment='left',
#                verticalalignment='top')


def plot_mollweid(matrix, star_directions, root_directory):
    child_dirname = os.path.split(root_directory)[-1] + '_24h'
    plt.clf()
    try:
        cmap_save = pd.DataFrame(matrix)
        cmap_save.to_csv(os.path.join(root_directory, 'data_GS_' + child_dirname + '.csv'), index=False)
        print('Matrix.csv saved!')
    except:
        print()

    ra = linspace(-math.pi, math.pi, len(matrix))
    dec = linspace(-math.pi / 2, math.pi / 2, len(matrix[0]))

    X, Y = meshgrid(ra, dec)
    Z = matrix.T
    plt.figure()
    ax = pl.subplot(111)  # , projection = 'mollweide')
    fig = ax.contourf(X, Y, Z, 100)

    for k, v in star_directions.items():
        add_star_annotated(v[0], v[1], k, ax)

    # ax.set_title('---$(24h)$', fontsize=15)  # , fontweight='bold')
    plt.xlabel(r'$\theta$', fontsize=15)  # Italic font method
    plt.ylabel(r'$\phi$', fontsize=15)  # Bold font method without fontweight parameters
    pl.colorbar(fig)
    ax.grid()
    # ax.contour(X,Y,Z,10,colors='k')
    # pl.show()

    fig_name1 = os.path.join(root_directory, 'GS_' + child_dirname + '.png')
    pl.savefig(fig_name1, bbox_inches='tight')


# try:
# 	pl.savefig(fig_name1, bbox_inches='tight')

# except:
# 	print()


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


def cartesian_to_galactic_alpha(system, vector):
    vector = unit(vector)
    s_phi = round(dot(vector, system[2]), 4)
    phi = arcsin(round(s_phi, 3))

    s_theta = -round(dot(vector, system[1]), 4) / (sqrt(1.0 - s_phi ** 2))
    theta = arcsin(round(s_theta, 3))
    return degrees(theta), degrees(phi)


def cartesian_to_galactic_beta(system, vector):
    vector = unit(vector)
    system = normvec(system)

    s_phi = round(dot(vector, system[2]), 4)
    phi = arcsin(round(s_phi, 3))

    tg_theta = vector[1] / vector[0]
    theta = math.atan(round(tg_theta, 3))

    if vector[1] < 0.0 and vector[0] < 0.0:
        theta = theta + radians(180.0)

    # elif vector[1] > 0.0 :
    # 	theta = theta + radians(180.0)

    # elif -vector[1] > 0.0 and vector[0] > 0.0:
    # 	theta = theta

    elif vector[1] < 0.0:
        theta = -theta  # + radians(180.0)
    # print(degrees(theta), degrees(phi))
    return degrees(theta), degrees(phi)


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


def cumulative_limits(user_mean, SD, usergrupped):
    l = len(SD)
    d1 = earth_DS(user_mean, SD)  # SD-n angles
    minange = max(abs(d1))
    d2, std_d2 = projUMEAN_users_mean(usergrupped, user_mean, l)
    # print('Mean projections: ', mean(absolute(d2)))
    # print('Mean projections: ', mean(d2))
    d1 = array(d1).astype('float64')
    d2 = array(d2).astype('float64')
    diff = p_m_byv_limits(d1, d2)
    return diff[0] - diff[1]


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


def posparser(file):
    a = []
    with open(file) as in_file:
        lineList = [line.rstrip('\n').split(" ") for line in in_file]
        for line in lineList[15:]:
            a.append(line[3], line[6], line[8])
    return np.array(a)


def prepare_data_for_CMBR_cmap(path, rtk_pos=None):
    if rtk_pos:
        user_ = posparser(path + '/user_pos_allsatellites.csv')
    else:
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


def calculate_using_directions_in_galactic_coor(user_mean, usergrupped, Nx, Ny, Nz, theta, phi):
    rDIRECTION = direction_3D_galactic_system(Nx, Ny, Nz, theta, phi)
    cum_lims = cumulative_limits(user_mean, rDIRECTION, usergrupped)
    # print(theta, " ", phi, "   :  ", cum_lims)
    return cum_lims


def vectors_dot_prods(v1, v2):
    v1 = normvec(v1)
    v2 = normvec(v2)
    dot_v1_v2 = []
    for i in range(len(v1)):
        dot_v1_v2.append(dot(v1[i], v2[i]))
    print('Local vs. Global directions: ', degrees(arccos(mean(dot_v1_v2))))
    return dot_v1_v2


def process_one_day(directory, star_dir, resolution):
    # S_galactic = array([-59.94133321, -42.06696326968])
    # D_galactic = array([109.37248674,  70.80087410094148])

    theta_max = math.pi  # * 2.0 #+ resolution * 0.1
    phi_max = math.pi / 2.0
    rot_theta = arange(-theta_max, theta_max, resolution)
    rot_phi = arange(-phi_max, phi_max, resolution)

    total_dir = os.path.join(directory, "allsatellites")
    print("Positional data: ", total_dir)
    user_mean, shifts, Nx_local, Ny_local, Nz_local, S_local, D_local = prepare_data_for_CMBR_cmap(total_dir)
    # print('S cos phi:', 90 - degrees(arccos(scal(Nz[0], S[0]))))
    # print('D cos phi:', 90 - degrees(arccos(scal(Nz[0], D[0]))))

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
                                                                     Nz_global, theta, phi)

    plot_mollweid(cmap, star_directions_in_GCS, directory)
    return cmap, star_directions_in_GCS


def process_calc_meancumulative_fastest(root_directory, resolution, star_dir, filter_24=False):
    list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]

    data = []
    for directory in list_subfolders_with_paths:
        try:
            cmap, star_direction = process_one_day(directory, star_dir, resolution)
            data.append(cmap)
            print("\n")
            print("\n")
        except:
            print("Data not found!")
            continue

    mean_data = mean(array(data), axis=0)
    print(star_direction)
    plot_mollweid(mean_data, star_direction, root_directory)


# =================================================================================================
# =================================================================================================


nr = radians(5.0)
# root_directories = [
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/januar",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/februar",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/marcius",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/aprilis",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/majus",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/junius",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/julius",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/augusztus",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/szeptember",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/oktober",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/november",

# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/process_NZLD/januar",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/process_NZLD/februar",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/process_NZLD/marcius",
# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/process_NZLD/oktober",

# 	r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/process_NASA/januar"
# 	]

# # 	, r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements/december2019",
# # ]


# root_directories = [
# 	r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/november",
# 	r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/aprilis",
# 	r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/majus"
# 	]

root_directories = [
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/majus"
]
root_directories = [
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/januar",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/februar",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/marcius",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/aprilis",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/majus",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/junius",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/julius",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/augusztus",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/oktober",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NZLD/november",

    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NASA/aprilis",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NASA/januar",
    r"/Users/kelemensz/Qsync/GPS/GPS_data/global_GCS_axis/process_NASA/marcius"
]

# root_directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/process_NASA/januar"

root_directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/iono_correction/with_to_send/CUTB20200104"

star_dir = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020"

# process_calc_meancumulative_fastest(root_directory, nr, star_dir)
process_one_day(root_directory, star_dir, nr)
# for directory in root_directories:
# 	# try:
# 	process_calc_meancumulative_fastest(directory, nr, star_dir)	
# 	# except:
# 		# pass


directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/process_NZLD/test/januar_PERTH/CUTB20200105"

# process_one_day(directory, star_dir, nr)
