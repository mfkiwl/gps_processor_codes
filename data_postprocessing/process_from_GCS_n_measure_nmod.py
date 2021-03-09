import re
import time

from numpy import *
import pandas as pd
import matplotlib.pyplot as plt

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
    D = normvec(pd.read_csv(path + '/Denebola_positions.csv', skiprows=0).values)
    S = normvec(pd.read_csv(path + '/Sadalmelik_positions.csv', skiprows=0).values)
    sun = pd.read_csv(path + '/SUN_positions.csv', skiprows=0).values
    galactic_n_p = pd.read_csv(path + '/GNP_positions.csv', skiprows=0).values
    galactic_center = pd.read_csv(path + '/GC_positions.csv', skiprows=0).values
    nunki = pd.read_csv(path + '/Nunki_positions.csv', skiprows=0).values
    capella = pd.read_csv(path + '/Capella_positions.csv', skiprows=0).values
    return sun, galactic_n_p, galactic_center, nunki, capella, S, D


def transform_matrix(f1, f2):  # transforms from f1 to f2
    R = array([
        [dot(f2[0], f1[0]), dot(f2[0], f1[1]), dot(f2[0], f1[2])],
        [dot(f2[1], f1[0]), dot(f2[1], f1[1]), dot(f2[1], f1[2])],
        [dot(f2[2], f1[0]), dot(f2[2], f1[1]), dot(f2[2], f1[2])]
    ])
    return R


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


def get_mean_direction_over_time(systems, directions):
    l = min(len(directions), len(systems[0]))
    phi = 0.
    theta = 0.
    for i in range(l):
        out = cartesian_to_galactic(array([systems[0][i], systems[1][i], systems[2][i]]), directions[i])
        theta += out[0]
        phi += out[1]
    return theta / float(l), phi / float(l)


def scal(v1, v2):
    return vdot(array(v1), array(v2))


def unit(v):
    return 1 / sqrt(scal(v, v)) * array(v)


def get_theta_phi(v):
    if 0.99 < v[0] and 1.01 > v[0] and -0.01 < v[1] and 0.01 > v[1] and -0.01 < v[2] and 0.01 > v[2]:
        return 0.0, 0.0
    if 0.99 < v[1] and 1.01 > v[1] and -0.01 < v[0] and 0.01 > v[0] and -0.01 < v[2] and 0.01 > v[2]:
        return 90.0, 0.0
    if 0.99 < v[2] and 1.01 > v[2] and -0.01 < v[1] and 0.01 > v[1] and -0.01 < v[0] and 0.01 > v[0]:
        return 0.0, 90.0
    return None, None


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
def normvec(vec):
    a = empty(shape(vec))
    for i in range(len(vec)):
        a[i] = unit(vec[i])
    return a


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
        # cmap_mod_n[i][j] += mod_n
        cmap_count[i][j] += 1
    cmap_count[cmap_count < 1] = 0
    return cmap_v, cmap_count, cmap_mod_n


def get_ij_on_map(direction, resolution):
    # theta_max = math.pi
    # phi_max = math.pi / 2.0
    # rot_theta = arange(-theta_max, theta_max, resolution)
    # rot_phi = arange(-phi_max, phi_max, resolution)
    theta, phi = cart2sph(direction[0], direction[1], direction[2])
    I_f = 0
    J_f = 0
    l_theta = 36  # int(len(rot_theta) / 2.0)
    l_phi = 18  # int(len(rot_phi) / 2.0)
    for i in range(-l_theta, l_theta):
        if i * resolution > theta:
            I_f = i + l_theta
            break
    for i in range(-l_phi, l_phi):
        if i * resolution > phi:
            J_f = i + l_phi
            break
    return I_f, J_f


def calc_correct_average(H, M, new_shape):
    shape = (new_shape[0], H.shape[0] // new_shape[0], new_shape[1], H.shape[1] // new_shape[1])
    H = nan_to_num(H, nan=0.0)
    M = nan_to_num(M, nan=0.0)
    MH = multiply(M, H)
    MH = nan_to_num(MH, nan=0.0)
    H_reshaped_summed = H.reshape(shape).sum(-1).sum(1)
    MH_reshaped_summed = MH.reshape(shape).sum(-1).sum(1)
    reshaped = divide(MH_reshaped_summed, H_reshaped_summed)
    return reshaped


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


def get_stars_in_GCS(star_dir, day_date):
    Nx, Ny, Nz, S, D = get_global_stars(star_dir, day_date)
    star_directions_in_GCS = {'GNP': get_mean_direction_over_time(array([Nx, Ny, Nz]), Nz),
                              'GC': get_mean_direction_over_time(array([Nx, Ny, Nz]), Nx),
                              'OY of GCS': get_mean_direction_over_time(array([Nx, Ny, Nz]), Ny),
                              'Sadalmelic': get_mean_direction_over_time(array([Nx, Ny, Nz]), S),
                              'Denebola': get_mean_direction_over_time(array([Nx, Ny, Nz]), D)}
    return star_directions_in_GCS


def func_string(lista, i):
    # lista[i] = list(map(float, list(lista[i].strip('[').strip(']').split())))
    lista[i] = list(map(float, lista[i].strip('[').strip(']').split()))


def get_raw_GCS_data(day_path, filename="GCS_Ndir_measure_Nmod_AB_and_BA.csv"):
    df = pd.read_csv(os.path.join(day_path, filename), index_col=0)
    rows = df.values
    ndirection = rows[:, 0].tolist()
    measure_nmod = rows[:, 1:].tolist()
    list(map(lambda i: func_string(ndirection, i), range(0, len(ndirection))))
    if len(ndirection) == len(measure_nmod):
        for i in range(len(ndirection)):
            # measure_nmod[i] = ndirection[i] + measure_nmod[i]
            measure_nmod[i].insert(0, ndirection[i])
        return array(measure_nmod, dtype=object)
    return 0


def operations_on_raw_data(data):
    # data = pd.DataFrame(data)
    # print(data)
    # data[1] = data[1] * data[2]
    # print(data)
    l = int(len(data) / 2)
    data = data[:l]
    return data  # .values


def process_one_day_rawGCS(day_path, resolution, fill_out=0.0):
    raw_results_GCS = get_raw_GCS_data(day_path)
    # raw_results_GCS = operations_on_raw_data(raw_results_GCS)
    if type(raw_results_GCS) != int:
        print("Raw results in GCS are in! (data size)", len(raw_results_GCS))
        day_data, day_count, day_cmap_n_mod = process_raw_GCS_data(raw_results_GCS, resolution)
        day_data = nan_to_num(day_data, nan=fill_out)
        day_cmap_n_mod = nan_to_num(day_cmap_n_mod, nan=fill_out)
        day_count = nan_to_num(day_count, nan=fill_out)
        return day_data, day_count, day_cmap_n_mod
    return 0, 0, 0


def create_averaged_plots_from_root(root_0, star_dir, months=None):
    execution_start = time.time()
    sum_all_cmap = []
    sum_all_hist = []
    sum_all_n_mod = []
    subfolders_with_paths_months = [f.path for f in os.scandir(root_0) if f.is_dir()]
    # star_directions_in_GCS = get_stars_in_GCS(star_dir, day_date="20200101")
    for month_root in subfolders_with_paths_months:
        month_name = str(month_root).split("/")[-1]
        if months and month_name in months:
            days_with_paths = [f.path for f in os.scandir(month_root) if f.is_dir()]
            print("Month name: ", month_name, "  nr days: ", len(days_with_paths))
            for day_root in days_with_paths:
                start = time.time()
                M, H, N = process_one_day_rawGCS(day_root, resolution)
                if type(M) != int:
                    sum_all_cmap.append(M)
                    pd.DataFrame(M).to_csv(
                        os.path.join(day_root, "sum_measure_r_inv_r_" + str(int(degrees(resolution))) + '.csv'),
                        index=True)
                    sum_all_hist.append(H)
                    sum_all_n_mod.append(N)

                print('Elapsed time of the current day: ', time.time() - start, day_root.split("/")[-1])
    print("Total number of days:  ", len(sum_all_cmap))
    sum_all_cmap = sum(array(sum_all_cmap), axis=0)
    sum_all_hist = sum(array(sum_all_hist), axis=0)
    sum_all_n_mod = sum(array(sum_all_n_mod), axis=0)
    print("\n Running time is {} minutes.\n".format((time.time() - execution_start) / 60.0))
    return sum_all_cmap, sum_all_hist, sum_all_n_mod


def handle_raw_not_averaged_matrices(M, H, N):
    ind_no_data = array(H < 1)

    M[ind_no_data] = 0
    N[ind_no_data] = 0.0

    H[H < 1] = 0
    N[N < 0] = 0

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
    M[M < 0] = -1
    M[M>0] = 1

    # M = nan_to_num(M, nan=0)
    # H = log(H)
    # plot_save_imshow_3_maps([H, M, N], ["Histogram", "(|1-r|/|n|)", "<n_mod>"], root_directory=None, resolution="5", logplot=False, show=True)

    plt.imshow(M)
    plt.colorbar()

    # plot_mollweid_simple(M[::-1].T)
    # plt.title("<|1-r|/|n|>")
    plt.show()


star_dir = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020"
resolution = radians(5.0)

all_months = ["julius", "szeptember", "februar", "marcius", "augusztus", "januar", "december2019", "oktober",
              "november", "majus", "aprilis", "junius", "december2020"]
months = ["januar"]
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_NASA/r_inv_r_over_Nmod_symmetrized"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_NZLD/r_inv_r_symmetrized"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/NASA_IIGC/r_inv_r_symmetrized"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_IIGC/r_inv_r_symmetrized"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/NZLD_IIGC/r_inv_r_symmetrized"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_IIGC/r_inv_r_symmetrized"

result_roots = [
    r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/NASA_IIGC/r_inv_r_symmetrized",
    r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_IIGC/r_inv_r_symmetrized",
    r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/NZLD_IIGC/r_inv_r_symmetrized",
    r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_IIGC/r_inv_r_symmetrized"]

# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/NZLD_TIDV/r_inv_r_symmetrized"
results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_test"
m, h, n = create_averaged_plots_from_root(results_root, star_dir, all_months)
handle_raw_not_averaged_matrices(m, h, n)

# for result_root in result_roots:
#     m, h, n = create_averaged_plots_from_root(result_root, star_dir, all_months)

# =================================================================================================
# =================================================================================================
