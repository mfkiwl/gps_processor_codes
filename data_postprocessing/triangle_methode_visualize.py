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


def rotateAntiClockwise(array):
    return rotate(array, 90)


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


def save_matrices(matrices, names, directory, resolution):
    for matrix, name in zip(matrices, names):
        cmap_save = pd.DataFrame(matrix)
        cmap_save.to_csv(os.path.join(directory, name + "_not_averaged_" + str(int(degrees(resolution))) + '.csv'),
                         index=True)
        print('Matrix.csv saved!: ', name)


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
        if "measure_not_averaged" in str(os.path.split(file)[-1]):  # or "measure" in str(os.path.split(file)[-1]):
        # if "divNmod" in str(os.path.split(file)[-1]) or "sum_measure_" in str(os.path.split(file)[-1]):

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
    # print(hist, "\n", list(H[-1]), "\n", shape(array(H)))
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
    M[ind_no_data] = 0.0
    N[ind_no_data] = 0.0

    H[H < 1] = 0

    # N[N<0]=0

    # nan_to_num(H, nan=0.0)
    nan_to_num(M, nan=0.0)
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
    M[M == 0.0] = nan
    # M = M * -1
    M[M < 0] = -1
    M[M > 0] = 1

    # M = nan_to_num(M, nan=0)
    # H = log(H)
    # plot_save_imshow_3_maps([H, M, N], ["Histogram", "(|1-r|/|n|)", "<n_mod>"], root_directory=None, resolution="5", logplot=False, show=True)

    plt.imshow(M)
    plt.colorbar()

    # plot_mollweid_simple(M[::-1].T)
    # plt.title("<|1-r|/|n|>")
    plt.show()


# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_NASA/r_inv_r_BA"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_NASA/r_inv_r_AB"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/NZLD_HKKS/r_inv_r_AB"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_NZLD/r_inv_r_symmetrized"


# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_PERTH/r_inv_r_BA"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_PERTH/r_inv_r_AB"


# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_NZLD/r_inv_r_AB"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_NZLD/r_inv_r_BA"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_NZLD/r_inv_r_symmetrized"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_NASA/r_inv_r_symmetrized"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_NASA/r_inv_r_over_Nmod_symmetrized"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/NASA_IIGC/r_inv_r_symmetrized"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_PERTH/r_inv_r_symmetrized"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_IIGC/r_inv_r_symmetrized"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_IIGC/r_inv_r_symmetrized"

results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/NZLD_TIDV/r_inv_r_symmetrized"

all_months = ["julius", "szeptember", "februar", "marcius", "augusztus", "januar", "december2019", "oktober",
              "november", "majus", "aprilis", "junius", "december2020"]
months1 = ["julius", "szeptember", "augusztus", "november", "junius", "december2020"]
months2 = ["majus", "februar", "marcius", "aprilis", "januar"]
months3 = ["februar"]


m, h, n = create_averaged_plots_from_root(results_root, all_months)
handle_raw_not_averaged_matrices(m, h, n)


def call_separatelly(results_root, all_months):
    for month in all_months:
        try:
            print("Month name:  ", month)
            m, h, n = create_averaged_plots_from_root(results_root, [month])
            handle_raw_not_averaged_matrices(m, h, n)
        except:
            pass


# call_separatelly(results_root, all_months)


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

# M = prepear_for_sphere(m, h)


# plot_on_sphere(M)

# =====================================================================================
