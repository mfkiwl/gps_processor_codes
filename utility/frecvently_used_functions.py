import matplotlib.pyplot as plt
from numpy import *
import pandas as pd

from itertools import chain
import pymap3d as pm
import os
from vpython import rotate
import pylab as pl
import math
from scipy.ndimage import rotate


def create_dir(root_path, dir_name):
    results_dir = os.path.join(root_path, dir_name)
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    return results_dir


def create_generalinfo_path(directory):
    generalinfo_path = os.path.join(directory, 'allsatellites')
    if not os.path.exists(generalinfo_path):
        os.makedirs(generalinfo_path)
    return generalinfo_path


def get_obs_files_simple(directory):
    obs_ext = ".obs"
    obs_files = []
    files = next(os.walk(directory))[2]
    for file in files:
        if os.path.splitext(file)[1] == obs_ext:
            obs_files.append(os.path.abspath(os.path.join(directory, file)))
    return obs_files


def get_pos_files_simple(directory):
    obs_ext = ".pos"
    pos_files = []
    files = next(os.walk(directory))[2]
    for file in files:
        if os.path.splitext(file)[1] == obs_ext:
            pos_files.append(os.path.abspath(os.path.join(directory, file)))
    return pos_files


def posparser(file):
    a = []
    with open(file) as in_file:
        lineList = [line.rstrip('\n').split(" ") for line in in_file]
        for line in lineList[15:]:
            # print(len(line))
            # print(line)
            # print(line[3], line[6], line[8])
            # out_file.write("{},{},{}\n".format(line[3], line[6], line[8]))
            a.append([line[3], line[6], line[8]])
    return array(a)


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


def too_bad_positions_std(path, filename, std_limit=10.0):
    path = os.path.join(path, "allsatellites")
    file = os.path.join(path, filename)
    if os.path.isfile(file):
        user_ = pd.read_csv(file, skiprows=1).values  # .transpose()
        std_pos = mean(std(user_, axis=0))
        # print("                             STD: ", std_pos)
        if std_pos > std_limit:
            return True
        return False
    return True


def is_reliable(A_day, B_day, needed_file):
    stdUA = too_bad_positions_std(A_day, needed_file)
    stdUB = too_bad_positions_std(B_day, needed_file)
    if stdUA or stdUB:
        return False
    return True


def are_reliable(A_day, B_day, needed_file, std_limits=None):
    if std_limits is None:
        std_limits = [10.0, 10.0]
    if is_all_data(A_day, [needed_file], add_allsatellites=True) and is_all_data(B_day, [needed_file],
                                                                                 add_allsatellites=True):
        stdUA = too_bad_positions_std(A_day, needed_file, std_limit=std_limits[0])
        stdUB = too_bad_positions_std(B_day, needed_file, std_limit=std_limits[1])
        if not stdUA and not stdUB:
            return True
    return False


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
        return False


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


def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


def get_matrix(path):
    return pd.read_csv(path, skiprows=0).values[:, 1:]  # .transpose()[0]#.astype('float64')


def get_csv_file(directory):
    csv_ext = ".csv"
    files = [f.path for f in os.scandir(directory) if f.is_file()]
    files_with_path = []
    for file in files:
        if os.path.splitext(file)[1] == csv_ext:
            files_with_path.append(os.path.abspath(os.path.join(directory, file)))
    return files_with_path


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


from matplotlib import cm

"""
https://stackoverflow.com/questions/22128909/plotting-the-temperature-distribution-on-a-sphere-with-python
"""


def plot_on_sphere(WW):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1, projection='3d')
    # WW = get_mean_matrix(root_dir)
    u = linspace(0, 2 * pi, len(WW))
    v = linspace(0, pi, len(WW[0]))

    # create the sphere surface
    XX = 10 * outer(cos(u), sin(v))
    YY = 10 * outer(sin(u), sin(v))
    ZZ = 10 * outer(ones(size(u)), cos(v))

    WW = WW + abs(amin(WW))
    myheatmap = WW / amax(WW)

    # ~ ax.scatter( *zip( *pointList ), color='#dd00dd' )
    ax.plot_surface(XX, YY, ZZ, cstride=1, rstride=1, facecolors=cm.jet(myheatmap))
    # plt.colorbar(cm.jet( myheatmap ))
    ax.set_xlabel('$X$')
    ax.set_ylabel('$Y$')
    ax.set_zlabel('$Z$')
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


#  -----------------------------------
# M = prepear_for_sphere(m, h)
# plot_on_sphere(M)
#  -----------------------------------


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


def handle_raw_not_averaged_matrices(M, H, N, fig_directory, name, nr_days, not_symmetrised=False, round=True):
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

    # M[M > nanmax(M) / 10] = 0
    # M[M < -1 * abs(nanmin(M)) / 10] = 0

    M[M == 0.0] = nan
    # M = M * -1
    if round:
        M[M < 0] = -1
        M[M > 0] = +1

    # M = nan_to_num(M, nan=0)
    # H = log(H)
    # plot_save_imshow_3_maps([H, M, N], ["Histogram", "(|1-r|/|n|)", "<n_mod>"], root_directory=None, resolution="5", logplot=False, show=True)

    # plot_mollweid_simple(M[::-1].T)
    id = 'symmetrized'
    id_png = 'symmetrized'
    if not_symmetrised:
        id = '-'
        id_png = 'not_symmetrized'

    plt.imshow(M)
    # plt.imshow(H)
    plt.colorbar()

    plt.title("<r-1/r> {} ({}_{})".format(id, name, nr_days))
    # plt.title("histogram {} ({}_{})".format(id, name, nr_days))
    fig_name1 = os.path.join(fig_directory, '{}_{}.png'.format(name, id_png))
    # fig_name1 = os.path.join(fig_directory, '{}_{}_histogram.png'.format(name, id_png))
    plt.savefig(fig_name1, bbox_inches='tight')
    plt.clf()


def select_cmap_hist_n_mod(file_list, histogram_string='', measure_string='', nmod_string=''):
    hist = None
    cmap = None
    n_mod = None
    for file in file_list:
        if histogram_string in str(os.path.split(file)[-1]):
            hist = file
        if measure_string in str(os.path.split(file)[-1]):  # or "measure" in str(os.path.split(file)[-1]):
            cmap = file
        if nmod_string in str(os.path.split(file)[-1]):
            n_mod = file
    return cmap, hist, n_mod
