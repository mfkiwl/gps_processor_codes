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
    D = normvec(pd.read_csv(path + '/Denebola_positions.csv', skiprows=0).values)
    S = normvec(pd.read_csv(path + '/Sadalmelik_positions.csv', skiprows=0).values)

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


def groupvec(vec, l):
    groups = []
    groups_std = []
    ll = int(len(vec) / l)
    for i in range(l):
        groups.append(vec[i * ll:(i + 1) * ll])
        groups_std.append(std(vec[i * ll:(i + 1) * ll]))
    groups = array(groups)
    groups_std = array(groups_std)
    return groups


def p_m_byv(d1, d2, v):
    pn = 0
    pp = 0
    for x in range(len(d1)):
        if float(d1[x]) < v and float(d1[x]) >= 0:
            pp += d2[x] / (abs(d2[x]) + 0.1)
        if float(d1[x]) > -v and float(d1[x]) <= 0:
            pn += d2[x] / (abs(d2[x]) + 0.1)
    return pn, pp


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
    # print(list(d1_n_idx),'\n',list(d1_p_idx))
    pn = sum(d2[d1_n_idx])
    pp = sum(d2[d1_p_idx])
    nr_p = len(d1_p_idx)
    nr_n = len(d1_n_idx)
    try:
        pp = pp / nr_p
    except:
        pass
    try:
        pn = pn / nr_n
    except:
        pass
    return pn, pp


def earth_DS(user_mean, SD):
    proj = []
    for x in SD:
        proj.append(dot(unit(user_mean), x))
    return array(proj)


def projSD_users(groupusers, SD):
    l = len(SD)
    proj2 = empty(l)
    for i in range(l):
        a = 0
        for j in groupusers[i]:
            # a += dot(norm(j),SD[i])/abs(dot(norm(j),SD[i]))
            a += dot(j, SD[i])  # /abs(dot(j,SD[i]))
        proj2[i] = a / len(groupusers[i])
        a = 0
    return proj2


def projUMEAN_users(groupusers, user_mean, l):
    proj2 = empty(l)
    std_shift_groups = empty(l)
    for i in range(l):
        a = 0
        for j in groupusers[i]:
            a += dot(j, unit(user_mean))  # /abs(dot(j,norm(user_mean)))
        proj2[i] = a
        std_shift_groups[i] = mean(std(groupusers[i], axis=0))
        a = 0
    return proj2, std_shift_groups


def projUMEAN_users_mean(groupusers, user_mean, l):
    proj2 = empty(l)
    std_shift_groups = empty(l)
    for i in range(l):
        proj2[i] = dot(mean(groupusers[i], axis=0), unit(user_mean))
        std_shift_groups[i] = mean(std(groupusers[i], axis=0))
        a = 0
    return proj2, std_shift_groups


def phi_theta_to_matrix(matrix, S_phi_theta, D_phi_theta, resolution, root_directory, val=100):
    phi_S, phi_D = S_phi_theta[1], D_phi_theta[1]
    theta_S, theta_D = S_phi_theta[0], D_phi_theta[0]
    l_p = len(matrix.T)
    l_t = len(matrix.T[0])
    print('Matrix: ')
    print(l_p, l_t)
    M = zeros(shape(matrix))
    theta_max = math.pi
    phi_max = math.pi / 2.0
    thetas = degrees(arange(-theta_max, theta_max, resolution))
    phis = degrees(arange(-phi_max, phi_max, resolution))
    # print(thetas, phis)
    for i in range(l_p - 1):
        for j in range(l_t - 1):
            # print('phis D: ', phi_D, phis[i], phis[i+1])
            # print('thetas D: ', theta_D, thetas[i], thetas[i+1])
            if ((phi_D > phis[i]) and (phi_D < phis[i + 1]) and (theta_D > thetas[j]) and (theta_D < thetas[j + 1])):
                M[j][i] = val
                print("D directions identified on the map: ", D_phi_theta)
            if ((phi_S > phis[i]) and (phi_S < phis[i + 1]) and (theta_S > thetas[j]) and (theta_S < thetas[j + 1])):
                M[j][i] = val
                print("S directions identified on the map: ", S_phi_theta)
    name = os.path.split(root_directory)[-1] + '_24h'
    plot_mollweid_SD(M.T, root_directory, name)
    return M.T


def plot_mollweid_SD(matrix, root_directory, name):
    plt.clf()
    try:
        cmap_save = pd.DataFrame(matrix)
        cmap_save.to_csv(os.path.join(root_directory, name + '_GS_SD_data.csv'), index=False)
        print('Matrix.csv saved!')
    except:
        print()

    ra = linspace(-math.pi, math.pi, len(matrix))
    dec = linspace(-math.pi / 2, math.pi / 2, len(matrix[0]))
    X, Y = meshgrid(ra, dec)
    Z = matrix.T
    fig = plt.figure()
    ax = pl.subplot(111, projection='mollweide')
    fig = ax.contourf(X, Y, Z, 100)
    ax.set_title('Perth$(24h)$', fontsize=20, fontweight='bold')

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


def plot_cmap(color_matrix, angle_lim1, angle_lim2, root_directory):
    plt.clf()
    try:
        cmap_save = pd.DataFrame(color_matrix)
        cmap_save.to_csv(os.path.join(root_directory, 'galactic_cmap.csv'), index=False)
        print('Matrix.csv saved!')
    except:
        print()
    fig = plt.imshow(color_matrix, extent=[0, degrees(angle_lim1), 0, degrees(angle_lim2)], aspect='auto')

    plt.title('24h_average, SD rotated')
    # plt.legend()
    # plt.tight_layout()
    plt.ylabel('Theta')
    plt.xlabel('Phi')
    plt.colorbar(fig)
    # plt.show()
    fig_name1 = os.path.join(root_directory, '_galactic_byXdeg.png')

    try:
        plt.savefig(fig_name1)
    except:
        print()


def plot_mollweid(matrix, root_directory):
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
    # ra = np.linspace(0, 2*np.pi, 72)
    # dec= np.linspace(0, np.pi, 36)
    X, Y = meshgrid(ra, dec)
    Z = matrix.T
    plt.figure()
    ax = pl.subplot(111, projection='mollweide')
    fig = ax.contourf(X, Y, Z, 100)
    ax.set_title('RO$(24h)$', fontsize=20, fontweight='bold')
    plt.xlabel(r'$\theta$', fontsize=20)  # Italic font method
    plt.ylabel(r'$\phi$', fontsize=20)  # Bold font method without fontweight parameters
    pl.colorbar(fig)
    # ax.grid()
    # ax.contour(X,Y,Z,10,colors='k')
    # pl.show()

    fig_name1 = os.path.join(root_directory, 'GS_' + child_dirname + '.png')
    pl.savefig(fig_name1, bbox_inches='tight')

    try:
        pl.savefig(fig_name1, bbox_inches='tight')

    except:
        print()


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
                       -math.cos(phi) * math.sin(theta),
                       math.sin(phi)]
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R2 = transform_matrix(system, ECEF)
    rotated_vector = R2.dot(v_new_in_system)
    return rotated_vector


def cartesian_to_galactic(system, vector):
    vector = unit(vector)
    s_phi = round(dot(vector, system[2]), 4)
    s_theta = -round(dot(vector, system[1]), 4) / (sqrt(1.0 - s_phi ** 2))
    phi = arcsin(round(s_phi, 3))
    theta = arcsin(round(s_theta, 3))
    # print(degrees(theta), degrees(phi))
    return degrees(theta), degrees(phi)


def direction_3D_galactic_system(Nx, Ny, Nz, S, D, theta, phi):
    galactic_direction = empty(shape(Nx))
    l = min(len(Nx), len(Ny), len(Nz), len(S), len(D))
    S_lat_lon = empty((l, 2))
    D_lat_lon = empty((l, 2))
    for i in range(l):
        system = array([Nx[i], Ny[i], Nz[i]])
        galactic_direction[i] = vector_3D_galactic_system(system, theta, phi)
        S_lat_lon[i] = cartesian_to_galactic(system, S[i])
        D_lat_lon[i] = cartesian_to_galactic(system, D[i])
    return galactic_direction, S_lat_lon, D_lat_lon


def cumulative_limits(user_mean, SD, usergrupped):
    l = len(SD)
    # print('Std. of the shifts:', std(shifts, axis=0))
    d1 = earth_DS(user_mean, SD)  # SD-n angles
    minange = max(abs(d1))
    # d2 = projSD_users(usergrupped,SD)
    # d2, std_d2 = projUMEAN_users(usergrupped,user_mean,l)	# sum(+-1) a csoporton belul
    d2, std_d2 = projUMEAN_users_mean(usergrupped, user_mean, l)
    # print('Mean projections: ', mean(absolute(d2)))
    # print('Mean projections: ', mean(d2))
    d1 = array(d1).astype('float64')
    d2 = array(d2).astype('float64')
    # diff = p_m_byv_beta(d1,d2,minange)
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


def prepare_data_for_CMBR_cmap(path):
    user_ = pd.read_csv(path + '/user_pos_allsatellites.csv', skiprows=1).values  # .transpose()
    user_mean = mean(user_, axis=0).astype('float64')
    shifts = user_ - user_mean
    sun, galactic_n_p, galactic_center, nunki, galactic_anti_center, S, D = getSTARS_for_galactic_system(path)
    S = normvec(S)
    D = normvec(D)
    # SD = normvec(S - D).astype('float64')

    sun = normvec(sun)
    galactic_n_p = normvec(galactic_n_p)
    galactic_center = normvec(galactic_center)
    galactic_anti_center = normvec(galactic_anti_center)
    nunki = normvec(nunki)

    # Nx = normvec(galactic_anti_center - galactic_center).astype('float64')
    # Nz = normvec(galactic_n_p - sun).astype('float64')

    Nx = normvec(galactic_center).astype('float64')
    Nz = normvec(galactic_n_p).astype('float64')

    Ny = get_third_dir_by_cross_product(Nz, Nx)

    # Nx = normvec(- sun + galactic_center).astype('float64')
    # Nz = normvec(galactic_n_p - sun).astype('float64')

    # for i in range(len(SD)-1):
    # 	print(degrees(arccos(dot(Nx[i], Nz[i]))))

    return user_mean, shifts, Nx, Ny, Nz, S, D


def cut24h(filepath):
    sec_per_day = 24 * 3600
    times = pd.read_csv(filepath + '/times_allsatellites.csv', sep=' ', skiprows=1).values.T[-1]
    print(times)
    tot_measurements = len(times)
    n = 0
    while times[n] - times[0] < sec_per_day:
        n += 1
    return n


def prepare_data_for_CMBR_cmap_filter24h(path, n_24h):
    user_ = pd.read_csv(path + '/user_pos_allsatellites.csv', skiprows=1).values
    n_tot = len(user_)

    user_mean = mean(user_, axis=0).astype('float64')
    shifts = user_ - user_mean
    sun, galactic_n_p, galactic_center, nunki, galactic_anti_center, S, D = getSTARS_for_galactic_system(path)

    S = normvec(S)
    D = normvec(D)
    # SD = normvec(S - D).astype('float64')

    sun = normvec(sun)
    galactic_n_p = normvec(galactic_n_p)
    galactic_center = normvec(galactic_center)
    galactic_anti_center = normvec(galactic_anti_center)
    nunki = normvec(nunki)

    Nx = normvec(galactic_center).astype('float64')
    Nz = normvec(galactic_n_p).astype('float64')

    Ny = get_third_dir_by_cross_product(Nz, Nx)
    l = len(Nx)

    l_24h = int(n_24h * l / float(n_tot))
    print('24h cut Star checks: {}; measurements: {}. Total: {}; {}'.format(l_24h, n_24h, l, n_tot))
    return user_mean, shifts[:n_24h], Nx[:l_24h], Ny[:l_24h], Nz[:l_24h], S[:l_24h], D[:l_24h]


def calculate_using_directions_in_galactic_coor(user_mean, usergrupped, Nx, Ny, Nz, S, D, theta, phi):
    rDIRECTION, S_lat_lon, D_lat_lon = direction_3D_galactic_system(Nx, Ny, Nz, S, D, theta, phi)
    cum_lims = cumulative_limits(user_mean, rDIRECTION, usergrupped)

    return cum_lims, S_lat_lon, D_lat_lon


def process_calc_meancumulative_fastest(root_directory, resolution, filter_24=False):
    list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]

    theta_max = math.pi  # * 2.0 #+ resolution * 0.1
    phi_max = math.pi / 2.0
    rot_theta = arange(-theta_max, theta_max, resolution)
    rot_phi = arange(-phi_max, phi_max, resolution)
    print(len(rot_theta), len(rot_phi))
    data = []
    S_Dir = []
    D_Dir = []
    for directory in list_subfolders_with_paths:
        try:
            S_dir = []
            D_dir = []
            total_dir = os.path.join(directory, "allsatellites")
            print(total_dir)
            if filter_24:
                n_24h = cut24h(total_dir)
                user_mean, shifts, Nx, Ny, Nz, S, D = prepare_data_for_CMBR_cmap_filter24h(total_dir, n_24h)
            else:
                user_mean, shifts, Nx, Ny, Nz, S, D = prepare_data_for_CMBR_cmap(total_dir)
            print('STD of positions: ', std(shifts, axis=0))
            usergrupped = groupvec(shifts, len(Nx))
            del shifts
            cmap = empty((len(rot_theta), len(rot_phi)))
            for j in range(len(rot_phi)):
                print('Phi index: {}/{}'.format(j + 1, len(rot_phi)))
                for i in range(len(rot_theta)):
                    theta = rot_theta[i]
                    phi = rot_phi[j]
                    cmap[i][j], s_dir, d_dir = calculate_using_directions_in_galactic_coor(user_mean, usergrupped, Nx,
                                                                                           Ny, Nz, S, D, theta, phi)
                    S_dir.append(mean(array(s_dir), axis=0))
                    D_dir.append(mean(array(d_dir), axis=0))

            S_Dir.append(mean(array(S_dir), axis=0))
            D_Dir.append(mean(array(D_dir), axis=0))
            print(S_Dir, '\n', D_Dir)

            data.append(cmap)
            plot_mollweid(cmap, directory)
            phi_theta_to_matrix(cmap, S_Dir[-1], D_Dir[-1], resolution, directory)
        except:
            continue

    mean_data = mean(array(data), axis=0)
    plot_mollweid(mean_data, root_directory)
    S_Dir = mean(array(S_Dir), axis=0)
    D_Dir = mean(array(D_Dir), axis=0)
    phi_theta_to_matrix(mean_data, S_Dir, D_Dir, resolution, root_directory)


# =================================================================================================
# =================================================================================================


nr = radians(5.0)
root_directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/januar_5_12"

root_directories = [
    r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/mean_projection/januar_5_12",
    r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/mean_projection/perth_julius_06_11",
    r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/mean_projection/perth_augusztus_2_8",
    r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/mean_projection/perth_october_4_10",
    r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/mean_projection/perth_sept_6_12",
    r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/mean_projection/perth_september_27_03",
]

root_directories = [
    r"/Users/kelemensz/Documents/Research/GPS/reciever_data/perth_data/NewFolder/febr_9_15",
    r"/Users/kelemensz/Documents/Research/GPS/reciever_data/perth_data/NewFolder/marc_8_14"

]

root_directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/RO/cmap"
# list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]

# process_calc_meancumulative_fastest(root_directory, nr, filter_24=True)

for directory in root_directories:
    try:
        process_calc_meancumulative_fastest(directory, nr)
    except:
        pass
