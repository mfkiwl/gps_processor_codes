from numpy import *
import pandas as pd
import matplotlib.pyplot as plt
from sympy import Point3D, Line3D, Segment3D, Point, Line, Plane, Poly
import pymap3d as pm
import os
from vpython import vector, rotate


def getSTARS(path, pair_index):
    if pair_index == 0:
        star_1 = pd.read_csv(path + '/star_d_positions.csv', skiprows=0).values  # .transpose()[0]#.astype('float64')
        star_2 = pd.read_csv(path + '/star_s_positions.csv', skiprows=0).values  # .transpose()[0]#.astype('float64')
    elif pair_index == 1:
        star_1 = pd.read_csv(path + '/star_Mirphak_positions.csv',
                             skiprows=0).values  # .transpose()[0]#.astype('float64')
        star_2 = pd.read_csv(path + '/star_Vega_positions.csv',
                             skiprows=0).values  # .transpose()[0]#.astype('float64')

    star_1 = normvec(star_1)
    star_2 = normvec(star_2)
    stars = array([star_1, star_2])
    return stars


def norm_(v):
    a = array(v)
    return a / (sqrt(sum(power(a, 2))))


def normvec(vec):
    a = empty(shape(vec))
    for i in range(len(vec)):
        a[i] = norm_(vec[i])
    return a


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
        pn = pn / nr_n
    except:
        pass
    return pn, pp


def earth_DS(user_mean, SD):
    proj = []
    for x in SD:
        proj.append(dot(norm_(user_mean), x))
    return array(proj)


def u_mean_in_stationarry_system(user_mean, SD, s2, s3):
    l = len(SD)
    # d1 = earth_DS(user_mean,SD)  # SD-n angles
    # print('Phi: \n', d1)
    n_in_stars = empty((l, 2))
    for i in range(l):
        system = flip(gram_schmidt(flip(array([s2[i], s3[i], SD[i]]))))
        n_in_stars[i] = get_n_in_stars(SD[i], system, user_mean)
    return n_in_stars


def plot_cmap(color_matrix, angle_lim1, angle_lim2, root_directory):
    plt.clf()
    try:
        cmap_save = pd.DataFrame(color_matrix)
        cmap_save.to_csv(os.path.join(root_directory, 'dat_cmap_by5degree.csv'), index=False)
        print('Matrix.csv saved!')
    except:
        print()
    fig = plt.imshow(color_matrix, extent=[0, degrees(angle_lim1), 0, degrees(angle_lim2)], aspect='auto')

    plt.title('24h_average, SD rotated')
    plt.ylabel('Theta')
    plt.xlabel('Phi')
    plt.colorbar(fig)
    # plt.show()
    fig_name1 = os.path.join(root_directory, 'cmap_by5degree.png')
    fig_name2 = os.path.join(root_directory, '_cmap_by5degree.png')

    try:
        plt.savefig(fig_name1)
    except:
        print()
    try:
        plt.imsave(fig_name2, fig)
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


def project_onto_plane(v, n_plane):
    d = dot(v, n_plane) / linalg.norm(n_plane)
    p = [d * unit(n_plane)[i] for i in range(len(n_plane))]
    return [v[i] - p[i] for i in range(len(v))]


def get_n_in_stars(SD, system, u_mean):
    # phi = degrees(arccos(dot(n, SD)))
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R2 = transform_matrix(ECEF, system)  # from ECEF to system
    u_in_system = unit(R2.dot(u_mean))
    phi = degrees(arccos(round(dot(u_mean, system[2]), 4)))

    u_in_plane = unit(array([dot(u_in_system, system[0]), dot(u_in_system, system[1]), 0.0]))

    theta = degrees(arctan(u_in_plane[0] / u_in_plane[1]))
    # if theta<0:
    # 	thet = theta + 180
    return phi, theta


# =================================================================================================
# =================================================================================================
def prepare_data(path):
    user_ = pd.read_csv(path + '/user_pos_allsatellites.csv', skiprows=1).values  # .transpose()
    user_mean = mean(user_, axis=0).astype('float64')
    shifts = user_ - user_mean
    stars_SD = getSTARS(path, 0)
    SD = -array(stars_SD[0] - stars_SD[1]).astype('float64')
    SD = normvec(SD)
    stars_secondary = getSTARS(path, 1)
    stars_secondary_A, stars_secondary_B = stars_secondary[0], stars_secondary[1]
    return user_mean, shifts, SD, stars_secondary_A, stars_secondary_B


def data_to_matrix(data, phis, thetas, val):
    l_p = len(phis)
    l_t = len(thetas)
    M = zeros((l_p, l_t))
    # print(data)
    # print(thetas)
    # print(phis)
    for i in range(l_t - 1):
        for j in range(l_p - 1):
            for phi, theta in data:
                if phi > phis[j] and phi < phis[j + 1] and theta > thetas[i] and theta < thetas[i + 1]:
                    M[j][i] = val
    return M.T


def plot_cmap(color_matrix, angle_lim1, angle_lim2, root_directory):
    plt.clf()
    try:
        cmap_save = pd.DataFrame(color_matrix)
        cmap_save.to_csv(os.path.join(root_directory, 'dat_cmap_by5degree.csv'), index=False)
        print('Matrix.csv saved!')
    except:
        print()
    fig = plt.imshow(color_matrix, extent=[0, degrees(angle_lim1), 0, degrees(angle_lim2)], aspect='auto')

    plt.title('24h, User direction')
    # plt.legend()
    # plt.tight_layout()
    plt.ylabel('Theta')
    plt.xlabel('Phi')
    plt.colorbar(fig)
    plt.show()
    try:
        fig_name1 = os.path.join(root_directory, 'cmap_by5degree.png')
        fig_name2 = os.path.join(root_directory, '_cmap_by5degree.png')
    except:
        pass
    try:
        plt.savefig(fig_name1)
    except:
        print()
    try:
        plt.imsave(fig_name2, fig)
    except:
        print()


def process_calc_meancumulative_faster(root_directory, resolution):
    list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]

    theta_max = math.pi * 2.0  # + resolution * 0.1
    phi_max = math.pi
    rot_theta = arange(0.0, theta_max, resolution)
    rot_phi = arange(0.0, phi_max, resolution)
    print(len(rot_theta), len(rot_phi))
    data = []
    for directory in list_subfolders_with_paths[:1]:
        total_dir = os.path.join(directory, "allsatellites")
        print(total_dir)
        user_mean, shifts, SD, stars_secondary_A, stars_secondary_B = prepare_data(total_dir)
        cmap = empty((len(rot_theta), len(rot_phi)))
        phi_theta_1 = u_mean_in_stationarry_system(unit(user_mean), SD, stars_secondary_A, stars_secondary_B)
        phi_theta_2 = u_mean_in_stationarry_system(unit(-user_mean), SD, stars_secondary_A, stars_secondary_B)
        for i in range(len(phi_theta_2)):
            phi_theta_2[i, 1] = phi_theta_2[i, 1] + 180
        print(phi_theta_2)
        M1 = data_to_matrix(phi_theta_1, degrees(rot_phi), degrees(rot_theta), 1)
        M2 = data_to_matrix(phi_theta_2, degrees(rot_phi), degrees(rot_theta), -1)

        plot_cmap(M1 + M2, phi_max, theta_max, None)


# mean_data = mean(array(data), axis=0)
# plot_cmap(mean_data, phi_max, theta_max, root_directory)


# =================================================================================================
# =================================================================================================


nr = radians(5.0)
root_directory = r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/mean_projection/perth_julius_06_11"
list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]
# print(list_subfolders_with_paths)

process_calc_meancumulative_faster(root_directory, nr)
