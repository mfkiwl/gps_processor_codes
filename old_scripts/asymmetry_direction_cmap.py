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


def scalevec(vec):
    a = empty(len(vec))

    for i in range(len(vec)):
        a[i] = sqrt(sum(power(vec[i], 2)))
    s = 2 * std(a)
    return vec / s


def groupmeanvec(vec, l):
    groups = []
    for i in range(int(len(vec) / l) + 1):
        groups.append(mean(vec[i * l:(i + 1) * l], axis=0))
    groups = array(groups)
    return groups


def groupmeanvec_std(vec, l):
    groups = []
    groups_std = []
    for i in range(int(len(vec) / l) + 1):
        groups.append(mean(vec[i * l:(i + 1) * l], axis=0))
        groups_std.append(std(vec[i * l:(i + 1) * l]))
    groups = array(groups)
    groups_std = array(groups_std)
    return groups, groups_std


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


def p_m_byv_abs(d1, d2, v):
    p = 0
    for x in range(len(d1)):
        if abs(float(d1[x])) < v:
            # print(float(d2[x]/abs(d2[x])))
            p += d2[x] / (abs(d2[x]) + 1.0)
    return p


def p_m_byv_normed(d1, d2, std_d2, v):
    pn = 0
    pp = 0
    nr_p = 0
    nr_m = 0
    for x in range(len(d1)):
        if float(d1[x]) < v and float(d1[x]) >= 0:
            pp += (d2[x] / (abs(d2[x]) + 0.1)) * (1 / std_d2[x])
            nr_p += 1
        if float(d1[x]) > -v and float(d1[x]) <= 0:
            pn += (d2[x] / (abs(d2[x]) + 0.1)) * (1 / std_d2[x])
            nr_m += 1
    return pn, pp


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
        pn = pn / nr_n
    except:
        pass
    return pn, pp


def p_m_byv_limits(d1, d2):
    pn = 0
    pp = 0
    for x in range(len(d1)):
        if float(d1[x]) >= 0:
            pp += d2[x] / (abs(d2[x]) + 0.1)
        if float(d1[x]) <= 0:
            pn += d2[x] / (abs(d2[x]) + 0.1)
    return float(pp - pn)


def earth_DS(user_mean, SD):
    proj = []
    for x in SD:
        proj.append(dot(norm_(user_mean), x))
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


def projUMEAN_users(groupusers, user_mean, SD):
    l = len(SD)
    proj2 = empty(l)
    std_shift_groups = empty(l)
    for i in range(l):
        a = 0
        # gr = normvec(groupusers[i])
        for j in groupusers[i]:
            a += dot(j, norm_(user_mean))  # /abs(dot(j,norm(user_mean)))
        # print(a)
        proj2[i] = a
        std_shift_groups[i] = mean(std(groupusers[i], axis=0))
        a = 0
    # print(std_shift_groups)
    # print('sum d2: ',sum(proj2))
    return proj2, std_shift_groups


def cumulative_limits(user_mean, SD, shifts, cond=2):
    l = len(SD)
    # print('Std. of the shifts:', std(shifts, axis=0))

    d1 = earth_DS(user_mean, SD)  # SD-n angles
    minange = max(abs(d1))
    # print(minange)  # min angle

    cond = 0
    if cond == 0:
        usergrupped = groupvec(shifts, l)
        # d2 = projSD_users(usergrupped,SD)
        d2, std_d2 = projUMEAN_users(usergrupped, user_mean, SD)  # sum(+-1) a csoporton belul
        d2 = array(d2).astype('float64')

    d1 = array(d1).astype('float64')
    d2 = array(d2).astype('float64')
    # print(len(d1),len(d2))

    # nr = 200
    # lims = linspace(0.0,minange,nr)
    # cumulative = empty((nr,2))
    # for i in range(nr):
    # 	v = lims[i]
    # 	cumulative[i] = p_m_byv(d1,d2,v)
    # 	#cumulative[i] = p_m_byv_normed(d1,d2,std_d2,v)
    # 	#cumulative[i] = p_m_byv_abs(d1,d2,v)
    # diff = p_m_byv_limits(d1, d2)
    # diff = p_m_byv(d1,d2,minange)
    diff = p_m_byv_beta(d1, d2, minange)
    # print('Anisotropy: ', diff)
    return diff[0] - diff[1]


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
    # plt.legend()
    # plt.tight_layout()
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


def rot_vector_3D(vector, system, theta, phi):
    v_new_in_system = [math.sin(phi) * math.cos(theta),
                       math.sin(phi) * math.sin(theta),
                       math.cos(phi)]
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R2 = transform_matrix(system, ECEF)
    rotated_vector = R2.dot(v_new_in_system)
    return rotated_vector


def rot_SD_3D(SD_init, s2, s3, rot_theta, rot_phi):
    SD_rotated = empty(shape(SD_init))
    for i in range(len(SD_init)):
        system = flip(gram_schmidt(flip(array([s2[i], s3[i], SD_init[i]]))))
        SD_rotated[i] = rot_vector_3D(SD_init[i], system, rot_theta, rot_phi)
    return SD_rotated


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


def calculate_using_rot_angles(user_mean, shifts, SD, stars_secondary_A, stars_secondary_B, rot_theta, rot_phi):
    rSD = rot_SD_3D(SD, stars_secondary_A, stars_secondary_B, rot_theta, rot_phi)
    i = random.randint(0, len(SD) - 1)
    # print('SD vs. rSD: ', degrees(arccos(dot(SD[i], rSD[i]))))
    return cumulative_limits(user_mean, rSD, shifts)


def process_calc_meancumulative_faster(root_directory, resolution):
    list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]

    theta_max = math.pi * 2.0  # + resolution * 0.1
    phi_max = math.pi
    rot_theta = arange(0.0, theta_max, resolution)
    rot_phi = arange(0.0, phi_max, resolution)
    print(len(rot_theta), len(rot_phi))
    data = []
    for directory in list_subfolders_with_paths:
        total_dir = os.path.join(directory, "allsatellites")
        print(total_dir)
        user_mean, shifts, SD, stars_secondary_A, stars_secondary_B = prepare_data(total_dir)
        cmap = empty((len(rot_theta), len(rot_phi)))
        for j in range(len(rot_phi)):
            print('Phi index: {}/{}'.format(j, len(rot_phi) - 1))
            for i in range(len(rot_theta)):
                theta = rot_theta[i]
                phi = rot_phi[j]
                cmap[i][j] = calculate_using_rot_angles(user_mean, shifts, SD, stars_secondary_A, stars_secondary_B,
                                                        theta, phi)
        data.append(cmap)
        plot_cmap(cmap, phi_max, theta_max, directory)

    mean_data = mean(array(data), axis=0)
    plot_cmap(mean_data, phi_max, theta_max, root_directory)


# =================================================================================================
# =================================================================================================


nr = radians(5.0)
root_directory = r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/recalculated_stars_perth_augusztus"
list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]
# print(list_subfolders_with_paths)

process_calc_meancumulative_faster(root_directory, nr)
