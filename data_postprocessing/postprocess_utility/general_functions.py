from numpy import *
import matplotlib.pyplot as plt
import pandas as pd
from vpython import rotate
import pymap3d as pm


def rotateAntiClockwise(array):
    return rotate(array, 90)


def cart2sph(x, y, z):
    hxy = hypot(x, y)
    r = hypot(hxy, z)
    el = arctan2(z, hxy)
    az = arctan2(y, x)
    return az, el


def filter_collected_triangles(all_data, sat_identifier=None):
    sat_data = {}
    for epoch, data in all_data.items():
        for triange in data:
            # if next(iter(triange)) == sat_identifier:
            tr = triange.get(sat_identifier, None)
            if tr:
                sat_data[epoch] = tr[:1]
                # print(tr[:1])
    return sat_data


def filter_collected_triangles_many(all_data, sat_identifiers=None):
    sat_data = {}
    for epoch, data in all_data.items():
        for triange in data:
            if list(triange.keys())[0] in sat_identifiers:
                sat_data[list(triange.keys())[0]][epoch] = list(triange.values())[0]
                # print(tr[:1])
    return sat_data


def sat_data_to_spherical(sat_data):
    for epoch, n in sat_data.items():
        sat_data[epoch] = cart2sph(n[0][0], n[0][1], n[0][2])
    return sat_data


def plot_phi_theta_n(sat_data_spherical, sat_id=None, day=None):
    directions = array(list(sat_data_spherical.values())).T
    epochs = [float(i) for i in sat_data_spherical.keys()]
    plt.scatter(epochs, directions[1], label="phi")
    plt.scatter(epochs, directions[0], label="theta")
    plt.legend()
    if sat_id and day:
        plt.title(day + " " + sat_id)
    plt.show()


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


def get_mean_direction_over_time(systems, directions):
    l = min(len(directions), len(systems[0]))
    phi = 0.
    theta = 0.
    for i in range(l):
        out = cartesian_to_galactic(array([systems[0][i], systems[1][i], systems[2][i]]), directions[i])
        theta += out[0]
        phi += out[1]
    return theta / float(l), phi / float(l)


def normvec(vec):
    a = empty(shape(vec))
    for i in range(len(vec)):
        a[i] = unit(vec[i])
    return a


def scal(v1, v2):
    return vdot(array(v1), array(v2))


def unit(v):
    return 1 / sqrt(scal(v, v)) * array(v)


def dr_lenght(v1, v2):
    return sqrt((v1[0] - v2[0]) ** 2 + (v1[1] - v2[1]) ** 2 + (v1[2] - v2[2]) ** 2)


def transform_matrix(f1, f2):  # transforms from f1 to f2
    R = array([
        [dot(f2[0], f1[0]), dot(f2[0], f1[1]), dot(f2[0], f1[2])],
        [dot(f2[1], f1[0]), dot(f2[1], f1[1]), dot(f2[1], f1[2])],
        [dot(f2[2], f1[0]), dot(f2[2], f1[1]), dot(f2[2], f1[2])]
    ])
    return R


def get_theta_phi(v):
    if 0.99 < v[0] and 1.01 > v[0] and -0.01 < v[1] and 0.01 > v[1] and -0.01 < v[2] and 0.01 > v[2]:
        return 0.0, 0.0
    if 0.99 < v[1] and 1.01 > v[1] and -0.01 < v[0] and 0.01 > v[0] and -0.01 < v[2] and 0.01 > v[2]:
        return 90.0, 0.0
    if 0.99 < v[2] and 1.01 > v[2] and -0.01 < v[1] and 0.01 > v[1] and -0.01 < v[0] and 0.01 > v[0]:
        return 0.0, 90.0
    return None, None


def ecef_to_gcs(system, vector):
    vector = unit(vector)
    system = normvec(system)
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    R = transform_matrix(ECEF, system)
    return around(R.dot(vector), decimals=4)


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


def cartesian_to_spherical(vector):
    theta, phi = get_theta_phi(vector)
    ECEF = array(([1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]))
    if not theta and not phi:
        phi, theta, _ = pm.ecef2geodetic(vector[0], vector[1], vector[2])
        phi = 90 - degrees(arccos(scal(ECEF[2], vector)))
    return radians(theta), radians(phi)


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


def get_ij_on_map(direction, resolution):
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
        pass
    return False


def get_mean_pos_from_root(root_path, positions_file, max_deviations=0.5):
    positions = []
    month_folders = [f.path for f in os.scandir(root_path) if f.is_dir()]
    for month_root in month_folders:
        day_folders = [f.path for f in os.scandir(month_root) if f.is_dir()]
        for day_root in day_folders:
            day_root = os.path.join(day_root, "allsatellites")
            if is_all_data(day_root, [positions_file], add_allsatellites=False):
                mean_pos = get_mean_position(day_root, positions_file)
                # if str(os.path.split(day_root)[-2]).split("/")[-1][:4] in ['NZLD', 'NZDL']:
                if str(os.path.split(day_root)[-2]).split("/")[-1][:4] not in ['BLUF', 'HOKI', 'MAVL', 'LKTA', 'MTJO']:
                    # print(str(os.path.split(day_root)[-2]).split("/")[-1], mean_pos)
                    positions.append(mean_pos)
    positions = array(positions)
    std_pos = std(positions, axis=0)
    print("Number of days with positions determined and std of the positions before filter:\n", len(positions), std_pos)

    std_pos_norm = sqrt(std_pos.dot(std_pos))
    mean_ = mean(positions, axis=0)
    distance_from_mean = nrm(positions - mean_, axis=1)
    not_outlier = distance_from_mean < max_deviations * std_pos_norm
    # print(not_outlier)
    no_outliers = positions[not_outlier]
    print("After filter: ", len(no_outliers), std(array(no_outliers), axis=0), "\n ")
    # print(mean(array(positions), axis=0), mean(array(no_outliers), axis=0))
    return mean(array(no_outliers), axis=0)


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


def create_dir(root_path, dir_name):
    results_dir = os.path.join(root_path, dir_name)
    if not os.path.isdir(results_dir):
        os.makedirs(results_dir)
    return results_dir

