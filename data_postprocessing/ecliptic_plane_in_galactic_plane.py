
from numpy import *
import pandas as pd
import matplotlib.pyplot as plt
import pymap3d as pm
import os
import pylab as pl
import math
import sys

sys.path.insert(1, '/Users/kelemensz/Documents/Research/GPS/gps_processor_codes/utility')
from frecvently_used_functions import plot_on_sphere

sys.path.insert(1, '/Users/kelemensz/Documents/Research/GPS/gps_processor_codes/data_postprocessing')
from data_locations_handler import AllGPSDataLocations
from postprocess_utility.general_functions import get_mean_direction_over_time
from triangle_method_rework.rework_lib import get_global_stars, Defaults, ecef_to_gcs, \
    get_ij_on_map, create_dir, rotate2darray, plot_mollweid, normvec

# from data_postprocessing.triangle_method_rework.triangle_with_two_satellites import process_raw_GCS_data

s_o_l = 1.0  # 3.0*10**8
# satellite_positions = "all_sats_pos_time.csv"
satellite_positions = "sats_pos_time_id.csv"


# =====================================================================================================================


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


# =================================================================================================
# =================================================================================================


def phi_theta_map_from_directions(vectors, resolution):
    theta_max = math.pi
    phi_max = math.pi / 2.0
    rot_theta = arange(-theta_max, theta_max, resolution)
    rot_phi = arange(-phi_max, phi_max, resolution)
    cmap = zeros((len(rot_theta), len(rot_phi)))
    value = 10
    for direction in vectors:
        i, j = get_ij_on_map(direction, resolution)
        cmap[i][j] += value

    cmap[cmap < 1] = 0
    return cmap

def plot_mollweid_(matrix, star_directions, anot=True):
    matrix = nan_to_num(matrix, nan=0.0)
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
    pl.show()

    plt.clf()
    pl.clf()


def add_star_annotated(theta, phi, name, ax):
    theta = radians(theta)
    phi = radians(phi)
    ax.text(theta, phi, name, fontsize=12)
    ax.scatter(theta, phi, marker='x', c='k', s=15)


def plot_phi_theta_n(sat_data_spherical, sat_id=None, day=None):
    directions = array(list(sat_data_spherical.values())).T
    epochs = [float(i) for i in sat_data_spherical.keys()]
    plt.scatter(epochs, directions[1], label="phi")
    plt.scatter(epochs, directions[0], label="theta")
    plt.legend()
    if sat_id and day:
        plt.title(day + " " + sat_id)
    plt.show()


def test_one_day_sun_dirrections(star_dir, resolution, day=None):
    Nx, Ny, Nz, S, D, sun = get_global_stars(star_dir, day)
    print(sun)
    l = len(Nx)
    star_directions_in_GCS = {'GNP': get_mean_direction_over_time(array([Nx, Ny, Nz]), Nz),
                              'GC': get_mean_direction_over_time(array([Nx, Ny, Nz]), Nx),
                              'OY of GCS': get_mean_direction_over_time(array([Nx, Ny, Nz]), Ny),
                              'Sadalmelic': get_mean_direction_over_time(array([Nx, Ny, Nz]), S),
                              'Denebola': get_mean_direction_over_time(array([Nx, Ny, Nz]), D)}
    # =================================================================================================================
    GCS_all = [array([Nx[i], Ny[i], Nz[i]]) for i in range(l)]
    sun = array(sun)
    sun_GCS = []
    for i in range(l):
        sun_GCS.append(ecef_to_gcs(GCS_all[i], sun[i]))
    # print(sun_GCS)

    cmap = phi_theta_map_from_directions(sun_GCS, resolution)
    # cmap = nan_to_num(cmap, nan=0.0)
    # plt.imshow(cmap)
    # plt.colorbar()
    # plt.show()

    plot_mollweid(cmap, star_directions_in_GCS)

    return 0, 0, 0


# =================================================================================================
# =================================================================================================


# star_dir = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020"
# resolution = radians(5.0)
# dirname_date = 'date20200102'
#
# test_one_day_sun_dirrections(star_dir, resolution, dirname_date)


# =================================================================================================
# =================================================================================================
# ======================================USER TO GCS OVER ONE DAY===================================
# =================================================================================================
# =================================================================================================

def transform_GCS_data_to_cmap(directions_GCS, resolution):
    theta_max = math.pi
    phi_max = math.pi / 2.0
    rot_theta = arange(-theta_max, theta_max, resolution)
    rot_phi = arange(-phi_max, phi_max, resolution)
    l_theta = int(len(rot_theta) / 2.0)
    l_phi = int(len(rot_phi) / 2.0)
    cmap_count = zeros((len(rot_theta), len(rot_phi)))
    for direction in directions_GCS:
        i, j = get_ij_on_map(direction, l_theta, l_phi, resolution)
        cmap_count[i][j] += 1
    cmap_count[cmap_count < 1] = 0
    return cmap_count


def test_one_day_user_dirrections(star_dir, resolution, day_directory):
    Nx, Ny, Nz, S, D = get_global_stars(star_dir, day_directory)
    l = len(Nx)
    star_directions_in_GCS = Defaults.get_star_directions_in_GCS(Nx, Ny, Nz, D, S)
    # =================================================================================================================
    GCS_all = [[Nx[i], Ny[i], Nz[i]] for i in range(l)]
    # sun = array(sun)
    # sun_GCS = []
    user_ = pd.read_csv(os.path.join(day_directory, 'allsatellites',
                                     Defaults.USER_POSITIONS_FILENAME.get('user')), skiprows=1).values

    user_ = normvec(user_)
    group_length = int(len(user_)/l)
    vectors_GCS = []
    kk = 0
    for i in range(0, l-1):
        for k in range(i*group_length, (i+1)*group_length):
            # print('i: ', i, ' k: ', kk, group_length)
            out = ecef_to_gcs(array([GCS_all[i][0], GCS_all[i][1], GCS_all[i][2]]), user_[kk])
            vectors_GCS.append(out)
            kk += 1
    cmap = transform_GCS_data_to_cmap(vectors_GCS, resolution)
    cmap = nan_to_num(cmap, nan=0.0)
    # plt.imshow(cmap)
    # plt.colorbar()
    # plt.show()

    # plot_mollweid(cmap, star_directions_in_GCS)

    return cmap


def find_days_and_process(user_sat_data_root, result_root, star_dir, resolution, month_names=AllGPSDataLocations.all_months):
    d = 0
    all_hist = []
    if os.path.isdir(user_sat_data_root) and os.path.isdir(result_root):

        months = [f.path for f in os.scandir(user_sat_data_root) if f.is_dir()]

        for month in months[:]:
            month_name = os.path.split(month)[-1]
            # condition = month_name in month_names
            condition = d < 1
            if condition:
                day_folders = [f.path for f in os.scandir(month) if f.is_dir()]
                print("Number of days: ", month_name, len(day_folders))
                for day_folder in day_folders[:2]:
                    try:

                        date = str(os.path.split(day_folder)[-1])[-8:]
                        hist = test_one_day_user_dirrections(star_dir, resolution, day_folder)
                        if hist is not None:
                            result_month = create_dir(result_root, month_name)
                            result_day = create_dir(result_month, date)
                            # plot_mollweid(hist, '', result_day, "histogram", str(int(degrees(resolution))), anot=False)
                            all_hist.append(hist)
                            print('Date: ', date)
                            d += 1
                    except:
                        print("Problem with: \n{}".format(day_folder))
        print("Total nr of days: ", d)
        all_hist = sum(array(all_hist), axis=0)
        try:
            all_hist = rotate2darray(nan_to_num(all_hist, nan=0.0))
        except:
            pass
        all_hist[all_hist < 1] = 0.0
        all_hist = nan_to_num(all_hist, nan=0)
        print('amax: ', amax(all_hist))
        plot_on_sphere(all_hist.T)
        # plt.imshow(all_hist)
        # plt.colorbar()
        # plt.savefig(os.path.join(result_root, 'total_hist_{}.png'.format(d)))
        # # plt.show()
        plt.clf()
        # plot_mollweid(all_hist, '', result_root, "histogram", str(int(degrees(resolution))), anot=False)




star_dir = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020"
nr = radians(5.0)

# directory = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements/januar/CUTB20200111"
# directory = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/IIGC"

# dest_dir = r"/Users/kelemensz/Documents/Research/GPS/process/show_user_inGCS/CUTB"
# test_one_day_user_dirrections(star_dir, nr, directory, dest_dir)
# find_days_and_process(directory, dest_dir, star_dir, nr)

for ID, location in AllGPSDataLocations.user_and_satellites.items():
    dest_dir = r"/Users/kelemensz/Documents/Research/GPS/process/show_user_inGCS"
    # if str(ID) not in ['CUTB', 'PERTH']:
    if str(ID) in ['IIGC']:
        print('--------------------------         ', ID, '           --------------------------\n')
        # try:
        dest_dir = create_dir(dest_dir, ID)
        find_days_and_process(location, dest_dir, star_dir, nr)
        # except:
        #     continue