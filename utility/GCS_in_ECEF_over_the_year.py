import numpy as np
import pandas as pd
from tqdm import tqdm
from skyfield.api import Star, load, Topos
from skyfield.data import hipparcos
from skyfield.positionlib import *
import pymap3d as pm
from vpython import *
import os
from datetime import datetime, timedelta, date

# S1=57632  # Denebola
# S2=109074  # Sadalmelik
# S1 = 15863  # Mirphak
# S2 = 91262  # Vega
# geographic_position = [46.7502933, 23.5835082, 6367740.701600294]  # Kolozsvar
# geographic_position = [46.083, 26.233, 6367740.701600294]  # Almas
# geographic_position = [-31.953512, 115.857048, 6367740.701600294]  # Perth
# geographic_position = [-38.50441, 175.99478, 6363911.3053]  # New Zeland
# geographic_position = [36.21598, 127.37449, 6364405.5371]  # Daejeon, Del Korea, NASA
geographic_position = [51.477928, -0.001545, 6367740.701600294]  # Greenwich
time_diff = 0.0


def format_time_auto_from_days_sec(date, sec_of_day, time_difference=time_diff):
    terrestrial_time = np.zeros(6)
    year = int(date[0])  # adjust acording to the current year
    month = int(date[1])  # adjust acording to the current month
    day = int(date[0])  # adjust acording to the first day of the curent week
    h_of_day = int(sec_of_day / 3600.0)
    sec_of_h = int(sec_of_day % 3600.0)
    min_of_h = int(sec_of_h / 60.0)
    sec_of_min = int(sec_of_h % 60.0)
    terrestrial_time[0] = year
    terrestrial_time[1] = month
    terrestrial_time[2] = day
    terrestrial_time[3] = h_of_day + time_difference  # (shift the time because of some unclear reasons)
    terrestrial_time[4] = min_of_h
    terrestrial_time[5] = sec_of_min
    return terrestrial_time


def galactic_Sun_north_center_SD_directions(measurement_date, sec_of_day, user_longlatdist=geographic_position,
                                            prints=False):
    user_longlatdist = np.array(user_longlatdist)
    with load.open(hipparcos.URL) as f:
        df = hipparcos.load_dataframe(f)

    planets = load('de421.bsp')
    earth = planets['earth']
    user_position = earth + Topos('{} N'.format(user_longlatdist[0]), '{} E'.format(
        user_longlatdist[1]))  # user_longlatdist      user_longlatdist[0],user_longlatdist[1]
    terrestrial_time = format_time_auto_from_days_sec(measurement_date, sec_of_day)
    ts = load.timescale()
    t = ts.tt(terrestrial_time[0], terrestrial_time[1], terrestrial_time[2], terrestrial_time[3], terrestrial_time[4],
              terrestrial_time[5])

    sun = planets['sun']
    astrometric_sun = user_position.at(t).observe(sun)
    alt_sun, az_sun, distance_sun = astrometric_sun.apparent().altaz()  # astrometric.apparent().altaz() = astrometric_denebola.apparent().altaz()
    sun_position = pm.aer2ecef(az_sun.degrees, alt_sun.degrees, distance_sun.m, user_longlatdist[0],
                               user_longlatdist[1], user_longlatdist[2])

    # Északi Galaktikus Pólus (RA: 12h 51.42m; D: 27° 07.8′, epocha J2000.0)
    coma_berenices = Star(ra_hours=(12, 51, 25.2), dec_degrees=(27, 7, 48.0))
    astrometric_berenice = user_position.at(t).observe(coma_berenices)
    alt_b, az_b, distance_b = astrometric_berenice.apparent().altaz()
    berenice_position = pm.aer2ecef(az_b.degrees, alt_b.degrees, distance_b.m, user_longlatdist[0], user_longlatdist[1],
                                    user_longlatdist[2])

    # Galaktikus Center (RA: 12h 51.42m; D: 27° 07.8′, epocha J2000.0) 17h 45.6m  −28.94°
    galactic_center = Star(ra_hours=(17, 45, 36.0), dec_degrees=(-28, 56, 24.0))
    astrometric_center = user_position.at(t).observe(galactic_center)
    alt_c, az_c, distance_c = astrometric_center.apparent().altaz()
    center_position = pm.aer2ecef(az_c.degrees, alt_c.degrees, distance_c.m, user_longlatdist[0], user_longlatdist[1],
                                  user_longlatdist[2])

    nunki = Star.from_dataframe(df.loc[92855])
    astrometric_nunki = user_position.at(t).observe(nunki)
    alt_n, az_n, distance_n = astrometric_nunki.apparent().altaz()  # astrometric.apparent().altaz() = astrometric_denebola.apparent().altaz()
    star_pos_n = pm.aer2ecef(az_n.degrees, alt_n.degrees, distance_n.m, user_longlatdist[0], user_longlatdist[1],
                             user_longlatdist[2])

    capella = Star.from_dataframe(df.loc[24608])
    astrometric_capella = user_position.at(t).observe(capella)
    alt_ca, az_ca, distance_ca = astrometric_capella.apparent().altaz()  # astrometric.apparent().altaz() = astrometric_denebola.apparent().altaz()
    star_pos_ca = pm.aer2ecef(az_ca.degrees, alt_ca.degrees, distance_ca.m, user_longlatdist[0], user_longlatdist[1],
                              user_longlatdist[2])

    denebola = Star.from_dataframe(df.loc[57632])
    astrometric_denebola = user_position.at(t).observe(denebola)
    alt_d, az_d, distance_d = astrometric_denebola.apparent().altaz()  # astrometric.apparent().altaz() = astrometric_denebola.apparent().altaz()
    star_pos_d = pm.aer2ecef(az_d.degrees, alt_d.degrees, distance_d.m, user_longlatdist[0], user_longlatdist[1],
                             user_longlatdist[2])

    sadalmelik = Star.from_dataframe(df.loc[109074])
    astrometric_sadalmelik = user_position.at(t).observe(sadalmelik)
    alt_s, az_s, distance_s = astrometric_sadalmelik.apparent().altaz()
    star_pos_s = pm.aer2ecef(az_s.degrees, alt_s.degrees, distance_s.m, user_longlatdist[0], user_longlatdist[1],
                             user_longlatdist[2])
    if prints:
        print('Time:  ', terrestrial_time)
        print('Location: ', '{} N'.format(user_longlatdist[0]), '{} E'.format(user_longlatdist[1]))
        print('Sun (alt/az):   ', alt_sun, az_sun)
        print('Coma Berenice (alt/az):   ', alt_b, az_b)
        print('Galactic Center (alt/az):   ', alt_c, az_c)
        print('Nunki (alt/az):   ', alt_n, az_n)
        print('Capella (alt/az):   ', alt_ca, az_ca)
        print('Denebola (alt/az):   ', alt_d, az_d)
        print('Sadalmelik (alt/az):   ', alt_s, az_s)

    stars_pos = [sun_position, berenice_position, center_position, star_pos_n, star_pos_ca, star_pos_d, star_pos_s]
    return stars_pos


def get_GCS_stars_in_time(day_path, date, star_check_per_day):
    n_epochs = 86400
    star_check_distance = int(n_epochs / float(star_check_per_day))

    sun_pos = []
    northpole_pos = []
    center_pos = []
    nunki_pos = []
    capella_pos = []
    denebola_pos = []
    sadalmelik_pos = []

    star_checks_sec = range(0, n_epochs, star_check_distance)
    print('Determining galactic system orientation!')
    for sec_of_day in tqdm(star_checks_sec):  # loop over the time/epochs
        star_positions = galactic_Sun_north_center_SD_directions(date, sec_of_day)  # OK, ECEF coordinates

        sun_pos.append(np.array(star_positions[0]))
        northpole_pos.append(np.array(star_positions[1]))

        center_pos.append(np.array(star_positions[2]))
        nunki_pos.append(np.array(star_positions[3]))
        capella_pos.append(np.array(star_positions[4]))

        denebola_pos.append(np.array(star_positions[5]))
        sadalmelik_pos.append(np.array(star_positions[6]))
    galactic_Sun_north_center_SD_directions(date, recv_timee, prints=True)
    try:
        pd.DataFrame(np.array(sun_pos)).to_csv(generalinfo_path + '/SUN_positions.csv', index=False)
    except:
        pass

    try:
        pd.DataFrame(np.array(northpole_pos)).to_csv(generalinfo_path + '/GNP_positions.csv', index=False)
    except:
        pass

    try:
        pd.DataFrame(np.array(center_pos)).to_csv(generalinfo_path + '/GC_positions.csv', index=False)
    except:
        pass

    try:
        pd.DataFrame(np.array(nunki_pos)).to_csv(generalinfo_path + '/Nunki_positions.csv', index=False)
    except:
        pass

    try:
        pd.DataFrame(np.array(capella_pos)).to_csv(generalinfo_path + '/Capella_positions.csv', index=False)
    except:
        pass

    try:
        pd.DataFrame(np.array(denebola_pos)).to_csv(generalinfo_path + '/Denebola_positions.csv', index=False)
    except:
        pass

    try:
        pd.DataFrame(np.array(sadalmelik_pos)).to_csv(generalinfo_path + '/Sadalmelik_positions.csv', index=False)
        print("\n Files Saved! \n")
    except:
        pass


def day_folder_date(root_dir, year, day_of_year, prefix):
    day_of_year = str(day_of_year)
    day_of_year = '{}'.format(day_of_year[1:] if day_of_year.startswith('0') else day_of_year)
    day_of_year = int('{}'.format(day_of_year[1:] if day_of_year.startswith('0') else day_of_year))
    date = datetime(year, 1, 1) + timedelta(day_of_year - 1)
    month = '0{}'.format(date.month) if len(str(date.month)) == 1 else date.month
    day = '0{}'.format(date.day) if len(str(date.day)) == 1 else date.day
    day_folder_name = prefix + str(date.year) + str(month) + str(day)
    day_folder = os.path.join(root_dir, day_folder_name)
    if not os.path.isdir(day_folder):
        os.makedirs(day_folder)
    return day_folder, [year, date.month, date.day]


def daterange(start_date, end_date):
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)


def create_root_dir(root_directory):
    if not os.path.isdir(root_directory):
        os.mkdir(root_directory)
    return root_directory


def create_day_path(directory, date):
    generalinfo_path = os.path.join(directory, date)
    if not os.path.exists(generalinfo_path):
        os.makedirs(generalinfo_path)
    return generalinfo_path


destination_path = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/"


def get_stars_ECEF(root_location, year):
    root_directory = create_root_dir(os.path.join(root_location, "STARS_{}".format(year)))
    # start_date = date(year, 1, 1)
    # end_date = date(year, 12, 31)
    # for single_date in daterange(start_date, end_date):
    # print(single_date.strftime("%Y-%m-%d"))
    for day_of_year in range(1, 3, 1):
        path_of_day, date = day_folder_date(root_directory, year, day_of_year, "STAR")
        print(date)

        get_GCS_stars_in_time(path_of_day, date, 87)


get_stars_ECEF(destination_path, 2020)
