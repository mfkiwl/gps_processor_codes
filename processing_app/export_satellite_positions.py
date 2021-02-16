import numpy as np
import pandas as pd
from tqdm import tqdm
from laika import constants, AstroDog
from laika.lib.coordinates import ecef2geodetic, geodetic2ecef
from laika.rinex_file import RINEXFile
import laika.raw_gnss as raw
from laika.raw_gnss import process_measurements, correct_measurements, calc_pos_fix
from skyfield.api import Star, load, Topos
from skyfield.data import hipparcos
from skyfield.positionlib import *
import os
from datetime import datetime


def get_file_name(file):
    f, ext = os.path.splitext(file)
    return f.split("/")[-1]


def get_sats_pos_and_pr(rinex_processed_grouped, generalinfo_path):
    signal = 'C1C'
    all_sats_pos = []
    i = 0
    # rinex_processed_grouped_B = rinex_processed_grouped#[::10]
    print('Get satellite positions')
    dat_from_if = 0
    dat_from_else = 0
    for corr in tqdm(rinex_processed_grouped):  # loop over the time/epochs
        recv_timee = corr[0].recv_time.tow
        flag = ['Epoch', 'index', i, recv_timee]
        all_sats_pos.append(flag)
        for meas in corr:  # loop over the satellites
            if signal in meas.observables_final and np.isfinite(meas.observables_final[signal]):
                pr = meas.observables_final[signal]
                sat_pos = meas.sat_pos_final
                dat_from_if += 1
            elif signal in meas.observables and np.isfinite(meas.observables[signal]) and meas.processed:
                pr = meas.observables[signal]
                pr += meas.sat_clock_err * constants.SPEED_OF_LIGHT
                sat_pos = meas.sat_pos
                dat_from_else += 1
            else:
                print('\n\n Satelite position cannot be defined \n\n')
                continue
            sat_data = np.array([sat_pos[0], sat_pos[1], sat_pos[2], pr])
            all_sats_pos.append(sat_data)
        i = i + 1
    print("If met {}, if did not meet: {}".format(dat_from_if, dat_from_else))
    all_sats_pos = pd.DataFrame(np.array(all_sats_pos))
    try:
        all_sats_pos.to_csv(generalinfo_path + '/all_sats_pos_time.csv', index=False)
    except:
        pass


def get_sats_pos_and_pr_beta(rinex_processed_grouped, generalinfo_path):
    signal = 'C1C'
    all_sats_pos = []
    i = 0
    # rinex_processed_grouped_B = rinex_processed_grouped#[::10]
    print('Get satellite positions')
    dat_from_if = 0
    dat_from_else = 0
    sat_ident = []
    for corr in tqdm(rinex_processed_grouped):  # loop over the time/epochs
        recv_timee = corr[0].recv_time.tow
        flag = ['Epoch', 'index', i, recv_timee, '-']
        all_sats_pos.append(flag)
        for meas in corr:  # loop over the satellites
            # if signal in meas.observables_final and np.isfinite(meas.observables_final[signal]):
            #     pr = meas.observables_final[signal]
            #     sat_pos = meas.sat_pos_final
            #     dat_from_if += 1
            # elif signal in meas.observables and np.isfinite(meas.observables[signal]) and meas.processed:
            #     pr = meas.observables[signal]
            #     pr += meas.sat_clock_err * constants.SPEED_OF_LIGHT
            #     sat_pos = meas.sat_pos
            #     dat_from_else += 1
            # else:
            #     print('\n\n Satelite position cannot be defined \n\n')
            #     continue
            try:
                prn = meas.prn
                pr = meas.observables[signal]
                pr += meas.sat_clock_err * constants.SPEED_OF_LIGHT
                sat_pos = meas.sat_pos
            except:
                print('\n\n Satelite position cannot be defined \n\n')
                continue

            # sat_ident.append(prn)
            sat_data = np.array([sat_pos[0], sat_pos[1], sat_pos[2], pr, prn])
            all_sats_pos.append(sat_data)
        i = i + 1
    # print("If met {}, if did not meet: {}".format(dat_from_if, dat_from_else))
    # print("Satellite identifiers: ", list(set(sat_ident)))
    all_sats_pos = pd.DataFrame(np.array(all_sats_pos))
    try:
        all_sats_pos.to_csv(generalinfo_path + '/sats_pos_time_id.csv', index=False)
    except:
        pass


# GET MEASUREMENT DATA
def get_processed_data(filepath):
    dog = AstroDog()
    print('Preprocessing: ', filepath)
    obs_data = RINEXFile(filepath)
    rinex_meas_grouped = raw.read_rinex_obs(obs_data)
    del obs_data
    rinex_processed_grouped = []
    for meas in tqdm(rinex_meas_grouped):
        proc = raw.process_measurements(meas, dog=dog)
        rinex_processed_grouped.append(proc)
    print('Data is IN!')
    del rinex_meas_grouped
    return rinex_processed_grouped


def get_obs_files(directory):
    obs_ext = ".obs"
    obs_files = []
    files = next(os.walk(directory))[2]
    for file in files:
        if os.path.splitext(file)[1] == obs_ext and len(os.path.splitext(file)[0]) == 14:
            obs_files.append(os.path.abspath(os.path.join(directory, file)))
    return obs_files


def create_dir_for_results(obs_file, new_root_directory):
    obs_filename = get_file_name(obs_file)
    new_dir_name = obs_filename[:-2]
    new_dir = os.path.join(new_root_directory, new_dir_name)
    if not os.path.isdir(new_dir):
        os.mkdir(new_dir)
    return new_dir


def create_generalinfo_path(directory):
    generalinfo_path = os.path.join(directory, 'allsatellites')
    if not os.path.exists(generalinfo_path):
        os.makedirs(generalinfo_path)
    return generalinfo_path


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


def process_many_from_obs_root(obs_location, root_directory, needed_files, month_names):
    obs_files = get_obs_files(obs_location)
    month_name = os.path.split(obs_location)[-1]
    root_directory = os.path.join(root_directory, month_name)
    total_processed = 0
    for obs_file in obs_files:  # [:1]:
        obs_file_name = os.path.splitext(os.path.split(obs_file)[-1])[0][:-2]
        if len(obs_file_name) == 12:  # and obs_file_name == "NASA20200102":
            print("Obs filename: ", obs_file_name)
            path_of_results = create_dir_for_results(obs_file, root_directory)
            generalinfo_path = create_generalinfo_path(path_of_results)
            rinex_processed_grouped = get_processed_data(obs_file)

            get_sats_pos_and_pr_beta(rinex_processed_grouped, generalinfo_path)
            total_processed += 1
            print("\n                          Processed: {}/{} \n".format(total_processed, len(obs_files)))


#                    NASA
# obs_path = r"/Volumes/KingstonSSD/GPS/NASA/obs_files/januar"
# destination_path = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_NASA"

#                   PERTH
# obs_path = r"/Volumes/ADATA SE800/GPS/raw_data/PERTH/rinex_obs/januar"
# destination_path = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/PERTH_daily_measurements"

#                   HKKS
# obs_path = r"/Volumes/KingstonSSD/GPS/HKG/HKKS/obs_files/januar"
# destination_path = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_HKKS"

#                   NZLD
# obs_path = r"/Volumes/ADATA SE800/GPS/raw_data/NZ/obs_files/ARTA_smaller_obs/januar"
# destination_path = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NZLD"

#                   IIGC
obs_path = r"/Volumes/KingstonSSD/GPS/INDIA/IIGC/obs_files_IIGC/januar"
destination_path = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_IIGC"


needed_files = ["user_pos_allsatellites.csv", "all_sats_pos_time.csv"]
month_names = ["julius", "szeptember", "februar", "marcius", "augusztus", "januar", "december2019", "oktober",
               "november", "majus", "aprilis", "junius", "december2020"]

process_many_from_obs_root(obs_path, destination_path, needed_files, month_names)
