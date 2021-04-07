import os
from numpy import *
import pandas as pd


def get_positions_std(path):
    path = os.path.join(path, "allsatellites")

    if os.path.isfile(path + '/user_pos_allsatellites.csv'):
        user_ = pd.read_csv(path + '/user_pos_allsatellites.csv', skiprows=1).values  # .transpose()
        std_pos = std(user_, axis=0)
        print(path.split('/')[-2])
        print("STD position: ", len(user_), std_pos)



def get_subdirs(root):
    return [os.path.join(root, o) for o in os.listdir(root)
            if os.path.isdir(os.path.join(root, o))]


def calc_for_all_month(root):
    days = get_subdirs(root)
    for day in days:
        get_positions_std(day)


CUTA = r'/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements/CUTA/julius'
CUTB = r'/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements/januar'
NZLD = r'/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_NZLD/januar'
IISC = r'/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_IIGC/januar'
NASA = r'/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_NASA/januar'
HKKS = r'/Volumes/BlueADATA S/GPS/processed_data/global_GCS_axis/process_HKKS/januar'
TIDV = r'/Volumes/BlueADATA S/GPS/processed_data/global_GCS_axis/process_TIDV/aprilis'


calc_for_all_month(CUTA)
