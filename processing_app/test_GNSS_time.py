import numpy as np
import pandas as pd
from tqdm import tqdm
from laika import constants, AstroDog
from laika.lib.coordinates import ecef2geodetic,geodetic2ecef
from laika.rinex_file import RINEXFile
import laika.raw_gnss as raw
from laika.raw_gnss import process_measurements, correct_measurements, calc_pos_fix
from skyfield.api import Star, load, Topos
from skyfield.data import hipparcos
from skyfield.positionlib import *
import pymap3d as pm
from vpython import *
import os
from datetime import datetime, timedelta


# S1=57632  # Denebola
# S2=109074  # Sadalmelik
# S1 = 15863  # Mirphak
# S2 = 91262  # Vega
geographic_position = [46.7502933, 23.5835082, 6367740.701600294]  # Kolozsvar
# geographic_position = [46.083, 26.233, 6367740.701600294]  # Almas
# geographic_position = [-31.953512, 115.857048, 6367740.701600294]  # Perth
# geographic_position = [-38.50441, 175.99478, 6363911.3053]  # New Zeland
# geographic_position = [36.21598, 127.37449, 6364405.5371]  # Daejeon, Del Korea, NASA
# geographic_position = [51.477928, -0.001545, 6367740.701600294]  # Greenwich
time_diff = 0




def get_measurements_time(rinex_processed_grouped):
  all_meas_time = []
  print('Save measurements time (tow)!')
  for corr in tqdm(rinex_processed_grouped):     #loop over the time/epochs
    recv_timee = corr[0].recv_time
    all_meas_time.append(recv_timee)
  all_meas_time = pd.DataFrame(np.array(all_meas_time))
  print(all_meas_time)



#GET MEASUREMENT DATA
def get_processed_data(filepath):
    dog = AstroDog()
    print('Preprocessing: ', filepath)
    obs_data = RINEXFile(filepath)
    rinex_meas_grouped = raw.read_rinex_obs(obs_data) 
    del obs_data 
    rinex_processed_grouped = []
    for meas in tqdm(rinex_meas_grouped[:100]):
        proc = raw.process_measurements(meas, dog=dog)
        rinex_processed_grouped.append(proc)
    print('Data is IN!')
    del rinex_meas_grouped
    return rinex_processed_grouped




# obs_file = r"/Users/kelemensz/Qsync/GPS/reciever_data/NZLD/obs/to_process/januar/NZLD2020010324.obs"
obs_file = r"/Users/kelemensz/Qsync/GPS/reciever_data/24h_F/CUT02020071124.obs"



# rinex_processed_grouped = get_processed_data(obs_file)

# get_sun_north_position_in_time(rinex_processed_grouped, generalinfo_path, 87)
# calc_user_pos(rinex_processed_grouped,'allsatellites',generalinfo_path)
# get_measurements_time(rinex_processed_grouped)

# date = array([2020, 1, 3])
# get_GCS_axis_in_time(rinex_processed_grouped, None, 10, date)

obs_files = [r"/Users/kelemensz/Qsync/GPS/reciever_data/NZLD/obs/to_process/MTJO/MTJO2020010224.obs",
r"/Users/kelemensz/Qsync/GPS/reciever_data/downloads_from_nasa/obs_files/done_positions/NASA2020010224.obs",
r"/Users/kelemensz/Qsync/GPS/reciever_data/NZLD/obs/done/januar/NZLD2020010224.obs",
r"/Users/kelemensz/Qsync/GPS/reciever_data/NZLD/obs/to_process/obs/julius/NZLD2020071124.obs",
r"/Users/kelemensz/Qsync/GPS/reciever_data/24h_F/CUT02020071124.obs",
r"/Users/kelemensz/Qsync/GPS/reciever_data/perth_data/perth/CUTA2020071124.obs"
]

for obs in obs_files:
  print(obs.split("/")[-1])
  rinex_processed_grouped = get_processed_data(obs)
  get_measurements_time(rinex_processed_grouped)
  del rinex_processed_grouped
  print("\n")
  print("\n")
  print("\n")





