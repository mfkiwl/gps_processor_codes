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




def get_file_name(file):
  f, ext = os.path.splitext(file)
  return f.split("/")[-1]


def get_sats_pos(rinex_processed_grouped):
  signal='C1C'
  all_sats_pos = []
  i=0
  rinex_processed_grouped_B = rinex_processed_grouped#[::10]
  print('Get satellite positions')
  for corr in tqdm(rinex_processed_grouped_B):     #loop over the time/epochs
    flag=['Epoch','index',i]
    all_sats_pos.append(flag)
    for meas in corr:                              #loop over the satellites      
      if signal in meas.observables_final and np.isfinite(meas.observables_final[signal]):
        sat_pos = meas.sat_pos_final
      elif signal in meas.observables and np.isfinite(meas.observables[signal]) and meas.processed:
        sat_pos = meas.sat_pos
      else:
        print('\n\n Satelite position cannot be defined \n\n')
        continue
      all_sats_pos.append(sat_pos)
    i=i+1

  all_sats_pos = pd.DataFrame(np.array(all_sats_pos))
  try:
    all_sats_pos.to_csv(generalinfo_path+'/all_sats_pos.csv', index=False)
  except:
    pass


def calc_user_pos(measurements_in, dataset, resultspath):
  print('Determine position: ', dataset)
  rinex_processed_grouped = measurements_in
  times = []
  ests = []
  for corr in tqdm(rinex_processed_grouped[:]):     #loop over the time/epochs
    try:
      fix, _ = raw.calc_pos_fix(corr)
      ests.append(fix)  
      recv_timee = corr[0].recv_time 
      times.append(recv_timee.tow)
    except:
      pass
    
  ests = np.array(ests)
  times = np.array(times)
    
  try:
    times = pd.DataFrame(times)
    times.to_csv(resultspath+'/times_'+dataset+'.csv', index=False)
    user_positions = pd.DataFrame(ests[:, :3])
    user_positions.to_csv(resultspath+'/user_pos_'+dataset+'.csv', index=False)

    mean_fix = np.mean(ests[:, :3], axis=0)
    mean_user_positions = pd.DataFrame(mean_fix[:3].transpose())
    mean_user_positions.to_csv(resultspath+'/mean_user_pos_'+dataset+'.csv', index=False)  
  except:
    pass


#GET MEASUREMENT DATA
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





obs_path = r"/Users/kelemensz/Qsync/GPS/reciever_data/perth_data/perth"
destination_path = r"/Users/kelemensz/Documents/Research/GPS/process/24h/auto_results/augusztus_"


