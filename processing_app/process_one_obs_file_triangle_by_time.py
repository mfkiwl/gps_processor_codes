import numpy as np
import pandas as pd
from tqdm import tqdm
from laika import constants, AstroDog
# from laika.lib.coordinates import ecef2geodetic,geodetic2ecef
from laika.rinex_file import RINEXFile
import laika.raw_gnss as raw
from laika.raw_gnss import process_measurements, correct_measurements, calc_pos_fix
# from skyfield.api import Star, load, Topos
# from skyfield.data import hipparcos
# from skyfield.positionlib import *
# import pymap3d as pm
# from vpython import *
import os
# from datetime import datetime, timedelta
import sys


def get_file_name(file):
  f, ext = os.path.splitext(file)
  return f.split("/")[-1]



def get_measurements_time(rinex_processed_grouped, generalinfo_path):
  all_meas_time = []
  print('Save measurements time (tow)!')
  for corr in tqdm(rinex_processed_grouped):     #loop over the time/epochs
    recv_timee = corr[0].recv_time
    all_meas_time.append(recv_timee)
  all_meas_time = pd.DataFrame(np.array(all_meas_time))
  try:
    all_meas_time.to_csv(generalinfo_path+'/times_allsatellites.csv', index=False)
  except:
    pass



def get_sats_pos_and_pr(rinex_processed_grouped, generalinfo_path):
  signal='C1C'
  all_sats_pos = []
  i=0
  # rinex_processed_grouped_B = rinex_processed_grouped#[::10]
  print('Get satellite positions')
  for corr in tqdm(rinex_processed_grouped):     #loop over the time/epochs
    recv_timee = corr[0].recv_time.tow
    flag=['Epoch','index',i, recv_timee]
    all_sats_pos.append(flag)
    for meas in corr:                              #loop over the satellites      
      if signal in meas.observables_final and np.isfinite(meas.observables_final[signal]):
        pr = meas.observables_final[signal]
        sat_pos = meas.sat_pos_final
      elif signal in meas.observables and np.isfinite(meas.observables[signal]) and meas.processed:
        pr = meas.observables[signal]
        pr += meas.sat_clock_err * constants.SPEED_OF_LIGHT
        sat_pos = meas.sat_pos
      else:
        print('\n\n Satelite position cannot be defined \n\n')
        continue
      sat_data = array([sat_pos[0], sat_pos[1], sat_pos[2], pr])
      all_sats_pos.append(sat_data)
    i=i+1

  all_sats_pos = pd.DataFrame(np.array(all_sats_pos))
  try:
    all_sats_pos.to_csv(generalinfo_path+'/all_sats_pos_time.csv', index=False)
  except:
    pass


def calc_user_pos(measurements_in, dataset, resultspath):
  print('Determine position: ', dataset)
  # rinex_processed_grouped = measurements_in
  times = []
  ests = []
  for corr in tqdm(measurements_in[:]):     #loop over the time/epochs
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

    ests = np.mean(ests[:, :3], axis=0)
    ests = pd.DataFrame(ests[:3].transpose())
    ests.to_csv(resultspath+'/mean_user_pos_'+dataset+'.csv', index=False)  

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


def get_obs_file(directory):
  obs_ext = ".obs"
  files = next(os.walk(directory))[2]
  for file in files:
    if os.path.splitext(file)[1] == obs_ext:
      return os.path.abspath(os.path.join(directory, file))
  return None


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
  generalinfo_path = os.path.join(directory,'allsatellites')
  if not os.path.exists(generalinfo_path):
    os.makedirs(generalinfo_path)
  return generalinfo_path


def process_many(root_directory):
  directories = next(os.walk(root_directory))[1]
  print(directories)
  for subdir in directories:
    subdir = os.path.join(root_directory, subdir)
    obs_file = get_obs_file(subdir)
    if obs_file:
      generalinfo_path = create_generalinfo_path(subdir)
      rinex_processed_grouped = get_processed_data(obs_file)
      
      get_sun_north_position_in_time(rinex_processed_grouped, generalinfo_path, 87)
      calc_user_pos(rinex_processed_grouped,'allsatellites',generalinfo_path)
      get_measurements_time(rinex_processed_grouped, generalinfo_path)





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



def filter_pats_by_child(paths, listB):
  filtered = []
  for path in paths:
    if str(os.path.split(path)[-1]) in listB:
      filtered.append(path)
  return filtered



def get_allready_processed_days(root_0, needed_files, month_names):
  m_d = 0
  m_d_s = 0
  months = {}
  subfolders_with_paths_months = [f.path for f in os.scandir(root_0) if f.is_dir()]
  subfolders_with_paths_months = filter_pats_by_child(subfolders_with_paths_months, month_names)
  days_with_positions = []
  days_with_sat_positions = []
  for month_root in subfolders_with_paths_months:
    month_name = os.path.split(month_root)[-1]
    days_with_paths = [f.path for f in os.scandir(month_root) if f.is_dir()]
    days = []
    with_sats = []
    d = 0   
    d_s = 0
    for day_root in days_with_paths:
      day_name = os.path.split(day_root)[-1]
      if len(day_name) == 12:
        d += 1
        days_with_positions.append(day_name)
        day_ind = day_name[-2:]
        days.append(day_ind)
        if is_all_data(day_root, needed_files, True):
          d_s += 1
          days_with_sat_positions.append(day_name)
          with_sats.append(day_ind)
    # print(month_name, ":  ", d, "  with_sats: ", d_s)
      
    days = list(set(days))
    with_sats = list(set(with_sats))
    print(month_name, ":  ", len(days), "  with_sats: ", len(with_sats))  
    months[month_name] = {"positions": days, "satelites_too": with_sats}
    m_d += len(days)
    m_d_s += len(with_sats)
  # print(months)
  # print(m_d, m_d_s)
  return days_with_positions, days_with_sat_positions, months


def decide_which_process_is_needed(filename, days_with_pos, days_with_sats):
  if filename in days_with_pos and filename in days_with_sats:
    return 0
  if filename in days_with_pos and filename not in days_with_sats:
    return 1
  return 2




def process_many_from_obs_root(obs_location, root_directory, needed_files, month_names):
  obs_files = get_obs_files(obs_location)
  # print(obs_files)
  days_with_positions, days_with_sat_positions, months = get_allready_processed_days(root_directory, needed_files, month_names)
  month_name = os.path.split(obs_location)[-1]
  root_directory = os.path.join(root_directory, month_name)
  total_processed = 0
  for obs_file in obs_files:
    obs_file_name = os.path.splitext(os.path.split(obs_file)[-1])[0][:-2]
    
    if len(obs_file_name) == 12:  
      print("Obs filename: ", obs_file_name)
      try:
        cond = decide_which_process_is_needed(obs_file_name, days_with_positions, days_with_sat_positions)
        # print(obs_file_name, "   ", cond)
        if cond != 0:
          path_of_results = create_dir_for_results(obs_file, root_directory)
          generalinfo_path = create_generalinfo_path(path_of_results)
          
          rinex_processed_grouped = get_processed_data(obs_file)
          if cond == 2:
            get_sats_pos_and_pr(rinex_processed_grouped, generalinfo_path)
            calc_user_pos(rinex_processed_grouped, 'allsatellites', generalinfo_path)
          if cond == 1:
            get_sats_pos_and_pr(rinex_processed_grouped, generalinfo_path)
          
          total_processed += 1
          print("\n                          Processed: {}/{} \n".format(total_processed, len(obs_files) - len(months.get(month_name, None)["satelites_too"] )))
          # sys.exit()
      except:
        pass



#                                             PERTH
obs_path = r"/Volumes/ADATA SE800/GPS/raw_data/PERTH/rinex_obs/julius"
destination_path = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/PERTH_daily_measurements"

#                                             NZLD
# obs_path = r"/Users/kelemensz/Documents/Research/GPS/NZLD_update/januar"
# destination_path = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NZLD"


needed_files = ["user_pos_allsatellites.csv", "all_sats_pos_time.csv"]
month_names = ["julius", "szeptember", "februar", "marcius", "augusztus", "januar", "december2019", "oktober", "november", "majus", "aprilis" , "junius"]
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/PERTH_daily_measurements"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NZLD"




process_many_from_obs_root(obs_path, destination_path, needed_files, month_names)









