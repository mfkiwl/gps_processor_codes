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
from sympy import Point3D, Line3D, Segment3D, Point, Line, Plane
import pymap3d as pm
from vpython import *
import os
import sys


# geographic_position = [46.7502933, 23.5835082, 6367740.701600294]  # Kolozsvar
# geographic_position = [46.083, 26.233, 6367740.701600294]  # Almas
geographic_position = [-31.953512, 115.857048, 6367740.701600294]  # Perth

#A GPS adatok szerint (RINEX) vasarnappal kezdodik a het: a hetek szama eggyel novekedik ha az idopont meghaladja/ELERI az 'aktualis het' vasarnapjanak 00:00 idopontjat
def format_time(recv_time):                       # csak 2020, 1, 12 utani adatokra mukodik
                                                 #datetime(2020, 1, 12) = week: 2088, tow: 0.0
  terrestrial_time = np.zeros(6)
  #----------------------------------------------------------------------------------------------=
  #----------------------------------------------------------------------------------------------=
  year = 2020                                     #adjust acording to the current year  
  month = 2                                         #adjust acording to the current month
  first_day_of_week = 9                           #adjust acording to the first day of the curent week
  #----------------------------------------------------------------------------------------------=
  #----------------------------------------------------------------------------------------------=
  dow = int(recv_time.tow/(24.0*3600.0))
  day = first_day_of_week + dow 
  sec_of_day = recv_time.tow%(24.0*3600.0)
  h_of_day = int(sec_of_day/3600.0) 
  sec_of_h = int(sec_of_day%3600.0)
  min_of_h = int(sec_of_h/60.0)
  sec_of_min = int(sec_of_h%60.0)
  terrestrial_time[0]=year
  terrestrial_time[1]=month
  terrestrial_time[2]=day
  terrestrial_time[3]=h_of_day - 2 #(shift the time because of some unclear reasons)
  terrestrial_time[4]=min_of_h
  terrestrial_time[5]=sec_of_min 
  print('terrestrial_time: ',terrestrial_time)
  return terrestrial_time

def get_star_position( recv_time, user_longlatdist=geographic_position, star_name=57632, star_name2=109074, star_names=None):
  if not star_name:
  	star_name = 57632             #Denebola
  if not star_name2: 
    star_name2 = 109074         #Sadalmelik
  
  user_longlatdist = np.array(user_longlatdist)#46.07983746062874, 26.232615661106784, 659.5100082775696
  with load.open(hipparcos.URL) as f:
    df = hipparcos.load_dataframe(f)

  planets = load('de421.bsp')
  earth = planets['earth']
  user_position = earth + Topos('{} N'.format(user_longlatdist[0]), '{} E'.format(user_longlatdist[1])) #user_longlatdist      user_longlatdist[0],user_longlatdist[1]
  terrestrial_time = format_time(recv_time)
  ts = load.timescale()
  t = ts.tt(terrestrial_time[0], terrestrial_time[1], terrestrial_time[2], terrestrial_time[3], terrestrial_time[4], terrestrial_time[5])
  
  denebola = Star.from_dataframe(df.loc[star_name])
  astrometric_denebola = user_position.at(t).observe(denebola)
  alt_d, az_d, distance_d = astrometric_denebola.apparent().altaz()#astrometric.apparent().altaz() = astrometric_denebola.apparent().altaz()
  print(star_names)
  print('Denebola (alt/az):   ',alt_d, az_d)
  star_pos_d = pm.aer2ecef( az_d.degrees,alt_d.degrees, distance_d.m,user_longlatdist[0],user_longlatdist[1],user_longlatdist[2])
  
  sadalmelik = Star.from_dataframe(df.loc[star_name2])
  astrometric_sadalmelik = user_position.at(t).observe(sadalmelik)
  alt_s, az_s, distance_s = astrometric_sadalmelik.apparent().altaz()
  print('Sadalmelik (alt/az):   ',alt_s, az_s)
  star_pos_s = pm.aer2ecef( az_s.degrees,alt_s.degrees, distance_s.m,user_longlatdist[0],user_longlatdist[1],user_longlatdist[2])
  
  stars_pos=[star_pos_d,star_pos_s]
  return stars_pos

def galactic_Sun_north_center_SD_directions( recv_time, user_longlatdist=geographic_position):
  user_longlatdist = np.array(user_longlatdist)#46.07983746062874, 26.232615661106784, 659.5100082775696
  with load.open(hipparcos.URL) as f:
    df = hipparcos.load_dataframe(f)

  planets = load('de421.bsp')
  earth = planets['earth']
  user_position = earth + Topos('{} N'.format(user_longlatdist[0]), '{} E'.format(user_longlatdist[1])) #user_longlatdist      user_longlatdist[0],user_longlatdist[1]
  terrestrial_time = format_time(recv_time)
  ts = load.timescale()
  t = ts.tt(terrestrial_time[0], terrestrial_time[1], terrestrial_time[2], terrestrial_time[3], terrestrial_time[4], terrestrial_time[5])
  
  sun = planets['sun']
  astrometric_sun = user_position.at(t).observe(sun)
  alt_sun, az_sun, distance_sun = astrometric_sun.apparent().altaz()#astrometric.apparent().altaz() = astrometric_denebola.apparent().altaz()
  print('Location: ', '{} N'.format(user_longlatdist[0]), '{} E'.format(user_longlatdist[1]))
  print('Sun (alt/az):   ',alt_sun, az_sun)
  sun_position = pm.aer2ecef( az_sun.degrees, alt_sun.degrees, distance_sun.m, user_longlatdist[0], user_longlatdist[1], user_longlatdist[2])
  
  # ??szaki Galaktikus P??lus (RA: 12h 51.42m; D: 27?? 07.8???, epocha J2000.0) 
  coma_berenices = Star(ra_hours=(12, 51, 25.2), dec_degrees=(27, 7, 48.0))
  astrometric_berenice = user_position.at(t).observe(coma_berenices)
  alt_b, az_b, distance_b = astrometric_berenice.apparent().altaz()
  print('Coma Berenice (alt/az):   ',alt_b, az_b)
  berenice_position = pm.aer2ecef( az_b.degrees, alt_b.degrees, distance_b.m, user_longlatdist[0], user_longlatdist[1], user_longlatdist[2])
  
  # Galaktikus Center (RA: 12h 51.42m; D: 27?? 07.8???, epocha J2000.0) 17h???45.6m  ???28.94??
  galactic_center = Star(ra_hours=(17, 45, 36.0), dec_degrees=(-28, 56, 24.0))
  astrometric_center = user_position.at(t).observe(galactic_center)
  alt_c, az_c, distance_c = astrometric_center.apparent().altaz()
  print('Galactic Center (alt/az):   ',alt_c, az_c)
  center_position = pm.aer2ecef( az_c.degrees, alt_c.degrees, distance_c.m, user_longlatdist[0], user_longlatdist[1], user_longlatdist[2])

  nunki = Star.from_dataframe(df.loc[92855])
  astrometric_nunki = user_position.at(t).observe(nunki)
  alt_n, az_n, distance_n = astrometric_nunki.apparent().altaz()#astrometric.apparent().altaz() = astrometric_denebola.apparent().altaz()
  print('Nunki (alt/az):   ',alt_n, az_n)
  star_pos_n = pm.aer2ecef( az_n.degrees, alt_n.degrees, distance_n.m,user_longlatdist[0],user_longlatdist[1],user_longlatdist[2])

  capella = Star.from_dataframe(df.loc[24608])
  astrometric_capella = user_position.at(t).observe(capella)
  alt_ca, az_ca, distance_ca = astrometric_capella.apparent().altaz()#astrometric.apparent().altaz() = astrometric_denebola.apparent().altaz()
  print('Capella (alt/az):   ',alt_ca, az_ca)
  star_pos_ca = pm.aer2ecef( az_ca.degrees, alt_ca.degrees, distance_ca.m,user_longlatdist[0],user_longlatdist[1],user_longlatdist[2])

  denebola = Star.from_dataframe(df.loc[57632])
  astrometric_denebola = user_position.at(t).observe(denebola)
  alt_d, az_d, distance_d = astrometric_denebola.apparent().altaz()#astrometric.apparent().altaz() = astrometric_denebola.apparent().altaz()
  print('Denebola (alt/az):   ',alt_d, az_d)
  star_pos_d = pm.aer2ecef( az_d.degrees,alt_d.degrees, distance_d.m,user_longlatdist[0],user_longlatdist[1],user_longlatdist[2])
  
  sadalmelik = Star.from_dataframe(df.loc[109074])
  astrometric_sadalmelik = user_position.at(t).observe(sadalmelik)
  alt_s, az_s, distance_s = astrometric_sadalmelik.apparent().altaz()
  print('Sadalmelik (alt/az):   ',alt_s, az_s)
  star_pos_s = pm.aer2ecef( az_s.degrees,alt_s.degrees, distance_s.m,user_longlatdist[0],user_longlatdist[1],user_longlatdist[2])
  

  stars_pos=[sun_position, berenice_position, center_position, star_pos_n, star_pos_ca, star_pos_d, star_pos_s]
  return stars_pos

def get_sun_north_position_in_time(rinex_processed_grouped, generalinfo_path, star_check_per_day):
  n_epochs = len(rinex_processed_grouped)
  star_check_distance = int(n_epochs/float(star_check_per_day))
  dog = AstroDog()
  signal='C1C'
  

  sun_pos = []
  northpole_pos = []
  center_pos = []
  nunki_pos = []
  capella_pos = []
  denebola_pos = []
  sadalmelik_pos = []

  rinex_processed_grouped_B = rinex_processed_grouped[::star_check_distance]
  print('Start processing!')
  for corr in tqdm(rinex_processed_grouped_B):     #loop over the time/epochs
    recv_timee = corr[0].recv_time
    star_positions = galactic_Sun_north_center_SD_directions(recv_timee)   # OK, ECEF coordinates
    
    sun_pos.append(np.array(star_positions[0]))
    northpole_pos.append(np.array(star_positions[1]))

    center_pos.append(np.array(star_positions[2]))
    nunki_pos.append(np.array(star_positions[3]))
    capella_pos.append(np.array(star_positions[4]))
    
    denebola_pos.append(np.array(star_positions[5]))
    sadalmelik_pos.append(np.array(star_positions[6]))


  try:
    pd.DataFrame(np.array(sun_pos)).to_csv(generalinfo_path+'/SUN_positions.csv', index=False)
  except:
    pass

  try:
    pd.DataFrame(np.array(northpole_pos)).to_csv(generalinfo_path+'/GNP_positions.csv', index=False)
  except:
    pass
  
  try:
    pd.DataFrame(np.array(center_pos)).to_csv(generalinfo_path+'/GC_positions.csv', index=False)
  except:
    pass

  try:
    pd.DataFrame(np.array(nunki_pos)).to_csv(generalinfo_path+'/Nunki_positions.csv', index=False)
  except:
    pass
  
  try:
    pd.DataFrame(np.array(capella_pos)).to_csv(generalinfo_path+'/Capella_positions.csv', index=False)
  except:
    pass

  try:
    pd.DataFrame(np.array(denebola_pos)).to_csv(generalinfo_path+'/Denebola_positions.csv', index=False)
  except:
    pass

  try:
    pd.DataFrame(np.array(sadalmelik_pos)).to_csv(generalinfo_path+'/Sadalmelik_positions.csv', index=False)
  except:
    pass



def get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, star_check_per_day, star1=None, star2=None, star_names=["d", "s"]):
  n_epochs = len(rinex_processed_grouped)
  star_check_distance = int(n_epochs/float(star_check_per_day))
  dog = AstroDog()
  signal='C1C'
  all_sats_pos = []
  star1_pos = []
  star2_pos = []
  rinex_processed_grouped_B = rinex_processed_grouped[::star_check_distance]
  print('Start processing!')
  for corr in tqdm(rinex_processed_grouped_B):     #loop over the time/epochs
    recv_timee = corr[0].recv_time
    print(star_names, star1, star2)
    star_positions = get_star_position(recv_timee, user_longlatdist=geographic_position, star_name=star1, star_name2=star2, star_names=star_names)   # OK, ECEF coordinates
    star1_pos.append(np.array(star_positions[0]))
    star2_pos.append(np.array(star_positions[1]))

  star1_pos = pd.DataFrame(np.array(star1_pos)/10**16)
  try:
    star1_pos.to_csv(generalinfo_path+'/star_'+star_names[0]+'_positions.csv', index=False)
  except:
    pass

  star2_pos = pd.DataFrame(np.array(star2_pos)/10**16)
  try:
    star2_pos.to_csv(generalinfo_path+'/star_'+star_names[1]+'_positions.csv', index=False)
  except:
    pass


def get_measurements_time(rinex_processed_grouped, generalinfo_path):
  dog = AstroDog()
  signal='C1C'
  all_meas_time = []

  print('Start processing!')
  for corr in tqdm(rinex_processed_grouped):     #loop over the time/epochs
    recv_timee = corr[0].recv_time
    
    all_meas_time.append(recv_timee)

  all_meas_time = pd.DataFrame(np.array(all_meas_time))
  try:
    all_meas_time.to_csv(generalinfo_path+'/times_allsatellites.csv', index=False)
  except:
    pass


def get_sats_pos(rinex_processed_grouped):
  dog = AstroDog()
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


def calc_user_pos(measurements_in,dataset,resultspath):
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
    print('Processing: ', filepath)
    obs_data = RINEXFile(filepath)
    rinex_meas_grouped = raw.read_rinex_obs(obs_data) 
    del obs_data 
    print('rinex_meas_grouped size: ', sys.getsizeof(rinex_meas_grouped))
    rinex_processed_grouped = []
    for meas in tqdm(rinex_meas_grouped):
        proc = raw.process_measurements(meas, dog=dog)
        rinex_processed_grouped.append(proc)
    print('Data is IN!')
    del rinex_meas_grouped
    print('rinex_processed_grouped size: ', sys.getsizeof(rinex_processed_grouped))
    return rinex_processed_grouped

def get_obs_file(directory):
  obs_ext = ".obs"
  files = next(os.walk(directory))[2]
  for file in files:
    if os.path.splitext(file)[1] == obs_ext:
      return os.path.abspath(os.path.join(directory, file))
  return None


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
    #print(subdir)
    obs_file = get_obs_file(subdir)
    if obs_file:
      generalinfo_path = create_generalinfo_path(subdir)
      #print(obs_file)
      #print(generalinfo_path)
      rinex_processed_grouped = get_processed_data(obs_file)
      # S1=57632  # Denebola
      # S2=109074  # Sadalmelik

      # get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, 87, star1=57632, star2=109074)
      # S1 = 15863  # Mirphak
      # S2 = 91262  # Vega
      # star_names=["Mirphak", "Vega"]
      # get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, 87, star1=15863, star2=91262, star_names=star_names)
      get_sun_north_position_in_time(rinex_processed_grouped, generalinfo_path, 87)
      calc_user_pos(rinex_processed_grouped,'allsatellites',generalinfo_path)
      get_measurements_time(rinex_processed_grouped, generalinfo_path)

root_directory = r"/Users/kelemensz/Documents/Research/GPS/reciever_data/perth_data/NewFolder/febr_9_15"
# root_directory = r"/Users/kelemensz/Documents/Research/GPS/reciever_data/perth_data/NewFolder/marc_8_14"
#list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]

#directories = [os.path.join(path, name) for path, subdirs, files in os.walk(root) for name in files]
#directories = [os.walk(root_directory).next()[1]]

process_many(root_directory)

# directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/RO/cmap/24h_G"

# obs_file = get_obs_file(directory)
# if obs_file:
#   generalinfo_path = create_generalinfo_path(directory)
#   #print(obs_file)
#   #print(generalinfo_path)
#   rinex_processed_grouped = get_processed_data(obs_file)
#   # S1=57632  # Denebola
#   # S2=109074  # Sadalmelik

#   # get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, 87, star1=57632, star2=109074)
#   # S1 = 15863  # Mirphak
#   # S2 = 91262  # Vega
#   # star_names=["Mirphak", "Vega"]
#   # get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, 87, star1=15863, star2=91262, star_names=star_names)
#   # calc_user_pos(rinex_processed_grouped,'allsatellites',generalinfo_path)
#   get_measurements_time(rinex_processed_grouped, generalinfo_path)




