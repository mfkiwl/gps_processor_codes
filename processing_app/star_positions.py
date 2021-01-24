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

# S1=57632  # Denebola
# S2=109074  # Sadalmelik

#A GPS adatok szerint (RINEX) vasarnappal kezdodik a het: a hetek szama eggyel novekedik ha az idopont meghaladja/ELERI az 'aktualis het' vasarnapjanak 00:00 idopontjat
def format_time(recv_time):                       # csak 2020, 1, 12 utani adatokra mukodik
                                                 #datetime(2020, 1, 12) = week: 2088, tow: 0.0
  terrestrial_time = np.zeros(6)
  #----------------------------------------------------------------------------------------------=
  #----------------------------------------------------------------------------------------------=
  year = 2020                                     #adjust acording to the current year  
  month = 10                                       #adjust acording to the current month
  first_day_of_week = 4                           #adjust acording to the first day of the curent week
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
  terrestrial_time[3]=h_of_day - 3 #(shift the time because of some unclear reasons)
  terrestrial_time[4]=min_of_h
  terrestrial_time[5]=sec_of_min 
  print('terrestrial_time: ',terrestrial_time)
  return terrestrial_time

def get_star_position( recv_time, star_name, star_name2, user_longlatdist=geographic_position):
  # if not star_name:  
  # 	star_name = 57632             #Denebola
  # if not star_name2: 
  #   star_name2 = 109074         #Sadalmelik
  
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
  print('Location: ', '{} N'.format(user_longlatdist[0]), '{} E'.format(user_longlatdist[1]))
  print('S1 (alt/az):   ',alt_d, az_d)
  star_pos_d = pm.aer2ecef( az_d.degrees,alt_d.degrees, distance_d.m,user_longlatdist[0],user_longlatdist[1],user_longlatdist[2])
  
  sadalmelik = Star.from_dataframe(df.loc[star_name2])
  astrometric_sadalmelik = user_position.at(t).observe(sadalmelik)
  alt_s, az_s, distance_s = astrometric_sadalmelik.apparent().altaz()
  print('S2 (alt/az):   ',alt_s, az_s)
  star_pos_s = pm.aer2ecef( az_s.degrees,alt_s.degrees, distance_s.m,user_longlatdist[0],user_longlatdist[1],user_longlatdist[2])
  
  stars_pos=[star_pos_d,star_pos_s]
  return stars_pos


def galactic_Sun_north_directions( recv_time, user_longlatdist=geographic_position):
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
  
  # Északi Galaktikus Pólus (RA: 12h 51.42m; D: 27° 07.8′, epocha J2000.0) 
  coma_berenices = Star(ra_hours=(12, 51, 25.2), dec_degrees=(27, 7, 48.0))
  astrometric_berenice = user_position.at(t).observe(coma_berenices)
  alt_b, az_b, distance_b = astrometric_berenice.apparent().altaz()
  print('Coma Berenice (alt/az):   ',alt_b, az_b)
  berenice_position = pm.aer2ecef( az_b.degrees, alt_b.degrees, distance_b.m, user_longlatdist[0], user_longlatdist[1], user_longlatdist[2])
  
  stars_pos=[sun_position, berenice_position]
  return stars_pos


def get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, star_check_per_day, S1_, S2_, star_names=["d", "s"]):
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
    print(star_names, S1_, S2_)
    star_positions = get_star_position(recv_timee, star_name=S1_,star_name2=S2_)   # OK, ECEF coordinates
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


def get_sun_north_position_in_time(rinex_processed_grouped, generalinfo_path, star_check_per_day, star_names=["SUN", "GNP"]):
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
    star_positions = galactic_Sun_north_directions(recv_timee)   # OK, ECEF coordinates
    star1_pos.append(np.array(star_positions[0]))
    star2_pos.append(np.array(star_positions[1]))

  star1_pos = pd.DataFrame(np.array(star1_pos))
  try:
    star1_pos.to_csv(generalinfo_path+'/'+star_names[0]+'_positions.csv', index=False)
  except:
    pass

  star2_pos = pd.DataFrame(np.array(star2_pos))
  try:
    star2_pos.to_csv(generalinfo_path+'/'+star_names[1]+'_positions.csv', index=False)
  except:
    pass



#GET MEASUREMENT DATA
def get_processed_data(filepath):
    dog = AstroDog()
    print('Processing: ', filepath)
    obs_data = RINEXFile(filepath)
    rinex_meas_grouped = raw.read_rinex_obs(obs_data) 
    del obs_data 
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


S1 = 15863  # Mirphak
S2 = 91262  # Vega
star_names=["Mirphak", "Vega"]

def process_many(root_directory):
  directories = next(os.walk(root_directory))[1]
  print(directories)
  for subdir in directories:
    subdir = os.path.join(root_directory, subdir)
    print(subdir)
    obs_file = get_obs_file(subdir)
    if obs_file:
      generalinfo_path = create_generalinfo_path(subdir)
      #print(obs_file)
      #print(generalinfo_path)
      rinex_processed_grouped = get_processed_data(obs_file)
      # S1 = 15863  # Mirphak
      # S2 = 91262  # Vega
      # star_names=["Mirphak", "Vega"]
      # get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, 86, S1_=S1, S2_=S2, star_names=star_names)
      # S1 = 57632  # Denebola
      # S2 = 109074  # Sadalmelik
      # star_names=["d", "s"]
      # get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, 86, S1_=S1, S2_=S2, star_names=star_names)
      # S1 = 24608  # Capella
      # S2 = 92855  # Nunki
      # star_names=["Capella", "Nunki"]
      # get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, 86, S1_=S1, S2_=S2, star_names=star_names)

      get_sun_north_position_in_time(rinex_processed_grouped, generalinfo_path, 86)
    

      #calc_user_pos(rinex_processed_grouped,'allsatellites',generalinfo_path)

def process_one(directory):

  obs_file = get_obs_file(directory)
  if obs_file:
    generalinfo_path = create_generalinfo_path(directory)
    #print(obs_file)
    #print(generalinfo_path)
    rinex_processed_grouped = get_processed_data(obs_file)
    S1 = 15863  # Mirphak
    S2 = 91262  # Vega
    star_names=["Mirphak", "Vega"]
    get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, 86, S1_=S1, S2_=S2, star_names=star_names)
    S1 = 57632  # Denebola
    S2 = 109074  # Sadalmelik
    star_names=["beta_d", "beta_s"]
    get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, 86, S1_=S1, S2_=S2, star_names=star_names)
    S1 = 24608  # Capella
    S2 = 92855  # Nunki
    star_names=["Capella", "Nunki"]
    get_stars_position_in_time(rinex_processed_grouped, generalinfo_path, 86, S1_=S1, S2_=S2, star_names=star_names)
    
    # calc_user_pos(rinex_processed_grouped,'allsatellites',generalinfo_path)

root_directory = "/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/mean_projection/perth_october_4_10"

#directories = [os.path.join(path, name) for path, subdirs, files in os.walk(root) for name in files]
#directories = [os.walk(root_directory).next()[1]]
process_many(root_directory)
# /Users/kelemensz/Documents/Research/GPS/process/24h/RO/kolozs_szeptember_2/sept2_3_cluj
directory = "/Users/kelemensz/Documents/Research/GPS/process/24h/RO/24h_B"
# process_one(directory)

