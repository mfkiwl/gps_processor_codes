
from numpy import *
import pandas as pd
from itertools import chain
import os
import json


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



def get_position_from_day_root(root):
	root = os.path.join(root, "allsatellites")
	user_ = pd.read_csv(root+'/user_pos_allsatellites.csv',skiprows=1).values	
	return mean(user_, axis=0).astype('float64')



def get_all_possitions_std(root_0, month_names):

	subfolders_with_paths_months = filter_pats_by_child([f.path for f in os.scandir(root_0) if f.is_dir()], month_names)
	all_positions = []
	days_with_positions = []
	for month_root in subfolders_with_paths_months:
		month_name = os.path.split(month_root)[-1]
		days_with_paths = [f.path for f in os.scandir(month_root) if f.is_dir()]
		months_positions = []

		for day_root in days_with_paths:
			day_name = os.path.split(day_root)[-1]
			if len(day_name) == 12:
				position = get_position_from_day_root(day_root)
				months_positions.append(position)
				all_positions.append(position)
		

		print(month_name, ":  ", std(array(months_positions), axis=0))
		
	print(std(array(months_positions), axis=0))




needed_files = ["user_pos_allsatellites.csv", "all_sats_pos_time.csv"]
month_names = ["julius", "szeptember", "februar", "marcius", "augusztus", "januar", "december2019", "oktober", "november", "majus", "aprilis" , "junius"]
results_root = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/PERTH_daily_measurements"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/process_NZLD"


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
	print(months)
	# print(m_d, m_d_s)
	return days_with_positions, days_with_sat_positions



get_allready_processed_days(results_root, needed_files, month_names)

# get_all_possitions_std(results_root, month_names)


