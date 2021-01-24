from __future__ import division
import scipy as sci
import scipy.special as sp
from numpy import *
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm, colors
from sympy import Point3D, Line3D, Segment3D, Point, Line, Plane,Poly
import pymap3d as pm
import os
from vpython import vector, rotate
import pylab as pl
import math


def getSTARS_for_galactic_system(path):
	# D = normvec(pd.read_csv(path+'/star_d_positions.csv',skiprows=0).values)	
	# S = normvec(pd.read_csv(path+'/star_s_positions.csv',skiprows=0).values)
	D = normvec(pd.read_csv(path+'/Denebola_positions.csv',skiprows=0).values)	
	S = normvec(pd.read_csv(path+'/Sadalmelik_positions.csv',skiprows=0).values)

	sun = pd.read_csv(path+'/SUN_positions.csv',skiprows=0).values
	galactic_n_p = pd.read_csv(path+'/GNP_positions.csv',skiprows=0).values
	galactic_center = pd.read_csv(path+'/GC_positions.csv',skiprows=0).values

	# star_1 = pd.read_csv(path+'/star_Mirphak_positions.csv',skiprows=0).values	
	# star_2 = pd.read_csv(path+'/star_Vega_positions.csv',skiprows=0).values

	# nunki = pd.read_csv(path+'/star_Nunki_positions.csv',skiprows=0).values
	# capella = pd.read_csv(path+'/star_Capella_positions.csv',skiprows=0).values
	nunki = pd.read_csv(path+'/Nunki_positions.csv',skiprows=0).values
	capella = pd.read_csv(path+'/Capella_positions.csv',skiprows=0).values
	return sun, galactic_n_p, galactic_center, nunki, capella, S, D


def normvec(vec):
	a = empty(shape(vec))
	for i in range(len(vec)):
		a[i] = unit(vec[i])
	return a


def groupvec(vec,l):
	groups = []
	groups_std = []
	ll = int(len(vec)/l)
	for i in range (l):
		groups.append(vec[i*ll:(i+1)*ll])
		# groups_std.append(std(vec[i*ll:(i+1)*ll]))
	groups = array(groups)
	# groups_std = array(groups_std)
	return groups


def p_m_byv_beta(d1,d2,v):
	pn = 0
	pp = 0
	nr_p = 0
	nr_m = 0
	for x in range(len(d1)):
		if float(d1[x]) < v and float(d1[x]) >= 0:
			pp += d2[x]
			nr_p += 1
		if float(d1[x]) > -v and float(d1[x]) <= 0:
			pn += d2[x]
			nr_m += 1
	try:
		pp = pp/nr_p
	except:
		pass
	try:
		pn = pn/nr_m
	except:
		pass
	return pn, pp



def p_m_byv_limits(d1,d2):
	d1 = array(d1)
	d2 = array(d2)
	d1_n_idx = where(d1<0.0)
	d1_p_idx = where(d1>0.0)
	# print('indices:		', d2[d1_n_idx], d2[d1_p_idx])
	d2_poz = d2[d1_p_idx]
	d2_neg = d2[d1_n_idx]
	pn = 0.0
	pp = 0.0
	if not d2_poz.size == 0:
		pp = mean(d2_poz)
		
	if not d2_neg.size == 0:
		pn = mean(d2_neg)
	
	# if math.isnan(pn):
	# 	pn = 0.0
	# if math.isnan(pp):
	# 	pp = 0.0
	return pn, pp


def earth_DS(user_mean,SD):
	proj = []
	for x in SD:
		proj.append(dot(unit(user_mean),x))
	return array(proj)


def projUMEAN_users_mean(groupusers,user_mean,l):
	proj2 = empty(l)
	std_shift_groups = empty(l)
	for i in range(l):
		proj2[i] = dot(mean(groupusers[i], axis=0), unit(user_mean))
		std_shift_groups[i] = mean(std(groupusers[i], axis=0))
	return proj2, std_shift_groups	

def scal(v1, v2):
	return vdot(v1, v2)

def unit(v):
	return 1/sqrt(scal(v,v))*v


def get_third_dir_by_cross_product(A, B):
	l = min(len(A), len(B))
	C = empty((l, 3))
	for i in range(l):
		C[i] = cross(A[i], B[i])
	return normvec(C)

def prepare_data_for_CMBR_cmap(path):
	user_ = pd.read_csv(path+'/user_pos_allsatellites.csv',skiprows=1).values#.transpose()
	user_mean = mean(user_, axis=0).astype('float64')
	shifts = user_ - user_mean 
	sun, galactic_n_p, galactic_center, nunki, galactic_anti_center, S, D = getSTARS_for_galactic_system(path)
	# S = normvec(S)
	# D = normvec(D)
	# SD = normvec(S - D).astype('float64')

	sun = normvec(sun)
	galactic_n_p = normvec(galactic_n_p)
	galactic_center = normvec(galactic_center)
	galactic_anti_center = normvec(galactic_anti_center)
	nunki = normvec(nunki)
	
	# Nx = normvec(galactic_anti_center - galactic_center).astype('float64')
	# Nz = normvec(galactic_n_p - sun).astype('float64')
	
	Nx = normvec(galactic_center).astype('float64')
	Nz = normvec(galactic_n_p).astype('float64')
	Ny = -get_third_dir_by_cross_product(Nz, Nx)

	# Nx = normvec(- sun + galactic_center).astype('float64')
	# Nz = normvec(galactic_n_p - sun).astype('float64')

	# for i in range(len(SD)-1):
	# 	print(degrees(arccos(dot(Nx[i], Nz[i]))))


	return user_mean, shifts, Nx, Ny, Nz, S, D

def visualize_positions_over_time(root_directory):
	list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]
	for directory in list_subfolders_with_paths:
		try:
			total_dir = os.path.join(directory, "allsatellites")
			
			user_mean, shifts, Nx, Ny, Nz, S, D = prepare_data_for_CMBR_cmap(total_dir)
			shifts = shifts.T
			plt.plot(shifts[0])
			plt.plot(shifts[1])
			plt.plot(shifts[2])

			plt.show()	
		except:
			print(total_dir)

# =================================================================================================
# =================================================================================================


root_directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/januar_5_12"

root_directory = r"/Users/kelemensz/Dropbox/GPS_data/Cumulative_plots/mean_projection/perth/daily_measurements/januar"
# root_directory = r"//Users/kelemensz/Documents/Research/GPS/process/24h/process_NZLD/januar"
root_directory = r"//Users/kelemensz/Documents/Research/GPS/process/24h/process_NASA/januar"

# visualize_positions_over_time(root_directory)

def vectors_dx_dy_dz(v1, v2):
	v1 = v1.T
	v2 = v2.T
	dv = v2-v1
	return dv

def vectors_dot_prods(v1, v2):
	v1 = normvec(v1)
	v2 = normvec(v2)
	dot_v1_v2 = []
	for i in range(len(v1)):
		dot_v1_v2.append(dot(v1[i], v2[i]))
	# v1 = v1.T
	# v2 = v2.T
	# dv = v2-v1
	# print(dot_v1_v2)
	print(degrees(arccos(dot_v1_v2)))
	return dot_v1_v2

def plot_dV(dv, title):
	plt.plot(dv)
	# plt.plot(dv[1])
	# plt.plot(dv[2])
	plt.title(title)
	# plt.show()


def compare_stars_over_time(directory_1, directory_2):

		total_dir_1 = directory_1 # os.path.join(directory_1, "allsatellites")
		total_dir_2 = directory_2 # os.path.join(directory_2, "allsatellites")
		
		# sun1, galactic_n_p1, galactic_center1, nunki1, galactic_anti_center1, S1, D1 = getSTARS_for_galactic_system(total_dir_1)
		# sun2, galactic_n_p2, galactic_center2, nunki2, galactic_anti_center2, S2, D2 = getSTARS_for_galactic_system(total_dir_2)
		star1 = getSTARS_for_galactic_system(total_dir_1)
		star2 = getSTARS_for_galactic_system(total_dir_2)
		
		# print(array(star1[0]))
		# print(array(star2[0]))

		for i in range(len(star1)):

		# 	dS = vectors_dx_dy_dz(star1[i], star2[i])
		# 	plot_dV(dS, str(i))
		# dS = vectors_dx_dy_dz(star1[0], star2[0])
		# plot_dV(dS, str("sun"))
			d = vectors_dot_prods(star1[i], star2[i])
			plot_dV(d, str("sun"))




# directory_1 = r"/Users/kelemensz/Documents/Research/GPS/process/24h/process_NZLD/januar/NZLD20200105/allsatellites"
# directory_2 = r"/Users/kelemensz/Dropbox/GPS_data/Cumulative_plots/mean_projection/perth/daily_measurements/januar/jan_5/allsatellites"
# directory_2 = r"/Users/kelemensz/Documents/Research/GPS/STARS_2020/STAR20200105"


directory_1 = r"/Users/kelemensz/Documents/Research/GPS/process/24h/process_NZLD/januar/NZLD20200107/allsatellites"
directory_2 = r"/Users/kelemensz/Dropbox/GPS_data/Cumulative_plots/mean_projection/perth/daily_measurements/januar/jan_7/allsatellites"


directory_1 = r"/Users/kelemensz/Documents/Research/GPS/stars_test/STARS_2020/STAR20200930"
# directory_2 = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020/STAR20200105"
directory_2 = r"/Users/kelemensz/Documents/Research/GPS/STARS_GREENWICH/STARS_2020/STAR20200930"

compare_stars_over_time(directory_1, directory_2)














