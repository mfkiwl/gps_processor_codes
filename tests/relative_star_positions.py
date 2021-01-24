import pandas as pd
from numpy import *
import matplotlib.pyplot as plt


def norm_(v):
  a = array(v)
  return a/(sqrt(sum(power(a,2))))	


def normvec(vec):
	a = empty(shape(vec))
	for i in range(len(vec)):
		a[i] = norm_(vec[i])
	return a


def getSTARS(path, pair_index):
	if pair_index == 0:
		star_1 = pd.read_csv(path+'/star_d_positions.csv',skiprows=0).values#.transpose()[0]#.astype('float64')	
		star_2 = pd.read_csv(path+'/star_s_positions.csv',skiprows=0).values#.transpose()[0]#.astype('float64')	
	elif pair_index == 1:
		star_1 = pd.read_csv(path+'/star_Mirphak_positions.csv',skiprows=0).values#.transpose()[0]#.astype('float64')	
		star_2 = pd.read_csv(path+'/star_Vega_positions.csv',skiprows=0).values#.transpose()[0]#.astype('float64')	
	
	star_1 = normvec(star_1)
	star_2 = normvec(star_2)
	stars = array([star_1,star_2])
	return stars


def relative_angles(s1, s2):
	#s2 = s2[2:]
	l = min(len(s1), len(s2))
	d = empty(l)
	for i in range(l):
		d[i] = degrees(arccos(dot(s1[i], s2[i])))
	time = [*range(0,l)]
	
	plt.plot(time, d)
	plt.show()


path1 = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/perth_before_august/24h_net/allsatellites'
d1, s1 = getSTARS(path1, 1)


path2 = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/perth_before_august/24h_D/allsatellites'
d2, s2 = getSTARS(path2, 0)

#relative_angles(s1, d2)
#relative_angles(d1, s2)
relative_angles(d1, s1)
relative_angles(d2, s2)
