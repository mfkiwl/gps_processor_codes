# this script filters the usable data from the measurements that contain inrealistic positions
# For example: the 24h_C and sept12_cluj 24 h measurements


import pandas as pd
from numpy import *


def filter_by_std(data):
	indices = []
	std_cumulative_mean = zeros(len(data))
	std_cumulative_mean[0] = data[0]
	for i in range(1, len(data)):
		if std(data[:i]) > 10*std_cumulative_mean[i-1]:
			indices.append(i)
			data = delete(data, i)
			print("!", i)
			std_cumulative_mean[i] = std_cumulative_mean[i-1]
			continue
		std_cumulative_mean[i] = std(data[:i])
	return indices


def filter_by_difference(data):
	limit = 50
	dat_DF = pd.DataFrame(data)
	dat_DF.columns = ['x', 'y', 'z']
	d2 = dat_DF.diff()
	l1 = d2.index[d2['x'] > limit].tolist()
	l2 = d2.index[d2['y'] > limit].tolist()
	l3 = d2.index[d2['z'] > limit].tolist()	
	l = l1 +l2
	l = l + l3
	# print(l)
	lf = []
	for i in l:
		lf.append(i)
		lf.append(i-1)
		lf.append(i+1)
	print (len(lf))
	user_filtered = delete(data, lf, 0)
	return user_filtered


def filter_false_measurements(path):
	user_ = pd.read_csv(path+'/user_pos_allsatellites.csv',skiprows=1).values#.transpose()
	user_std = std(user_, axis=0).astype('float64')
	print(user_std)
	print(len(user_))
	user_ = filter_by_difference(user_)
	user_std = std(user_, axis=0).astype('float64')
	print(user_std)
	print(len(user_))
	
	df = pd.DataFrame(user_)
	df.to_csv(path+'/user_pos_allsatellites_filtered.csv', index=False)



# path = r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/sept12_cluj/allsatellites"
path = r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/24h_C/allsatellites"

filter_false_measurements(path)







