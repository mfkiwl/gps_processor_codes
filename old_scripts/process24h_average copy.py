from numpy import *
import pandas as pd
import matplotlib.pyplot as plt
from sympy import Point3D, Line3D, Segment3D, Point, Line, Plane,Poly
import pymap3d as pm
import os
from vpython import vector, rotate



def create_directory(path, directory_name):
	directory = os.path.join(path, directory_name)
	if not os.path.isdir(directory):
		os.makedirs(directory)
	return directory

def getSD(path):
	star_D = pd.read_csv(path+'/star_d_positions.csv',skiprows=0).values#.transpose()[0]#.astype('float64')	
	star_S = pd.read_csv(path+'/star_s_positions.csv',skiprows=0).values#.transpose()[0]#.astype('float64')	
	star_S = normvec(star_S)
	star_D = normvec(star_D)
	stars = array([star_D,star_S])
	return stars


def norm_(v):
  a = array(v)
  return a/(sqrt(sum(power(a,2))))	


def normvec(vec):
	a = empty(shape(vec))
	for i in range(len(vec)):
		a[i] = norm_(vec[i])
	return a


def scalevec(vec):
	a = empty(len(vec))
	
	for i in range(len(vec)):
		a[i] = sqrt(sum(power(vec[i],2)))
	s = 2*std(a)
	#s = max(a)
	return vec/s


def groupmeanvec(vec,l):
	groups = []
	for i in range (int(len(vec)/l)+1):
		groups.append(mean(vec[i*l:(i+1)*l],axis=0))
	groups = array(groups)
	return groups


def groupmeanvec_std(vec,l):
	groups = []
	groups_std = []
	for i in range (int(len(vec)/l)+1):
		groups.append(mean(vec[i*l:(i+1)*l],axis=0))
		groups_std.append(std(vec[i*l:(i+1)*l]))
	groups = array(groups)
	groups_std = array(groups_std)
	return groups, groups_std


def groupvec(vec,l):
	groups = []
	groups_std = []
	ll = int(len(vec)/l)
	for i in range (l):
		groups.append(vec[i*ll:(i+1)*ll])
		groups_std.append(std(vec[i*ll:(i+1)*ll]))
	groups = array(groups)
	groups_std = array(groups_std)
	return groups


def projonmeanuser(user_mean,usergrupped):
	dd = []
	for x in usergrupped:
		dd.append(dot(user_mean,x))
	dd = array(dd)
	plt.plot(dd,label='proj. of the shifts on user_mean\n of the mean shifts/SD update')
	plt.legend()
	plt.show()
	plt.clf()


def pearson(d1,d2):
	d1=array(d1)
	d2=array(d2)
	c = ( mean(d2*d1) - mean(d2)*mean(d1) )/( (mean(pow(d1,2)) - pow(mean(d1),2))  * (mean(pow(d2,2)) - pow(mean(d2),2)) )
	print ('Pearson correlation coefficient: ',c)
	return (c)


def plus_minus(d2):
	p = len(d2[d2>0])
	m = len(d2[d2<0])
	print('plus:_ , minus:_ ',p,m)


def plus_minus_condition(d1,d2,v):
	p = 0
	n = 0
	#print('gg ',len(d1),len(d2))
	print('max(d1): ',max(d1),'\n')

	for x in range(len(d1)):
		if abs(d1[x]) < v:
			n += 1
			#print(degrees(arccos(d1[x])))
			p += d2[x]/abs(float(d2[x]))
			#print(x,' - ',d2[x])
			#print(d2[x]/abs(d2[x]))
	#print(len(d2),n)
	print('With condition: ',p, '	v = ',v)	
	return 0


def p_m_byv_abs(d1,d2,v):
	p = 0
	for x in range(len(d1)):
		if abs(float(d1[x])) < v:
			#print(float(d2[x]/abs(d2[x])))
			p += d2[x]/(abs(d2[x])+1.0)
	return p 


def p_m_byv_normed(d1, d2, std_d2, v):
	pn = 0
	pp = 0
	nr_p = 0
	nr_m = 0
	for x in range(len(d1)):
		if float(d1[x]) < v and float(d1[x]) >= 0:
			#print(float(d2[x]/abs(d2[x])))
			pp += (d2[x]/(abs(d2[x])+0.1))*(1/std_d2[x])
			#print(d1[x])
			nr_p += 1
	#for x in range(len(d1)):
		if float(d1[x]) > -v and float(d1[x]) <= 0:
			#print(float(d2[x]/abs(d2[x])))
			pn += (d2[x]/(abs(d2[x])+0.1))*(1/std_d2[x])
			#print(d1[x])
			nr_m += 1

	return pn,pp


def p_m_byv(d1,d2,v):
	pn = 0
	pp = 0
	for x in range(len(d1)):
		if float(d1[x]) < v and float(d1[x]) >= 0:
			#print(float(d2[x]/abs(d2[x])))
			pp += d2[x]/(abs(d2[x])+0.1)
			#print(d1[x])
	#for x in range(len(d1)):
		if float(d1[x]) > -v and float(d1[x]) <= 0:
			#print(float(d2[x]/abs(d2[x])))
			pn += d2[x]/(abs(d2[x])+0.1)
			#print(d1[x])

	return pn,pp

def plotbyvalue(d1,d2,v):
	a=[]
	b=[]
	for i in range(len(d1)):
		if abs(d1[i])>v:
			a.append(d1[i])
			b.append(d2[i])
	a = array(a)	
	b = array(b)
	plt.clf()
	plt.plot(a)
	plt.plot(b,'o',ms = 3)
	plt.show()	


def projction_1(SD,user_mean,user_,l0,l):
	
	groups = groupmeanvec(user_,l)
	groups = array(groups)-user_mean
	user_mean = norm_(user_mean)
	SD = normvec(SD)
	d1 = []
	d2 = []
	d2_mean = []
	for i in range(0,len(SD)-1):
		cosSDn = dot(user_mean,SD[i])
		d1.append(cosSDn)
		k = i*int(l0/l)
		for j in range(0, int(l0/l),1):
			cosSDm = dot(groups[k+j],SD[i])
			#cosSDm = dot(user_mean,groups[i+j])
			d2.append(cosSDm)
		cosSDm = dot(mean(groups[k:k+j],axis=0),SD[i])
		#cosSDm = dot(user_mean,mean(groups[k:k+j],axis=0))
		d2_mean.append(cosSDm)
	d1 = array(d1)
	d2 = array(d2)
	d2_mean = array(d2_mean)
	print(mean(d2))
	print(len(d2))
	plus_minus(d2)
	#plus_minus_condition(d1,d2)
	plt.plot(d1)
	plt.plot(d2_mean,'o',ms = 3)
	plt.show()	
	plt.clf()
	plt.hist(d2,bins=200)
	plt.show()
	#pearson(d1,d2)


def earth_DS(user_mean,SD):
	proj = []
	for x in SD:
		proj.append(dot(norm_(user_mean),x))
	return array(proj)


def projSD_users_simple(user_,SD):
	l = len(SD)
	proj = empty(l)
	for i in range(l):
		proj[i] = dot(norm_(user_[i]),SD[i])
		#proj[i] = dot(norm(user_[i]),norm(user_mean))
	return array(proj)


def projUMEAN_users_simple(user_,user_mean,SD):
	l = len(SD)
	if len(user_) in arange(l-5,l+5,1):
		print('egyforma hosszuak!	: ',l,len(user_))
	proj = empty(l)
	for i in range(l):
		#proj[i] = dot(norm(user_[i]),SD[i])
		proj[i] = dot(norm_(user_[i]),norm_(user_mean))
	return array(proj)


def projSD_users(groupusers,SD):
	l = len(SD)
	proj2 = empty(l)
	for i in range(l):
		a=0
		#gr = normvec(groupusers[i])
		for j in groupusers[i]:
			#a += dot(norm(j),SD[i])/abs(dot(norm(j),SD[i]))
			a += dot(j,SD[i])#/abs(dot(j,SD[i]))
		#print(a)
		proj2[i] = a/len(groupusers[i])
		a=0

	#print('sum d2: ',sum(proj2))
	return proj2


def projUMEAN_users(groupusers,user_mean,SD):
	l = len(SD)
	proj2 = empty(l)
	std_shift_groups = empty(l)
	for i in range(l):
		a=0
		#gr = normvec(groupusers[i])
		for j in groupusers[i]:
			a += dot(j,norm_(user_mean))#/abs(dot(j,norm(user_mean)))
		#print(a)
		proj2[i] = a
		std_shift_groups[i] = mean(std(groupusers[i], axis=0))
		a=0
	#print(std_shift_groups)
	#print('sum d2: ',sum(proj2))
	return proj2, std_shift_groups


def cumulativebyvalue(user_mean,SD,user_,cond = 2):
	SD = normvec(SD)
	l=len(SD)
	shifts = user_ - user_mean
	print('Std. of the shifts:',std(shifts,axis=0))
		
	d1 = earth_DS(user_mean,SD)  # SD-n angles
	minange = max(abs(d1))
	print(minange)  # min angle

	
	cond = 0
	if cond == 0:
		usergrupped = groupvec(shifts,l)
		#d2 = projSD_users(usergrupped,SD)	
		
		d2, std_d2 = projUMEAN_users(usergrupped,user_mean,SD)	# sum(+-1) a csoporton belul
		d2 = array(d2).astype('float64')

	elif cond == 1:
		k = int(len(shifts)/l)
		#k=200
		#print('k: ',k)
		usergrupped = groupmeanvec(shifts,k)
		print('cond = 1: ',len(usergrupped))
		#d2 = projSD_users_simple(usergrupped,SD)
		d2 = projUMEAN_users_simple(usergrupped,user_mean,SD)
		projonmeanuser(user_mean, usergrupped)

	else:
		#k = int(len(shifts)/l)
		#b = int(len(shifts)/l)
		k=200
		usergruppedmean = groupmeanvec(shifts,k)
		print('cond > 1: ',len(usergruppedmean))

		usergrupped = groupvec(usergruppedmean,l)		
		#print(len(usergrupped))

		d2 = projSD_users(usergrupped,SD)
	

	#random array for testing: mean value=0.1 ; std = 1
	#d2 = random.normal(0.1,size=len(d2))

	# plt.plot(d2,'o',ms=2)
	# plt.show()
	d1 = array(d1).astype('float64')
	d2 = array(d2).astype('float64')
	print(len(d1),len(d2))
	
	# plt.plot(d1, d2)
	# plt.show()
	
	nr = 200
	lims = linspace(0.0,minange,nr)
	cumulative = empty((nr,2))
	for i in range(nr):
		v = lims[i]
		cumulative[i] = p_m_byv(d1,d2,v)
		#cumulative[i] = p_m_byv_normed(d1,d2,std_d2,v)
		#cumulative[i] = p_m_byv_abs(d1,d2,v)

	# plot_cumulative(lims, cumulative, str(random.randint(0,100)), '_')


	return lims, cumulative


def calc_mean_cumulative(data, rot_axe, rot_angle, directory):
	lims = []#empty(len(data))
	cumulative = []#empty(len(data))
	i=0

	for lims_cumulative in data:
		lims.append(lims_cumulative[0])
		cumulative.append(lims_cumulative[1])
		i += 1
	lims = mean(lims, axis=0)
	cumulative_plus = []#empty(len(data))
	cumulative_minus = []#empty(len(data))
	i=0
	for cumulative_i in cumulative:
		cum = array(cumulative_i).T
		cumulative_plus.append(array(cum[1]))
		cumulative_minus.append(array(cum[0]))

	cumulative_plus = mean(array(cumulative_plus), axis=0)
	cumulative_minus = mean(array(cumulative_minus), axis=0)
	cumulative_final = []#empty(2)
	cumulative_final.append(cumulative_minus)
	cumulative_final.append(cumulative_plus)
	rot_axe = "SD_X_n"
	plot_cumulative(lims, array(cumulative_final).T, rot_axe, rot_angle, directory)



def plot_cumulative(lims, cumulative, rot_axe, rot_angle, directory):
	#print(cumulative)
	plt.clf()
	plt.plot(-lims,cumulative[:,0],label = 'Cumulative -')
	plt.plot(lims,cumulative[:,1],label = 'Cumulative +')
	plt.title('24h_average, SD rotated: {}: {}'.format(rot_axe, str(rot_angle)))
	plt.legend()
	# plt.show()
	
	# plt.savefig(os.path.join('/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/24h_C', '{}_{}_not_weighted'.format(rot_axe, str(rot_angle))))
	plt.savefig(os.path.join(directory, '{}_{}_not_weighted.png'.format(rot_axe, str(rot_angle))))
	
	return 0
	

def proj_on_SD(user_mean,SD,user_):
	SD = normvec(SD)
	l=len(SD)
	shifts = user_ - user_mean
	print('Std. of the shifts:',std(shifts,axis=0))
		
	d1 = earth_DS(user_mean,SD)
	minange = max(abs(d1))
	print(minange)
	
	usergrupped = groupvec(shifts,l)
	d2 = projSD_users(usergrupped,SD)	
	d2 = array(d2).astype('float64')
	print('Value of the mean shift\'s projction on SD: ' , mean(d2))
	plt.plot(d2, 'o')
	plt.show()
	

def rot_SD(SD_init, rot_axe, rot_angle):
	if rot_axe == "X":
		rot_axe = [1.0,0.0,0.0]
	elif rot_axe == "Y":
		rot_axe = [0.0,1.0,0.0]
	elif rot_axe == "Z":
		rot_axe = [0.0,0.0,1.0]
	
	rot_axe = vector(rot_axe[0],rot_axe[1],rot_axe[2])
	angle = radians(float(rot_angle))
	SD_rotated = empty(shape(SD_init))
	i = 0
	for sd in SD_init:
		sd_v = vector(sd[0], sd[1], sd[2])
		#print(sd, '-sd\n', sd_v.x, sd_v.y, sd_v.z)
		#=====    rotate
		sd_r = rotate(sd_v, angle=angle, axis=rot_axe)
		#print( sd_r.x, sd_r.y, sd_r.z)
		SD_rotated[i] = array([sd_r.x, sd_r.y, sd_r.z])
		
		#print(dot(SD_rotated[i], sd))
		i += 1
	#print('dot() of SDs: ', sum(dot(SD_init[0], SD_rotated[0])))
	return SD_rotated


def rot_SD_u_mean(SD_init, user, rot_angle):
	user = array(user)
	angle = radians(float(rot_angle))
	#rot_axe = vector(rot_axe[0],rot_axe[1],rot_axe[2])
	SD_rotated = empty(shape(SD_init))
	i = 0
	for sd in SD_init:
		sd_v = vector(sd[0], sd[1], sd[2])
		p0 = Point3D(0.0, 0.0, 0.0)
		p1 = Point3D(sd[0], sd[1], sd[2])
		p2 = Point3D(user[0], user[1], user[2])
		plane_actual = Plane(p0, p1, p2)
		rot_axe = plane_actual.normal_vector
		rot_axe = [rot_axe[0].n(5),rot_axe[1].n(5),rot_axe[2].n(5)]
		#print(rot_axe)
		rot_axe = vector(float(rot_axe[0]), float(rot_axe[1]), float(rot_axe[2]))
		#print(angle)
		sd_r = rotate(sd_v, angle=angle, axis=rot_axe)
		SD_rotated[i] = array([sd_r.x, sd_r.y, sd_r.z])
		i += 1
	return SD_rotated

def dotvectors(pathA, pathB, rot_axe, rot_angle):
	path0 = pathA
	print(path0)

	user_ = pd.read_csv(path0+'/user_pos_allsatellites.csv',skiprows=1).values#.transpose()
	user_mean = mean(user_, axis=0).astype('float64')
	path0 = pathB
	stars = getSD(path0)
	SD = -array(stars[0]-stars[1]).astype('float64')
	SD = rot_SD_u_mean(SD, user_mean, rot_angle)
	return cumulativebyvalue(user_mean,SD,user_)
	




def getpaths(casenr):
	if casenr==0:
		#==================================================================================================== 24h_A
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/24h_delelott_all'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/starsby4000'
		pathB = path + '/allsatellites'
	if casenr==1:
		#==================================================================================================== 24h_B
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/24h_B'#/allsatellites'
		pathA = path + '/allsatellites_starsby20000'
		pathB = path + '/starsby5000'
	if casenr==2:
		#==================================================================================================== 24h_C
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/24h_C'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	if casenr==3:
		#==================================================================================================== 24h_D
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_D'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	if casenr==4:
		#==================================================================================================== 24h_E
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_E'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	if casenr==5:
		#==================================================================================================== 24h_F
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_F'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	
	if casenr==6:
		#==================================================================================================== 24h_G
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/24h_G'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/starsby5000'	
	if casenr==7:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_net'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'

	if casenr==8:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_net_augusztus/24h_H'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	if casenr==9:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_net_augusztus/24h_I'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	if casenr==10:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_net_augusztus/24h_J'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	if casenr==11:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_net_augusztus/24h_K'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	if casenr==12:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_net_augusztus/24h_L'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	if casenr==13:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_net_augusztus/24h_M'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	if casenr==14:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_net_augusztus/24h_N'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'

	if casenr==15:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/24h_B'#/allsatellites'
		pathA = path + '/allsatellites_45'
		pathB = path + '/allsatellites_45'

	if casenr==16:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/sept2_3_cluj/'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'
	
	if casenr==17:
		#==================================================================================================== 24h_net
		path = '/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/sept12_cluj'#/allsatellites'
		pathA = path + '/allsatellites'
		pathB = path + '/allsatellites'

	return pathA, pathB


def getpaths_2(root_directory):
	list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]
	data_dirs = []
	for data_dir in list_subfolders_with_paths:
		data_dirs.append(data_dir + "/allsatellites")
		#data_dirs.append(os.path(data_dir, "allsatellites"))
	return data_dirs

def time_hist(filepath):
	times = pd.read_csv(filepath + '/times_allsatellites.csv', skiprows=1).values/3600
	times = times - times.min()
	plt.hist(times, bins = 100)
	plt.show()


perth = [3, 4, 5, 7, 8, 9, 10, 11, 12, 13, 14]
cluj = [1, 6, 16]

rot_dicts = [
{"X": 0, "Y": 0, "Z": 0},
{"X": 30, "Y": 30, "Z": 30},
{"X": -30, "Y": -30, "Z": -30},
{"X": 90, "Y": 90, "Z": 90},
{"X": -90, "Y": -90, "Z": -90},
{"X": 15, "Y": 15, "Z": 15},
{"X": -10, "Y": -10, "Z": -10},
]

# paths = getpaths(16)
# dotvectors(paths[0],paths[1])

def process_calc_meancumulative(index_list, rot_dict):
	for rot_axe, rot_angle in rot_dict.items():
		data = []
		print (rot_axe, rot_angle)
		for x in index_list:
			paths = getpaths(x)
			data.append(dotvectors(paths[0],paths[1], rot_axe, rot_angle))

		# rot_axe = 'u_mean'
		calc_mean_cumulative(data, rot_axe, rot_angle)




def process_calc_meancumulative_from_rootdir(root_directory, rot_dict):
	print(rot_dict)
	for rot_axe, rot_angle in rot_dict.items():
		data = []
		print (rot_axe, rot_angle)
		dir_list = getpaths_2(root_directory)
		# dir_list = ["/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/24h_C/allsatellites"]
		for x in dir_list:
			lims, cumulative = dotvectors(x,x, rot_axe, rot_angle)
			data.append([lims, cumulative])
			plots_directory = create_directory(os.path.dirname(x), "cumulative_plots_rotSD")
			rot_axe = "SD_X_n"
			plot_cumulative(lims, cumulative, rot_axe, rot_angle, plots_directory)		# rot_axe = 'u_mean'
		
		calc_mean_cumulative(data, rot_axe, rot_angle, root_directory)


rot_list_meanuser = [*range(-90, +90, 5)]

rot_dicts_meanuser = {}
for i in rot_list_meanuser:
	rot_dicts_meanuser[str('SD_vs_SD2_'+str(i))] = i

#print(rot_dicts_meanuser)

root_directory = r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/24h_net_augusztus"
#/Users/kelemensz/Documents/Research/GPS/processed/24_h/RO/24h_C

process_calc_meancumulative_from_rootdir(root_directory, rot_dicts_meanuser)









