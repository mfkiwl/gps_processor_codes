import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import math
import pylab as pl








def get_matrix(path):
	# print('Path: ', path)
	M = pd.read_csv(path, skiprows=0).values#.transpose()[0]#.astype('float64')	
	return M


def get_csv_file(directory):
  csv_ext = ".csv"
  files = next(os.walk(directory))[2]
  # print(files)
  files_with_path = []
  for file in files:
  	if os.path.splitext(file)[1] == csv_ext:
  		files_with_path.append(os.path.abspath(os.path.join(directory, file)))
  return files_with_path



def plot_cmap(color_matrix, angle_lim1, angle_lim2, root_directory):
	plt.clf()
	try:
		cmap_save = pd.DataFrame(color_matrix)
		cmap_save.to_csv(os.path.join(root_directory, 'dat_cmap_by5degree.csv'), index=False)
		print('Matrix.csv saved!')
	except:
		pass
	fig = plt.imshow(color_matrix, extent=[0,angle_lim1,0,angle_lim2], aspect='auto')
	
	plt.title('24h_average, 32 days')
	
	plt.ylabel('Theta')
	plt.xlabel('Phi')
	plt.colorbar(fig)
	plt.show()
	# fig_name1 = os.path.join(root_directory, 'cmap_by5degree.png')
	# fig_name2 = os.path.join(root_directory, '_cmap_by5degree.png')
	
	# try:
	# 	plt.savefig(fig_name1)
	# except:
	# 	pass
	# try:
	# 	plt.imsave(fig_name2, fig)
	# except:
	# 	pass


def plot_mollweid(matrix, S, D, root_directory=None):
	# try:
	# 	cmap_save = pd.DataFrame(color_matrix)
	# 	cmap_save.to_csv(os.path.join(root_directory, 'data_galactic_system_byXdeg.csv'), index=False)
	# 	print('Matrix.csv saved!')
	# except:
	# 	print()

	ra = np.linspace(-math.pi, math.pi, len(matrix))
	dec= np.linspace(-math.pi/2, math.pi/2, len(matrix[0]))

	X,Y = np.meshgrid(ra, dec)
	Z = matrix.T
	plt.figure()
	ax = pl.subplot(111)#, projection = 'mollweide')
	fig = ax.contourf(X,Y,Z,100)
	print(S[0],S[1])
	print(D[0],D[1])
	x = np.radians(np.array([S[0], D[0]]))
	y = np.radians(np.array([S[1], D[1]]))
	print(x, y)

	ax.scatter(x,y, marker='x',c='k',s=15)
	
	# ax.annotate('Sadalmelic direction',
 #            xy=np.radians(S),
 #            xycoords='data',
 #            xytext=np.radians(S - np.array([-10, 20])),
 #            arrowprops=
 #                dict(facecolor='black', shrink=0.05),
 #                horizontalalignment='left',
 #                verticalalignment='top')
	# ax.annotate('Denebola direction',
 #            xy=np.radians(D),
 #            xycoords='data',
 #            xytext=np.radians(D+20),
 #            arrowprops=
 #                dict(facecolor='black', shrink=0.05),
 #                horizontalalignment='left',
 #                verticalalignment='top')


	# ax.set_title('---$(24h)$', fontsize=15)  # , fontweight='bold')
	# plt.xlabel(r'$\theta$',fontsize=15)#Italic font method
	# plt.ylabel(r'$\phi$',fontsize=15)#Bold font method without fontweight parameters
	# pl.colorbar(fig)
	ax.grid()
	# ax.contour(X,Y,Z,10,colors='k')
	pl.show()

	# try:
	#	fig_name1 = os.path.join(root_directory, 'galactic_system_byXdeg.png')
	# 	pl.savefig(fig_name1, bbox_inches='tight')
	# except:
	# 	print()


def plot_mollweid_save(matrix, S, D, root_directory=None):
	try:
		cmap_save = pd.DataFrame(matrix)
		cmap_save.to_csv(os.path.join(root_directory, 'data_2d.csv'), index=False)
		print('Matrix.csv saved!')
	except:
		print()

	plt.clf()
	ra = np.degrees(np.linspace(-math.pi, math.pi, len(matrix)))
	dec= np.degrees(np.linspace(-math.pi/2, math.pi/2, len(matrix[0])))

	X,Y = np.meshgrid(ra, dec)
	Z = matrix.T
	plt.figure()
	ax = pl.subplot(111)#, projection = 'mollweide')
	fig = ax.contourf(X,Y,Z,100)
	print(S[0],S[1])
	print(D[0],D[1])
	# x = np.radians(np.array([S[0], D[0]]))
	# y = np.radians(np.array([S[1], D[1]]))
	x = np.array([S[0], D[0]])
	y = np.array([S[1], D[1]])
	print(x, y)

	ax.scatter(x,y, marker='o',c='k',s=200)
	
	ax.annotate('Sadalmelic', fontsize=15,
            xy=np.radians(S),
            xycoords='data',
            xytext=np.radians(S + np.array([20, 20])),
            arrowprops=
                dict(facecolor='black', shrink=0.05),
                horizontalalignment='left',
                verticalalignment='top')
	ax.annotate('Denebola', fontsize=15,
            xy=np.radians(D),
            xycoords='data',
            xytext=np.radians(D-20),
            arrowprops=
                dict(facecolor='black', shrink=0.05),
                horizontalalignment='left',
                verticalalignment='top')


	# ax.set_title('---$(24h)$', fontsize=15)  # , fontweight='bold')
	plt.xlabel(r'$\theta$',fontsize=15)#Italic font method
	plt.ylabel(r'$\phi$',fontsize=15)#Bold font method without fontweight parameters
	pl.colorbar(fig)
	ax.grid()
	# ax.contour(X,Y,Z,10,colors='k')
	pl.show()
	# child_dirname = os.path.split(root_directory)[-1]+'_24h'
	# fig_name1 = os.path.join(root_directory, 'GS_'+child_dirname+'.png')
	# pl.savefig(fig_name1, bbox_inches='tight')
	# try:
	# 	fig_name1 = os.path.join(root_directory, 'GCS_.png')
	# 	pl.savefig(fig_name1, bbox_inches='tight')
	# except:
	# 	print()


def process_matrices(root_directory):
	csv_files = get_csv_file(root_directory)
	M_all = np.zeros(np.shape(get_matrix(csv_files[0])))
	
	n=0
	for file in csv_files:

		try:
			n += 1
			M_all = M_all + np.array(get_matrix(file))
		except:
			pass
	
	M = M_all/float(n)
	# M = np.array(get_matrix(file))
	# print(M)
	# plot_cmap(M,180,360,None)
	# D = np.array((70.73453543, 70.73453543))
	# S = np.array((-59.88270525, -42.0670648))
	S = np.array(np.flip([-59.94133321, -89.99984419]))
	D = np.array(np.flip([109.37248674,  89.99993096]))

	plot_mollweid(M, S, D)



def rebin(arr, new_shape):
    shape = (new_shape[0], arr.shape[0] // new_shape[0],
             new_shape[1], arr.shape[1] // new_shape[1])
    return arr.reshape(shape).mean(-1).mean(1)


# root_directory = r"/Users/kelemensz/Documents/Research/GPS/process/24h/RO/cmap/GS_matrix"
root_directory = r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/galactic_sys_matrices"

# root_directory = "/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/mean_projection/febr_9_15/febr_9"
# process_matrices(root_directory)


file = r"/Users/kelemensz/Documents/Research/GPS/processed/24_h/Perth/daily_measurements/data_GS_daily_measurements_24h.csv"
# m=np.array(get_matrix(file))
# S = np.array(np.flip([-59.94133321, -89.99984419]))
# D = np.array(np.flip([109.37248674,  89.99993096]))

# plot_mollweid(m, S, D)

def get_mean_matrix(root_directory):
	print(root_directory)
	list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]
	M_all = [] #np.zeros(np.shape(get_matrix(csv_files[0])))

	# for months in list_subfolders_with_paths_root:
	# 	list_subfolders_with_paths = [f.path for f in os.scandir(months) if f.is_dir()]

	for directory in list_subfolders_with_paths:
		csv_files = get_csv_file(directory)
		folder_name = directory.split("/")[-1]
		# print('------',folder_name)
		for file in csv_files:
			# a =file.split("/")[-1]
			# print ("matrix file: ", file) 
			if file.split("/")[-1].find(folder_name) != -1:#folder_name in file.split("/")[-1]
				if file.split("/")[-1].find("SD") == -1:#folder_name in file.split("/")[-1]				
					M_all.append(get_matrix(file))	
					# print(folder_name)
		
	# M = M_all/float(len(M))
	print("Day number: ", len(M_all))
	M = np.mean(np.array(M_all), axis=0)
	# print(M)
	# a = 8
	# b = a * 2
	# M = rebin(M, (int(len(M)/b), int(len(M[0])/a)))
	# imshow_bins(M)
	return M

def imshow_bins(M, S, D):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	plt.imshow(M.T, origin = 'lower',  extent = [0, len(M), 0, len(M[0])], aspect='auto')

	

	# ax.annotate('Sadalmelic', fontsize=20, xy=(1.5, 0.5),
	#             xycoords='data', xytext=(1, -1),
	#             textcoords='offset points',
	#             arrowprops=dict(arrowstyle="->",
	#                             linewidth = 5.,
	#                             color = 'red')
	#             )
	# ax.annotate('Denebola', fontsize=20, xy=(3.5, 3.5),
	#             xycoords='data', xytext=(-1, -0.5),
	#             textcoords='offset points',
	#             arrowprops=dict(width = 5.,
	#                             headwidth = 15.,
	#                             frac = 0.2,
	#                             shrink = 0.05,
	#                             linewidth = 2,
	#                             color = 'red')
	#             )
	
	plt.show()
	# close()


def process_matrices(root_directory):
	print(root_directory)
	list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]
	M_all = [] #np.zeros(np.shape(get_matrix(csv_files[0])))

	# for months in list_subfolders_with_paths_root:
	# 	list_subfolders_with_paths = [f.path for f in os.scandir(months) if f.is_dir()]

	for directory in list_subfolders_with_paths:
		csv_files = get_csv_file(directory)
		folder_name = directory.split("/")[-1]
		# print('------',folder_name)
		for file in csv_files:
			# a =file.split("/")[-1]
			# print ("matrix file: ", file) 
			if file.split("/")[-1].find(folder_name) != -1:#folder_name in file.split("/")[-1]
				if file.split("/")[-1].find("SD") == -1:#folder_name in file.split("/")[-1]				
					M_all.append(get_matrix(file))	
					# print(folder_name)
					# print(np.shape(M_all[-1]))
		
	# M = M_all/float(len(M))
	print("Day number: ", len(M_all))
	M = np.mean(np.array(M_all), axis=0)
	
	# S = np.array(np.flip([-59.94133321, -89.99984419]))
	# D = np.array(np.flip([109.37248674,  89.99993096]))

	S = np.array([-59.94133321, -42.06696326968])
	D = np.array([109.37248674,  70.80087410094148])
	# print(M)
	# a = 8
	# b = a * 2
	# M = rebin(M, (int(len(M)/b), int(len(M[0])/a)))
	# print(M)
	# plt.imshow(M.T  , origin = 'lower',  extent = [0, len(M), 0, len(M[0])], aspect = 'auto')
	# plt.colorbar()

	# plt.show()
	# imshow_bins(M, S, D)
	# cmap_save = pd.DataFrame(M)
	# cmap_save.to_csv(os.path.join(root_directory, 'data_2d.csv'), index=False)
	plot_mollweid_save(M, S, D)
	return M


root_dir = r"/Users/kelemensz/Dropbox/GPS_data/Cumulative_plots/mean_projection/perth/daily_measurements"
root_dir = r"/Users/kelemensz/Documents/Research/GPS/process/24h/auto_results"
root_dir = r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/PERTH_daily_measurements"
root_dir = r"/Users/kelemensz/Documents/Research/GPS/process/24h/global_GCS_axis/process_NZLD"
# root_dir = r"/Users/kelemensz/Qsync/GPS/GPS_data/local_GCS_axis/process_NZLD"
process_matrices(root_dir)
# get_mean_matrix(root_dir)


def process_matrices_2(root_directory):
	print(root_directory)
	list_subfolders_with_paths = [f.path for f in os.scandir(root_directory) if f.is_dir()]
	M_all = [] #np.zeros(np.shape(get_matrix(csv_files[0])))

	# for months in list_subfolders_with_paths_root:
	# 	list_subfolders_with_paths = [f.path for f in os.scandir(months) if f.is_dir()]

	for directory in list_subfolders_with_paths:
		csv_files = get_csv_file(directory)
		folder_name = directory.split("/")[-1]
		# print('------',folder_name)
		for file in csv_files:
			# a =file.split("/")[-1]
			# print ("matrix file: ", file) 
			if file.split("/")[-1].find(folder_name) != -1:#folder_name in file.split("/")[-1]
				if file.split("/")[-1].find("SD") == -1:#folder_name in file.split("/")[-1]				
					M_all.append(get_matrix(file))	
					# print(folder_name)
		
	# M = M_all/float(len(M))
	print("Day number: ", len(M_all))
	M = np.mean(np.array(M_all), axis=0)
	print(np.shape(M))

	plt.plot(np.mean(M, axis=0))
	# plt.set_title("".format(len(M)))
	print("".format(len(M)))
	plt.show()

	plt.clf()
	plt.plot(np.mean(M, axis=1))
	# plt.set_title("".format(len(M[0])))
	print("".format(len(M[0])))
	plt.show()


root_dir = r"/Users/kelemensz/Documents/Research/GPS/process/24h/auto_results"

# process_matrices_2(root_dir)


def count_days(root_directory):
	print(root_directory)
	list_subfolders_with_paths_root = [f.path for f in os.scandir(root_directory) if f.is_dir()]
	M_all = [] #np.zeros(np.shape(get_matrix(csv_files[0])))

	for months in list_subfolders_with_paths_root:
		list_subfolders_with_paths = [f.path for f in os.scandir(months) if f.is_dir()]

		for directory in list_subfolders_with_paths:
			csv_files = get_csv_file(directory)
			folder_name = directory.split("/")[-1]
			# print('------',folder_name)
			for file in csv_files:
				# a =file.split("/")[-1]
				# print ("matrix file: ", file) 
				if file.split("/")[-1].find(folder_name) != -1:#folder_name in file.split("/")[-1]
					if file.split("/")[-1].find("SD") == -1:#folder_name in file.split("/")[-1]				
						M_all.append(get_matrix(file))	
						# print(folder_name)
		
	# M = M_all/float(len(M))
	print("Day number: ", len(M_all))

# count_days(root_dir)





# =====================================================================================
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm


def random_point( r=1 ):
    ct = 2*np.random.rand() - 1
    st = np.sqrt( 1 - ct**2 )
    phi = 2* np.pi *  np.random.rand()
    x = r * st * np.cos( phi)
    y = r * st * np.sin( phi)
    z = r * ct
    return np.array( [x, y, z ] )

def near( p, pntList, d0 ):
    cnt=0
    for pj in pntList:
        dist=np.linalg.norm( p - pj )
        if dist < d0:
            cnt += 1 - dist/d0
    return cnt


"""
https://stackoverflow.com/questions/22128909/plotting-the-temperature-distribution-on-a-sphere-with-python
"""
def plot_on_sphere():
	pointList = np.array([ random_point( 10.05 ) for i in range( 65 ) ] )

	fig = plt.figure()
	ax = fig.add_subplot( 1, 1, 1, projection='3d')
	WW = get_mean_matrix(root_dir)

	u = np.linspace( 0, 2 * np.pi, len(WW))
	v = np.linspace( 0, np.pi, len(WW[0]) )

	# create the sphere surface
	XX = 10 * np.outer( np.cos( u ), np.sin( v ) )
	YY = 10 * np.outer( np.sin( u ), np.sin( v ) )
	ZZ = 10 * np.outer( np.ones( np.size( u ) ), np.cos( v ) )

	# WW = XX.copy()
	# for i in range( len( XX ) ):
	#     for j in range( len( XX[0] ) ):
	#         x = XX[ i, j ]
	#         y = YY[ i, j ]
	#         z = ZZ[ i, j ]
	#         WW[ i, j ] = near(np.array( [x, y, z ] ), pointList, 3)
	# WW = WW / np.amax( WW )
	# print (np.shape(WW))
	WW = WW + abs(np.amin( WW ))

	myheatmap = WW/np.amax( WW ) 

	# ~ ax.scatter( *zip( *pointList ), color='#dd00dd' )
	ax.plot_surface( XX, YY,  ZZ, cstride=1, rstride=1, facecolors=cm.jet( myheatmap ) )
	# plt.colorbar(cm.jet( myheatmap ))
	plt.show()


# plot_on_sphere()


# =====================================================================================
