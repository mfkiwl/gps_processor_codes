import os

root_directory = "/Users/kelemensz/Documents/Research/GPS/process/24h/perth/perth_september_october"



directories = ['sept27', 'sept28', 'sept29', 'sept30', 'oct01', 'oct02', 'oct03', 'oct04', 'oct05', 'oct06', 'oct07', 'oct08', 'oct09', 'oct10']

for subdir in directories:
	sub = os.path.join(root_directory, subdir)
	if not os.path.isdir(sub):
		os.makedirs(sub)