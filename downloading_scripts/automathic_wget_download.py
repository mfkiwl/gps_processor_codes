import os




def get_folder_names(destination_path):
	return [f for f in os.listdir(destination_path) if os.path.isdir(os.path.join(destination_path, f))]
    
def get_day_folder_names(n = 365):
	a = []
	for i in range(1, n):
		istr = str(i)
		if len(istr) == 1:
			istr = "00" + istr
		if len(istr) == 2:
			istr = "0" + istr
		a.append(istr)
	return a


def download_with_wget(folders_to_be_downloaded):
	for folder in folders_to_be_downloaded:
		command = "wget --ftp-user anonymous -r -A 'ARTA*' ftp://ftp.geonet.org.nz/rtgps/rinex1Hz/PositioNZ/2020/{}/".format(folder)
		os.system(command)




# get_day_folder_names()
destination_path = r"/Users/kelemensz/Documents/Research/GPS/keep_from_qsync/ARTA/ftp.geonet.org.nz/rtgps/rinex1Hz/PositioNZ/2020"


days_to_be_downloaded = list(set(get_day_folder_names()).difference(set(get_folder_names(destination_path))))
print(len(days_to_be_downloaded))
print(days_to_be_downloaded)

download_with_wget(days_to_be_downloaded)












