
import numpy as np


# posfile = r"/Users/kelemensz/Documents/Research/GPS/process/24h/iono_correction/without/nocor_CUTB2020120724.pos"
# outfile = r"/Users/kelemensz/Documents/Research/GPS/process/24h/iono_correction/without/CUTB2020120724.csv"

posfile = r"/Users/kelemensz/Documents/Research/GPS/process/24h/iono_correction/with_to_send/CUTB20200104/CUTB2020010424.pos"
outfile = r"/Users/kelemensz/Documents/Research/GPS/process/24h/iono_correction/with_to_send/CUTB20200104/user_pos_allsatellites.csv"


def posparser(file, outfile):
	a = []
	with open(file) as in_file, open(outfile, 'w') as out_file:
		lineList = [line.rstrip('\n').split(" ") for line in in_file]
		for line in lineList[15:]:
			# print(len(line))
			# print(line)
			# print(line[3], line[6], line[8])
			# out_file.write("{},{},{}\n".format(line[3], line[6], line[8]))
			a.append(line[3], line[6], line[8])
		
	return np.array(a)

posparser(posfile, outfile)




# with open(in_filename) as in_file, open(out_filename, 'w') as out_file:
#    for line in in_file:
#      ...
#      ... 
#      out_file.write(parsed_line)
