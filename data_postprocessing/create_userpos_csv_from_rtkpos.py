import os
import pandas as pd

from utility.frecvently_used_functions import create_dir, create_generalinfo_path, is_all_data, \
    get_pos_files_simple, posparser


def extract_user_position_from_rtklib_pos_file(obs_location, root_directory, needed_files):
    pos_files = get_pos_files_simple(obs_location)
    month_name = os.path.split(obs_location)[-1]
    # root_directory = os.path.join(root_directory, month_name)
    root_directory = create_dir(root_directory, month_name)
    total_processed = 0
    # print(obs_files)
    for pos_file in pos_files[:]:
        pos_file_name = os.path.splitext(pos_file)[-2].split('/')[-1][:-2]
        # print("Pos filename: ", pos_file_name)
        path_of_results = create_dir(root_directory, pos_file_name)
        generalinfo_path = create_generalinfo_path(path_of_results)
        print(pos_file, generalinfo_path)
        if os.path.isdir(generalinfo_path) and not is_all_data(generalinfo_path, needed_files[:1]):
            # rinex_processed_grouped = get_processed_data(pos_file)

            # calc_user_pos(rinex_processed_grouped, 'allsatellites', generalinfo_path)
            outfile = os.path.join(generalinfo_path, needed_files[0])
            user_positions = pd.DataFrame(posparser(pos_file))
            user_positions.to_csv(outfile, index=False)

            total_processed += 1
        print("\n                          Processed: {}/{} \n".format(total_processed, len(pos_files)))


needed_files = ["user_pos_allsatellites.csv", '']

#                   PERTH_CUTA
# obs_path = r"/Volumes/BlueADATA S/GPS/raw_data/PERTH/CUTA_obs/aprilis"
obs_path = r"/Volumes/ADATA SE800/GPS/raw_data/PERTH/CUTA_obs/aprilis"
destination_path = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements/CUTA"


extract_user_position_from_rtklib_pos_file(obs_path, destination_path, needed_files)
