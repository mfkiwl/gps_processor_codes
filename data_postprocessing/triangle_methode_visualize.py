from numpy import *
import os
from data_postprocessing.data_locations_handler import AllGPSDataLocations
from utility.frecvently_used_functions import find_corresponding_dirs_in_different_roots, are_reliable, is_all_data, \
    get_csv_file, select_cmap_hist_n_mod, get_matrices_from_paths, handle_raw_not_averaged_matrices

# read_matrix(path)

f_h = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/test/20200803/data_GS_20200803_24h_histogram_5.csv"
f_d = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/test/20200803/data_GS_20200803_24h_cmap_5.csv"

# read_matrix(f_d, f_h)


#  ====================================================================================================================
#  ====================================================================================================================
#  ====================================================================================================================


def create_averaged_plots_from_root(root_0, months=None):
    sum_all_cmap = []
    sum_all_hist = []
    sum_all_n_mod = []
    subfolders_with_paths_months = [f.path for f in os.scandir(root_0) if f.is_dir()]
    for month_root in subfolders_with_paths_months:
        month_name = str(month_root).split("/")[-1]
        if months and month_name in months:

            days_with_paths = [f.path for f in os.scandir(month_root) if f.is_dir()]
            print("Month name: ", month_name, "  nr days: ", len(days_with_paths))
            for day_root in days_with_paths:
                csv_files = get_csv_file(day_root)
                cmap, hist, n_mod = select_cmap_hist_n_mod(csv_files)
                try:
                    if (cmap and hist and n_mod) and (
                            os.path.isfile(cmap) and os.path.isfile(hist) and os.path.isfile(n_mod)):
                        M, H, N = get_matrices_from_paths([cmap, hist, n_mod])

                        sum_all_cmap.append(M)
                        sum_all_hist.append(H)
                        sum_all_n_mod.append(N)
                except:
                    pass
    # print(hist, "\n", list(H[-1]), "\n", shape(array(H)))
    nr_days = len(sum_all_cmap)
    print("Total number of days:  ", nr_days)
    sum_all_cmap = sum(array(sum_all_cmap), axis=0)
    sum_all_hist = sum(array(sum_all_hist), axis=0)
    sum_all_n_mod = sum(array(sum_all_n_mod), axis=0)
    # plot_save_imshow_3_maps([sum_all_cmap, sum_all_hist, sum_all_n_mod], ["|1-1/r|", "Histogram", "<n_mod>"],
    #                         root_0, resolution="5")
    # plot_save_imshow(sum_all_cmap, root_0, "|1-1/r|")
    # print(sum_all_hist)
    return sum_all_cmap, sum_all_hist, sum_all_n_mod, nr_days


def create_averaged_plots_from_root_and_filter_by_positionsSTD(root_0, needed_files, months=None):
    sum_all_cmap = []
    sum_all_hist = []
    sum_all_n_mod = []
    locationA, locationB = root_0.split("/")[-2].split('_')
    path_A = AllGPSDataLocations.user_and_satellites.get(locationA)
    path_B = AllGPSDataLocations.user_and_satellites.get(locationB)
    print(path_A, path_B)
    month_pairs = find_corresponding_dirs_in_different_roots(path_A, path_B)
    subfolders_with_paths_months = [f.path for f in os.scandir(root_0) if f.is_dir()]
    for A_month, B_month in month_pairs:
        print(A_month)
        for month_root in subfolders_with_paths_months:
            month_name_of_position_data = os.path.split(A_month)[-1]
            month_name_triangle_data = str(month_root).split("/")[-1]

            if months and month_name_of_position_data == month_name_triangle_data:
                days_with_paths = [f.path for f in os.scandir(month_root) if f.is_dir()]
                day_pairs = find_corresponding_dirs_in_different_roots(A_month, B_month)
                # print("Month name: ", month_name_triangle_data, "  nr days (triangle): ", len(days_with_paths),
                #       'position data:', len(day_pairs))
                for A_day, B_day in day_pairs:
                    for day_root in days_with_paths:
                        date_from_positions = str(os.path.split(B_day)[-1])[-8:]
                        date_from_triangle_results = str(os.path.split(day_root)[-1])

                        if date_from_positions == date_from_triangle_results:
                            # print('Curent day: ', A_day, is_all_data(A_day, needed_files, add_allsatellites=True))

                            if are_reliable(A_day, B_day, needed_files):
                                # print('Reliable: ', date_from_positions)

                                csv_files = get_csv_file(day_root)
                                cmap, hist, n_mod = select_cmap_hist_n_mod(csv_files)
                                try:
                                    if (cmap and hist and n_mod) and (
                                            os.path.isfile(cmap) and os.path.isfile(hist) and os.path.isfile(n_mod)):
                                        M, H, N = get_matrices_from_paths([cmap, hist, n_mod])

                                        sum_all_cmap.append(M)
                                        sum_all_hist.append(H)
                                        sum_all_n_mod.append(N)
                                except:
                                    pass
    nr_days = len(sum_all_cmap)
    print("Total number of days:  ", nr_days)
    sum_all_cmap = sum(array(sum_all_cmap), axis=0)
    sum_all_hist = sum(array(sum_all_hist), axis=0)
    sum_all_n_mod = sum(array(sum_all_n_mod), axis=0)
    return sum_all_cmap, sum_all_hist, sum_all_n_mod, nr_days




# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/PERTH_NZLD/r_inv_r_symmetrized"

# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/PERTH_NASA/r_inv_r_symmetrized"

# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/HKKS_NASA/r_inv_r_symmetrized"

# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/HKKS_PERTH/r_inv_r_symmetrized"

# results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data/HKKS_IIGC/r_inv_r_symmetrized"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_method/PERTH_IIGC/r_inv_r_symmetrized"

# results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data/NZLD_TIDV/r_inv_r_symmetrized"

# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/CUTA_NZLD/r_inv_r_symmetrized"

# results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data/NASA_IIGC/r_inv_r_symmetrized"

# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/NZLD_NASA/r_inv_r_symmetrized"

# results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data/NZLD_IIGC/r_inv_r_symmetrized"

# results_root = r"/Volumes/KingstonSSD/GPS/processed_data/triangular_method/processed_data/NZLD_HKKS/r_inv_r_symmetrized"

# results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data/PERTH_TIDV/r_inv_r_symmetrized"

# results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data/HKKS_TIDV/r_inv_r_symmetrized"

# results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data/NASA_TIDV/r_inv_r_symmetrized"

# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_method/CUTB30s_NZLD/r_inv_r_symmetrized"

# results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data/CUTA_CUTB/r_inv_r_symmetrized"

# results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data/PERTH_IIGC/r_inv_r_symmetrized"

# results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data/IIGC_TIDV/r_inv_r_symmetrized"


# pair = results_root.split("/")[-2]
# print(pair)
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_test"

all_months = ["julius", "szeptember", "februar", "marcius", "augusztus", "januar", "december2019", "oktober",
              "november", "majus", "aprilis", "junius", "december2020"]
months1 = ["julius", "szeptember", "augusztus", "november", "junius", "december2020"]

# fig_dir = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/figures_filtered_histograms"
results_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_divided_by_n/CUTA_CUTB/r_inv_r_symmetrized"

pair = results_root.split("/")[-2]
print(pair)
# fig_dir = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_divided_by_n"
# m, h, n, n_days = create_averaged_plots_from_root_and_filter_by_positionsSTD(results_root, needed_files=[
#     "user_pos_allsatellites.csv"], months=all_months)
# handle_raw_not_averaged_matrices(m, h, n, fig_directory=fig_dir, name=pair + "_int", nr_days=n_days, round=True)
# handle_raw_not_averaged_matrices(m, h, n, fig_directory=fig_dir, name=pair, nr_days=n_days, round=False)

result_roots = [
    r"CUTA_CUTB/r_inv_r_symmetrized", r"PERTH_NZLD/r_inv_r_symmetrized", r"PERTH_NASA/r_inv_r_symmetrized",
    r"HKKS_PERTH/r_inv_r_symmetrized", r"HKKS_IIGC/r_inv_r_symmetrized", r"PERTH_IIGC/r_inv_r_symmetrized",
    r"NZLD_TIDV/r_inv_r_symmetrized", r"CUTA_NZLD/r_inv_r_symmetrized", r"NASA_IIGC/r_inv_r_symmetrized",
    r"NZLD_NASA/r_inv_r_symmetrized", r"NZLD_IIGC/r_inv_r_symmetrized", r"NZLD_HKKS/r_inv_r_symmetrized",
    r"PERTH_TIDV/r_inv_r_symmetrized", r"HKKS_TIDV/r_inv_r_symmetrized", r"NASA_TIDV/r_inv_r_symmetrized",
    r"CUTB30s_NZLD/r_inv_r_symmetrized", r"PERTH_IIGC/r_inv_r_symmetrized", r"HKKS_NASA/r_inv_r_symmetrized",
    r"IIGC_TIDV/r_inv_r_symmetrized"]

main_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data"
# main_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_divided_by_n"
fig_dir = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_divided_by_n"


def save_imshows_for_all(main_root, result_roots, fig_dir):
    needed_files = ["user_pos_allsatellites.csv"]

    for result_root in result_roots[:]:
        root = os.path.join(main_root, result_root)
        if not os.path.isdir(root):
            print(root)
            continue
        print(root, '\n', '\n')
        pair = result_root.split("/")[-2]
        print(pair)
        # m, h, n, n_days = create_averaged_plots_from_root(root, all_months)
        m, h, n, n_days = create_averaged_plots_from_root_and_filter_by_positionsSTD(root, needed_files, all_months)
        try:
            handle_raw_not_averaged_matrices(m, h, n, fig_directory=fig_dir, name=pair + "_int", nr_days=n_days,
                                             round=True)
            handle_raw_not_averaged_matrices(m, h, n, fig_directory=fig_dir, name=pair, nr_days=n_days, round=False)
        except:
            pass


# save_imshows_for_all(main_root, result_roots, fig_dir)



# M = rotate(M, -90)
# plot_mollweid_simple(M)


# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/HKKS_PERTH/r_inv_r"
# results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangular_method/processed_data/PERTH_NZLD/r_inv_r"
# plot_the_three_raw_matrix_from_path(results_root)

# =================================================================================================================================================================
# =================================================================================================================================================================
# =================================================================================================================================================================

