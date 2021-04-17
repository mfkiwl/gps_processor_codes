import os


def setup_triangle_results_locations(main_root, result_roots):
    results = [os.path.join(main_root, result_root) for result_root in result_roots]
    identifiers = [result_root.split("/")[-2] for result_root in result_roots]
    return dict(zip(identifiers, results))


class AllGPSDataLocations:
    OBS_file_locations = {
        'CUTB': '',
        'CUTA': '',
        'NZLD': '',
        'HKKS': '',
        'IISC': '',
        'NASA': '',
        'TIDV': ''
    }

    user_and_satellites = {
        'CUTB': r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements",
        'PERTH': r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements",
        'CUTA': r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements/CUTA",
        'NZLD': r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_NZLD",
        'HKKS': r"/Volumes/BlueADATA S/GPS/processed_data/global_GCS_axis/process_HKKS",
        'IIGC': r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_IIGC",
        'IISC': r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_IISC",
        'NASA': r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_NASA",
        'TIDV': r"/Volumes/BlueADATA S/GPS/processed_data/global_GCS_axis/process_TIDV",
        'CUTB30s': r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/CUTB2020_30s_rinex2"
    }

    # # --------------------------------------------CUTB30s_rinexbol-NZLD-symmetrized-------------------------------------------
    # place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/CUTB2020_30s_rinex2"
    # place_B = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/process_NZLD"
    # results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_method/CUTB30s_NZLD/tests"

    # # --------------------------------------------CUTB30s_rinexbol-CUTB-symmetrized-------------------------------------------
    # place_A = r"/Users/kelemensz/Documents/Research/GPS/process/global_GCS_axis/CUTB2020_30s_rinex2"
    # place_B = r"/Volumes/KingstonSSD/GPS/processed_data/user_and_sat_positions_and_ionospheric_effects/PERTH_daily_measurements"
    # results_root = r"/Users/kelemensz/Documents/Research/GPS/process/triangle_method/CUTB30s_CUTB/r_inv_r_symmetrized"

    result_roots = [
        r"CUTA_CUTB/r_inv_r_symmetrized", r"PERTH_NZLD/r_inv_r_symmetrized", r"PERTH_NASA/r_inv_r_symmetrized",
        r"HKKS_PERTH/r_inv_r_symmetrized", r"PERTH_IISC/r_inv_r_symmetrized", r"PERTH_IIGC/r_inv_r_symmetrized",
        r"NZLD_TIDV/r_inv_r_symmetrized", r"CUTA_NZLD/r_inv_r_symmetrized", r"NASA_IIGC/r_inv_r_symmetrized",
        r"NZLD_NASA/r_inv_r_symmetrized", r"NZLD_IIGC/r_inv_r_symmetrized", r"NZLD_HKKS/r_inv_r_symmetrized",
        r"PERTH_TIDV/r_inv_r_symmetrized", r"HKKS_TIDV/r_inv_r_symmetrized", r"NASA_TIDV/r_inv_r_symmetrized",
        r"CUTB30s_NZLD/r_inv_r_symmetrized", r"HKKS_NASA/r_inv_r_symmetrized", r"HKKS_IIGC/r_inv_r_symmetrized",
        r"IIGC_TIDV/r_inv_r_symmetrized"]

    main_root = r"/Volumes/BlueADATA S/GPS/processed_data/triangular_method/processed_data"

    triangle_results_locations = setup_triangle_results_locations(main_root, result_roots)

    all_months = ["julius", "szeptember", "februar", "marcius", "augusztus", "januar", "december2019", "oktober",
                  "november", "majus", "aprilis", "junius", "december2020"]


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
