from data_postprocessing.postprocess_utility.general_functions import *


gcs = array([[-1, 0, 0], [0, 1, 0], [0, 0, 1]])

vector_ecef = array([-1, 0, 0])

vector_gcf = ecef_to_gcs(gcs, vector_ecef)

print(vector_gcf)