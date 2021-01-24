import os



def run_convbin(day_dir, output_dir):

    command1 = "python process_one_obs_file_triangle_by_time.py"
    os.system(command1)


def run():
    max_days = 31
    for day in range(max_days):
    	run_convbin()


run()






