import os
from time import sleep


def run_convbin():
    command1 = "python process_one_obs_file_triangle_by_time.py"
    os.system(command1)


def run():
    max_days = 31
    for day in range(max_days):
        run_convbin()
        print("\n\n\nWaiting 30 seconds before the next execution!\nExecution status: {}/{}\n\n".format(day, max_days))
        sleep(30)


run()
