import os
from datetime import datetime


def get_subfolders_with_source(root_path):
    return [os.path.join(root_path, f) for f in os.listdir(root_path) if os.path.isdir(os.path.join(root_path, f))]


def convert_number_to_date(num_str, year="2020"):
    if num_str[0] == "0":
        num_str = num_str[1:]
        if num_str[0] == "0":
            num_str = num_str[1:]
    date = datetime.strptime(year + "-" + num_str, "%Y-%j").strftime("%Y%m%d")
    return date


def run_convbin(day_dir, output_dir):
    day_index = os.path.split(day_dir)[-1]
    day_dir = os.path.join(day_dir, "hkks/1s")

    file_name = "HKKS" + convert_number_to_date(day_index) + "24"
    input_files = os.path.normpath(str(day_dir) + "/*.20o")
    # print(input_files)
    nav_file = os.path.join(output_dir, file_name + ".nav")
    obs_file = os.path.join(output_dir, file_name + ".obs")
    pos_file = os.path.join(output_dir, file_name + ".pos")
    tool = os.path.normpath(r"C:\SzabolcsKelemen\PhD\GPS\RTKLIB\RTKLIB_bin-rtklib_2.4.3\RTKLIB_bin-rtklib_2.4.3\bin\convbin")
    command1 = "{} -r rinex -f 3 -v 2.10 -od -os -oi -ot -ol -halfc -o {} -n {} {}".format(tool, obs_file, nav_file,
                                                                                                input_files)
    os.system(command1)

    command2 = "rnx2rtkp -o {} -p 0 -f 3 -e {} {}".format(pos_file, obs_file, nav_file)
    # os.system(command2)


def run_converting(root_path, output_dir):
    binexes = [os.path.join(root_path, f) for f in os.listdir(root_path) if os.path.isdir(os.path.join(root_path, f))]
    for folder in binexes[2:]:
        print(folder)
        run_convbin(folder, output_dir)


input_ = r"C:\SzabolcsKelemen\PhD\GPS\HKKS\downloaded"
output = r"C:\SzabolcsKelemen\PhD\GPS\HKKS\obs_files"

run_converting(input_, output)
