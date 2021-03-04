import os


def run_convbin(day_binex, output_dir):
    file_name = str(os.path.split(day_binex)[-1]).split(".")[0]

    if str(day_binex).split(".")[-1] != "gz":

        nav_file = os.path.join(output_dir, file_name + ".nav")
        obs_file = os.path.join(output_dir, file_name + ".obs")
        pos_file = os.path.join(output_dir, file_name + ".pos")

        tool = os.path.normpath(
            r"C:\SzabolcsKelemen\PhD\GPS\RTKLIB\RTKLIB_bin-rtklib_2.4.3\RTKLIB_bin-rtklib_2.4.3\bin\convbin")
        tool2 = os.path.normpath(
            r"C:\SzabolcsKelemen\PhD\GPS\RTKLIB\RTKLIB_bin-rtklib_2.4.3\RTKLIB_bin-rtklib_2.4.3\bin\rnx2rtkp")

        command1 = "{} -r binex -f 3 -v 2.10 -od -os -oi -ot -ol -halfc -o {} -n {} {}".format(tool, obs_file, nav_file,
                                                                                               day_binex)
        command2 = "{} -o {} -p 0 -f 3 -e {} {}".format(tool2, pos_file, obs_file, nav_file)

        if not os.path.isfile(obs_file):
            print(day_binex)
            os.system(command1)
            os.system(command2)


def run_converting(root_path, output_dir):
    source_folders = [os.path.join(root_path, f) for f in os.listdir(root_path) if
                      os.path.isfile(os.path.join(root_path, f))]
    for folder in source_folders:
        run_convbin(folder, output_dir)


input_ = r"C:\SzabolcsKelemen\PhD\GPS\PERTH\CUTA"
output = r"C:\SzabolcsKelemen\PhD\GPS\PERTH\CUTA_obs"

run_converting(input_, output)
