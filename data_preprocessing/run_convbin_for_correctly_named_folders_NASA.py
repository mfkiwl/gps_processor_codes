import os


def run_convbin(day_dir, output_dir):
    file_name = os.path.split(day_dir)[-1]
    print(day_dir)
    input_files = os.path.join(day_dir, "*.20o")

    nav_file = os.path.join(output_dir, file_name + ".nav")
    obs_file = os.path.join(output_dir, file_name + ".obs")
    pos_file = os.path.join(output_dir, file_name + ".pos")

    tool = os.path.normpath(
        r"C:\SzabolcsKelemen\PhD\GPS\RTKLIB\RTKLIB_bin-rtklib_2.4.3\RTKLIB_bin-rtklib_2.4.3\bin\convbin")
    command1 = "{} -r rinex -f 3 -v 2.10 -od -os -oi -ot -ol -halfc -o {} -n {} {}".format(tool, obs_file, nav_file,
                                                                                           input_files)
    # print(command1)
    os.system(command1)

    command2 = "rnx2rtkp -o {} -p 0 -f 3 -e {} {}".format(pos_file,  obs_file , nav_file)
    # os.system(command2)


def run_converting(root_path, output_dir):
    source_folders = [os.path.join(root_path, f) for f in os.listdir(root_path) if
                      os.path.isdir(os.path.join(root_path, f))]
    for folder in source_folders:
        run_convbin(folder, output_dir)


input_ = r"C:\SzabolcsKelemen\PhD\GPS\NASA\2020_extracted_collected"
output = r"C:\SzabolcsKelemen\PhD\GPS\NASA\obs_files"

run_converting(input_, output)
