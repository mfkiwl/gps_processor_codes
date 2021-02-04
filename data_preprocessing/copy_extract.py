import os
import zipfile
import datetime
from shutil import copyfile


def copy_from_dir_to_dir(src, dest):
    files = os.listdir(src)  # [os.path.abspath(name) for name in os.listdir(src) if os.path.isfile(name)]
    for file in files:
        # print(file)
        abs_file = os.path.join(src, file)
        dest_file = os.path.join(dest, file)
        copyfile(abs_file, dest_file)


def copy_data_of_day(src_day, dest_day):
    rinex_src = os.path.join(src_day, "20o")
    # hod_dirs = [os.path.abspath(name) for name in os.listdir(rinex_src) if os.path.isdir(name)]
    hod_dirs = [f.name for f in os.scandir(rinex_src) if f.is_dir()]
    for hour in hod_dirs:
        src_hour = os.path.join(rinex_src, hour)
        copy_from_dir_to_dir(src_hour, dest_day)


def copy_data_of_day_by_condition(src_day, dest_day, receiver_name, nr_days):
    files = os.listdir(src_day)
    for file in files:
        if receiver_name in file[::nr_days]:
            abs_file = os.path.join(src_day, file)
            dest_file = os.path.join(dest_day, file)
            copyfile(abs_file, dest_file)


def download_ftp_data_of_day_by_condition(src_day, dest_day, receiver_name, nr_days):
    files = os.listdir(src_day)
    for file in files:
        if receiver_name in file[::nr_days]:
            abs_file = os.path.join(src_day, file)
            dest_file = os.path.join(dest_day, file)
            copyfile(abs_file, dest_file)


def download_from_ftp(dir_source, dir_dest, pattern):
    import sys
    import os
    import ftplib
    import ftputil
    import fnmatch
    import time
    from time import mktime
    import datetime
    import os.path, time
    from ftplib import FTP

    print("logging into GSP FTP")  # print

    with ftputil.FTPHost(dir_source, 'anonymous', '') as host:  # ftp host info
        recursive = host.walk(dir_source, topdown=True, onerror=None)  # recursive search
        for root, dirs, files in recursive:
            for name in files:
                seeked_list = fnmatch.filter(files, pattern)
                for name in seeked_list:
                    fpath = host.path.join(root, name)
                    print(fpath)
                    if host.path.isfile(fpath):
                        host.download_if_newer(fpath, os.path.join(dir_dest, name), 'b')

    host.close()


def day_folder_name(year, dest_year, day_of_year, prefix):
    day_of_year = '{}'.format(day_of_year[1:] if day_of_year.startswith('0') else day_of_year)
    day_of_year = int('{}'.format(day_of_year[1:] if day_of_year.startswith('0') else day_of_year))
    date = datetime.datetime(year, 1, 1) + datetime.timedelta(day_of_year - 1)
    month = '0{}'.format(date.month) if len(str(date.month)) == 1 else date.month
    day = '0{}'.format(date.day) if len(str(date.day)) == 1 else date.day
    day_folder_name = prefix + str(date.year) + str(month) + str(day) + str(24)
    day_folder = os.path.join(dest_year, day_folder_name)
    if not os.path.isdir(day_folder):
        os.makedirs(day_folder)
    return day_folder


def extract_content(path):
    copy_from_dir_to_dir(os.path.normpath("C:\SzabolcsKelemen\PhD\GPS\Data_Processing\gps_processor_codes\shel_for_7z"),
                         path)
    files = os.listdir(path)
    print(len(files))
    os.chdir(path)
    command = ".\extract_with_7z.sh"
    os.system(command)

    # files = [f.path for f in os.scandir(path) if (f.is_file() and os.path.splitext(f.path)[1] != "20o")]
    print(len(files))

    # for file in files:
    #     file_name = os.path.abspath(file)
    #     zip_ref = zipfile.ZipFile(file_name)
    #     zip_ref.extractall(path)
    #     zip_ref.close()

    for file in files:
        os.remove(file)


def create_daydir_names():
    days = []
    for i in range(1, 300, 1):
        if len(str(i)) == 3:
            day_str = str(i)
        elif len(str(i)) == 2:
            day_str = "0{}".format(str(i))
        elif len(str(i)) == 1:
            day_str = "00{}".format(str(i))
        days.append(day_str)
    return days


def process_all_year(src_year, dest_year, prefix):
    day_dirs = [f.name for f in os.scandir(src_year) if f.is_dir()]
    # day_dirs = create_daydir_names()
    for day in day_dirs:
        try:
            day_src_path = os.path.join(src_year, day)
            # print(day_src_path)
            dest_day = day_folder_name(2020, dest_year, day, prefix)
            print(dest_day)
            copy_data_of_day(day_src_path, dest_day)
            # copy_data_of_day_by_condition(day_src_path, dest_day, "ARTA", 100)
            # download_from_ftp(day_src_path, dest_day, "ARTA")
            extract_content(dest_day)
        except:
            pass


src_year = r"C:\SzabolcsKelemen\PhD\GPS\NASA\2020"
dest_year = r"C:\SzabolcsKelemen\PhD\GPS\NASA\2020_extracted_collected"
prefix = "NASA"

# src_year = r"ftp://ftp.geonet.org.nz/rtgps/rinex1Hz/PositioNZ/2020"
# dest_year = r"/Users/kelemensz/Qsync/GPS/reciever_data/NZLD"
# prefix = "NZLD"


process_all_year(src_year, dest_year, prefix)
