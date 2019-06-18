import os
from datetime import datetime
from shutil import rmtree

TIME_STAMP = str(datetime.now()).split()[0]


def make_default_output_dir(dirpath=None, overwrite=False):
    """ pass it any dirpath to make a directory """
    if not dirpath:
        dirpath = "_".join(['default_out', str(datetime.now())])

    if os.path.exists(dirpath):
        if overwrite:
            rmtree(dirpath)
            os.mkdir(dirpath)
        else:
            raise IOError("Output dir exists.")
    else:
        os.mkdir(dirpath)

    return dirpath


# def validate_fies_exist(data_paths):
#     """Validate all files in config exist"""
#     for file_name in data_paths:
#         assert os.path.isfile(file), "One or more files is not correct in the config."


# def check_for_even_number_of_files(self, data_paths):
#     """
#     Assert that there are an even number of file paths in the config
#     This should run AFTER the check that the the paths exist.
#     """
#     assert len(data_paths) % 2 == 0, "Even number of input files required. Use empty file if necessary."
#     assert len(self.args.time) % 2 == 0


# def file_is_empty(file_name):
#     """Determine if a given file is empty"""
#     return file_name if int(os.stat(file_name).st_size) == 0 else False  # if file is empty, record the list position

