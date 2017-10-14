import os


class PrepareOutputDirectory:
    """
    Must instantiate after args are set
    """

    def __init__(self):
        pass

    @staticmethod
    def make_folder(path):
        outdir = os.path.join(path)

        if not os.path.isdir(outdir):
            os.makedirs(outdir)


def validate_fies_exist(data_paths):
    """Validate all files in config exist"""
    for file_name in data_paths:
        assert os.path.isfile(file), "One or more files is not correct in the config."


def check_for_even_number_of_files(self, data_paths):
    """
    Assert that there are an even number of file paths in the config
    This should run AFTER the check that the the paths exist.
    """
    assert len(data_paths) % 2 == 0, "Even number of input files required. Use empty file if necessary."
    assert len(self.args.time) % 2 == 0


def file_is_empty(self, file_name):
    """Determine if a given file is empty"""
    if int(os.stat(file_name).st_size) == 0:  # if file is empty, record the list position
        return file_name
    else:
        return False