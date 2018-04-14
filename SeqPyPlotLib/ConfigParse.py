from ConfigParser import ConfigParser
from ast import literal_eval


def config_parser(config_path):
    """
    Read the provided config.ini file and parse the file paths to the data.
    """

    config = ConfigParser()
    config.read(config_path)

    control_data = literal_eval(config.get('data', 'controls'))
    treated_data = literal_eval(config.get('data', 'treated'))

    control_names = literal_eval(config.get('names', 'controls'))
    treated_names = literal_eval(config.get('names', 'treated'))

    dir_path = config.get('data_directory', 'dir')

    data = control_data + treated_data
    names = control_names + treated_names

    return (dir_path, data, names)
