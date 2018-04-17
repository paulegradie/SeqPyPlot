from ConfigParser import ConfigParser
from ast import literal_eval


def config_parser(config_path):
    """
    Read the provided config.ini file and parse the file paths to the data.
    """

    config = ConfigParser()
    config.read(config_path)

    dir_path, data, names, num_file_pairs = configure_input_data(config)
    config.set('data', 'paths', value=data)
    config.set('names', 'sample_names', value=names)
    config.set('misc', 'num_file_pairs', value=num_file_pairs)

    return config


def configure_input_data(config_object):

    control_data = literal_eval(config_object.get('data', 'controls'))
    treated_data = literal_eval(config_object.get('data', 'treated'))

    control_names = literal_eval(config_object.get('names', 'controls'))
    treated_names = literal_eval(config_object.get('names', 'treated'))

    dir_path = config_object.get('data_directory', 'dir')

    data = control_data + treated_data
    names = control_names + treated_names

    assert len(control_data) == len(treated_data), "Use empty file to fill gaps in data."
    return dir_path, data, names, len(control_data)
