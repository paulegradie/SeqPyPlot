import os
from ast import literal_eval
from configparser import ConfigParser
from pathlib import Path


# Method to add to configparser
def getlist(self, section, value):
    raw_list = self.get(section, value)
    return literal_eval(raw_list)


def config_parser(config_path):
    """
    Read the provided config.ini file and parse the file paths to the data.
    """
    ConfigParser.getlist = getlist

    config_obj = ConfigParser()
    config_obj.read(config_path)

    paths, names, pairs, num_file_pairs = _configure_input_data_(config_obj)

    config_obj.set('data', 'paths', value=str(paths))

    config_obj.set('names', 'sample_names', value=str(names))
    config_obj.set('names', 'file_pairs', value=str(pairs))

    config_obj.add_section('misc')
    config_obj.set('misc', 'num_file_pairs', value=str(num_file_pairs))

    return config_obj


def _configure_input_data_(config_object):

    dir_path = Path(config_object.get('data_directory', 'dir'))
    control_data = config_object.getlist('data', 'controls')
    treated_data = config_object.getlist('data', 'treated')
    paths = [str(dir_path / x) for x in control_data] + [str(dir_path / x) for x in treated_data]

    control_names = config_object.getlist('names', 'controls')
    treated_names = config_object.getlist('names', 'treated')
    names = control_names + treated_names
    control_names = config_object.getlist('names', 'controls')
    treated_names = config_object.getlist('names', 'treated')

    pairs = [(x, y) for x, y in zip(control_names, treated_names)]

    assert len(control_data) == len(treated_data), "Use empty file to fill gaps in data."
    return paths, names, pairs, len(control_data)
