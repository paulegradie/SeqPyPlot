from ConfigParser import ConfigParser
from ast import literal_eval
import os


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

    data, names, pairs, num_file_pairs = configure_input_data(config_obj)
    
    config_obj.set('data', 'paths', value=data)

    config_obj.set('names', 'sample_names', value=names)
    config_obj.set('names', 'file_pairs', value=pairs)  

    config_obj.add_section('misc')
    config_obj.set('misc', 'num_file_pairs', value=str(num_file_pairs))

    return config_obj


def configure_input_data(config_object):

    dir_path = config_object.get('data_directory', 'dir')
    control_data = config_object.getlist('data', 'controls')
    treated_data = config_object.getlist('data', 'treated')
    data = [os.path.join(dir_path, x) for x in control_data] + [os.path.join(dir_path, x) for x in treated_data]
    
    control_names = config_object.getlist('names', 'controls')
    treated_names = config_object.getlist('names', 'treated')
    names = control_names + treated_names
    control_names = config_object.getlist('names', 'controls')
    treated_names = config_object.getlist('names', 'treated')

    pairs = zip(control_names, treated_names)

    assert len(control_data) == len(treated_data), "Use empty file to fill gaps in data."
    return data, names, pairs, len(control_data)
