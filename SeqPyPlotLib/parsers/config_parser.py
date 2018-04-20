from ConfigParser import ConfigParser
from ast import literal_eval
import os

def config_parser(config_path):
    """
    Read the provided config.ini file and parse the file paths to the data.
    """

    config = ConfigParser()
    config.read(config_path)

    data, names, pairs, num_file_pairs, times, conds = configure_input_data(config)
    
    config.set('data', 'paths', value=data)
    
    config.set('names', 'sample_names', value=names)
    config.set('names', 'file_pairs', value=pairs)  
    config.set('names', 'times', value=times)
    config.st('names', 'conditions', values=conds)
    config.add_section('misc')
    config.set('misc', 'num_file_pairs', value=str(num_file_pairs))

    return config


def configure_input_data(config_object):

    control_data = literal_eval(config_object.get('data', 'controls'))
    treated_data = literal_eval(config_object.get('data', 'treated'))

    control_names = literal_eval(config_object.get('names', 'controls'))
    treated_names = literal_eval(config_object.get('names', 'treated'))

    times = literal_eval(config_object.get('names', 'times'))
    conds = literal_eval(config_obj.get('names', 'conditions'))

    dir_path = config_object.get('data_directory', 'dir')

    data = [os.path.join(dir_path, x) for x in control_data] + [os.path.join(dir_path, x) for x in treated_data]
    names = control_names + treated_names

    pairs = zip(control_names, treated_names)

    assert len(control_data) == len(treated_data), "Use empty file to fill gaps in data."
    return data, names, pairs, len(control_data), times, conds
