import os

def config_parser(config_path):
    """
    Read the file paths in the order given in the config file.
    
    Config should be tab deliminted with the file first, and its name second:
    ~/path/to/file <TAB> name

    Returns:
        :list of lists: [[file]]
    """
    with open(config_path, 'r') as conf: 
        paths = list()
        names = list()
        for line in conf.readlines():
            line = line.split('\t')
            paths.append(line[0].strip())
            names.append(line[1].strip())
    
    return paths, names

