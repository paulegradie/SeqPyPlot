
# coding=utf-8
import csv
import os
import codecs
import subprocess
import sys
import numpy as np
import pandas as pd

from opertor import reduce
from normalizer import norm_tmm as TMM

# from utils import validate_fies_exist, check_for_even_number_of_files

# TODO: parallelize the file read
class Htseqparser(object):
    """ 
    Read a directory of expression counts in ht-seq format. Each sample should be an
    individual file in the directory. File names and sample order are specified in
    the config file (order is determined by order IN the config.)

    This class is intended to return the raw dataframe of samples with missinlg sample columns
    as NaN.
    """

    def __init__(self, number_of_conditions=2):

        self.number_of_conditions = number_of_conditions


    def load_htseq_dataframe_single(self, data_paths, sample_names, empty_files):
        """
        Read the input files from the config file and load in to a pandas dataframe.

        params
            data_paths: list of file paths specified in the config. Returned from config parse
            sample_names: list of sample names specified in the config returned from config parse
        
        """

        # Load the data
        dfs = [pd.read_csv(file, sep='\t', names=['gene', sample_names[idx]]) for idx, file in enumerate(data_paths)]
        df = reduce(lambda x: pd.merge(x, on='gene', how='outer'), dfs)
        df = df[df['gene'].str.startswith('__') == False]
        df.set_index('gene')
        df.fillna(value=0, inplace=True)
        
        #write out raw data for access later
        df.to_csv('{}_raw_count_data.csv'.format(self.args.out))

        return df
