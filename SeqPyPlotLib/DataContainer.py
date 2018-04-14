"""
Initialize data container object. Reads in the directory specified in the config and 

Strategy for imputing missing sampels:
1. The data is read in, missing samples are added as NaN columns. 
2. The matrix for normalization is created by:
    1. dropping missing sample columns (imputed later)
    2. dropping zeroed rows (TMM doesn't deal with these...?)
3. The matrix is normalized
4. The zeroed rows are added back in.
5. Imputation is performed on the missing columns by averaging the flanking data points

The config should specify an empty file for missing samples. Parsers will automatically
fill empty files with the 

:param args:
:param Optimize: Set to true when using the Go-term add on script

The initialization of the DataContainer should automatically call the correct parser, 
create the normalization matrix, and then 
"""

import codecs
import csv
import os
import subprocess
import sys

import numpy as np
import pandas as pd

from ConfigParse import config_parser
from normalizer import norm_tmm as TMM
from parsers import CuffNormParser, HtSeqParser, PlotDataParser
from utils import file_is_empty, write_to_csv

PARSERS = {'htseq': HtSeqParser,
           'cuffnorm': CuffNormParser,
           'plot': PlotDataParser}


class DataContainer(object):
    """"
    Class for holding normalized data in a standard format. Calls parsers to collect data in to 
    pandas data frames. Matrix generation and normalization happens here
    """

    def __init__(self, args):
        self.args = args
        self.data_directory, self.paths, self.names = self._unpack_config_()
        self.raw_df = self._parse_input_()


    @property
    def is_parsed(self):
        pass

    def contain_data(self):
        pass

    def _unpack_config_(self):
        """Unpack config contents
        Returns:
            tuple:
                data_dir: A string
                paths: A list of tuples -- [(control_data1, test_data1), etc]
                names: A list of tuples naming the paths
        """
        data_dir, paths, names = config_parser(self.args.config_path)
        paths = [os.path.join(data_dir, x) for x in paths]
        return data_dir, paths, names

    def _parse_input_(self):
        # Instantiante parser
        parser = PARSERS[self.args.datatype]()

        # Execute parser given the data paths and the sample names
        raw_df = parser.parse_data(self.paths, self.names)

        return raw_df












    def normalize_data(self, raw_df, data_pairs):
        pass


    def create_unnormalized_matrix(self, df):
        # TODO This function won't work. We absolutely must treat each pair of samples per time
        # point separately. We should only remove zeroed genes right before we perform normalization,
        # and then reinsert them immediately aftewards. This is because
        """Create the raw matrix for TMM normalization."""

        empty_files = [sample_names[idx] for ix, _file_ in enumerate(data_paths)
                       if file_is_empty(_file_)]
        keep_cols = set(df.columns) - set(empty_files)
        df = df[keep_cols]

        # drop zeroed rows
        nonzero_df = df[df.values.sum(axis=1) != 0]
        zero_df = df[df.values.sum(axis=1) == 0]

        # write out matrix csv
        matrix_path = os.path.join('.', self.args.out, self.args.prefix)
        write_to_csv(matrix_path, suffix='_count_matrix.txt')

        return nonzero_df, zero_df

    def execute_normalization(self, unnormalized_matrix):
        """
        This function should take in an unnormalized matrix of only two time points:
        1. The control time point
        2. The test time point
        """
        return TMM(unnormalized_matrix)

    def recombine_genes(self, nonzero_df, zero_df):
        return pd.concatenate([nonzero_df, zero_df])



        df = self.recombine_genes(normalized_matrix, zeroed_genes)
        if True:
            pass
        # else if there is analyzed plot data
        elif self.args.datatype == 'plot':
            # TODO add check to plotdataparser that checks input
            parser = PlotDataParser()
            self.df = parser.parse_data(datafile)

        print "Data Container initialized Successfully....\n"

        return df



    def __average_flanking__(self, value):
        """
        :param value: a list with missing values to be filled in
        """
        flanked_averaged = []

        for data in enumerate(value):
            if data[1] is not None:
                flanked_averaged.append(data[1])
            else:
                if data[0] == 0 or data[0] == len(value) - 1:  # if its the first or last position -> none
                    flanked_averaged.append(None)
                else:  # if the value is internal to the series list
                    if value[data[0]-1] is not None and value[data[0] + 1] is not None:  # if there is flanking data
                        flanking = [value[data[0] + 1], value[data[0] - 1]]
                        flanked_averaged.append(np.mean(flanking))   # average the data
                    else:
                        flanked_averaged.append(None)
        return flanked_averaged
