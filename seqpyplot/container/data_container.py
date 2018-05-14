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
fill empty files with the average of the flanking time points.

:param args:
:param Optimize: Set to true when using the Go-term add on script

The initialization of the DataContainer should automatically call the correct parser, 
create the normalization matrix, and then 
"""
import os

import numpy as np
import pandas as pd

from ..parsers import CuffNormParser, HtSeqParser, PlotDataParser
from ..utils.utils import write_to_csv
from .normalizer import norm_tmm as TMM

try:
    from functools import reduce # python 3 compatibility
except ImportError:
    pass


PARSERS = {'htseq': HtSeqParser,
           'cuffnorm': CuffNormParser,
           'plot': PlotDataParser}


class DataContainer(object):
    """"
    Class for holding normalized data in a standard format. Calls parsers to collect data in to 
    pandas data frames. Matrix generation and normalization happens here
    """

    def __init__(self, config_obj):

        self.config_obj = config_obj

        self.data_directory = self.config_obj.get('data_directory', 'dir')
        self.paths = self.config_obj.getlist('data', 'paths')
        self.names = self.config_obj.getlist('names', 'sample_names')
        self.times = self.config_obj.getlist('names', 'times')
        self.file_pairs = self.config_obj.get('names', 'file_pairs')
        self.num_file_pairs = self.config_obj.getint('misc', 'num_file_pairs')
    
    def split(self, normalized_df):
        return [normalized_df[[control_col, treated_col]] for (control_col, treated_col) in self.file_pairs]

    def parse_input(self):
        # Instantiante parser
        parser = PARSERS[self.config_obj.get('data', 'data_type')]()

        # Execute parser given the data paths and the sample names
        raw_df, ercc_data = parser.parse_data(self.paths, self.names)

        # self.is_parsed = True
        return raw_df, ercc_data

    def make_col_pairs(self, df):
        """
        Input: df where a1, a2, b1, b2
        Return: col names as [(a1, b1), (a2, b2)]
        """
        control_cols = df.columns[:self.num_file_pairs]
        treated_cols = df.columns[self.num_file_pairs:]

        return zip(control_cols, treated_cols)

    def normalize_file_pairs(self, raw_df):
        """
        TMM shouldn't be used to normalize across developmental time
        series data. This func splits off pairs at the same step
        and normalized them togther
        """
        normalized_pairs = list()
        for control, treat in self.make_col_pairs(raw_df):

            sub_df = raw_df[[control, treat]]

            normalized_sub_df = self.execute_normalization(sub_df)
            normalized_pairs.append(normalized_sub_df)

        remerged_df = self.merge_dfs(normalized_pairs)
        return remerged_df

    def execute_normalization(self, unnormalized_matrix):
        return TMM(unnormalized_matrix)

    def merge_dfs(self, dfs):
        return reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True, how='outer'), dfs)

    def reorder_cols(self, df):
        controls = [x[1] for x in enumerate(df.columns) if x[0] % 2 == 0]
        treated = [x[1] for x in enumerate(df.columns) if x[0] % 2 != 0]
        return df[controls + treated]

    #TODO implement support for missing data (data imputation)
    def _average_flanking_(self, value):
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
