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

        self.is_parsed = False

        self.args = args
        self.data_directory, self.paths, self.names, self.num_file_pairs = self._unpack_config_()
        self.raw_df, self.ercc_data = self._parse_input_()
        self.normalized_df = self.normalize_file_pairs()

        if impute_by_nieghbors:
            self.normalized_imputed_df = self._average_flanking_()

    @property
    def is_parsed(self):
        pass

    def _unpack_config_(self):
        """Unpack config contents
        Returns:
            tuple:
                data_dir: A string
                paths: A list of tuples -- [(control_data1, test_data1), etc]
                names: A list of tuples naming the paths
        """
        data_dir, paths, names, num_file_pairs = config_parser(self.args.config_path)
        paths = [os.path.join(data_dir, x) for x in paths]
        return data_dir, paths, names, num_file_pairs

    def _parse_input_(self):
        # Instantiante parser
        parser = PARSERS[self.args.datatype]()

        # Execute parser given the data paths and the sample names
        raw_df, ercc_data = parser.parse_data(self.paths, self.names)

        self.is_parsed = True
        return raw_df, ercc_data

    def normalize_data(self, raw_df, data_pairs):
        pass

    def make_col_pairs(self):
        control_cols = self.raw_df.columns[:self.num_file_pairs]
        treated_cols = self.raw_df.columns[self.num_file_pairs:]

        return zip(control_cols, treated_cols)

    def normalize_file_pairs(self):

        normalized_pairs = list()
        for control, treat in self.make_col_pairs():

            sub_df = self.raw_df[[control, treat]]
            nonzero_df, zero_df = self.extract_usable_data(sub_df)

            normalized_nonzero = self.execute_normalization(nonzero_df)
            normalized_pairs.append(pd.concat([normalized_nonzero, zero_df]))

        remerged_df = self.merge_dfs(normalized_pairs)

        return remerged_df

    def extract_usable_data(self, df):

        nonzero_df = df[df.values.sum(axis=1) != 0]
        zero_df = df[df.values.sum(axis=1) == 0]

        return nonzero_df, zero_df

    def execute_normalization(self, unnormalized_matrix):
        """
        This function should take in an unnormalized matrix of only two time points:
        1. The control time point
        2. The test time point
        """
        return TMM(unnormalized_matrix, unnormalized_matrix.columns[0])

    def merge_dfs(self, dfs):
        return reduce(lambda x, y: pd.merge(x, y, on='gene', how='outer'), dfs)

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
