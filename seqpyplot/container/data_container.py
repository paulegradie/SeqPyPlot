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
import math

import numpy as np
import pandas as pd

from pathlib import Path
from scipy import linalg
from scipy.stats import boxcox
from sklearn.linear_model import LinearRegression   

from seqpyplot.parsers import CuffNormParser, HtSeqParser, PlotDataParser
from seqpyplot.container.normalizer import norm_tmm as TMM

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
        self.data_directory = Path(self.config_obj.get('data_directory', 'dir'))
        self.paths = [str(Path(x)) for x in self.config_obj.getlist('data', 'paths')]
        self.names = self.config_obj.getlist('names', 'sample_names')
        self.times = self.config_obj.getlist('names', 'times')
        self.file_pairs = self.config_obj.getlist('names', 'file_pairs')
        self.num_file_pairs = self.config_obj.getint('misc', 'num_file_pairs')
    
    def split(self, normalized_df):
        self.normalized_df = normalized_df
        self.split_normalized_dfs = [normalized_df[[control_col, treated_col]] for (control_col, treated_col) in self.file_pairs]
        return self.split_normalized_dfs

    def parse_input(self):
        # Instantiante parser
        parser = PARSERS[self.config_obj.get('data', 'data_type')]()

        # Execute parser given the data paths and the sample names
        data_df, ercc_df = parser.parse_data(self.paths, self.names)
        
        self.data_df = data_df
        self.ercc_df = ercc_df

        self.complete_gene_list = data_df.index.tolist()
        
        return data_df, ercc_df

    def make_col_pairs(self, df):
        """
        Input: df where a1, a2, b1, b2
        Return: col names as [(a1, b1), (a2, b2)]
        """
        control_cols = df.columns[:self.num_file_pairs]
        treated_cols = df.columns[self.num_file_pairs:]

        return zip(control_cols, treated_cols)

    def normalize_file_pairs(self, data_df):
        """
        TMM shouldn't be used to normalize across developmental time
        series data. This func splits off pairs at the same step
        and normalized them togther
        """
        normalized_pairs = list()
        for control, treat in self.make_col_pairs(data_df):

            sub_df = data_df[[control, treat]]

            normalized_sub_df = self.execute_normalization(sub_df)
            normalized_pairs.append(normalized_sub_df)

        remerged_df = self.merge_dfs(normalized_pairs)
        reordered_df = self.reorder_cols(remerged_df)
        self.normalized_data = reordered_df
        return reordered_df

    def execute_normalization(self, unnormalized_matrix):
        return TMM(unnormalized_matrix)

    def merge_dfs(self, dfs):
        return reduce(lambda x, y: pd.merge(x, y, left_index=True, right_index=True, how='outer'), dfs)

    def reorder_cols(self, df):
        controls = [x[1] for x in enumerate(df.columns) if x[0] % 2 == 0]
        treated = [x[1] for x in enumerate(df.columns) if x[0] % 2 != 0]
        return df[controls + treated]

    def rotate(self, point_tuple, theta):
        """
        Rotate a point_tuple counterclockwise by a given angle around a given origin.
        The angle should be given in radians.
        """

        rotation_matrix = np.array([[np.cos(theta), np.sin(theta)],
                                    [-np.sin(theta), np.cos(theta)]])

        # # angle = np.radians(angle)
        x, y = point_tuple
        new_x = (x * np.cos(theta)) + (x * np.sin(theta))
        new_y = (y * -np.sin(theta)) + (y * np.cos(theta))

        return new_x, new_y

    def calc_theta(self, coef1, coef2):
       return np.abs(
           np.arctan(np.abs(coef1 - coef2) / (1. + (coef1 * coef2)))
       )

    def correct_heteroskedacity(self, dfs):
        """
        Employ Box-cox transformation (power tranfsormation) followed by rotation
        to correct heteroskedacitiy in the data.
        Use sklearn.linearregression to compute slopes. Use scipy to comput boxcox
        
        """
        result = list()
        for df in dfs:
            df = df.copy()

            cols = df.columns.tolist()

            df.loc[:, 'mean'] = df.mean(axis=1)
            control = df[cols[0]].values.reshape(-1, 1)
            treated = df[cols[1]].values.reshape(-1, 1)
            mean = df['mean'].values.reshape(-1, 1)

            control_regression = LinearRegression(fit_intercept=True, n_jobs=2)
            control_regression.fit(control, mean)
            treated_regression = LinearRegression(fit_intercept=True, n_jobs=2)
            treated_regression.fit(treated, mean)

            # First, ajust all points to the origin (remove bias, zero out bias, etc)
            intercepts = [control_regression.intercept_, treated_regression.intercept_]
            max_idx = np.argmax(intercepts)
            min_idx = np.argmin(intercepts)
            df.loc[:, col[max_idx]] = df[col[max_idx]] - intercepts[max_idx]
            df.loc[:, col[min_idx]] = df[col[min_idx]] - intercepts[min_idx]



            # Transform
            theata = self.calc_theta(coef1, coef2)
            # First, rotate the lower line (lower y intercept) around its y intercept

            df[cols[min_idx]] = self.rotate(point_tuple=(df['mean'], df[cols[min_idx]]), 
                                            angle=theta)
            df[cols[min_idx]] = df[cols[min_idx]].apply(lambda x: x + vertical_shift)


            result.append(df.drop('mean', axis=1))
        
        return result

    def apply_box_cox(self, col_data):
        return boxcox(col_data)


    def remove_variance(self, df, num_components=1):
        "Use Single Value Decomposition to remove first component of variation"
        u, s, v = linalg.svd(df, full_matrices=False)
        s = np.diag(s)
        s[:num_components, :] = 0.0
        
        reconstructed = np.dot(u, np.dot(s, v))
        return reconstructed

    # #TODO reimplement support for missing data (data imputation)
    # def _average_flanking_(self, value):
    #     """
    #     :param value: a list with missing values to be filled in
    #     """
    #     flanked_averaged = []

    #     for data in enumerate(value):
    #         if data[1] is not None:
    #             flanked_averaged.append(data[1])
    #         else:
    #             if data[0] == 0 or data[0] == len(value) - 1:  # if its the first or last position -> none
    #                 flanked_averaged.append(None)
    #             else:  # if the value is internal to the series list
    #                 if value[data[0]-1] is not None and value[data[0] + 1] is not None:  # if there is flanking data
    #                     flanking = [value[data[0] + 1], value[data[0] - 1]]
    #                     flanked_averaged.append(np.mean(flanking))   # average the data
    #                 else:
    #                     flanked_averaged.append(None)
    #     return flanked_averaged
