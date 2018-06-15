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
        self.num_components = self.config_obj.getint('params', 'num_components')

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

    def calc_theta(self, coef1, coef2):
       return np.abs(
           np.arctan(np.abs(coef1 - coef2) / (1. + (coef1 * coef2)))
       )

    def compute_rot_mat(self, rad):
        " Compute a rotation matrix using rad for a given regression coefficient "
        rotation_matrix = np.array([[np.cos(rad), -np.sin(rad)],
                                    [np.sin(rad), np.cos(rad)]])
        return rotation_matrix

    def correct_via_rotation(self, dfs, min_itvl=100, max_intvl=500):
        """
        Correct heteroscedacitiy using linear regression.

        """
        result = list()
        for frame in dfs:

            dataframe = frame.copy()
            # dataframe is the sample pair - add mean
            dataframe.loc[:, 'mean'] = dataframe.mean(axis=1)

            # make filtered copy for linear regression
            df = dataframe[(dataframe['mean'] > min_itvl) & (dataframe['mean'] < max_intvl)].copy()

            # make ordered list of col names
            cols = df.columns.tolist()

            # prep the data for sklearn linear regression
            control = df[cols[0]].values.reshape(-1, 1)
            treated = df[cols[1]].values.reshape(-1, 1)
            mean = df['mean'].values.reshape(-1, 1)

            # perform linear regression
            control_regression = LinearRegression(fit_intercept=True, n_jobs=2)
            control_regression.fit(control, mean)
            treated_regression = LinearRegression(fit_intercept=True, n_jobs=2)
            treated_regression.fit(treated, mean)

            # remove bias
            dataframe.loc[:, cols[0]] = dataframe[cols[0]].sub(control_regression.intercept_[0])
            dataframe.loc[:, cols[1]] = dataframe[cols[1]].sub(treated_regression.intercept_[0])

            # Compute angle between two data groups
            coefficients = list(map(float, [np.squeeze(control_regression.coef_[0]), np.squeeze(treated_regression.coef_[0])]))
            theta = self.calc_theta(*coefficients)

            # compute rotation matrix
            rotation_matrix = self.compute_rot_mat(theta)

            # keep track of lower and upper col
            lower_col = cols[np.argmin(coefficients)]
            upper_col = cols[np.argmax(coefficients)]

            old_points = dataframe[[lower_col, 'mean']].values
            new_points = np.dot(old_points, rotation_matrix)

            corrected_df = pd.DataFrame(new_points[:, 0], columns=[lower_col], index=dataframe.index)
            corrected_df.loc[:, upper_col] = dataframe[upper_col]
            corrected_df = corrected_df[[cols[0], cols[1]]]
            corrected_df[corrected_df < 0] = 0

            result.append(corrected_df)

        return result

    def _compute_svd_(self, df, num_components):
        "Use Single Value Decomposition to remove first component of variation"
        u, s, v = linalg.svd(df, full_matrices=False)
        s = np.diag(s)
        s[:num_components, :] = 0.0
        
        reconstructed = np.dot(u, np.dot(s, v))
        return reconstructed

    def remove_variance(self, dfs, num_components=1):
        result = list()
        for dataframe in dfs:
            df = dataframe.copy()
            
            cols = df.columns
            index = df.index
            
            thinned = self._compute_svd_(df, self.num_components)

            result.append(pd.DataFrame(thinned, columns=cols, index=index))
        return result

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
