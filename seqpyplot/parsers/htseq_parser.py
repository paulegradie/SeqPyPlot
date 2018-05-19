"""
Read a directory of expression counts in ht-seq format. Each sample
should be an individual file in the directory. File names and
sample order are specified in the config file (order is determined
by order IN the config.)

This class is intended to return the raw dataframe of samples with
missing sample columns as NaN.
"""

import pandas as pd
from pathos.multiprocessing import ProcessPool
try:
    from functools import reduce # for py3 compatibility
except ImportError:
    pass

class HtSeqParser(object):

    def __init__(self, nodes=2):

        self.nodes = nodes

    def parse_data(self, data_paths, sample_names):
        """
        Read the input files from the config file and load in to a
        pandas dataframe.

        params
            data_paths: list of file paths specified in the config. Returned
            from config parse sample_names: list of sample names specified in
            the config returned from config parse

        """
        output = self.load_data(data_paths, sample_names)
        data, ercc_df = (self.merge_dfs(output)
                         .pipe(self.df_cleanup)
                         .pipe(self.split_on_ercc))

        return data, ercc_df

    def load_data(self, data_paths, sample_names):
        " Multiprocess load of files in to a list of dfs "
        pool = ProcessPool(nodes=self.nodes)
        dfs = pool.map(self.load_func, zip(data_paths, sample_names))
        return dfs

    @staticmethod
    def load_func(data_tuple):
        path, sample_name = data_tuple
        return pd.read_csv(path, sep='\t', names=['gene', sample_name])

    def merge_dfs(self, dfs):
        return reduce(lambda x, y: pd.merge(x, y, on='gene', how='outer'), dfs)

    def df_cleanup(self, df_old):
        " Clean away unwanted columns, reset index, and fillna "
        df = df_old.copy()
        df = df[df['gene'].str.startswith('__') == False]
        df.set_index('gene', inplace=True)
        df.fillna(value='Nan', inplace=True)
        return df

    def split_on_ercc(self, df):
        " Extract the ERCC data "
        ercc_cols = df.index.str.startswith('ERCC-')

        ercc_df = df[ercc_cols]
        data = df[~ercc_cols]

        return data, ercc_df
