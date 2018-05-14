import pandas as pd
import os
import numpy as np

try:
    from functools import reduce
except ImportError:
    pass

class PairedSampleFilter(object):
    """
    This class handles filtering paired normalized samples.
    It inherits from the DataContainer object and only needs a valid config.ini

    Operations work on file pairs. Paired columns are extracted from the
    normalized df and filtering is applied. Operations that act on the
    entire extracted df should be executed before ops that split the extracted df.
    i.e. fold change and diff ops act on the entire data frame, where as hi and low
    ops split the extracted df.
    
    
    Arguments:
        DataContainer {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """
    

    def __init__(self, config_obj, log2fold=None, low=None, hi=None, diff=None):
        
        self.config_obj = config_obj

        self.log2fold = log2fold or self.config_obj.getfloat('params', 'log2fold')
        self.low = low or self.config_obj.getint('params', 'low')
        self.hi = hi or self.config_obj.getint('params', 'hi')
        self.diff = diff or self.config_obj.getlist('params', 'diff')
        assert isinstance(self.diff, list), 'diff must be iterable - 2 values'

        self.times = self.config_obj.getlist('names', 'times')
        self.file_pairs = self.config_obj.get('names', 'file_pairs')

    def main_filter_process(self, input_df_list):
        
        fold_change_dfs = self.apply_fold_change(input_df_list)
        diff_dfs = self.apply_diff(input_df_list)
        merged_dfs = list()
        for fcd, dd in zip(fold_change_dfs, diff_dfs):
            merged_dfs.append(fcd.reset_index().merge(dd,how='inner').set_index('gene'))

        result = self.apply_low(merged_dfs)
        result = self.apply_hi(result)

        # Set Attributes
        self.set_attributes(result)
        return result  # A list of filtered dataframes

    def set_attributes(self, result):
        """Set attributes
        Arguments:
            result {list of pandas dfs} -- filter results
        """
        self.de_count_by_stage = {time: len(df) for time, df in zip(self.times, result)}
        self.de_count_by_gene = self.count_by_gene(result)
        self.de_gene_list_by_stage = {time: df.index for time, df in zip(times, result)}
        self.complete_de_gene_list = set(sorted(reduce(lambda x, y: pd.concat([x, y], axis=0), result).index.tolist()))
        
        return None

    def count_by_gene(self, result):
        "Result: list of data frames"
        gene_count = dict()
        for df in result:

            for gene in df.index:
                if gene in gene_count.keys():
                    gene_count[gene] += 1  
                else:
                    gene_count[gene] = 1

        return gene_count

    def apply_fold_change(self, input_df_list):

        fold_change_dfs = list()
        filtered_results = list()
        for (control_col, treated_col), df in zip(self.file_pairs, input_df_list):

            result = np.log2(df[treated_col].div(df[control_col])).abs()

            fold_change_dfs.append(result)
            filtered_results.append(df[result > self.log2fold])            

        self.fold_change_filtered_dfs = pd.concat(fold_change_dfs, axis=1)
        return filtered_results

    def apply_diff(self, input_df_list):
        """
        
        Arguments:
            input_df_list {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """

        diff_dfs = list()
        filtered_results = list()
        for (control_col, treated_col), df in zip(self.file_pairs, input_df_list):

            result = df[control_col].sub(df[treated_col]).abs()
            diff_dfs.append(result)

            df = df[(result >= self.diff[0]) & (result <= self.diff[1])]
            filtered_results.append(df)            

        self.diff_filtered_dfs = pd.concat(diff_dfs, axis=1)
        return filtered_results     

    def apply_low(self, input_df_list):
        """
        This function is applied, like apply_hi, differently than apply_diff or apply_fold_change.
        It takes the output from diff and fold_change functions

        Returns:
            [type] -- [description]
        """
        filtered_results = list()
        state_change = list()
        for (control_col, treated_col), df in zip(self.file_pairs, input_df_list):
            df.columns = [control_col, treated_col]
            passing = df[(df[control_col] >= self.low) | (df[treated_col] >= self.low)]
            filtered_results.append(passing)
            
            deactivated = df[(df[control_col] > self.low) & (df[treated_col] < self.low)]
            activated = df[(df[control_col] < self.low) & (df[treated_col] > self.low)]    
            undetected = df[(df[control_col] < self.low) & (df[treated_col] < self.low)]
            state_change.append((deactivated[deactivated[control_col].sub(deactivated[treated_col]).abs() > self.diff[0]],
                                 activated[activated[control_col].sub(activated[treated_col]).abs() > self.diff[0]],
                                 undetected))

        self.state_change = state_change

        return filtered_results

    def apply_hi(self, input_df_list):

        filtered_results = list()
        saturated_dfs = list()
        for (control_col, treated_col), df in zip(self.file_pairs, input_df_list):

            passing = df[
                (df[control_col] < self.hi) | 
                (df[treated_col] < self.hi)
                ]
            filtered_results.append(passing)

            saturated = df[(df[control_col] > self.hi) & (df[treated_col] > self.hi)]
            saturated_dfs.append(saturated)

        self.saturated_dfs = saturated_dfs

        return filtered_results

    def calculate_dispersion_estimate(self, normalized_df):

        dispersion_df = list()
        for cont, treat in self.file_pairs:
            dispersion_df.append((normalized_df[cont] - normalized_df[treat]).abs())

        disp_df = pd.concat(dispersion_df, axis=1)
        row_means = disp_df.mean(axis=1)
        dispersion_estimate = row_means.mean()

        return dispersion_estimate

    # TODO implement this idea
    def apply_min_dispersion_estimate(self, input_df_list):
        pass
    #     filtered_results = list()
    #     for (control_col, treated_col), df in zip(self.file_pairs, input_df_list):
    #         result = df[(df[control_col])]
    #         filtered_results.append(s)



if __name__ == "__main__":
    pass
    # self.filtered_genes = self.main_filter_process(self.split_normalized_dfs)
