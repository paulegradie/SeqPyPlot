
from ..container.DataContainer import DataContainer
import pandas as pd
import os

class PairedSampleFilter(DataContainer):
    """
    This class handles filtering paired normalized samples.
    It inherits from the DataContainer object and only needs a valid config.ini

    Operations word on file pairs. Paired columns are extracted from the
    normalized df and filtering is applied. Operations that act on the
    entire extracted df should be executed before ops that split the extracted df.
    i.e. fold change and diff ops act on the entire data frame, where as hi and low
    ops split the extracted df.
    
    
    Arguments:
        DataContainer {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """

    def __init__(self, config):
        """[summary]
        
        Arguments:
            config {[type]} -- [description]
        
        Other properties available after __init__:
        self.saturated_dfs   -- a list of dfs
        self.filtered_genes  -- a list of dfs 
        self.state_change    -- a list of dfs

        """
        super(PairedSampleFilter, self).__init__(config=config)    
        self.log2fold = self.config.getfloat('params', 'log2fold')
        self.low = self.config.getint('params', 'low')
        self.hi = self.config.getint('params', 'hi')
        self.diff = self.config.getint('params', 'diff')

        self.split_dfs = self.split()
        self.main_filter_process(self.split_dfs)

        # Collect other information concerning the data
        times = self.config.get('names', 'times')
        self.de_count_by_stage = {time: len(df) for time, df in zip(times, self.filtered_genes)}

        self.de_count_by_gene = self.count_by_gene()
        self.de_gene_list_by_stage = {time: df.index for time, df in zip(times, self.filtered_genes)}
        self.complete_de_gene_list = set(sorted(reduce(lambda x, y: pd.concat([x, y], axis=1), self.filtered_genes).index.tolist()))

        
    def count_by_gene(self):
        
        gene_count = dict()
        for df in self.filtered_genes:

            for gene in df.index:
                if gene in gene_count.keys():
                    gene_count[gene] += 1  
                else:
                    gene_count[gene] = 0

        return gene_count

    def main_filter_process(self, input_df_list):
        
        fold_change_dfs = self.apply_fold_change(input_df_list)
        diff_dfs = self.apply_diff(input_df_list)

        merged_dfs = list()
        for fcd, dd in zip(fold_change_dfs, diff_dfs):
            merged_dfs.append(pd.concat([fcd, dd], axis=1).dropna())

        result = self.apply_low(merged_dfs)
        result = self.apply_hi(result)
        self.filtered_genes = result  # A list of dataframes
        
    def split(self):
        return [self.normalized_df[[control_col, treated_col]] for (control_col, treated_col) in self.file_pairs]

    def apply_fold_change(self, input_df_list):

        fold_change_dfs = list()
        filtered_results = list()
        for (control_col, treated_col), df in zip(self.file_pairs, input_df_list):

            result = df[control_col].div(df[treated_col])

            fold_change_dfs.append(result)
            filtered_results.append(result[result > self.log2fold])            

        self.fold_change_filtered_dfs = pd.concat(fold_change_dfs, axis=1)
        return filtered_results

    def apply_diff(self, input_df_list):

        diff_dfs = list()
        filtered_results = list()
        for (control_col, treated_col), df in zip(self.file_pairs, input_df_list):

            result = df[control_col].sub(df[treated_col]).abs()

            diff_dfs.append(result)
            filtered_results.append(result[result > self.diff])            

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
            passing = df[(df[control_col] > self.low) | (df[treated_col] > self.low)]
            filtered_results.append(passing)
            
            deactivated = df[(df[control_col] > self.low) & (df[treated_col] < self.low)]
            activated = df[(df[control_col] < self.low) & (df[treated_col] > self.low)]    
            undetected = df[(df[control_col] < self.low) & (df[treated_col] < self.low)]
            state_change.append((deactivated[deactivated[control_col].sub(deactivated[treated_col]).abs() > self.diff],
                                 activated[activated[control_col].sub(activated[treated_col]).abs() > self.diff],
                                 undetected))

        self.state_change = state_change

        return filtered_results

    def apply_hi(self, input_df_list):

        filtered_results = list()
        saturated_dfs = list()
        for (control_col, treated_col), df in zip(self.file_pairs, input_df_list):

            passing = df[(df[control_col] < self.hi) | (df[control_col] < self.hi)]
            filtered_results.append(passing)

            saturated = df[(df[control_col] > self.hi) & (df[control_col] > self.hi)]
            saturated_dfs.append(saturated)

        self.saturated_dfs = saturated_dfs

        return filtered_results