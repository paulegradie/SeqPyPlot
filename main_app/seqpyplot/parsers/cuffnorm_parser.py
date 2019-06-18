import pandas as pd


class CuffNormParser(object):
    """
    A class for parsing cufflnorm outputs and reading them in to the standard
    pandas df datafram structure used by the DataContainer
    """
    def __init__(self):
        pass

    def parse_data(self, datafile):
        """
        This function expects output from Cuffnorm. If the ERCC option is set, it will also return any ERCCs.
            Returns: dictionaries with Gene:[ExpressionData]
        """

        df = pd.read_csv(datafile)
        df.set_index('gene')

        ercc_cols = df.index.str.startswith('ERCC-')
        ercc_df = df[ercc_cols]
        data = df[~ercc_cols]

        return data, ercc_df
