import pandas as pd


class PlotDataParser(object):

    def __init_(self):
        pass

    def parse_data(self, datafile, *args, **kwargs):
        """This reads a preformatted data file output from SeqPyPlot."""
        df = pd.read_csv(datafile, sep='\t')
        df.set_index('Gene', inplace=True)
        ercc_df = df[df.index.str.startswith('ERCC-')]
        return df, ercc_df
