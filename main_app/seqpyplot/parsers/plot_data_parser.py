import pandas as pd


class PlotDataParser(object):

    def __init_(self):
        pass

    def parse_data(self, datafile):
        """This reads a preformatted data file output from SeqPyPlot."""
        df = pd.read_csv(datafile)
        df.set_index('gene')
        ercc_df = df[df.index.str.startswith('ERCC-')]
        return df, ercc_df
