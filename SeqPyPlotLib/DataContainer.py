# coding=utf-8
import csv
import os
import codecs
import subprocess
import sys
import numpy as np
import pandas as pd

from opertor import reduce
from normalizer import norm_tmm as TMM

from parsers import PlotDataParser, Htseqparser
from ConfigParse import config_parser

from utils, import file_is_empty


PARSERS = {'htseq': htSeqParser,
           'cuffnorm': cuffnormParser}

class DataContainer(object):
    """"
    Class for holding normalized data in a standard format. Calls parsers to collect data in to 
    pandas data frames. Matrix generation and normalization happens here
    
    """

    def __init__(self, args, Optimize=False):
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
        self.analyzed = self.analyzed()  # Function call to datermine if contain contains flitered data

        self.args = args
 
        # if Optimize is True:
        #     print("Engaging optimization mode.")
        #     self.gene_map, self.ercc_map, self.data_frame_header = self.parse_plot_data(self.args.plot_data)

        # First need to perform config parsing
        paths, names = config_parser(self.args.config_path)

        # we need to call a parser, and then create the matrix for normalization
        # WOo! Good design - if the raw_data arg is set, this gets bypassed (doesn't matter what its set as)

        # if the raw_data argument is set to true (by default)   
        if self.args.raw_data:

            # Instantiante parser
            parser = PARSERS[self.args.datatype](self.args.num)
            # Execute parser given the data paths and the sample names
            raw_df = parser.parse_data(paths, names)
            matrix = self._create_unnormalized_matrix(raw_df)
            # TODO test this and make sure it actually... returns a df...
            normalized_matrix = TMM(matrix)
            # TODO recombine the zeroed genes in to the normalized matrix
            self.df = self.recombine_genes(normalized_matrix, zeroed_genes)

        # if the raw_data argument is false, an analyzed data file is passed 
        else:
            # TODO add check to plotdataparser that checks input
            parser = PlotDataParser()
            self.df = parser.parse_data(datafile)
        
        print "Data Container initialized Successfully....\n"

        return df

    def _create_unnormalized_matrix(self, df):
        """Create the raw matrix for TMM normalization."""

        matrix_path = os.path.join('.', self.args.out, self.args.prefix + '_count_matrix.txt')
        empty_files = [sample_names[idx] for ix, file in enumerate(data_paths) 
                        if file_is_empty(file)]
        keep_cols = set(df.columns) - set(empty_files)
        df = df[keep_cols]

        #drop zeroed rows
        df = df[df.values.sum(axis=1) != 0]
        
        #write out matrix csv
        df.to_csv(matrix_path, sep='\t')

        return df


    #Sanity Check
    def analyzed(self):
        if self.args.plot_data is None:
            return False
        else:
            return True



    def __average_flanking__(self, value):
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



