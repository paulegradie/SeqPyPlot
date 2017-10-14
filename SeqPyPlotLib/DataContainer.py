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
        parser_dict = {'htseq': htSeqParser,
                       'cuffnorm': cuffnormParser}

        # if the raw_data argument is set to true (by default)   
        if self.args.raw_data:

            # Generalize the parser and feed in the same arguments (paths and sample names) to each parser
            # Instantiante parser
            parser = parser_dict[self.args.datatype](self.args.num)

            # Execute parser given the data paths and the sample names
            df = parser.parse_data(paths, names)
            matrix = self._create_unnormalized_matrix(df)
            # TODO test this and make sure it actually... returns a df...
            normalized_matrix = TMM(matrix)
            # TODO recombine the zeroed genes in to the normalized matrix
            refilled_matrix = self.recombine_genes(normalized_matrix, zeroed_genes)



        # if the raw_data argument is false, an analyzed data file is passed 
        else:
            # TODO add check to plotdataparser that checks input
            parser = PlotDataParser()
            df = parser.




        print "Data Container initialized Successfully....\n"


    def _create_unnormalized_matrix(self, df):
        """Create the raw matrix for TMM normalization."""

        matrix_path = os.path.join('.', self.args.out, self.args.prefix + '_count_matrix.txt')

        empty_files = [sample_names[idx] for ix, file in enumerate(data_paths) 
                        if self.file_is_empty(file)]
        keep_cols = set(df.columns) - set(empty_files)
        df = df[keep_cols]

        #drop zeroed rows
        df = df[df.values.sum(axis=1) != 0]
        
        #write out matrix csv
        df.to_csv(matrix_path, sep='\t')

        return df







        #     # If any data points are missing:::
        #     # reinsert missing data points with zeros
        #     if len(insert_zero_list) > 0:  # If any data is missing...
        #         insertion_list = zip(insert_zero_list, missing_file_name)
        #         for position in insertion_list:
        #             self.data_frame_header["Gene"].insert(position[0], position[1])
        #             for key, value in self.gene_map.items():
        #                 value.insert(position[0], None)

        #     if self.args.num == 1:
        #         pass

        #     elif self.args.num == 2:
        #         # Reorder the data so that series1 comes before series2
        #         self.reorder(self.data_frame_header)
        #         self.reorder(self.gene_map)
        #         self.reorder(self.ercc_map)

        #         # produce flaking averages for each series if any stage is missing.
        #         if len(insert_zero_list) > 0:  # If any data is missing...average flanking data
        #             for key, value in self.gene_map.items():
        #                 series1 = value[:len(value) / 2]  # split the data
        #                 series2 = value[len(value) / 2:]
        #                 self.gene_map[key] = self.__average_flanking__(series1) + self.__average_flanking__(series2)

        #         # finaly, add in zeroed out genes
        #         for key in self.raw_counts.keys():
        #             if key not in self.gene_map.keys():
        #                 self.gene_map[key] = [0] * int(len(self.data_frame_header["Gene"]))
        #             else:
        #                 pass
        #     elif self.args.num > 2:
        #         print "Data Container line 229 - does not support more than three series yet!"

        # return self.gene_map, self.ercc_map, self.data_frame_header



    #Sanity Check
    def analyzed(self):
        if self.args.plot_data is None:
            return False
        else:
            return True

    #Load in premade results from elsewhere
    def load_results(self):
        if self.args.filter_results is None:
            self.de_gene_list = []
            return self.de_gene_list
        else:
            with open(self.args.filter_results, 'r') as de_results:
                self.de_gene_list = [result.rstrip() for result in de_results.readlines()]
                # print self.de_gene_list
                return self.de_gene_list



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





class PrepareOutputDirectory:
    """
    Must instantiate after args are set
    """

    def __init__(self):
        pass

    @staticmethod
    def make_folder(path):
        outdir = os.path.join(path)

        if not os.path.isdir(outdir):
            os.makedirs(outdir)


