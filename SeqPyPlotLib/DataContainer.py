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

from htSeqParser import (_config_parser_,
                         file_is_empty,
                         load_htseq_dataframe,
                         create_unnormalized_matrix)


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
        """

        self.args = args
        self.analyzed = self.analyzed()  # Function call to datermine if contain contains flitered data
        datatype = self.args.datatype
        self.gene_map = dict()
        self.ercc_map = dict()
        self.data_frame_header = dict()
        self.raw_counts = dict()
        self.args.unformatted_plot_data = dict()

        if Optimize is True:
            print("Engaging optimization mode.")
            self.gene_map, self.ercc_map, self.data_frame_header = self.parse_plot_data(self.args.plot_data)


        elif self.args.raw_data is not None:

            
            if datatype == 'cuffnorm':
                print "Reading data as cuffnorm. (This method is not recommended.)"
                self.gene_map, self.ercc_map, self.data_frame_header = self.parse_cuffnorm(self.args.raw_data)

            elif datatype == 'htseq':

                print "Reading data from FOLDER as HTSeq counts and normalizing using edgeR Methodology..."
                print "See: Loraine, A.E. et al., Analysis and Visualization of RNA-Seq Expression Data Using RStudio, Bioconductor, and Integrated Genome Browser.(2015)\n"
                
                data_paths, sample_names = self._config_parser_(config_path)
                number_of_conditions = self.args.num
                self.gene_map, self.ercc_map, self.data_frame_header = self.load_htseq_dataframe(data_paths, 
                                                                                                 sample_names,
                                                                                                 number_of_conditions)

            else:
                print "Raw_data type selected: ", self.args.datatype
                print "Software under development. Other raw data types will be added for parsing in the future."
                sys.exit()

        elif self.args.plot_data is not None:
            if self.args.gene_list is not None:
                self.gene_map, self.ercc_map, self.data_frame_header = self.parse_plot_data(self.args.plot_data)
            else:
                self.gene_map, self.ercc_map, self.data_frame_header = self.parse_plot_data(self.args.plot_data)
        else:
            "Anomaly. Exiting."
            sys.exit()
        print "Data Container initialized Successfully....\n"


    def create_unnormalized_matrix(self, df):
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




    def parse_plot_data(self, datafile):
        """This reads a preformatted data file output from SeqPyPlot."""
        # type: (input normalized data file) -> dictionary, key=gene name, value=normalized expression data
        # convert plotter data to a dictionary for quick access
        gene_map = dict()
        data_frame_header = dict()
        ercc_map = dict()

        if os.path.exists(datafile):
            with open(datafile, 'rb') as datafile:
                data_reader = csv.reader(datafile, delimiter='\t')

                header = True
                for row in data_reader:
                    if header:
                        data_frame_header[row[0]] = row[1:]
                        header = False
                    else:
                        temprow = row[1:]
                        finalrow = []
                        for i in temprow:
                            if i == '':
                               finalrow.append(None)
                            elif i < 1.0:
                                finalrow.append(0.0)
                            else:
                                finalrow.append(i)
                        gene_map[row[0].capitalize()] = finalrow
                if self.args.num == 1:
                    pass
                elif self.args.num == 2:
                    if self.args.unformatted_plot_data:
                        self.reorder(data_frame_header)
                        # print self.data_frame_header
                        self.reorder(gene_map)
                        self.reorder(ercc_map)
                        for key, value in gene_map.items():
                            series1 = value[:len(value) / 2]  # split the data
                            series2 = value[len(value) / 2:]
                            gene_map[key] = self.__average_flanking__(series1) + self.__average_flanking__(series2)
                else:
                    print "num == more than 2 - does not support. DataConta. line 347"
                    sys.exit()
        else:
            print "Couldn't open *_plotter_file.txt"
            sys.exit()

        return gene_map, ercc_map, data_frame_header


class MakeFigureList(object):
    """
    Function to parse the input gene list.
    """

    def __init__(self, args):
        """
        :param args: Args object
        """
        self.args = args
        self.gene_list = []  # The full gene list from the input file
        self.figure_list, self.gene_list = self.input_list_parser()  # nested list - gene_list broken in to gruops of 6

    def input_list_parser(self):

        """
        # type: (input_file object)
        :return: 2D nested list object

        """
        # Handle file encodings when you open the input file
        for e in ["utf-8", "ascii", "ansi"]:
            try:
                # create list of gene names ['str1', 'str2', etc]
                genes_in = codecs.open(self.args.gene_list, 'r', encoding=e)
                for row in genes_in:
                    if row.rstrip() != '':
                        self.gene_list.append(str(row.rstrip().capitalize()))
            except UnicodeDecodeError:
                print("File is encoded in a format other than {}.".format(e))

            else:
                print("Parsing file using {} encoding.".format(e))
                break

        # Handle the case where ultimately the file can't be opened
        if len(self.gene_list) == 0:
            print("File list not parsed properly. Save file as either utf-8 or ascii encoded text.")
            raise RuntimeError

        # subdivide in put files with more than 6 gene names.
        if len(self.gene_list) >= 6:

            sub_list = []

            for i in range((len(self.gene_list) / 6) + 1):
                sub_list.append(self.gene_list[:6])
                del (self.gene_list[:6])
                self.figure_list = [x for x in sub_list if x != []]

            return self.figure_list, self.gene_list

        else:
            self.figure_list = [[x for x in self.gene_list]]

            return self.figure_list, self.gene_list

    def figure_list_length(self):
        if len(self.figure_list) == 0:
            print("Your gene list is currently empty or couldn't be opened.")
        else:
            return len(self.figure_list)

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


