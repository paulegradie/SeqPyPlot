# coding=utf-8
import csv
import os
import codecs
import subprocess
import sys
import numpy as np

class DataContainer(object):
    """"Functions for parsing input files from various programs."""

    def __init__(self, args, Optimize=False):

        self.args = args
        self.analyzed = self.analyzed()  # Function call to datermine if contain contains flitered data
        datatype = self.args.datatype
        self.gene_map = dict()
        self.ercc_map = dict()
        self.data_frame_header = dict()
        self.raw_counts = dict()
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

                self.gene_map, self.ercc_map, self.data_frame_header = self.parse_htseq()

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
        print "Data Container initialized Successfully...."
    def __getitem__(self, key):
        return self.data_frame_header[key]

    def analyzed(self):
        if self.args.plot_data is None:
            return False
        else:
            return True

    def load_results(self):
        if self.args.filter_results is None:
            self.de_gene_list = []
            return self.de_gene_list
        else:
            with open(self.args.filter_results, 'r') as de_results:
                self.de_gene_list = [result.rstrip() for result in de_results.readlines()]
                # print self.de_gene_list
                return self.de_gene_list

    @staticmethod
    def reorder(dictionary):
        for key, value in sorted(dictionary.items()):
            new_order = [y for y in range(len(value)) if y % 2 == 0] + [y for y in range(len(value)) if y % 2 != 0]
            dictionary[key] = [value[t] for t in new_order]
            key.capitalize()

        return dictionary

    @staticmethod
    def __average_flanking__(value):

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

    def parse_htseq(self):

        insert_zero_list = []
        missing_file = 0  # counter for recording empty files for later  fill with zeros
        missing_file_name = []
        self.data_frame_header = dict()

        if os.path.isdir(self.args.raw_data):

            for datafolder in os.walk(os.path.join('.', self.args.raw_data)):
                if self.args.num == 2:
                    if len(datafolder[2]) % 2 != 0:
                        print "Gotta have even number of input files for num = 2, You can use an empty file."
                        sys.exit()
                else:
                    print("Still working on support for single time series data...")
                current_file_count = 1
                tot = len(datafolder[2])
                self.data_frame_header["Gene"] = []

                if len(self.args.time) != int(tot/2):
                    print "Check your -time option, make sure it matches the no. of input files."
                    sys.exit()

                for _file in datafolder[2]:
                    file_string = os.path.join('.', self.args.raw_data, _file)

                    if int(os.stat(file_string).st_size) == 0:  # if file is empty, record the list position
                        print(str(_file) + '...  ' + str(current_file_count) + '/' + str(tot) + ' --empty file')
                        insert_zero_list.append(missing_file)
                        missing_file_name.append('.'.join([str(_file[0:2]), str(_file[5:8])]))
                        missing_file += 1
                        current_file_count += 1

                        continue

                    self.data_frame_header["Gene"] += ['E' + '.'.join([str(_file[0:2]), str(_file[5:8])])]
                    print(str(_file) + '...  ' + str(current_file_count) + '/' + str(tot))
                    with open(file_string, 'r') as countfile:
                        # ht_reader = csv.reader(countfile, delimiter='\t')
                        # for row in ht_reader:
                        try:
                            for row in countfile.readlines():
                                if row.split()[0].rstrip() not in self.raw_counts.keys():
                                    self.raw_counts[str(row.split()[0].rstrip())] = [int(row.split()[1].rstrip())]
                                else:
                                    self.raw_counts[str(row.split()[0].rstrip())].append(row.split()[1].rstrip())
                        except IndexError:
                            "Sorry for that! Make sure there are no unrelated files in the counts directory."
                    current_file_count += 1
                    missing_file += 1
        else:
            print("ERROR...")
            print("Directory listed not available. Please pass the path to a folder with ordered htseq data.")
            print("Alternatively, you may need to reset the -data_type parameter. Use -h for usage.")
            sys.exit()

        matrix_path = os.path.join('.', self.args.out, self.args.prefix + '_count_matrix.txt')
        print 'matrix path: ', matrix_path

        with open(matrix_path, 'wb+') as matrix:
            matrix_writer = csv.writer(matrix, delimiter='\t')
            for key, value in self.data_frame_header.items():
                matrix_writer.writerow([key] + value)
            for key, value in self.raw_counts.items():
                if '_' in key:
                    pass
                else:
                    matrix_writer.writerow([key] + value)
        # Run Rscript to normalize the data via EdgeR
                # The following normalization method was modified from\n""")
                # Loraine, A.E., Blakley, I.C., Jagadeesan, S. Harper, J., Miller, G., and Firon, N. (2015).
                # Analysis and Visualization of RNA-Seq Expression Data Using
                # RStudio, Bioconductor, and Integrated Genome Browser. pp. 481–501.

        print "Attempting data normalization by edgeR...\n"
        try:
            norm_path = os.path.join('.', self.args.out, self.args.prefix + '_normalized_count_data.txt')
            print 'norm_path: ', norm_path
            subprocess.call('Rscript '
                            + os.path.join('.',
                                            'SeqPyPlotLib',
                                            'Normalization_Method.R ')
                            + matrix_path
                            + ' '
                            + norm_path)

            print "Data Normalized Successfully"

        except WindowsError:
            print(""" You need R installed and in your path to use this option.

            Do you have Rscript added in your system path?

            Google 'Windows Environment Variables' for help.

            In windows for example its location is
            C:\Program Files\R\R-3.2.1\bin

            You have to add the path of Rscript.exe in you system path in environment variables.""")

        # normalized_data_path = os.path.join('.', self.args.out,  'normalized_count_data.txt')
        # normalized_data_path = os.path.join('normalized_count_data.txt')

        with open(norm_path, 'r') as normalized:
            normalized_reader = csv.reader(normalized, delimiter='\t')

            # add normalized data to the gene_map
            for row in normalized_reader:
                self.gene_map[row[0]] = [int(x) for x in row[1:]]

            # If any data points are missing:::
            # reinsert missing data points with zeros
            if len(insert_zero_list) > 0:  # If any data is missing...
                insertion_list = zip(insert_zero_list, missing_file_name)
                for position in insertion_list:
                    self.data_frame_header["Gene"].insert(position[0], position[1])
                    for key, value in self.gene_map.items():
                        value.insert(position[0], None)

            # Reorder the data so that series1 comes before series2
            self.reorder(self.data_frame_header)
            self.reorder(self.gene_map)
            self.reorder(self.ercc_map)

            # produce flaking averages for each series if any stage is missing.
            if len(insert_zero_list) > 0:  # If any data is missing...average flanking data
                for key, value in self.gene_map.items():
                    series1 = value[:len(value) / 2]  # split the data
                    series2 = value[len(value) / 2:]
                    self.gene_map[key] = self.__average_flanking__(series1) + self.__average_flanking__(series2)

            # finaly, add in zeroed out genes
            for key in self.raw_counts.keys():
                if key not in self.gene_map.keys():
                    self.gene_map[key] = [0] * int(len(self.data_frame_header["Gene"]))
                else:
                    pass

        return self.gene_map, self.ercc_map, self.data_frame_header

    def parse_cuffnorm(self, infile):
        """This function expects output from Cuffnorm. If the ERCC option is set, it will also return any ERCCs.
            Returns: dictionaries with Gene:[ExpressionData]"""
        self.data_frame_header = dict()
        self.gene_map = dict()
        self.ercc_map = dict()
        self.gene_map_const = dict()
        # open the files one at a time from within the loop
        try:
            with open(infile, 'rb') as dataFile:

                data_frame = csv.reader(dataFile, delimiter='\t')

                # for each row, ignore the data_frame_header row, then execute the following code
                for gene in data_frame:

                    zero = False
                    count = 9
                    if "tracking_id" in str(gene[0]):  # generate data_frame_header
                        self.data_frame_header["Gene"] = []
                        for column_header in enumerate(gene):
                            if column_header[0] >= count:
                                self.data_frame_header['Gene'] += [column_header[1]]
                                count += 4

                    elif "ERCC-" in gene[3]:  # collect ERCC data if ERCC option set
                        self.ercc_map[gene[3]] = []
                        for column in enumerate(gene):
                            if column[0] == count:
                                self.ercc_map[gene[3]] += [column[1]]
                                count += 4

                    else:  # Collect all data

                        self.gene_map[gene[3]] = []
                        for column in enumerate(gene):
                            if int(column[0]) == count:
                                self.gene_map[gene[3]] += [column[1]]
                                if '0' == str(column[1]):
                                    zero = True
                                count += 4

                        if zero is False:
                            self.gene_map_const[gene[3]] = self.gene_map[gene[3]]
        except IOError:
            print "/nAre you loading cuffnorm data? - Try resetting data_type argument."

        self.reorder(self.data_frame_header)
        self.reorder(self.gene_map)
        self.reorder(self.ercc_map)

        return self.gene_map, self.ercc_map, self.data_frame_header

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
                rowcount = 0

                for row in data_reader:
                    if rowcount == 0:
                        data_frame_header[row[0]] = row[1:]
                    else:
                        gene_map[row[0].capitalize()] = row[1:]
                    rowcount += 1

        else:
            print "Couldn't open file."
            sys.exit()

        return gene_map, ercc_map, data_frame_header

    def parse_edgeR(self, infile, filelist=None):
        pass

    def parse_deseq2(self, infile):
        pass

    def parse_cuffdiff(self, infile):
        pass

    def generic_parser(self, inputfile):
        """This Generic filter will simply split your input in to two files - a data file and ERCC file."""
        # TODO copy over code to fill this out
        pass


class MakeFigureList(object):
    """Function to parse the input gene list."""

    def __init__(self, args):

        self.args = args
        self.gene_list = []  # The full gene list from the input file
        self.figure_list, self.gene_list = self.input_list_parser()  # nested list - gene_list broken in to gruops of 6

    def input_list_parser(self):  # type: (input_file object) -> 2D nested list object

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


class PrepareOutputDirectory:  # must instantiate after args are set

    def __init__(self):
        pass

    @staticmethod
    def make_folder(path):
        outdir = os.path.join(path)

        if not os.path.isdir(outdir):
            os.makedirs(outdir)


class DataPrinter:

    """
    This should have a data container and analyzer instance passed upon instantiation.

    Available upon passing container object:
        self.args = args
        self.datatype = self.args.datatype
        self.data_frame_header = dict()
        self.gene_map = dict()
        self.ercc_map = dict()
        self.gene_map_const = dict()
        self.de_gene_list = []
        self.de_count_by_stage = None
        self.de_gene_list = None
        self.analyzed = self.analyzed()
        self.filtered_data = dict()
    """
    def __init__(self, args, container, analyzer):

        self.container = container
        self.analyzer = analyzer
        self.filtered_data = analyzer.filtered_data
        self.args = args

        if self.args.out is not None:
            self.path = os.path.join(self.args.out, self.args.prefix)
        else:
            self.path = os.path.join('.', self.args.prefix)

    def write_ercc_data(self):

        with open((self.path + ".all_ercc.txt"), "wb+") as ercc:
            ercc_data_writer = csv.writer(ercc, delimiter='\t')
            for key, value in self.container.data_frame_header.items():
                ercc_data_writer.writerow(["ERCC"] + value)
            for key, value in self.container.ercc_map.items():
                ercc_data_writer.writerow([key] + value)

    def write_plot_data(self):
        with open(self.path + '_plotter_data.txt', 'wb+') as plot:
            plot_data_writer = csv.writer(plot, delimiter='\t')

            for key, value in self.container.data_frame_header.items():
                plot_data_writer.writerow([key] + value)

            for key, value in sorted(self.container.gene_map.items()):
                plot_data_writer.writerow([key] + value)

    def write_filtered_data(self):
        with open((self.path + '_filtered.txt'), "wb+") as filt:

            filtered_data_writer = csv.writer(filt, delimiter='\t')
            for key, value in self.container.data_frame_header.items():
                filtered_data_writer.writerow([key] + value)

            for key, value in self.filtered_data.items():
                filtered_data_writer.writerow([key] + value)

    def write_de_results(self):
        with open(self.path + "_de_gene_list.txt", 'wb+') as results_out:
            result_writer = csv.writer(results_out, delimiter='\n')
            for de_gene in sorted(self.analyzer.filtered_data):
                result_writer.writerow([de_gene])
