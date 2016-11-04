__author__="Paul Gradie"

import argparse
import csv
import os
import platform
import sys
import codecs

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np


class PrepareOutputDirectory:
    def __init__(self, out_folder):

        if platform.system() == "Windows":
            results_dir = str(out_folder) + '\\'
            print "Windows OS"
        else:
            results_dir = str(out_folder) + '/'
            print "Unix OS"

        if not os.path.isdir(results_dir):
            os.makedirs(results_dir)


class ReadNormalizedData:
    # Parse the data sets.

    def __init__(self):
        pass

    @staticmethod
    def read_data(datafile):
        # type: (input normalized data file) -> dictionary, key=gene name, value=normalized expression data
        gene_map = dict()

        # convert plotter data to a dictionary for quick access
        with open(datafile, 'rb') as datafile:

            data_reader = csv.reader(datafile, delimiter='\t')

            for gene in data_reader:
                gene_map[gene[0].title()] = gene[1:]

        return gene_map

    ##TODO Add functions to autoread more data from other normalization sources

    @staticmethod
    def data_collector(gene_map, gene_name, args):

        data_length = len(gene_map.values()[0])
        found = True

        try:
            gene_row = gene_map[gene_name]

        except KeyError:
            # print "The current gene: -- {} -- was not found in the plot data.".format(gene_name)
            found = False
            gene_row = [0] * data_length

            # TODO Rename this argument to match
        if args.num == '1':
            if found is True:
                series1 = np.array(
                    [float(x) if x is not '' else None for x in gene_row[:(len(gene_row))]]).astype(np.double)
            else:
                series1 = np.array(
                    [float(0) for x in range(data_length)]).astype(np.double)
            s1mask = None
            s2mask = None
            series2 = None

        elif args.num == '2':
            if found is True:
                try:
                    series1 = np.array(
                        [float(x) if x is not '' else None for x in gene_row[:(len(gene_row) / 2)]]).astype(np.double)
                    series2 = np.array(
                        [float(x) if x is not '' else None for x in gene_row[(len(gene_row) / 2):]]).astype(np.double)
                    s1mask = np.isfinite(series1)
                    s2mask = np.isfinite(series2)
                except IndexError:
                    print "The program expects even number of data columns with either data or '' for missing time points."
                    print "Reformat your data."
                    sys.exit()
            else:
                series1 = np.array([float(0) for x in range(data_length / 2)]).astype(np.double)
                series2 = np.array([float(0) for x in range(data_length / 2)]).astype(np.double)
                s1mask = np.isfinite(series1)
                s2mask = np.isfinite(series2)

        elif args.num > 2:
            print "You can only plot up to two experiments at the moment. Future updates will be more flexible."
            sys.exit()

        # TODO-use found  as: if TRUE - keep gene name, if FALSE - gene name=gene_name + "-Not Found"
        else:
            print "Cannot use negative number for -num argument."
            sys.exit()

        return found, series1, s1mask, series2, s2mask


class ParseGeneList:
    # Parse and format the input file gene list.

    def __init__(self):

        self.gene_list = []  # The full gene list from the input file
        self.extensions = ["utf-8", "ascii", "ansi"]
        self.figure_list = []

    def input_list_parser(self, gene_list_file):  # type: (input_file object) -> 2D nested list object

        # Handle file encodings when you open the input file
        for e in self.extensions:
            try:
                # create list of gene names ['str1', 'str2', etc]
                genes_in = codecs.open(gene_list_file, 'r', encoding=e)
                self.gene_list = [row.rstrip().title() for row in genes_in]

            except UnicodeDecodeError:
                print "File is encoded in a format other than {}.".format(e)

            else:
                print "Parsing file using {} encoding.".format(e)
                break

        # Handle the case where ultimately the file can't be opened
        if len(self.gene_list) == 0:
            print "File list not parsed properly. Save file as either utf-8 or ascii encoded text."
            raise RuntimeError

        # subdivide in put files with more than 6 gene names.
        if len(self.gene_list) >= 6:

            sub_list = []

            for i in range((len(self.gene_list) / 6) + 1):
                sub_list.append(self.gene_list[:6])
                del (self.gene_list[:6])
                self.figure_list = [x for x in sub_list if x != []]

            return self.figure_list

        else:
            self.figure_list = [x for x in self.gene_list]

            return self.figure_list

    def figure_list_length(self):
        if len(self.figure_list) == 0:
            print "Your gene list is currently empty or couldn't be opened."
        else:
            return len(self.figure_list)


class MainDataPlotter:
    # For  building a all of the figures

    def __init__(self):
        pass

    def plot_figures(self, plot_data, time_label, condition_label, figure_list, figure_name, args):
        print
        # Figure list is a 2D nested list.
        # Outerlist = Figure list.
        # Innerlist = gene lists for each figure.

        figure_list_count = 1
        for figure in figure_list:

            if len(figure) == 1:
                fig_size = (7, 7)
            else:
                fig_size = (10, 10)

            fig = plt.figure(num=1,
                             dpi=300,
                             figsize=fig_size,
                             edgecolor='black',
                             frameon=False,
                             )
            fig.suptitle(figure_name,
                         verticalalignment='top',
                         horizontalalignment='center',
                         fontsize=12,
                         x=0.415
                         )
            series1_line = mlines.Line2D([], [], color="blue",
                                         marker='o',
                                         linewidth=1.4,
                                         linestyle="-",
                                         dash_capstyle='round',
                                         dash_joinstyle='bevel',
                                         fillstyle='full',
                                         markeredgecolor='blue',
                                         markerfacecolor='white',
                                         markeredgewidth=.95,
                                         markersize=6)

            series2_line = mlines.Line2D([], [], color="red",
                                         marker='o',
                                         linewidth=1.4,
                                         linestyle="-",
                                         dash_capstyle='round',
                                         dash_joinstyle='bevel',
                                         fillstyle='full',
                                         markeredgecolor='black',
                                         markerfacecolor='red',
                                         markeredgewidth=0.95,
                                         markersize=6)
            gene_count = 1
            for gene in figure:
                if len(figure) == 1:
                    fig.add_subplot(1, 1, gene_count).set_title('- ' + str(gene) + ' -', fontsize=20)
                else:
                    fig.add_subplot(3, 2, gene_count).set_title('- ' + str(gene) + ' -', fontsize=20)

                self.__subplot(plot_data, gene, time_label, args)
                gene_count += 1

            # generate condition labels
            try:
                condition_labels = (x for x in condition_label.split(','))
            except AttributeError:
                print "The condition labels have to be comma separated."
                sys.exit()

            fig.legend(handles=[series1_line, series2_line],
                       labels=condition_labels,
                       loc='upper right')

            plt.savefig("{0}\\{1}_{2}.png".format(figure_name,
                                                  str(figure_list_count),
                                                  str(figure)),
                        format='png',
                        bbox_inches='tight')

            plt.savefig("{0}\\{1}_{2}.svg".format(figure_name,
                                                  figure_list_count,
                                                  str(figure)),
                        format='svg',
                        bbox_inches='tight')

            # plt.show()
            plt.clf()
            plt.cla()
            plt.close()

            figure_list_count += 1

    def __subplot(self, plot_data, gene, time_label, args):

        data_collector = ReadNormalizedData()
        found, series1_data, series1_mask, series2_data, series2_mask = data_collector.data_collector(plot_data, gene, args)

        x_axis = np.arange(float(max(len(series1_data), len(series2_data))))

        try:
            # x axis parameters for each subplot
            plt.tick_params(axis='both',
                            direction='out',
                            which='major',
                            labelsize=14)

            # xlim is the length of the label list
            plt.xlim(-0.5, (len(time_label)))

            # Plot Details
            ax = plt.gca()
            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.spines['bottom'].set_position(('data', 0))
            plt.xticks([i for i in range(len(time_label))], time_label)

            # fill the plot with data points and graph them
            self.__series1_plot(series1_data, series1_mask, x_axis, label='Series1')
            self.__series2_plot(series2_data, series2_mask, x_axis, label='Series2')

            try:
                y_maximum = float(max(max(series1_data), max(series2_data))) * 1.3
            except ValueError:
                print "Hmm, {0} wasn't found in the plot-able dataset.".format(gene)
                y_maximum = 10

            # Extend subplot frame for high fpkm values
            plt.ylim(0, y_maximum)
            plt.tight_layout()

            # adjust subplot location
            plt.subplots_adjust(top=0.85)

            return plt

        except UnboundLocalError:
            print "Don't worry - this just means the gene is probably off."

    @staticmethod
    def __series1_plot(series1, s1mask, xs, label):
        return plt.plot(xs[s1mask],
                        series1[s1mask],
                        'bo',
                        color="blue",
                        marker='o',
                        linewidth=1.4,
                        linestyle="-",
                        dash_capstyle='round',
                        dash_joinstyle='bevel',
                        label=label,
                        fillstyle='full',
                        markeredgecolor='blue',
                        markerfacecolor='white',
                        markeredgewidth=.95,
                        markersize=6)

    @staticmethod
    def __series2_plot(series2, s2mask, xs, label):
        # xs = [float(img_folder) for img_folder in range(len(series2))]

        return plt.plot(xs[s2mask],
                        series2[s2mask],
                        'ro',
                        color="red",
                        marker='o',
                        linewidth=1.4,
                        linestyle="-",
                        dash_capstyle='round',
                        dash_joinstyle='bevel',
                        label=label,
                        fillstyle='full',
                        markeredgecolor='black',
                        markerfacecolor='red',
                        markeredgewidth=.95,
                        markersize=6)


class ArgumentParser:
    def __init__(self):
        pass

    @staticmethod
    def argument_parser():
        # type: () -> NameSpace()

        parser = argparse.ArgumentParser(description='Plot time series expression data without replicates. \n')

        parser.add_argument('-t', '--time',
                            metavar='\b',
                            default='1,2,3,4,5',
                            type=str,
                            help='\tA comma separated list of time points.')
        parser.add_argument('-o', '--out',
                            metavar='\b',
                            default='Default_out',
                            type=str,
                            help='\tOutput Folder Name')

        parser.add_argument('-c', '--condition',
                            metavar='\b',
                            default=('Series1,Series2'),
                            type=str,
                            help='\tA comma separated list of conditions (max 2)')
        parser.add_argument('-n', '--num',
                            metavar='\b',
                            default='2',
                            type=str,
                            help='\tNumber of samples per plot.')
        parser.add_argument('gene_list',
                            nargs='?',
                            type=str,
                            default='genefile.txt',
                            help='\tSingle Column Gene list in txt file.')

        parser.add_argument('plotter_data',
                            nargs='?',
                            type=str,
                            default='plot_data.txt',
                            help='\t_Plotter_data.txt output from GetCuffData_v2.0.py')
        return parser.parse_args()

    @staticmethod
    def __label_parser(argument):

        # type: (character_string) -> list of strings
        try:
            parsed_list = [x for x in argument.split(',')]
            return parsed_list

        except AttributeError:
            print "The group labels have to be comma separated."
            sys.exit()

    def time_parser(self, time_label_arg):
        return self.__label_parser(str(time_label_arg))

    def condition_label_parser(self, condition_label_arg):
        return self.__label_parser(str(condition_label_arg))
