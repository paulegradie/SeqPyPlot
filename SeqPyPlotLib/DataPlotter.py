from __future__ import division
# Test
import matplotlib.lines as mlines
import matplotlib.pyplot as plt

import numpy as np
import os
import sys


class MainDataPlotter(object):
    # For  building a all of the figures

    def __init__(self, args, analyzer, figurelist):

        self.args = args

        if figurelist is not None:
            self.gene_list = figurelist.gene_list
            self.figure_list = figurelist.figure_list
        else:
            self.gene_list = []
            self.figure_list = []

        self.de_genes = analyzer.de_gene_list
        self.de_count_by_time = analyzer.de_count_by_stage
        self.analyzer = analyzer
        self.gene_map = analyzer.gene_map
        self.data_frame_header = analyzer.data_frame_header

        if self.args.out is not None:
            self.path = os.path.join(self.args.out, self.args.prefix)
        else:
            self.path = os.path.join('.', self.args.prefix)

    def __check_for_de(self, gene):
        if gene in self.de_genes:
            return True
        else:
            return False

    def plot_figures(self):

        # Figure list is a 2D nested list.
        # Outerlist = Figure list.
        # Innerlist = gene lists for each figure.

        figure_list_count = 1
        for figure in self.figure_list:

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
            fig.suptitle(self.args.out,
                         verticalalignment='top',
                         horizontalalignment='center',
                         fontsize=24,
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

            series2_line = mlines.Line2D([], [],
                                         color="red",
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

            log_line = mlines.Line2D([], [], color='white')

            gene_count = 1
            for gene in figure:
                print "Plotting: ", gene
                if len(figure) == 1:
                    fig.add_subplot(1, 1, gene_count).set_title('- ' + str(gene) + ' -', fontsize=20)
                else:
                    fig.add_subplot(3, 2, gene_count).set_title('- ' + str(gene) + ' -', fontsize=20)

                self.__subplot(gene)

                gene_count += 1

            try:
                condition_labels = [x for x in self.args.condition]
            except AttributeError:
                print("The condition labels have to be comma separated.")
                sys.exit()

            fig.legend(handles=[series1_line, series2_line, log_line],
                       labels=(condition_labels + ["Log2: " + str(self.args.log)]),
                       loc='upper right')

            path = os.path.join(self.args.out,
                                self.args.prefix)

            plt.savefig("{}{}_{}.png".format(path,
                                             str(figure_list_count).replace(" ", ""),
                                             str([str(fig.replace(" ", "")) for fig in figure])),
                        format='png', bbox_inches='tight')

            # TODO implement svg argument
            # if self.args.svg:
            # plt.savefig("{}_{}_{}.png".format(path,
            #                                   str(figure_list_count).replace(" ", ""),
            #                                   str(figure).replace(" ", "")),
            #                 format='svg',
            #                 bbox_inches='tight')

            # plt.show()
            plt.clf()
            plt.cla()
            plt.close()

            figure_list_count += 1

    def __subplot(self, gene):

        found, series1_data, series1_mask, series2_data, series2_mask = self.data_collector(gene)

        np.asarray(series1_mask, dtype=bool)
        np.asarray(series2_mask, dtype=bool)

        x_axis = np.arange(float(max(len(series1_data), len(series2_data))))

        # ---------------
        def __calculate_mean__(series1, series2):
            """Return a numpy array and mask that contains the average of the series1 and series2 data.
                This array is built for each gene to be plotted."""

            temp_series_mean = []

            for itr in range(max(len(series1), len(series2))):
                try:
                    current_mean = np.mean([series1[itr], series2[itr]])
                except TypeError:
                    current_mean = None
                temp_series_mean.append(current_mean)

            series_mean_array = np.array([float(x) if x is not None else None for x in temp_series_mean]).astype(
                np.double)

            mean_mask = np.isfinite(series_mean_array)

            return series_mean_array, mean_mask

        try:
            # x axis parameters for each subplot
            plt.tick_params(axis='both',
                            direction='out',
                            which='major',
                            labelsize=14)

            # xlim is the length of the label list
            plt.xlim(-0.5, (len(series1_data)))

            # Plot Details
            ax = plt.gca()

            is_de = self.__check_for_de(gene)
            if is_de:
                ax.set_axis_bgcolor('yellow')
            else:
                pass

            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.spines['bottom'].set_position(('data', 0))

            if self.args.time[0] == 'None':
                # print "FALSE"
                labels = [x for x in self.analyzer.data_frame_header['Gene'] if
                          self.analyzer.data_frame_header['Gene'].index(x) % 2 != 0]
            else:
                # print "TRUE"
                labels = [i for i in self.args.time]

            plt.xticks([i for i in range(len(series1_data))], labels)

            # fill the plot with data points and graph them
            self.__series1_plot(series1_data, series1_mask, x_axis, label='Series1')
            self.__series2_plot(series2_data, series2_mask, x_axis, label='Series2')

            # error bars
            # return mean, and mask
            series_mean, mask = __calculate_mean__(series1_data, series2_data)

            err, upper = self.__err_range_plot(series_mean, x_axis, label="Range around the mean = log2(range)=1 ")

            try:
                # y_maximum = float(max(max(series1_data), max(series2_data))) * 1.3
                y_maximum = float(max(max(series1_data), max(series2_data), max(upper))) * 1.4
            except ValueError:
                print("Hmm, {0} wasn't found in the plot-able data set.".format(gene))
                y_maximum = 10

            # Extend subplot frame for high fpkm values
            plt.ylim(0, y_maximum)
            plt.tight_layout()

            # adjust subplot location
            plt.subplots_adjust(top=0.85)

            return plt

        except UnboundLocalError:
            print("Don't worry - this just means the gene is probably turned off.")

    def data_collector(self, gene_name):

        data_length = len(self.gene_map.values()[0])
        found = True

        try:
            gene_row = self.gene_map[gene_name]

        except KeyError:
            print "The current gene: -- {} -- was not found in the plot data.".format(gene_name)
            found = False
            gene_row = [0] * data_length

        if self.args.num == 1:
            if found is True:
                series1 = np.array(
                    [float(x) if x is not None else None for x in gene_row[:(len(gene_row))]]).astype(np.double)
            else:
                series1 = np.array(
                    [float(0) for x in range(data_length)]).astype(np.double)
            s1mask = None
            s2mask = None
            series2 = None

        elif self.args.num == 2:
            if found is True:
                try:
                    series1 = np.array(
                        [float(x) if x is not None else None for x in gene_row[:(len(gene_row) // 2)]]).astype(
                        np.double)
                    series2 = np.array(
                        [float(x) if x is not None else None for x in gene_row[(len(gene_row) // 2):]]).astype(
                        np.double)
                    s1mask = np.isfinite(series1)
                    s2mask = np.isfinite(series2)
                except IndexError:
                    print(
                        "The program expects even number of data columns with either data or '' for missing time points.")
                    print("Reformat your data.")
                    sys.exit()
            else:
                series1 = np.array([float(0) for x in range(data_length // 2)]).astype(np.double)
                series2 = np.array([float(0) for x in range(data_length // 2)]).astype(np.double)
                s1mask = np.isfinite(series1)
                s2mask = np.isfinite(series2)

        elif self.args.num > 2:
            print("You can only plot up to two experiments at the moment. Future updates will be more flexible.")
            sys.exit()

        # TODO-use found  as: if TRUE - keep gene name, if FALSE - gene name=gene_name + "-Not Found"
        else:
            print("Cannot use negative number for -num argument.")
            sys.exit()

        return found, series1, s1mask, series2, s2mask  # series masks are a True/False list for use in _calculate_mean

    # TODO: consolidate the series plot and make the differnces parameters to add in at calltime
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

    def __err_range_plot(self, series_mean, xs, label):

        upper, lower = self._calc_error(series_mean)

        # xs is x_axis
        return plt.errorbar(xs,
                            series_mean,
                            linestyle='None',
                            mfc='black',
                            yerr=[upper, lower],
                            ecolor='black',
                            elinewidth='1',
                            label=label), upper

    def _calc_error(self, series_mean):

        var = self.args.log

        upper_list = []
        lower_list = []

        for i in series_mean:
            b = (2.0 * i) / ((2.0 ** var) + 1)
            a = (2.0 ** var) * b

            upper_list.append(a - i)
            lower_list.append(i - b)
        upper = np.asarray([float(c) for c in upper_list])
        lower = np.asarray([float(d) for d in lower_list])

        if min(upper) < 0 or min(lower) < 0:
            print("Error bars out of range - but still continuing.")

        return upper, lower

    def de_bar(self, colour):
        plt.close()
        dfh = self.data_frame_header
        # assert len(self.de_count_by_stage) == len(self.args.time)

        y_values = np.asarray([int(v) for k, v in sorted(self.de_count_by_time.items())])
        x_axis = range(len(y_values))
        bar_width = float(0.6)  # the width of the bars
        try:
            ymax = max(y_values) * 1.3
        except ValueError:
            print("Didn't find any DE genes.")
            print("These plots will be empty")

        fig, ax = plt.subplots()
        ax.bar(x_axis, y_values, bar_width, color=colour, align="center")

        if self.args.time[0] == 'None':
            ax.set_xticklabels([''] + [dfh["Gene"][x] for x in range(int(len(dfh["Gene"]) / 2))], rotation='vertical')
        else:
            # print "TRUE"
            ax.set_xticklabels([''] + [i for i in self.args.time])
        # ax.set_xticklabels([''] + [str(x) for x in range(len(y_values))])
        try:
            ax.set_ylim([0, ymax])
        except UnboundLocalError:
            pass

        # add some text for labels, title and axes ticks
        ax.set_ylabel('No. of Flagged Genes')
        ax.set_xlabel('Experimental Stage')
        ax.set_title('Number of flagged genes per stage', loc='left')
        ax.yaxis.grid()

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)

        plt.tick_params(
            axis='both',  # changes apply to the x-axis
            which='both',  # both major and minor ticks are affected
            bottom='on',  # ticks along the bottom edge are off
            top='off',  # ticks along the top edge are off
            left='on',
            right='off',
            labelbottom='on')  # labels along the bottom edge are off

        log_line = mlines.Line2D([], [], color='yellow')
        expression_upper = mlines.Line2D([], [], color='orange')
        expression_lower = mlines.Line2D([], [], color='blue')
        difference = mlines.Line2D([], [], color='black')
        if self.args.hi > 99999:
            hi = 'inf'
        else:
            hi = self.args.hi

        # range =

        fig.legend(handles=[expression_upper, expression_lower, log_line, difference],
                   labels=(["Log2: " + str(self.args.log), "Range: " + str(self.args.low) + "-" + str(self.args.hi),
                            "Diff: " + str(self.args.dif)]),
                   loc='upper right')

        plt.savefig(self.path + "{}_DE_Gene_by_time.png".format(self.args.prefix), format='png', bbox_inches='tight')

        return None

    def plot_tally(self):

        cutoffs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9,
                   2.0, 3, 4, 5]

        line_plot_y_list = []
        iteration = len(cutoffs)
        current = 1
        print "Iterating over log2fold cutoff values... "

        for iterator in cutoffs:
            # print iterator
            self.analyzer.seqpyfilter(iterator)
            temp_de_count = self.analyzer.de_gene_list_length
            # print temp_de_count
            line_plot_y_list.append(int(temp_de_count))
            print "{} of {}.".format(current, iteration)

            self.analyzer.de_gene_list_length = 0
            current += 1

        # open a new line plot
        self.de_line_plot(cutoffs, line_plot_y_list, "DE cutoffs by log2fold")

    def de_line_plot(self, cutoffs, y_value, label):
        plt.close()

        # create a figure for the current figure
        fig = plt.figure(num=1,
                         figsize=(7, 7),
                         dpi=300,
                         edgecolor='black',
                         frameon=False,
                         )

        # set figure tit
        fig.suptitle("DE genes Detected vs log2Fold cutoff",
                     verticalalignment='top',
                     horizontalalignment='center',
                     fontsize=12,
                     x=0.415
                     )

        expression_upper = mlines.Line2D([], [], color='white')
        expression_lower = mlines.Line2D([], [], color='white')
        expression_dif = mlines.Line2D([], [], color='white')

        fig.legend(handles=[expression_upper, expression_lower, expression_dif],
                   labels=(["Upper: " + str(self.args.hi),
                            "Lower: " + str(self.args.low),
                            "Dif: " + str(self.args.dif)]),
                   loc='upper right')

        plt.plot(cutoffs,
                 y_value,
                 'bo',
                 # ymax=ymax,
                 color="black",
                 marker='o',
                 linewidth=1.4,
                 linestyle="-",
                 dash_capstyle='round',
                 dash_joinstyle='bevel',
                 label=label,
                 fillstyle='full',
                 markeredgecolor='black',
                 markerfacecolor='white',
                 markeredgewidth=.95,
                 markersize=6)

        # xlim is the length of the label list
        plt.xlim(0, 5.5)

        # Plot Details
        ax = plt.gca()
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data', 0))
        ax.set_ylabel('Number of DE genes Detected')
        ax.set_xlabel('log2Fold Cutoff')

        plt.savefig(self.path + "_DE_cutoff_by_log2fold.png",
                    format='png',
                    bbox_inches='tight')

    def plot_histograms(self):

        plt.close()
        gene_map = self.gene_map
        histogram_list = []
        setup = True
        for gene in gene_map.values():
            for index, value in enumerate(gene):
                if setup:
                    histogram_list.append([value])
                else:
                    histogram_list[index].append(value)
            setup = False

        figurelist = []
        sublist = []
        counter = 0
        for sample in histogram_list:
            counter += 1
            sublist.append(sample)
            if counter == 4:
                figurelist.append(sublist)
                counter = 0
                sublist = []

        counter = 0
        sublist = []
        figure_labels = []
        for name in self.data_frame_header.values()[0]:
            counter += 1
            sublist.append(name)
            if counter == 4:
                figure_labels.append(sublist)
                counter = 0
                sublist = []

        filecnt = 1
        fig_pos = 0
        for figure in figurelist:

            n_bins = 100
            # color = 'black'
            rang = tuple([float(x) for x in self.args.hist_range.split(',')])

            fig, axes = plt.subplots(nrows=2, ncols=2)
            ax0, ax1, ax2, ax3 = axes.flatten()

            ax0.hist([float(x) for x in figure[0]], n_bins, normed=100, range=rang)
            ax0.set_title(figure_labels[fig_pos][0])

            if len(figure) > 1:
                ax1.hist([float(x) for x in figure[1]], n_bins, normed=100, range=rang)
                ax1.set_title(figure_labels[fig_pos][1])

            if len(figure) > 2:
                ax2.hist([float(x) for x in figure[2]], n_bins, normed=100, range=rang)
                ax2.set_title(figure_labels[fig_pos][2])

            if len(figure) > 3:
                ax3.hist([float(x) for x in figure[3]], n_bins, normed=100, range=rang)
                ax3.set_title(figure_labels[fig_pos][3])

            fig.tight_layout()
            # plt.show()
            fig_pos += 1

            path = os.path.join('.', self.args.out, self.args.prefix)

            plt.savefig("{}_{}_{}.png".format(path,
                                              str(filecnt),
                                              'histogram_plots'),
                        format='png',
                        bbox_inches='tight')
            filecnt += 1
