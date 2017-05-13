from __future__ import division
import os
import sys

try:
    import matplotlib.lines as mlines
    import matplotlib.pyplot as plt
    import numpy as np
except ImportError:
    print("Check Dependencies. You need numpy and MatPlotLib")
    sys.exit()
try:
    import seaborn as sns
    sns.set_palette('Blues', n_colors=256)
except ImportError:
    print("You should install Seaborn! For now these plots are gonna look awful!")
    sys.exit()

class MainDataPlotter(object):
    """
        For  building a all of the figures

    """
    def __init__(self, args, analyzer, figurelist):
        """
        :param args: Args Object
        :param analyzer: Analyzer Object
        :param figurelist: Figurelist var from DataContainer.MakeFigureList
        """
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

        # num == 1
        self.sing_comp_header = analyzer.sing_comp_header

        if self.args.out is not None:
            self.path = os.path.join(self.args.out, self.args.prefix)
        else:
            self.path = os.path.join('.', self.args.prefix)

    def plot_figures(self):
        """
        Figure list is a 2D nested list.
        Outerlist = Figure list.
        Innerlist = gene lists for each figure.

        :return: Line plots
        """

        de_set = set(self.de_genes)
        figure_list_count = 1
        for figure in self.figure_list:

            if len(figure) == 1:
                fig_size = (10, 10)
            else:
                fig_size = (10, 10)

            fig = plt.figure(num=1,
                             dpi=600,
                             figsize=fig_size,
                             edgecolor='black',
                             frameon=False,
                             )
            fig.suptitle(self.args.prefix,
                         verticalalignment='top',
                         horizontalalignment='right',
                         fontsize=24
                         # x=0.415
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
            fold_line = mlines.Line2D([], [], color='white')
            dif_line = mlines.Line2D([], [], color='white')


            gene_count = 1
            for gene in figure:
                print "Plotting: ", gene
                if len(figure) == 1:
                    fig.add_subplot(1, 1, gene_count).set_title('- ' + str(gene) + ' -', fontsize=20)
                else:
                    fig.add_subplot(3, 2, gene_count).set_title('- ' + str(gene) + ' -', fontsize=20)

                self.__subplot(gene, de_set)

                gene_count += 1

            try:
                condition_labels = [x for x in self.args.condition]
            except AttributeError:
                print("The condition labels have to be comma separated.")
                sys.exit()

            if self.args.num == 1:
                handles = [series1_line]
                labels = condition_labels
            elif self.args.num == 2:
                handles = [series1_line, series2_line, log_line, fold_line, dif_line]
                labels = (condition_labels
                          + ["Log2: " + str(self.args.log)]
                          + ["Fold: " + str(round(2.0**self.args.log, ndigits=1))]
                          + ["Diff: " + str(int(self.args.dif_range[0])) + ', ' + str(int(self.args.dif_range[1]))])
            else:
                print "reset -num argument"
                sys.exit()
            fig.legend(handles=handles,
                       labels=labels,
                       loc='upper right')

            path = os.path.join(self.args.out,
                                self.args.prefix)

            plt.savefig("{}_{}_{}.png".format(path,
                                             str(figure_list_count).replace(" ", ""),
                                             str([str(fig.replace(" ", "")) for fig in figure])),
                        format='png', bbox_inches='tight')

            if self.args.svg:
                plt.savefig("{}_{}_{}.png".format(path,
                                                  str(figure_list_count).replace(" ", ""),
                                                  str(figure).replace(" ", "")),
                            format='svg',
                            bbox_inches='tight')

            # TODO add argument to show plot windows
            # plt.show()
            plt.clf()
            plt.cla()
            plt.close()

            figure_list_count += 1

    def __subplot(self, gene, deset):
        """
        HIDDEN - don't worry this is handled internally
        :param gene:
        :param deset:
        :return:
        """
        found, series1_data, series1_mask, series2_data, series2_mask = self.data_collector(gene)

        if series2_data is None:
            np.asarray(series2_mask, dtype=bool)
            x_axis = np.arange(float(len(series1_data)))
        else:
            np.asarray(series1_mask, dtype=bool)
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
            plt.xlim(-0.85, (len(series1_data))-0.15)

            # Plot Details
            ax = plt.gca()

            is_de = self.__check_for_de(gene, deset)
            if is_de:
                ax.set_axis_bgcolor('#FFFC97')
            else:
                pass

            ax.spines['right'].set_color('none')
            ax.spines['top'].set_color('none')
            ax.xaxis.set_ticks_position('bottom')
            ax.yaxis.set_ticks_position('left')
            ax.spines['bottom'].set_position(('data', 0))

            labels = self.args.time

            plt.xticks([i for i in range(len(series1_data))], labels)

            # fill the plot with data points and graph them
            self.__series1_plot(series1_data, series1_mask, x_axis, self.args.condition[0], self.args.num)

            if self.args.num == 2:
                self.__series2_plot(series2_data, series2_mask, x_axis, self.args.condition[1])

                # error bars
                # return mean, and mask
                series_mean, mask = __calculate_mean__(series1_data, series2_data)

                err, upper = self.__err_range_plot(series_mean, x_axis, label="Range around the mean = log2(range)=1 ")

            try:
                if self.args.num == 2:
                    y_maximum = float(max(max(series1_data), max(series2_data), max(upper))) * 1.4
                else:
                    y_maximum = float(max(series1_data)) * 1.3

            except ValueError:
                print("Hmm, {0} wasn't found in the plot-able data set.".format(gene))
                y_maximum = 10

            # Extend subplot frame for high fpkm values
            if y_maximum < self.args.low:
                plt.ylim(0, self.args.low + 10)
            else:
                plt.ylim(0, y_maximum)

            plt.axhline(self.args.low, color='gray', linestyle='--', label="-low")
            # plt.legend()
            plt.tight_layout()

            # adjust subplot location
            plt.subplots_adjust(top=0.85)

            return plt

        except UnboundLocalError:
            print("Don't worry - this just means the gene is probably turned off.")

    def data_collector(self, gene_name):
        """
        Another interally used method (should be hidden...)
        :param gene_name: titled gene name (str.Title())
        :return: formatted list of expression data
        """
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

    def __check_for_de(self, gene, de__set):
        if str(gene) in de__set:
            return True
        else:
            return False

    def _calc_error(self, series_mean):
        """
        **HIDDEN**
        Use some linear algebra to calculate log2fold range around a given value based on user logfold parameter
        :param series_mean:
        :return:
        """
        var = self.args.log
        upper_list = []
        lower_list = []
        low_plot = []
        for i in series_mean:
            b = (2.0 * i) / ((2.0 ** var) + 1)
            dif = i - b
            a = i + dif

            #matplotlib requires the DIFFERENCE between the error value and the point
            upper_list.append(dif)
            lower_list.append(dif)
            low_plot.append(b)

        upper = np.asarray([float(c) for c in upper_list])
        lower = np.asarray([float(d) for d in lower_list])
        low_plot = np.asarray([float(e) for e in low_plot])

        if min(upper) < 0 or min(lower) < 0 or min(low_plot) < 0:
            print("Error bars out of range - but still continuing.")

        return upper, lower, low_plot

    def __err_range_plot(self, series_mean, xs, label):

        upper, lower, tempx = self._calc_error(series_mean)

        # xs is x_axis
        return plt.errorbar(xs,
                            series_mean,
                            linestyle='None',
                            mfc='black',
                            yerr=[upper, lower],
                            ecolor='black',
                            elinewidth='1',
                            label=label), upper

    @staticmethod
    def __series1_plot(series1, s1mask, xs, label, num):
        if num == 1:
            width = 1.8
            return plt.plot(xs,
                            series1,
                            'bo',
                            color="blue",
                            marker='o',
                            linewidth=width,
                            linestyle="-",
                            dash_capstyle='round',
                            dash_joinstyle='bevel',
                            label=label,
                            fillstyle='full',
                            markeredgecolor='blue',
                            markerfacecolor='white',
                            markeredgewidth=.95,
                            markersize=4)

        else:
            width = 1.4
            return plt.plot(xs[s1mask],
                            series1[s1mask],
                            'bo',
                            color="blue",
                            marker='o',
                            linewidth=width,
                            linestyle="-",
                            dash_capstyle='round',
                            dash_joinstyle='bevel',
                            label=label,
                            fillstyle='full',
                            markeredgecolor='blue',
                            markerfacecolor='white',
                            markeredgewidth=.95,
                            markersize=4)

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
                        markersize=4)

    def de_bar(self, colour):

        plt.close()
        fig = plt.figure(1, figsize=(7, 7))
        xlabel = []
        y_values = []
        bar_width = 0.5  # the width of the bars

        ax = plt.subplot()

        if self.args.num == 1:
            for k, v in sorted(self.de_count_by_time.items()):
                y_values.append(int(v))
                xlabel.append(k)
            x_axis = xrange(len(y_values))
            ax.bar(x_axis, y_values, bar_width, color=colour, align="center")
            ax.set_xticks(x_axis)
            ax.set_xticklabels(xlabel, rotation=45)

        elif self.args.num == 2:
            for k in self.args.time:
                y_values.append(int(self.de_count_by_time[k]))
                xlabel.append(k)

            x_axis = range(len(y_values))
            ax.set_xticks(x_axis)
            ax.bar(x_axis, y_values, bar_width, color=colour, align="center")
            ax.set_xticklabels(xlabel, rotation='vertical')
        else:
            sys.exit()
        try:
            ymax = max(y_values) * 1.3
            plt.ylim([0, ymax])

        except ValueError:
            print("Didn't find any DE genes.")
            print("These plots will be empty")

        # add some text for labels, title and axes ticks
        ax.set_ylabel('No. of Flagged Genes')
        ax.set_xlabel('Experimental Stage')
        ax.set_title('Number of flagged genes per stage', loc='left')
        ax.yaxis.grid()

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)

        if len(xlabel) > 6:
            ax.tick_params(axis='both', which='major', labelsize=8)
        else:
            plt.tick_params(
                axis='both',  # changes apply to the x-axis
                which='both',  # both major and minor ticks are affected
                bottom='on',  # ticks along the bottom edge are off
                top='off',  # ticks along the top edge are off
                left='on',
                right='off',
                labelbottom='on')  # labels along the bottom edge are off

        log_line = mlines.Line2D([], [], color='white')
        expression_upper = mlines.Line2D([], [], color='white')
        expression_lower = mlines.Line2D([], [], color='white')
        difference = mlines.Line2D([], [], color='white')
        if self.args.hi > 99999:
            hi = 'inf'
        else:
            hi = self.args.hi

        fig.legend(handles=[expression_upper, expression_lower, log_line, difference],
                   labels=(["Log2: " + str(self.args.log), "Range: " + str(self.args.low) + "-" + str(hi),
                            "Diff: {}, {}".format(str(int(self.args.dif_range[0])), str(int(self.args.dif_range[1])))]),
                   loc='upper right')

        plt.savefig(self.path + "{}_DE_Gene_by_time.png".format(self.args.prefix), format='png', bbox_inches='tight')
        plt.close()
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
            text = "{} of {}.".format(current, iteration)
            print '{:^43}'.format(text)
            self.analyzer.de_gene_list_length = 0
            current += 1

        # open a new line plot
        self.de_line_plot(cutoffs, line_plot_y_list, "DE cutoffs by log2fold")

    def de_line_plot(self, cutoffs, y_value, label):
        plt.close()

        # create a figure for the current figure
        fig = plt.figure(num=1,
                         figsize=(7, 7),
                         dpi=600,
                         edgecolor='black',
                         frameon=False,
                         )

        # set figure tit
        fig.suptitle("DE genes Detected vs log2Fold cutoff",
                     verticalalignment='top',
                     horizontalalignment='center',
                     fontsize=12,
                     x=0.315
                     )

        expression_upper = mlines.Line2D([], [], color='white')
        expression_lower = mlines.Line2D([], [], color='white')
        expression_dif = mlines.Line2D([], [], color='white')

        fig.legend(handles=[expression_upper, expression_lower, expression_dif],
                   labels=(["Upper: " + str(self.args.hi),
                            "Lower: " + str(self.args.low),
                            "Dif: " + str(self.args.dif_range[0]) + ', ' + str(self.args.dif_range[1])]),
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
        plt.close()

    def plot_histograms(self):

        # plt.close()
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
        if len(sublist) != 0:
            figurelist.append(sublist)

        counter = 0
        sublist = []
        figure_labels = []

        for name in self.data_frame_header["Gene"]:
            counter += 1
            sublist.append(name)
            if counter == 4:
                figure_labels.append(sublist)
                counter = 0
                sublist = []
        if len(sublist) != 0:
            figure_labels.append(sublist)

        filecnt = 1
        fig_pos = 0

        for figure in figurelist:

            n_bins = 100
            color = 'black'
            rang = tuple([float(x) for x in self.args.hist_range.split(',')])

            fig, axes = plt.subplots(nrows=2, ncols=2)
            fig.suptitle("Gene Density Per Count",
                         verticalalignment='top',
                         horizontalalignment='right',
                         fontsize=12,
                         y=1.05
                         )

            expression_upper = mlines.Line2D([], [], color='white')
            fig.legend(handles=[expression_upper],
                       labels=(["Range: " + str(self.args.hist_range)]),
                       loc='upper right')

            ax0, ax1, ax2, ax3 = axes.flatten()

            ax0.hist([float(x) for x in figure[0] if x is not None], n_bins, color=color, range=rang)
            ax0.set_title(figure_labels[fig_pos][0])

            if len(figure) > 1:
                ax1.hist([float(x) for x in figure[1] if x is not None], n_bins, color=color, range=rang)
                ax1.set_title(figure_labels[fig_pos][1])

            if len(figure) > 2:
                ax2.hist([float(x) for x in figure[2] if x is not None], n_bins, color=color, range=rang)
                ax2.set_title(figure_labels[fig_pos][2])

            if len(figure) > 3:
                ax3.hist([float(x) for x in figure[3] if x is not None], n_bins, color=color, range=rang)
                ax3.set_title(figure_labels[fig_pos][3])

            lower = int(rang[0])
            upper = int(rang[1])
            difference = int(upper - lower)
            inc = int(difference // 6)
            xlabels = [lower // 100,
                       lower + inc * 1,
                       lower + inc * 2,
                       lower + inc * 3,
                       lower + inc * 4,
                       lower + inc * 5,
                       upper]

            for ax in axes.flatten():
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
                ax.set_xticklabels(xlabels)
                ax.set_xlabel("Count", fontsize=8)
                ax.set_ylabel("No. of Geens", fontsize=8)

            fig.tight_layout()
            # plt.show()

            path = os.path.join('.', self.args.out, self.args.prefix)

            plt.savefig("{}_{}_{}.png".format(path,
                                              str(filecnt),
                                              'histogram_plots'),
                        format='png',
                        bbox_inches='tight')
            filecnt += 1
            fig_pos += 1
            plt.close()
            if fig_pos == len(figure_labels):
                return

    def make_scatter_plots(self, flagged=False):

        if flagged is False:
            if self.args.num == 1:
                gene_map = self.analyzer.sing_time_series_data

            elif self.args.num == 2:
                gene_map = self.analyzer.unflagged_genes
            else:
                sys.exit()
            plot_name = 'Scatter_Plots_unflagged_genes'

        else:
            if self.args.num == 1:
                gene_map = self.analyzer.sing_time_series_data

            elif self.args.num == 2:
                gene_map = self.analyzer.filtered_data
            else:
                sys.exit()
            plot_name = 'Scatter_Plots_flagged_genes'

        scatter_dict = dict()
        # print gene_map.items()[0][1]
        width = len(gene_map.items()[0][1])
        spacer = width // 2

        sublist = []
        counter = 0

        plot_labels = []
        figure_labels = []

        if self.args.num == 1:
            for name in self.sing_comp_header:
                counter += 1
                sublist.append(name)
                if counter == 4:
                    figure_labels.append(sublist)
                    counter = 0
                    sublist = []

            if len(sublist) != 0:
                figure_labels.append(sublist)

            for i in self.sing_comp_header:
                sublist = i.split('/')
                plot_labels.append(sublist)

            times = self.sing_comp_header

        if self.args.num == 2:

            for name in self.args.time:
                counter += 1
                sublist.append(name)
                if counter == 4:
                    figure_labels.append(sublist)
                    counter = 0
                    sublist = []
            if len(sublist) != 0:
                figure_labels.append(sublist)

            sublist = []
            name = self.data_frame_header["Gene"]

            for i in range(spacer):
                sublist.append(name[i])
                sublist.append(name[i + spacer])
                plot_labels.append(sublist)
                sublist = []
            times = self.args.time

        # make scatter dictionary
        setup1 = True
        while setup1:
            for gene in gene_map.values():
                for i in range(spacer):
                    scatter_dict[times[i]] = [[], []]
                setup1 = False
                break

        scatrang = tuple([float(x) for x in self.args.scatt_range.split(',')])
        if flagged:
            for plot in enumerate(times):

                timeidx = plot[0]
                timekey = plot[1]
                try:
                    for gene in self.analyzer.de_gene_list_by_stage[timekey]:
                        scatter_dict[timekey][0].append(self.analyzer.gene_map[gene][timeidx])
                        if self.args.num == 2:
                            scatter_dict[timekey][1].append(self.analyzer.gene_map[gene][timeidx+spacer])
                except KeyError:
                    pass
        else:
            for plot in enumerate(times):

                timeidx = plot[0]
                timekey = plot[1]

                try:
                    for gene in self.analyzer.gene_map.keys():

                        if gene not in self.analyzer.de_gene_list_by_stage[timekey]:
                            scatter_dict[timekey][0].append(self.analyzer.gene_map[gene][timeidx])
                            if self.args.num == 2:
                                scatter_dict[timekey][1].append(self.analyzer.gene_map[gene][timeidx + spacer])
                except KeyError:
                    pass

        filecnt = 1
        fig_pos = 0

        # bounds
        y_range = range(int(scatrang[1]))
        # print y_range
        up, temp, lowerbound = self._calc_error(y_range)
        upperbound = []
        for i in range(len(y_range)):
            upperbound.append(y_range[i]+up[i])

        namecount = 0
        for figure in figure_labels:

            fig, axes = plt.subplots(nrows=2, ncols=2)
            fig = plt.figure(num=1,
                             dpi=1200,
                             figsize=(70, 70),
                             edgecolor='black',
                             frameon=False,
                             )
            if flagged:
                suptitle = "Flagged Genes"
            else:
                suptitle = "Unflagged Genes"

            fig.suptitle(suptitle,
                         verticalalignment='top',
                         horizontalalignment='right',
                         fontsize=14,
                         y=1.09
                         # x=0.4
                         )

            L1 = mlines.Line2D([], [], color='white')
            L2 = mlines.Line2D([], [], color='white')
            L3 = mlines.Line2D([], [], color='white')

            fig.legend(handles=[L1,L2, L3],
                       labels=(["ScatRange: " + str(self.args.scatt_range),
                                "DiffRange: {}, {} ".format(str(self.args.dif_range[0]), str(self.args.dif_range[1])),
                                "Log2: " + str(self.args.log)]),
                       loc='upper right')

            i = 0
            for ax in axes.flatten():
                try:
                    x = scatter_dict[figure[i]][0]
                    y = scatter_dict[figure[i]][1]
                    z = [np.mean([float(x[index]), float(y[index])]) for index in range(len(x))]
                except ValueError:
                    print "Double check the -time argument. Does it match?"
                    sys.exit()

                ax.scatter(z, x,
                           s=8,

                           alpha=0.5
                           )
                ax.scatter(z, y,
                           s=8,
                           # c=colors,
                           alpha=0.5)

                ax.plot(y_range,
                        lowerbound,
                        color='black',
                        linestyle='--')

                ax.plot(y_range,
                        upperbound,
                        color='black',
                        linestyle='--')

                ax.plot(range(int(scatrang[1])),
                        range(int(scatrang[1])))

                ax.axvline(self.args.low, color='gray', linestyle='--')
                ax.axvline(self.args.hi, color='gray', linestyle='--')
                ax.axvline(self.args.low, color='gray', linestyle='--')
                ax.axvline(self.args.hi, color='gray', linestyle='--')


                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
                ax.set_xlim(scatrang)
                ax.set_ylim(scatrang)
                ax.set_xlabel('Mean', fontsize=10)
                ax.set_ylabel('Count', fontsize=10)
                ax.set_title(figure[i], fontsize=12)
                ax.tick_params(axis='both', which='major', labelsize=8)
                i += 1
                namecount += 1
                if i == len(figure):
                    break

            fig.tight_layout()
            path = os.path.join('.', self.args.out, self.args.prefix)


            plt.savefig("{}_{}_{}.png".format(path,
                                              plot_name,
                                              str(filecnt)),
                        format='png',
                        bbox_inches='tight')

            filecnt += 1
            fig_pos += 1

            plt.close()
            del fig

    def bland_altman_plot(self):

        gene_map = self.analyzer.unflagged_genes
        if len(gene_map) == 0:
            print "No unflagged genes available!"
            return

        ba_dict = dict()
        width = len(gene_map.items()[0][1])
        spacer = width // 2

        sublist = []
        counter = 0

        plot_labels = []
        figure_labels = []

        for name in self.args.time:
            counter += 1
            sublist.append(name)
            if counter == 4:
                figure_labels.append(sublist)
                counter = 0
                sublist = []
        if len(sublist) != 0:
            figure_labels.append(sublist)

        sublist = []
        name = self.data_frame_header["Gene"]

        for i in range(spacer):
            sublist.append(name[i])
            sublist.append(name[i + spacer])
            plot_labels.append(sublist)
            sublist = []
        times = self.args.time

        # make BA dictionary
        setup1 = True
        while setup1:
            for _ in gene_map.values():
                for i in range(spacer):
                    ba_dict[times[i]] = [[], []]
                setup1 = False
                break

        # barange = tuple([float(x) for x in self.args.ba_range.split(',')])
        barange = float(self.args.ba_range)
        for plot in enumerate(times):

            timeidx = plot[0]
            timekey = plot[1]
            try:

                for gene in self.analyzer.gene_map.keys():

                    if gene not in self.analyzer.de_gene_list_by_stage[timekey]:
                        ba_dict[timekey][0].append(self.analyzer.gene_map[gene][timeidx])
                        ba_dict[timekey][1].append(self.analyzer.gene_map[gene][timeidx + spacer])
            except KeyError:
                pass

        filecnt = 1
        fig_pos = 0

        namecount = 0
        for figure in figure_labels:

            fig, axes = plt.subplots(nrows=2, ncols=2)
            fig = plt.figure(num=1,
                             dpi=1200,
                             figsize=(70, 70),
                             edgecolor='black',
                             frameon=False,
                             )

            fig.suptitle("Bland-Altman Plots - Unflagged Genes",
                         verticalalignment='top',
                         horizontalalignment='center',
                         fontsize=14,
                         y=1.05
                         )

            i = 0
            for ax in axes.flatten():
                try:
                    x = ba_dict[figure[i]][0]
                    y = ba_dict[figure[i]][1]
                except ValueError:
                    print "Double check the -time argument. Does it match?"
                    sys.exit()

                self.BA_plot(x, y, ax)

                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
                ax.set_xlim(0, barange)
                ax.set_ylim(barange*-1.0, barange)
                ax.set_title(figure[i])
                ax.tick_params(axis='both', which='major', labelsize=8)
                i += 1
                namecount += 1
                if i == len(figure):
                    break

            fig.tight_layout()
            path = os.path.join('.', self.args.out, self.args.prefix)

            plt.savefig("{}_{}_{}.png".format(path,
                                              'Bland_Altman_plots',
                                              str(filecnt)),
                        format='png',
                        bbox_inches='tight')

            filecnt += 1
            fig_pos += 1

            plt.close()
            del fig

    def BA_plot(self, data1, data2, ax):
        """ Generate a Bland-Altman plot.
        Arguments:
            :type s1: numpy.array
            :param s1: An array of sample1 data.
            :type s2: numpy.array
            :param s2: An array of sample2 data.
            :param ax:
            :type ax: matplotlib.axes.AxesSubplot
        Returns:
            :returns: If avaiable returns a matplotlib.figure.Figure else adds plot
                to current axis.
            :rtype: matplotlib.figure.Figure
        """

        # Make sure s1 and s2 are numpy arrays
        s1 = np.asarray(data1).astype(np.double)
        s2 = np.asarray(data2).astype(np.double)
        mean = []
        diff = []
        md = np.mean(diff)
        std = np.std(diff, axis=0)

        for counter in range(len(data1)):
            mean.append((float(data1[counter]) + float(data2[counter]))/2.0)
            diff.append(float(data1[counter]) - float(data2[counter]))

        ax.scatter(mean, diff)
        ax.set_xlabel('Mean')
        ax.set_ylabel('Difference')

        ax.axhline(md, color='r', ls='--', lw=2)
        ax.axhline(md + 1.96 * std, color='gray', linestyle='--')
        ax.axhline(md - 1.96 * std, color='gray', linestyle='--')

        return ax

    def collective_log_plot(self):

        # sort data from decreasing to increaseing
        # Get range len(log_list)
        # plot all of the log list against the range. Boom
        logs = []
        for l in self.analyzer.foldchange_map.values():
            for o in l:
                if isinstance(o, int) or isinstance(o, float):
                    if abs(float(o)) < self.args.log:
                        pass
                    else:
                        logs += (l)
        loglist = []
        for val in logs:
            if isinstance(val, str):
                pass
            else:
                if val != 0.0:
                    loglist.append(float(val))

        n_bins = 100
        color = 'black'

        fig, ax = plt.subplots(nrows=1, ncols=1)

        ax.hist(loglist, n_bins, color=color)
        ax.set_title("General Log2Fold Change Range")

        ax.set_xlim(min(loglist), max(loglist))

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set_xlabel("Log2Fold Value", fontsize=8)
        ax.set_ylabel("No. of Genes", fontsize=8)
        ax.set_xlim([-6, 6])

        fig.tight_layout()
        path = os.path.join('.', self.args.out, self.args.prefix)
        plt.savefig("{}_{}.png".format(path,
                                          'log2_histogram'),
                    format='png',
                    bbox_inches='tight')

        plt.close()
        del fig

    def single_log_plots(self):

        data_list = self.analyzer.foldchange_map
        histogram_list = dict()
        setup = True
        for time in enumerate(self.args.time):
            timeidx = time[0]
            timekey = time[1]

            for gene in data_list.values():
                if setup is True:
                    histogram_list[timekey] = []
                    setup = False
                else:
                    if isinstance(gene[timeidx], int) or isinstance(gene[timeidx], float):
                        if abs(float(gene[timeidx])) < self.args.log:
                            pass
                        else:
                            histogram_list[timekey].append(gene[timeidx])
            setup = True

        counter = 0
        sublist = []
        figure_labels = []
        for name in self.args.time:
            counter += 1
            sublist.append(name)
            if counter == 4:
                figure_labels.append(sublist)
                counter = 0
                sublist = []
        if len(sublist) != 0:
            figure_labels.append(sublist)

        filecnt = 1
        fig_pos = 0

        for figure in figure_labels:

            n_bins = 100
            color = 'black'

            fig, axes = plt.subplots(nrows=2, ncols=2)
            ax0, ax1, ax2, ax3 = axes.flatten()
            try:
                ax0.hist(sorted(histogram_list[figure[0]]),
                         n_bins,
                         color=color,
                         range=(min(histogram_list[figure[0]]),
                                max(histogram_list[figure[0]])))
            except ValueError:
                pass

            ax0.set_title(figure_labels[fig_pos][0])

            if len(figure) > 1:
                try:
                    ax1.hist(sorted(histogram_list[figure[1]]),
                             n_bins,
                             color=color,
                             range=(min(histogram_list[figure[1]]),
                                    max(histogram_list[figure[1]])))
                except ValueError:
                    pass
                ax1.set_title(figure_labels[fig_pos][1])

            if len(figure) > 2:
                try:
                    ax2.hist(sorted(histogram_list[figure[2]]),
                             n_bins,
                             color=color,
                             range=(min(histogram_list[figure[2]]),
                                    max(histogram_list[figure[2]])))
                except ValueError:
                    pass
                ax2.set_title(figure_labels[fig_pos][2])

            if len(figure) > 3:
                try:
                    ax3.hist(sorted(histogram_list[figure[3]]),
                             n_bins,
                             color=color,
                             range=(min(histogram_list[figure[3]]),
                                    max(histogram_list[figure[3]])))
                except ValueError:
                    pass

                ax3.set_title(figure_labels[fig_pos][3])

            for ax in axes.flatten():
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
                ax.set_xlabel("Log2Fold", fontsize=8)
                ax.set_ylabel("No. of Geens", fontsize=8)
                ax.set_xlim([-6,6])
            fig.tight_layout()
            path = os.path.join('.', self.args.out, self.args.prefix)
            plt.savefig("{}_{}_{}.png".format(path,
                                              str(filecnt),
                                              'sample_log2fold_histograms'),
                        format='png',
                        bbox_inches='tight')
            filecnt += 1
            fig_pos += 1
            plt.close()
            if fig_pos == len(figure_labels):
                return

    def bland_alt_log2_plot(self, flagged=False):
        if flagged:
            suptitle = "Bland-Altman-Log2 Plots - Flagged Genes"
            gene_map = self.analyzer.filtered_data
        else:
            suptitle = "Bland-Altman-Log2 Plots - Unflagged Genes"
            gene_map = self.gene_map

        bl_dict = dict()
        width = len(gene_map.items()[0][1])
        spacer = width // 2

        sublist = []
        counter = 0

        plot_labels = []
        figure_labels = []

        for name in self.args.time:
            counter += 1
            sublist.append(name)
            if counter == 4:
                figure_labels.append(sublist)
                counter = 0
                sublist = []
        if len(sublist) != 0:
            figure_labels.append(sublist)

        sublist = []
        name = self.data_frame_header["Gene"]

        for i in range(spacer):
            sublist.append(name[i])
            sublist.append(name[i + spacer])
            plot_labels.append(sublist)
            sublist = []
        times = self.args.time

        # make BA dictionary
        setup1 = True
        while setup1:
            for _ in gene_map.values():
                for i in range(spacer):
                    bl_dict[times[i]] = [[], []]
                setup1 = False
                break

        blrange = tuple([float(x) for x in self.args.bl_range.split(',')])

        for plot in enumerate(times):

            timeidx = plot[0]
            timekey = plot[1]

            if flagged:
                try:
                    for gene in self.gene_map.keys():
                        if gene in self.analyzer.de_gene_list_by_stage[timekey]:
                            bl_dict[timekey][0].append(self.gene_map[gene][timeidx])
                            bl_dict[timekey][1].append(self.gene_map[gene][timeidx + spacer])
                except KeyError:
                    pass
            else:
                try:
                    for gene in self.gene_map.keys():
                        if gene not in self.analyzer.de_gene_list_by_stage[timekey]:
                            bl_dict[timekey][0].append(self.gene_map[gene][timeidx])
                            bl_dict[timekey][1].append(self.gene_map[gene][timeidx + spacer])
                except KeyError:
                    pass

        filecnt = 1
        fig_pos = 0

        namecount = 0
        for figure in figure_labels:

            fig, axes = plt.subplots(nrows=2, ncols=2)
            fig = plt.figure(num=1,
                             dpi=1200,
                             figsize=(70, 70),
                             edgecolor='black',
                             frameon=False,
                             )

            fig.suptitle(suptitle,
                         verticalalignment='top',
                         horizontalalignment='center',
                         fontsize=14,
                         y=1.05
                         )

            i = 0
            for ax in axes.flatten():
                try:
                    x = bl_dict[figure[i]][0]
                    y = bl_dict[figure[i]][1]
                except ValueError:
                    print "Double check the -time argument. Does it match?"
                    sys.exit()

                self.BLOG_plot(x, y, ax)

                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.get_xaxis().tick_bottom()
                ax.get_yaxis().tick_left()
                ax.set_ylim(blrange)
                ax.set_xlim([-4, 4])
                ax.set_title(figure[i])

                ax.tick_params(axis='both', which='major', labelsize=8)
                i += 1
                namecount += 1
                if i == len(figure):
                    break

            fig.tight_layout()
            path = os.path.join('.', self.args.out, self.args.prefix)

            if flagged:
                plt.savefig("{}_{}_{}.png".format(path,
                                                  'Bland_Alt_log2_plots_flagged_genes',
                                                  str(filecnt)),
                            format='png',
                            bbox_inches='tight')
            else:
                plt.savefig("{}_{}_{}.png".format(path,
                                                  'Bland_Alt_log2_plots_unflagged_genes',
                                                  str(filecnt)),
                            format='png',
                            bbox_inches='tight')


            filecnt += 1
            fig_pos += 1
            plt.close()

    def BLOG_plot(self, data1, data2, ax):
        """ Generate a Bland-Altman plot.
        Arguments:
            :type s1: numpy.array
            :param s1: An array of sample1 data.
            :type s2: numpy.array
            :param s2: An array of sample2 data.
            :param ax:
            :type ax: matplotlib.axes.AxesSubplot
        Returns:
            :returns: If avaiable returns a matplotlib.figure.Figure else adds plot
                to current axis.
            :rtype: matplotlib.figure.Figure
        """

        # Make sure s1 and s2 are numpy arrays
        s1 = np.asarray(data1).astype(np.double)
        s2 = np.asarray(data2).astype(np.double)
        mean = []
        log2fold = []

        for counter in range(len(data1)):
            if np.mean([float(data1[counter]), float(data2[counter])]) < self.args.low:
                pass
            else:

                mean.append(np.log2(np.mean([float(data1[counter]), float(data2[counter])])))
                log2fold.append(np.log2(float(data1[counter]) / float(data2[counter])))

        sd = np.std(np.asarray(log2fold), axis=0)

        # Calculate mean and difference
        # mean = (s1 + s2) / 2
        # diff = abs(s1 - s2)
        assert len(mean) == len(log2fold)

        ax.scatter(log2fold, mean)
        ax.axvline(0, color='r', ls='--', lw=2)
        # ax.axvline(0.0 + 1.96 * sd, color='gray', linestyle='--')
        # ax.axvline(0.0 - 1.96 * sd, color='gray', linestyle='--')
        ax.set_ylabel('Log2(Mean)')
        ax.set_xlabel('Log2(Fold)')



        return ax

if __name__ == '__main__':
    pass