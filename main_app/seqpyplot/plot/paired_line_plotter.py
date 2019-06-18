from __future__ import division

import logging
import os
import sys
from ast import literal_eval
from operator import concat

import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from seqpyplot.plot.base.plot_base import PlotBase
from tqdm import tqdm

if sys.platform == 'darwin':
    import matplotlib
    matplotlib.use('Agg')

try:
    from functools import reduce
except ImportError:
    pass


class PairedDataLinePlotter(PlotBase):

    def __init__(
        self,
        output_dir,
        normalized_df,
        complete_de_gene_list,
        log2fold,
        expression_min,
        time_point_names,
        condition_labels,
        experiment_name,
        min_diff,
        max_diff
    ):

        super(PairedDataLinePlotter, self).__init__()
        plt.close()

        self.output_dir = output_dir

        self.normalized_df = normalized_df
        self.normalized_df_genes = [x.lower() for x in self.normalized_df.index.tolist()]
        self.complete_de_gene_list = [x.lower() for x in complete_de_gene_list]

        self.log2fold = log2fold
        self.expression_min = expression_min
        self.time_point_names = time_point_names
        self.condition_labels = condition_labels
        self.experiment_name = experiment_name
        self.diffrange = [min_diff, max_diff]

    def plot_figure(self, figure_list, plottable_data):
        """
        Arguments:
            figure_list {object containing figure lists} -- 2d array
            plottable_data {[type]} -- df of gene expression values
        """

        line_plot_kwargs = {'marker': 'o',
                            'linewidth': 1.4,
                            'linestyle': "-",
                            'dash_capstyle': 'round',
                            'dash_joinstyle': 'bevel',
                            'fillstyle': 'full',
                            'markeredgewidth': 0.95,
                            'markersize': 6}

        pbar = tqdm(total=sum([len(x) for x in figure_list]))

        image_files = list()
        for fig_idx, figure in enumerate(figure_list, start=1):

            fig = self.set_figure(figure_prefix=self.experiment_name)

            line_plot_kwargs.update({'color': 'blue', 'markeredgecolor': 'blue', 'markerfacecolor': 'white'})
            series1 = self.set_line(kwargs=line_plot_kwargs)

            line_plot_kwargs.update({'color': 'red', 'markeredgecolor': 'red', 'markerfacecolor': 'white'})
            series2 = self.set_line(kwargs=line_plot_kwargs)

            legend_line_kwargs = {'color': 'white'}
            log_line = self.set_line(legend_line_kwargs)
            fold_line = self.set_line(legend_line_kwargs)
            diff_line = self.set_line(legend_line_kwargs)

            for idx, gene in enumerate(figure, start=1):

                fig.add_subplot(3, 2, idx).set_title('- ' + str(gene) + ' -', fontsize=20)

                if not self.gene_exists(gene):
                    logging.warning(
                        'Data for -- {} -- not found. This just means the gene is probably turned off.'.format(gene))
                    pbar.update(1)

                else:
                    is_de = self.is_de(gene)  # check if the gene is DE for bg highlighting
                    data = self.retrieve_data(gene)
                    self.create_subplot(is_de, *data)
                    pbar.update(1)

            handles = [series1, series2, log_line, fold_line, diff_line]
            fig = self.tidy_up_figure(fig, handles)

            file_name = self.save_figure(fig, fig_idx, figure_list)

            image_files.append(file_name)
            plt.clf()
            plt.cla()
            plt.close()

        pbar.close()

    def tidy_up_figure(self, fig, handles):

        diffstart = str(int(self.diffrange[0]))
        diffend = str(int(self.diffrange[1]))

        figlabels = [
            self.condition_labels,
            ["Log2: " + str(self.log2fold)],
            ["Fold: " + str(round(2.0 ** self.log2fold, ndigits=1))],
            ["Diff: " + " - ".join([diffstart, diffend])]
        ]
        labels = (reduce(concat, figlabels))
        fig.legend(handles=handles,
                   labels=labels,
                   loc='upper right')

        return fig

    def set_figure(self, figure_prefix, **args):
        fig = plt.figure(num=1,
                         dpi=600,
                         figsize=(10, 10),
                         edgecolor='black',
                         frameon=False
                         )
        fig.suptitle(figure_prefix,
                     verticalalignment='top',
                     horizontalalignment='right',
                     fontsize=24,
                     x=0.415,
                     y=1.09
                     )
        return fig

    def save_figure(self, fig, fig_idx, figure_list):

        path_ = os.path.join(self.output_dir, self.experiment_name)
        genes = "_".join([fi.strip() for fi in figure_list[fig_idx-1]])
        file_name = "{}_{}_{}.png".format(path_, str(fig_idx), genes)
        plt.savefig(str(file_name), format='png', bbox_inches='tight')
        print(file_name)
        return file_name

    def gene_exists(self, gene):
        return gene.lower() in self.normalized_df_genes

    def is_de(self, gene):
        return gene.lower() in self.complete_de_gene_list

    def retrieve_data(self, gene):
        """
        Return plot data and associated masks

        Arguments:
            plottable_data {[type]} -- [description]

        Returns:
            [type] -- [description]
        """
        data_length = self.normalized_df.shape[1]
        try:
            data = self.normalized_df.loc[gene.title()]

        except KeyError:
            print("The current gene: -- {} -- was not found in the plot data.".format(gene))
            data = [0] * data_length

        # These complicated comprehension is to handle cases where the user includes an
        # empty data file (to keep the data paired)
        series1 = np.array(
            [float(x) if x is not None else None for x in data[:(len(data) // 2)]]).astype(
             np.double)
        series2 = np.array(
            [float(x) if x is not None else None for x in data[(len(data) // 2):]]).astype(
             np.double)
        s1mask = np.isfinite(series1)
        s2mask = np.isfinite(series2)

        return series1, s1mask, series2, s2mask

    def create_subplot_template(self, is_de, y_max):
        ax = plt.gca()  # origin of the plot

        ax.tick_params(axis='both',
                       direction='out',
                       which='major',
                       labelsize=14)

        num_stages = len(self.time_point_names)

        # xlim is the length of the label list
        ax.set_xlim(-0.85, (num_stages)-0.15)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data', 0))

        ax.set_xticks([i for i in range(num_stages)], self.time_point_names)

        if is_de:  # set background color
            ax.set_facecolor('#FFFC97')

        ax.axhline(self.expression_min, color='gray', linestyle='--', label="-low")

        if y_max < self.expression_min:
            ax.set_ylim(0, self.expression_min + 10)
        else:
            ax.set_ylim(0, y_max)

        return ax

    def create_subplot(self, is_de, series1_data, series1_mask, series2_data, series2_mask):
        """
        Function for plotting individual plots in the figure.

        Arguments:
            gene {[type]} -- [description]
            is_de {bool} -- [description]
        """

        x_axis = np.arange(float(max(len(series1_data), len(series2_data))))

        data_mean, _ = self.compute_mean(series1_data, series2_data)

        diffs = self.compute_bounds(data_mean)
        y_max = self.compute_max_yval(series1_data, series2_data, diffs)
        ax = self.create_subplot_template(is_de, y_max)

        plot_kwargs = {'marker': 'o',
                       'linewidth': 1.4,
                       'linestyle': "-",
                       'dash_capstyle': 'round',
                       'dash_joinstyle': 'bevel',
                       'fillstyle': 'full',
                       'markeredgewidth': 0.95,
                       'markersize': 4}

        # Layer up plot elements
        plot_kwargs.update({'color': 'blue', 'markeredgecolor': 'blue', 'markerfacecolor': 'white', 'label': self.condition_labels[0]}) # .update returns NoneType
        ax = self.plot_series(series1_data, series1_mask, ax, x_axis, plot_kwargs)

        plot_kwargs.update({'color': 'red', 'markeredgecolor': 'red', 'markerfacecolor': 'white', 'label': self.condition_labels[1]})
        ax = self.plot_series(series2_data, series2_mask, ax, x_axis, plot_kwargs)

        _ = self.plot_log_fold_bounds(ax, data_mean, diffs, x_axis)  # May need the mean mask here

        self.tidy_up_plot()

        return ax

    def compute_mean(self, series1, series2):
        """Return a numpy array and mask that contains the average of the series1 and series2 data.
            This array is built for each gene to be plotted."""

        data_temp_mean = []

        for itr in range(max(len(series1), len(series2))):
            try:
                current_mean = np.mean([series1[itr], series2[itr]])
            except TypeError:
                current_mean = None
            data_temp_mean.append(current_mean)

        data_mean = np.array([float(x) if x is not None else None for x in data_temp_mean]).astype(
                             np.double)

        mean_mask = np.isfinite(data_mean)

        return data_mean, mean_mask

    def compute_bounds(self, series_mean):

        diffs = list()
        for i in series_mean:
            b = (2.0 * i) / ((2.0 ** self.log2fold) + 1)
            dif = i - b
            diffs.append(float(dif))

        diffs = np.asarray(diffs)

        if min(diffs) < 0 or min(series_mean - diffs) < 0:
            print("Error bars out of range - but still continuing.")

        return diffs

    def plot_log_fold_bounds(self, ax, series_mean, diffs, x_axis):

        ax = ax.errorbar(x_axis,
                         series_mean,
                         linestyle='None',
                         mfc='black',
                         yerr=[diffs, diffs],
                         ecolor='black',
                         elinewidth=1,
                         label='Range around the mean = log2(range)=1')
        return ax

    def compute_max_yval(self, series1_data, series2_data, diffs, scale=1.4):
        return max(np.mean(np.mean(np.asarray([series1_data, series2_data]), axis=0)) + diffs) * float(scale)

    def plot_series(self, series_data, series_mask, ax, xs, kwargs):
        ax.plot(xs[series_mask],
                series_data[series_mask],
                'bo',
                **kwargs)
        return ax

    def tidy_up_plot(self):

        plt.subplots_adjust(top=0.85)
        plt.tight_layout()
        return plt
