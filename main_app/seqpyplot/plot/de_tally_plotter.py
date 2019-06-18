import os

import matplotlib.lines as mlines
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import rcParams
import seaborn as sns

from seqpyplot.plot.base.plot_base import PlotBase
from tqdm import tqdm

from ..analyzer.paired_sample_filter import PairedSampleFilter



class TallyDe(PlotBase):

    def __init__(
        self,
        output_dir,
        container_obj,
        experiment_name,
        log2fold,
        expression_max,
        expression_min,
        min_diff,
        max_diff,
        file_name_pairs,
        time_point_names
    ):

        plt.close()
        super()

        self.cutoffs = [x/100. for x in range(200) if x % 5 == 0][1:]
        self.container_obj = container_obj

        self.output_dir = output_dir
        self.experiment_name = experiment_name

        self.log2fold = log2fold
        self.expression_min = expression_min
        self.expression_max = expression_max
        self.diff = [min_diff, max_diff]
        self.file_name_pairs = file_name_pairs
        self.time_point_names = time_point_names

    def compute_tally(self, input_df_list):
        y_values = []
        print("\nIterating over log2fold cutoff values... ")
        for idx, cutoff in tqdm(enumerate(self.cutoffs), total=len(self.cutoffs)):
            analyzer = PairedSampleFilter(
                log2fold=cutoff,
                expression_min=self.expression_min,
                expression_max=self.expression_max,
                min_diff=self.diff[0],
                max_diff=self.diff[1],
                time_point_names=self.time_point_names,
                file_name_pairs=self.file_name_pairs
                )
            _ = analyzer.main_filter_process(input_df_list)

            # print temp_de_count
            y_values.append(len(analyzer.complete_de_gene_list))
            text = "{} of {}.".format(idx + 1, len(self.cutoffs)+1)
            # print('{:^43}'.format(text))

        return y_values

    def set_figure(self, handles):
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
        labels = [
            " ".join(['Upper:',  str(self.expression_max)]),
            " ".join(['Lower:', str(self.expression_min)]),
            " ".join(['Dif:', str(self.diff[0]), '-', str(self.diff[1])])
                  ]

        fig.legend(handles=handles,
                    labels=(labels),
                    loc='upper right')
        rcParams['legend.frameon'] = 'False'

        return fig

    def create_subplot(self, ax, y_values):
        ax.plot(self.cutoffs,
                 y_values,
                 'bo',
                 color="black",
                 marker='o',
                 linewidth=1.4,
                 linestyle="-",
                 dash_capstyle='round',
                 dash_joinstyle='bevel',
                 label="DE cutoffs by log2fold",
                 fillstyle='full',
                 markeredgecolor='black',
                 markerfacecolor='white',
                 markeredgewidth=.95,
                 markersize=6)
        return ax

    def format_tally_plot(self, ax):

        # xlim is the length of the label list
        ax.set_xlim(0, 3)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data', 0))
        ax.set_ylabel('Number of DE genes Detected')
        ax.set_xlabel('log2Fold Cutoff')
        return ax

    def save_fig(self, fig):

        path_ = os.path.join(self.output_dir, self.experiment_name)
        file_name = "_".join([path_, "_DE_Tally_Plot.png"])
        fig.savefig(file_name, format='png', bbox_inches='tight')
        fig.close()

    def save_plot(self, fig):

        path_ = os.path.join(self.output_dir, self.experiment_name)
        file_name = "_".join([path_,"DE_Tally_Plot.png"])

        fig.savefig(file_name, format='png', bbox_inches='tight')

        return fig

    def create_tally_plot(self, input_df_list):
        handles = [self.set_line() for _ in range(3)]
        fig = self.set_figure(handles)
        ax = plt.subplot()

        y_values = self.compute_tally(input_df_list)
        ax = self.create_subplot(ax, y_values)
        ax = self.format_tally_plot(ax)
        self.save_plot(fig)
        plt.close()
