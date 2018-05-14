from seqpyplot.plot.base.plot_base import PlotBase
from ..analyzer.paired_sample_filter import PairedSampleFilter
import matplotlib.pyplot as plt
from matplotlib import rcParams
import matplotlib.lines as mlines
import os
import numpy as np


class TallyDe(PlotBase):

    def __init__(self, config_obj, container_obj):
        plt.close()
        super(TallyDe, self).__init__()

        self.cutoffs = [x/100. for x in range(200) if x % 5 == 0]
        
        self.container_obj = container_obj
        self.config_obj = config_obj

        self.output_dir = self.create_output_directory()
        self.prefix = self.config_obj.get('names', 'experiment_name')

        self.log2fold = self.config_obj.getfloat('params', 'log2fold')
        self.low = self.config_obj.getint('params', 'low')
        self.hi = self.config_obj.getint('params', 'hi')
        self.diff = self.config_obj.getlist('params', 'diff')

    def compute_tally(self, input_df_list):

        y_values = []
        iteration = len(self.cutoffs)
        print("Iterating over log2fold cutoff values... ")
        for idx, cutoff in enumerate(self.cutoffs):
            print("current cutoff: ", cutoff)
            analyzer = PairedSampleFilter(self.config_obj,
                                          log2fold=cutoff)
            _ = analyzer.main_filter_process(input_df_list)

            # print temp_de_count
            y_values.append(len(analyzer.complete_de_gene_list))
            text = "{} of {}.".format(idx + 1, len(self.cutoffs)+1)
            print('{:^43}'.format(text))

        return y_values

    def set_figure(self, handles):
        # create a figure for the current figure
        fig = plt.figure(num=1,
                         figsize=(7, 7),
                         dpi=1200,
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
            " ".join(['Upper:',  str(self.hi)]),
            " ".join(['Lower:', str(self.low)]),
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

        path_ = os.path.join(self.output_dir, self.prefix)
        file_name = "_".join([path_, "_DE_Tally_Plot.png"])
        fig.savefig(file_name, format='png', bbox_inches='tight')
        fig.close()

    def save_plot(self, fig):

        path_ = os.path.join(self.output_dir, self.prefix)
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