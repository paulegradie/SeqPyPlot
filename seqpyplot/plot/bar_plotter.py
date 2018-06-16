import os

import matplotlib.pyplot as plt
from seqpyplot.plot.base.plot_base import PlotBase

plt.style.use('bmh')


class PairedBarPlot(PlotBase):

    def __init__(self, config_obj):
        
        super(PairedBarPlot, self).__init__()

        plt.close()

        self.config_obj = config_obj

        self.labels = self.config_obj.getlist('names', 'times')
        self.output_dir = self.create_output_directory()
        self.prefix = self.config_obj.get('names', 'experiment_name')

        self.log = self.config_obj.getfloat('params', 'log2fold')
        self.hi = self.config_obj.getint('params', 'hi')
        self.low = self.config_obj.getint('params', 'low')
        self.diff = self.config_obj.getlist('params', 'diff')

    def save_plot(self, fig):

        path_ = os.path.join(self.output_dir, self.prefix)
        file_name = "_".join([path_,"DE_Gene_by_time.png"])

        fig.savefig(file_name, format='png', bbox_inches='tight')

        return fig

    def tidy_figure_up(self, fig, handles):

        hi = 'inf' if self.hi > 99999 else self.hi

        diffstart = str(int(self.diff[0]))
        diffend = str(int(self.diff[1]))
        fig.legend(handles=[self.set_line({'color': 'white'}) for _ in range(3)],
                   labels=(
                       ["Log2: " + str(self.log),
                        "Range: " + " - ".join([str(self.low), str(hi)]),
                        "Diff: " + " - ".join([diffstart, diffend])]
                        ),
                   loc='upper right')
        return fig


    def set_figure(self, figure_prefix, **args):
        fig = plt.figure(num=1,
                         figsize=(7, 7),
                         dpi=800,
                         edgecolor='black',
                         frameon=False,
                         )
        fig.suptitle(figure_prefix,
                     verticalalignment='top',
                     horizontalalignment='left',
                     fontsize=16,
                     x=0.108,
                     y=1.08
                     )
        return fig


    def create_subplot(self, x_axis, y_values, x_label, ymax, bar_width=0.5, color='black'):

        ax = plt.subplot()
        ax.bar(x_axis, y_values, bar_width, color=color, align="center")
        ax = self.format_subplot(ax, x_axis, x_label, ymax, num_stages=len(x_label))
        return ax

    def format_subplot(self, ax, x_axis, x_label, ymax, num_stages):

        ax.set_xticks(x_axis)
        ax.set_xticklabels(x_label, rotation='vertical')
        ax.set_ylim([0, ymax])

        # add some text for labels, title and axes ticks
        ax.set_ylabel('No. of Flagged Genes')
        ax.set_xlabel('Experimental Stage')
        ax.yaxis.grid()

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(True)
        ax.spines['left'].set_visible(True)

        if num_stages > 6:
            ax.tick_params(axis='both', which='major', labelsize=8)
        else:
            ax.tick_params(axis='both',  # changes apply to the x-axis
                           which='both',  # both major and minor ticks are affected
                           bottom=True,  # ticks along the bottom edge are off
                           top=False,  # ticks along the top edge are off
                           left=True,
                           right=False,
                           labelbottom=True)  # labels along the bottom edge are off

        return ax

    def create_bar_plot(self, de_count_by_time_dict):

        y_values = [int(de_count_by_time_dict[x]) for x in self.labels]
        x_label = self.labels
        ymax = max(y_values) * 1.2
  
        x_axis = range(len(y_values))

        fig = self.set_figure(figure_prefix='Paired Sample DE Bar Plot')
        self.create_subplot(x_axis, y_values, x_label, ymax)

        line_handles = [self.set_line() for _ in range(4)]
        fig = self.tidy_figure_up(fig, line_handles)
        plt.tight_layout()

        self.save_plot(fig)
