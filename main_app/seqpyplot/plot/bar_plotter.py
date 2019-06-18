import os

import matplotlib.pyplot as plt
from seqpyplot.plot.base.plot_base import PlotBase

plt.style.use('bmh')


class PairedBarPlot(PlotBase):

    def __init__(
        self,
        time_point_names,
        output_dir,
        experiment_name,
        log2fold,
        expression_max,
        expression_min,
        min_diff,
        max_diff
        ):

        super()

        plt.close()

        self.time_point_names = time_point_names
        self.output_dir = output_dir
        self.experiment_name = experiment_name
        self.log2fold = log2fold
        self.expression_max = expression_max
        self.expression_min = expression_min
        self.diff = [min_diff, max_diff]

    def save_plot(self, fig):

        path_ = os.path.join(self.output_dir, self.experiment_name)
        file_name = "_".join([path_,"DE_Gene_by_time.png"])

        fig.savefig(file_name, format='png', bbox_inches='tight')

        return fig

    def tidy_figure_up(self, fig, handles):

        # leading underscore to avoid confusion
        _hi = 'inf' if self.expression_max > 99999 else self.expression_max

        diffstart = str(int(self.diff[0]))
        diffend = str(int(self.diff[1]))
        fig.legend(
            handles=[self.set_line({'color': 'white'}) for _ in range(3)],
            labels=(
                [
                    "Log2: " + str(self.log2fold),
                    "Range: " + " - ".join([str(self.expression_min), str(_hi)]),
                    "Diff: " + " - ".join([diffstart, diffend])
                ]
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

        y_values = [int(de_count_by_time_dict[x]) for x in self.time_point_names]
        x_label = self.time_point_names
        ymax = max(y_values) * 1.2

        x_axis = range(len(y_values))

        fig = self.set_figure(figure_prefix='Paired Sample DE Bar Plot')
        self.create_subplot(x_axis, y_values, x_label, ymax)

        line_handles = [self.set_line() for _ in range(4)]
        fig = self.tidy_figure_up(fig, line_handles)
        plt.tight_layout()

        self.save_plot(fig)
