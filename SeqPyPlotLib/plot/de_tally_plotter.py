from base.plot_base import PlotBase
from ..analyzer.paired_sample_filter import PairedSampleFilter
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

class TallyDe(PlotBase):

    def __init__(self, config_obj, container_obj):

        self.cutoffs = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0,
                        1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0,
                        3, 4, 5]
        self.container_obj = container_obj
        self.config_obj = config_obj

        self.output_dir = self.create_output_directory()
        self.prefix = self.config_obj.get('names', 'experiment_name')

        self.log2fold = self.config_obj.getfloat('params', 'log2fold')
        self.low = self.config_obj.getint('params', 'low')
        self.hi = self.config_obj.getint('params', 'hi')
        self.diff = self.config_obj.getlist('params', 'diff')

    def compute_tally(self)

        line_plot_y_list = []
        iteration = len(cutoffs)
        print("Iterating over log2fold cutoff values... ")

        for idx, cutoff in enumerate(self.cutoffs):
            
            analyzer = PairedSampleFilter(self.config_obj, self.container_obj)
        
            # print temp_de_count
            line_plot_y_list.append(len(analyzer.complete_de_gene_list))
            text = "{} of {}.".format(idx + 1, len(self.cutoffs)+1)
            print('{:^43}'.format(text))

        return line_plot_y_list

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
            " ".join(['Upper:',  str(self.hi)]),
            " ".join(['Lower:', str(self.low)]),
            " ".join(['Dif:', str(self.diff[0]), '-', str(self.diff[1])])
                  ]

        fig.legend(handles=handles,
                    labels=(labels),
                    loc='upper right')
        return fig

    def set_line(self, color='white'):
        return mlines.Line2D([], [], color=color)

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
                 label=label,
                 fillstyle='full',
                 markeredgecolor='black',
                 markerfacecolor='white',
                 markeredgewidth=.95,
                 markersize=6)
        return ax

    def format_tally_plot(self, ax):

        # xlim is the length of the label list
        ax.set_xlim(0, 5.5)
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
        file_name = "_".join([path_, "_DE_cutoff_by_log2fold.png"])
        fig.savefig(file_name, format='png', bbox_inches='tight')
        fig.close()

    def save_plot(self, fig):

        path_ = os.path.join(self.output_dir, self.prefix)
        file_name = "_".join([path_,"DE_Gene_by_time.png"])

        fig.savefig(file_name, format='png', bbox_inches='tight')

        return fig


    def create_tall_plot(self):

        handles = [self.set_line() for _ in range(3)]
        fig = self.set_figure(handles)
        ax = plt.subplot()

        y_values = self.compute_tally()
        ax = self.create_subplot(ax, y_values)
        ax = format_tally_plot(ax)

        self.save_plot(self, fig)