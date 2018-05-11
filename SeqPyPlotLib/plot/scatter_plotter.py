from base.plot_base import PlotBase
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os
import numpy as np
import seaborn as sns

plt.style.use('bmh')


class ScatterPlots(PlotBase):
    """Every plot gets its own figure and thus its own save file.
    It is up to the user to present these plots in a combined format
    if they wish to.

    Line plots are the only exception.
    
    Arguments:
        PlotBase {[type]} -- [description]
    
    Returns:
        [type] -- [description]
    """

    def __init__(self, config_obj, container_obj, analyzer_obj):
        super(ScatterPlots, self).__init__()
        plt.close()

        self.config_obj = config_obj
        self.container_obj = container_obj
        self.analyzer_obj = analyzer_obj

        self.output_dir = self.create_output_directory()
        self.prefix = self.config_obj.get('names', 'experiment_name')

        self.log = self.config_obj.getfloat('params', 'log2fold')
        self.low = self.config_obj.getint('params', 'low')
        self.hi = self.config_obj.getint('params', 'hi')
        self.diff = self.config_obj.getlist('params', 'diff')
        self.times = self.config_obj.getlist('names', 'times')

        self.scatrange = self.config_obj.getlist('plot_options', 'scatrange')

        # organize plot data
        self.flagged_data = self.analyzer_obj.filtered_genes  # a list of dfs
        self.unflagged_data = self.collect_unflagged_data()  # a list of dfsff
        import pdb;pdb.set_trace()

    def collect_unflagged_data(self):

        complete_gene_set = set(self.container_obj.complete_gene_list)

        unflagged_data = list()
        for (filtered_df,
             normalized_df) in zip(self.analyzer_obj.filtered_genes,
                                   self.container_obj.split_normalized_dfs):
            
            flagged_genes = set(filtered_df.index.tolist())
            unflagged =  (complete_gene_set - flagged_genes)
            unflagged_data.append(normalized_df.loc[unflagged])
        return unflagged_data

    def create_scatter_plots(self):
        """
        This need only take in the result from analyzer_obj:
            analyzer.filterd_genes
        which is a list of data frames
        
        """

        for flagged_df, unflagged_df, time in zip(self.flagged_data,
                                                   self.unflagged_data,
                                                   self.times):

            suptitle = " ".join([time, 'Expression Scatter Plot'])
            fig, axes = self.create_figure(suptitle)
            
            titles = ['Flagged Genes', 'Unflagged Genes']
            lims = (self.scatrange[0], self.scatrange[1])

            dfs = [flagged_df, unflagged_df]
            for ax, df, title in zip(axes.flatten(), dfs, titles):
                
                df['mean'] = df.mean(axis=1)
                # upper, lower = self.calc_bounds(df['mean'].tolist())

                y_range=range(lims[1])
                up, _, lowerbound = self.calc_bounds(y_range)

                upperbound = []
                for i in range(len(y_range)):
                    upperbound.append(y_range[i]+up[i])

                ax = self.format_plot(ax, title, lims)
                cols = df.columns.tolist()

                ax.scatter(df['mean'], df[cols[0]], s=2, color='blue')
                ax.scatter(df['mean'], df[cols[1]], s=2, color='red')
                ax.plot(range(len(upperbound)), upperbound, color='black', linestyle='--')
                ax.plot(range(len(lowerbound)), lowerbound, color='black', linestyle='--')
                ax.plot(range(lims[1]), range(lims[1]))

            self.save_fig(time)
            plt.close()


    def create_figure(self, suptitle):

        fig, axes = plt.subplots(nrows=1, ncols=2)
        fig = plt.figure(num=1,
                         dpi=1200,
                         figsize=(14, 16),
                         edgecolor='black',
                         frameon=False,
                         )
        fig.suptitle(suptitle,
                     verticalalignment='top',
                     horizontalalignment='right',
                     fontsize=14,
                     y=1.09
                     )

        handles = [self.set_line('white') for _ in range(3)]
        labels = ['ScatRange:' + " - ".join([str(x) for x in self.scatrange]),
                  "DiffRange:" + " - ".join([str(x) for x in self.diff]),
                  "Log2: " + str(self.log)]
        
        fig.legend(handles=handles,
                   labels=(labels),
                   loc='upper right',
                   prop={'size': 5})
        fig.tight_layout()

        return fig, axes

    def format_plot(self, ax, time, lims):

        ax.axvline(self.low, color='gray', linestyle='--')
        ax.axvline(self.hi, color='gray', linestyle='--')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set_ylim(lims)
        ax.set_xlim(lims)
        ax.set_xlabel('Mean', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title(time, fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=8)

        return ax

    def set_line(self, color):
        return mlines.Line2D([], [], color=color)


    def save_fig(self, time):
      
        path_ = os.path.join(self.output_dir, self.prefix)
        file_name = "_".join([path_, time, 'ScatterPlot.png'])

        plt.savefig(file_name, format='png', bbox_inches='tight')


    def calc_bounds(self, series_mean):
        """
        **HIDDEN**
        Use some linear algebra to calculate log2fold range around a given value based on user logfold parameter
        :param series_mean:
        :return:
        """
        var = self.log
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
