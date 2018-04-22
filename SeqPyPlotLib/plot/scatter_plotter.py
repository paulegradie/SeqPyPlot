from base.plot_base import PlotBase
import matplotlib.pyplot as plt
import matplotlib.lines as mlines


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

        self.config_obj = config_obj
        self.container_obj = container_obj
        self.analyzer_obj = analyzer_obj

        self.output_dir = self.create_output_directory()
        self.prefix = self.config_obj.get('names', 'experiment_name')

        self.log2fold = self.config_obj.getfloat('params', 'log2fold')
        self.low = self.config_obj.getint('params', 'low')
        self.hi = self.config_obj.getint('params', 'hi')
        self.diff = self.config_obj.getlist('params', 'diff')
        self.times = self.config_obj.get('names', 'times')


        self.scatrange = self.config_obj.getlist('plot_options', 'scatrange')

        # organize plot data
        self.flagged_data = self.analyzer_obj.filtered_genes  # a list of dfs
        self.unflagged_data = self.collect_unflagged_data()  # a list of dfsff

    def collect_unflagged_data(self):
        complete_gene_set = set(self.container_obj.complete_gene_list)

        unflagged_data = list()
        for (filtered_df,
             normalized_df) in zip(self.analyzer_obj.filtered_genes,
                                   self.container_obj.split_dfs):
            
            flagged_genes = set(filtered_df_list.index.tolist())
            
            unflagged = list(complete_gene_set - flagged_genes)
            unflagged_data.append(normalized_df.iloc[unflagged])
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
            dfs = [flagged_df, unflagged_df]
            for ax, df, title in zip(axes.flatten(), dfs, titles):
                ax = self.format_plot(ax, title)
                cols = df.columns.tolist()
                df.plot(kind='scatter',
                        ax=ax,
                        x=cols[0],
                        y=cols[1])
            
            self.save_fig(fig)
            plt.close()


    def create_figure(self, suptitle):

        fig, axes = plt.subplots(nrows=1, ncols=2)
        fig = plt.figure(num=1,
                         dpi=1200,
                         figsize=(16, 12),
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
        labels = [" ".join(['ScatRange:', str(self.scatt_range)]),
                  "DiffRange: {}, {} ".format(str(self.diff[0]), str(self.diff[1])),
                  "Log2: " + str(self.log)]
        
        fig.legend(handles=handles,
                   labels=(labels),
                   loc='upper right')
        fig.tight_layout()

        return fig, axes

    def format_plot(self, ax, time):

        ax.axvline(self.low, color='gray', linestyle='--')
        ax.axvline(self.hi, color='gray', linestyle='--')

        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.get_xaxis().tick_bottom()
        ax.get_yaxis().tick_left()
        ax.set_xlim(self.scatrang[0])
        ax.set_ylim(self.scatrang[1])
        ax.set_xlabel('Mean', fontsize=10)
        ax.set_ylabel('Count', fontsize=10)
        ax.set_title(title, fontsize=12)
        ax.tick_params(axis='both', which='major', labelsize=8)

        return ax

    def set_line(self, color):
        return mlines.Line2D([], [], color=color)


    def save_fig(self, fig, time):
      
        path_ = os.path.join(self.output_dir, self.prefix)
        file_name = "_".join([path_, time, 'ScatterPlot.png'])

        fig.savefig(file_name, format='png', bbox_inches='tight')

        return fig
