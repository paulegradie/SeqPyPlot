from __future__ import division
import os
import sys
import numpy as np
import logging
import seaborn as sns
sns.set_palette('Blues', n_colors=256)


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



class PairedDataLinePlotter(object):

    def __init__(self, config_obj, analyzer_obj, de_gene_list):

        self.config_obj = config_obj
        self.output_dir = self.create_output_directory()
        
        self.analyzer_obj = analyzer_obj
        self.normalized_df = analyzer_obj.normalized_df
        self.de_gene_list = de_gene_list

        self.log = self.analyzer_obj.log2fold
        self.low = self.analyzer_obj.low
        self.times = self.config_obj.get('names', 'times')
        self.labels = self.config_obj.get('names', 'conditions')
        self.prefix = self.config_obj.get('names', 'experiment_name')

    def create_output_directory(self):
        #TODO Fix this to make it work
        dir_name = 'default_dir'
        if not self.config_obj.get('data_directory', 'output'):
            os.mkdir(dir_name)
        else:
            dir_name = self.config_obj.get('data_directory', 'output')
            os.mkdir(dir_name)
        return dir_name

    def set_figure(self, figure_prefix, **args):
        fig = plt.figure(num=1,
                         dpi=600,
                         figsize=(10, 10),
                         edgecolor='black',
                         frameon=False,
                         )
        fig.suptitle(figure_prefix,
                     verticalalignment='top',
                     horizontalalignment='right',
                     fontsize=24
                     )
        return fig

    def set_line(self, line_color, edge_color, face_color, kwargs={}):
        series_line = mlines.Line2D([], [],
                                    color=line_color,
                                    markeredgecolor=edge_color,
                                    markerfacecolor=face_color,
                                    **kwargs)
        return series_line
            
    def plot_figure(self, figure_list, plottable_data):
        """[summary]
        
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
        for fig_idx, figure in enumerate(figure_list, start=1):

            fig = self.set_figure(figure_prefix=self.config_obj.get('file_names', 'prefix'))

            series1 = self.set_line(line_color='blue',
                                    edge_color='blue',
                                    face_color='white',
                                    kwargs=line_plot_kwargs)
            series2 = self.set_line(line_color='red',
                                    edge_color='red',
                                    face_color='black',
                                    kwargs=line_plot_kwargs)
            
            log_line = self.set_line(line_color='white',
                                     edge_color='white',
                                     face_color='white')
            fold_line = self.set_line(line_color='white',
                                      edge_color='white',
                                      face_color='white')
            diff_line = self.set_line(line_color='white',
                                      edge_color='white',
                                      face_color='white')

            for idx, gene in enumerate(figure, start=1):

                fig.add_subplot(3, 2, idx).set_title('- ' + str(gene) + ' -', fontsize=20)

                if not self.gene_exists(gene):
                    logging.warning('Don"t worry - this just means the gene is probably turned off.')
                    continue
                
                is_de = self.is_de(gene)  # check if the gene is DE for bg highlighting
                self.create_subplot(gene, is_de)

            handles = [series1, series2, log_line, fold_line, diff_line]
            fig = self.tidy_up_figure(fig, handles)

            # self.save_figure(fig, fig_idx, figure_list)
            # plt.show()
            plt.clf()
            plt.cla()
            plt.close()

    def save_figure(self, fig, fig_idx, figure_list):
        path_ = os.path.join(self.output_dir, self.prefix)
        genes = str([fi.strip() for fi in figure_list[fig_idx-1]])
        file_name = "{}_{}_{}.png".format(path_, str(fig_idx), genes)

        fig.savefig(file_name, format='png', bbox_inches='tight')
        

    def tidy_up_figure(self, fig, handles):

        labels = (self.labels
                  + ["Log2: " + str(self.log)]
                  + ["Fold: " + str(round(2.0**self.log, ndigits=1))])
                #   + ["Diff: " + str(int(self.args.dif_range[0])) + ', ' + str(int(self.args.dif_range[1]))])
          
        fig.legend(handles=handles,
                   labels=labels,
                   loc='upper right')

        return fig

    def gene_exists(self, gene):
        return gene in self.normalized_df.index.tolist()

    def is_de(self, gene):
        return gene in self.de_gene_list

    def retrieve_data(self, gene):
        """
        Return plot data and associated masks
        
        Arguments:
            plottable_data {[type]} -- [description]
        
        Returns:
            [type] -- [description]
        """
        data_length = self.analyzer_obj.normalized_df.shape[1]

        try:
            data = self.analyzer_obj.normalized_df.loc[gene]
            
        except KeyError:
            print("The current gene: -- {} -- was not found in the plot data.".format(gene_name))
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

        num_stages = len(self.times)

        # xlim is the length of the label list
        ax.set_xlim(-0.85, (num_stages)-0.15)
        ax.spines['right'].set_color('none')
        ax.spines['top'].set_color('none')
        ax.xaxis.set_ticks_position('bottom')
        ax.yaxis.set_ticks_position('left')
        ax.spines['bottom'].set_position(('data', 0))

        ax.set_xticks([i for i in range(num_stages)], self.times)

        if is_de:  # set background color
            ax.set_facecolor('#FFFC97')

        ax.axhline(self.low, color='gray', linestyle='--', label="-low")

        if y_max < self.low:
            ax.set_ylim(0, self.low + 10)
        else:
            ax.set_ylim(0, y_max)

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
            b = (2.0 * i) / ((2.0 ** self.log) + 1)
            dif = i - b
            diffs.append(float(dif))
       
        diffs = np.asarray(diffs)

        if min(diffs) < 0 or min(series_mean - diffs) < 0:
            print("Error bars out of range - but still continuing.")

        return diffs

    def plot_log_fold_bounds(self, series_mean, diffs, x_axis):
        return plt.errorbar(x_axis,
                            series_mean,
                            linestyle='None',
                            mfc='black',
                            yerr=[diffs, diffs],
                            ecolor='black',
                            elinewidth='1',
                            label="Range around the mean = log2(range)=1 ")

    def tidy_up_plot(self):

        plt.subplots_adjust(top=0.85)
        plt.tight_layout()
        return plt

    def compute_max_yval(self, series1_data, series2_data, diffs, scale=1.4):
        return max(np.mean(np.mean(np.asarray([series1_data, series2_data]), axis=0)) + diffs) * float(scale)

    def plot_series(self, series_data, series_mask, ax, xs, kwargs):
            return ax.plot(xs[series_mask],
                           series_data[series_mask],
                           'bo',
                           **kwargs)

    def create_subplot(self, gene, is_de):
        """Function for plotting individual plots in the figure.
        
        Arguments:
            gene {[type]} -- [description]
            is_de {bool} -- [description]
        """

        (series1_data, series1_mask,
         series2_data, series2_mask) = self.retrieve_data(gene)    

        x_axis = np.arange(float(max(len(series1_data), len(series2_data))))

        data_mean, mean_mask = self.compute_mean(series1_data, series2_data)

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
        plot_kwargs.update({'color': 'blue', 'label': self.labels[0]}) # .update returns NoneType
        self.plot_series(series1_data, series1_mask, ax, x_axis,
                         plot_kwargs
                         )
        plot_kwargs.update({'color': 'red', 'label': self.labels[1]})
        self.plot_series(series2_data, series2_mask, ax, x_axis,
                         plot_kwargs
                         )
        _ = self.plot_log_fold_bounds(data_mean, diffs, x_axis)  # May need the mean mask here

        self.tidy_up_plot()

        return ax  # may need to return plt

