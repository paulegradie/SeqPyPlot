from sklearn.decomposition import PCA
from seqpyplot.plot.base.plot_base import PlotBase
import matplotlib.pyplot as plt
import os

class PCADecomposition(PlotBase):

    def __init__(self, output_dir, container_obj, experiment_name, plottable_data, col_names):
        """perform PCA on input data columnwise

        Arguments:
            container_obj {obj} -- data container object
            experiment_name {str} -- name of experiment (AKA PREFIX)
            plottable_data {data frame} -- data frame of normalized expression data (column normalized)
            col_names {list} -- list of sample names
        """
        super()

        self.plottable_data = plottable_data  # container_obj.normalized_df

        self.names = col_names  #self.plottable_data.columns.tolist()
        self.output_dir = output_dir
        self.experiment_name = experiment_name  #self.config_obj.get('names', 'experiment_name')

    def execute_pca(self, n_components):

        pca = PCA(n_components=n_components, svd_solver='full')
        pca.fit(self.plottable_data)

        x = pca.components_[0, :]
        y = pca.components_[1, :]

        return x, y

    def set_figure(self):
        # create a figure for the current figure
        fig = plt.figure(num=1,
                         figsize=(7, 7),
                         dpi=1200,
                         edgecolor='black',
                         frameon=False,
                         )

        # set figure tit
        fig.suptitle("PCA Decomposition",
                     verticalalignment='top',
                     horizontalalignment='center',
                     fontsize=14,
                     )

        # rcParams['legend.frameon'] = 'False'

        return fig

    def generate_subplot(self, x, y):

        ax = plt.subplot()
        ax.scatter(x, y)
        return ax

    def tidy_plot(self, ax):

        ax.set_xlabel("First principle component")
        ax.set_ylabel("Second principle component")

        return ax

    def create_pca_plot(self, n_components=2):

        fig = self.set_figure()
        x, y = self.execute_pca(n_components)
        ax = self.generate_subplot(x, y)
        ax = self.tidy_plot(ax)

        for idx, name in enumerate(self.names):
            ax.annotate(name, (x[idx], y[idx]))

        self.save_plot(fig)

    def save_plot(self, fig):

        path_ = os.path.join(self.output_dir, self.experiment_name)
        file_name = "_".join([path_,"PCA_plot.png"])

        fig.savefig(file_name, format='png', bbox_inches='tight')

        return fig

