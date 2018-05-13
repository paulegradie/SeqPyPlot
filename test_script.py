from SeqPyPlot.seqpyplot.analyzer.paired_sample_filter import PairedSampleFilter
from SeqPyPlot.seqpyplot.container.data_container import DataContainer

from SeqPyPlot.seqpyplot.parsers.config_parser import config_parser
from SeqPyPlot.seqpyplot.parsers.gene_list_parser import MakeFigureList

from SeqPyPlot.seqpyplot.plot.paired_line_plotter import PairedDataLinePlotter
from SeqPyPlot.seqpyplot.plot.bar_plotter import PairedBarPlot
from SeqPyPlot.seqpyplot.plot.de_tally_plotter import TallyDe
from SeqPyPlot.seqpyplot.plot.scatter_plotter import ScatterPlots
from SeqPyPlot.seqpyplot.plot.PCA import PCADecomposition


from shutil import rmtree
import pandas as pd
from sys import exit

try:
    rmtree('./test_output')
except:
    pass


if __name__ == "__main__":

    config = 'config.ini'
    config_obj = config_parser(config)

    container = DataContainer(config_obj=config_obj)
    analyzer = PairedSampleFilter(config_obj=config_obj, container_obj=container)

    de_gene_list = analyzer.complete_de_gene_list
    plottable_data = analyzer.normalized_df
    fig_list = MakeFigureList(config_obj)


    # line_plotter = PairedDataLinePlotter(config_obj, analyzer, de_gene_list)
    # line_plotter.plot_figure(figure_list=fig_list.plot_groups, plottable_data=plottable_data)

    # bar_plotter = PairedBarPlot(config_obj=config_obj)
    # bar_plotter.create_bar_plot(analyzer.de_count_by_stage)

    scatter_plotter = ScatterPlots(config_obj=config_obj, container_obj=container, analyzer_obj=analyzer)
    scatter_plotter.create_scatter_plots()

    # tally_plotter = TallyDe(config_obj, container)
    # tally_plotter.create_tally_plot()

    # pca_decomp = PCADecomposition(config_obj, container)
    # pca_decomp.create_pca_plot()

    print("Script completed no errors")

    