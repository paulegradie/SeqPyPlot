
from argparse import ArgumentParser

from seqpyplot.parsers.config_parser import config_parser
from seqpyplot.parsers.gene_list_parser import MakeFigureList

from seqpyplot.analyzer.paired_sample_filter import \
    PairedSampleFilter
from seqpyplot.container.data_container import DataContainer


from seqpyplot.plot.paired_line_plotter import PairedDataLinePlotter
from seqpyplot.plot.bar_plotter import PairedBarPlot
from seqpyplot.plot.de_tally_plotter import TallyDe
from seqpyplot.plot.scatter_plotter import ScatterPlots
from seqpyplot.plot.PCA import PCADecomposition



if __name__ == "__main__":
    """
    Script to normalize data frame independantly
    """
    usage = """
    This is how you use this program
    """
    name = "SeqpyPlot v0.4"

    
    parser = ArgumentParser(description=name,
                            prog='SeqPyPlot v0.4')#,
                            # usage=usage,
                            # epilog=epilog)


    parser.add_argument('-c', '--config', type=str, required=True)
    parser.add_argument('-o', '--output', type=str, default='default')
    parser.add_argument('-i', '--impute', action='store_true')
    args = parser.parse_args()

    config_obj = config_parser(args.config)
    
    # load the data container
    container = DataContainer(config_obj)
    raw_df, ercc_data = container.parse_input()
    
    normalized_df = container.reorder_cols(container.normalize_file_pairs(raw_df))
    split_normalized_dfs = container.split(normalized_df)

    complete_gene_list = normalized_df.index.tolist()

    Filter = PairedSampleFilter(config_obj)
    filter_result = Filter(split_normalized_dfs)
    
    fig_list = MakeFigureList(config_obj)

    # line_plotter = PairedDataLinePlotter(config_obj, analyzer, de_gene_list)
    # line_plotter.plot_figure(figure_list=fig_list.plot_groups, plottable_data=plottable_data)

    # bar_plotter = PairedBarPlot(config_obj=config_obj)
    # bar_plotter.create_bar_plot(analyzer.de_count_by_stage)

    scatter_plotter = ScatterPlots(config_obj=config_obj, container_obj=container, analyzer_obj=Filter)
    scatter_plotter.create_scatter_plots()

    # tally_plotter = TallyDe(config_obj, container)
    # tally_plotter.create_tally_plot()

    # pca_decomp = PCADecomposition(config_obj, container)
    # pca_decomp.create_pca_plot()

    print("Script completed no errors")