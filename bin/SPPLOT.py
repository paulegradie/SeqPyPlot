import os
from argparse import ArgumentParser

from seqpyplot.analyzer.paired_sample_filter import PairedSampleFilter
from seqpyplot.container.data_container import DataContainer
from seqpyplot.printers.data_printer import DataPrinter
from seqpyplot.parsers.config_parser import config_parser
from seqpyplot.parsers.gene_list_parser import MakeFigureList
from seqpyplot.plot.bar_plotter import PairedBarPlot
from seqpyplot.plot.de_tally_plotter import TallyDe
from seqpyplot.plot.paired_line_plotter import PairedDataLinePlotter
from seqpyplot.plot.PCA import PCADecomposition
from seqpyplot.plot.scatter_plotter import ScatterPlots
from seqpyplot.utils import make_default_output_dir


if __name__ == "__main__":
#--------------------------------------------------------------------
# Arguments

    usage = """
    SPPLOT only requires a path to a configuration file. The config file
    is expected to be in .ini format. See the examples/ directory for an
    example config file.

    NYI indicates not-yet-implemented.
    """

    epilog = "\n\n"
    
    parser = ArgumentParser(description=usage,
                            prog='SPPLOT v0.4',
                            epilog=epilog)


    parser.add_argument('-c', '--config', type=str, required=True)

    parser.add_argument('-o', '--overwrite', action='store_true', help='Overwrite existing output directory.')
    parser.add_argument('-i', '--impute', action='store_true', help='Impute missing data using neighbors: NYI')
    parser.add_argument('-p', '--plot',   action='store_true', help='Create all plots.')
    parser.add_argument('-t', '--correct', action='store_true', help='Attempt to correct for heteroskedasticity: NYI')
    parser.add_argument('-u', '--unnorm', action='store_true', help='Do not normalize data.')

    args = parser.parse_args()

#--------------------------------------------------------------------
# Parse and correct incoming data

    config_obj = config_parser(args.config)

    # load the data container_obj
    container_obj = DataContainer(config_obj)
    data, ercc_data = container_obj.parse_input()

    if not args.unnorm:
        data = container_obj.normalize_file_pairs(data) # Single df of normalized data

    split_data = container_obj.split(data)  # List of normalized dfs

    if args.correct:
        split_data = container_obj.correct_heteroskedacity(split_data)

#--------------------------------------------------------------------
#  Filter data

    filter_obj = PairedSampleFilter(config_obj)
    filter_result = filter_obj.main_filter_process(split_data)


#--------------------------------------------------------------------
# Save filter results

    output_path = config_obj.get('data_directory', 'output')
    output_path = make_default_output_dir(output_path or None, args.overwrite)

    data_printer = DataPrinter(config_obj, container_obj or None, filter_obj or None)
    data_printer()

#--------------------------------------------------------------------
# Generate Plots

    if args.plot:
        line_plotter = PairedDataLinePlotter(config_obj, filter_obj, data)
        fig_list = MakeFigureList(config_obj)
        line_plotter.plot_figure(figure_list=fig_list.plot_groups, plottable_data=data)

        bar_plotter = PairedBarPlot(config_obj=config_obj)
        bar_plotter.create_bar_plot(filter_obj.de_count_by_stage)

        scatter_plotter = ScatterPlots(config_obj=config_obj, container_obj=container_obj, filter_obj=filter_obj)
        scatter_plotter.create_scatter_plots()

        tally_plotter = TallyDe(config_obj, container_obj)
        tally_plotter.create_tally_plot(split_data)

        pca_decomp = PCADecomposition(config_obj, container_obj)
        pca_decomp.create_pca_plot()



#--------------------------------------------------------------------
# Tidy up

    print("Script completed no errors")
