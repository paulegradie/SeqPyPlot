import time
import os
import sys

import SeqPyPlotLib.DataContainer as DataContainer

from SeqPyPlotLib.DataPlotter import MainDataPlotter
from SeqPyPlotLib.DataAnalyzer import DataAnalyzer
from SeqPyPlotLib.ArgumentCollector import Args

from SeqPyPlotLib import timing


print '\n' + '{:^68}'.format(
    '***SeqPyPlot v0.2***'), '\nA tool for helping you analyze pilot data (data without replicates).\n'

if len(sys.argv) == 0:
    print "Set options."
    sys.exit()
start = float(time.time())
print '\n\n'
# Collect args
argz = Args()
argz.make_logs()
args = argz.args


# if neither raw data and plot data is given... exit
if args.raw_data is None and args.plot_data is None:
    print("\nProvide either a (raw_data file) or (plot_data file)\n")
    os.system('python SeqPyPlot.py -h')
    sys.exit()

# if both raw data and plot data is given... exit
elif args.raw_data is not None and args.plot_data is not None:
    print("\nProvide either a (raw_data file) or (plot_data file)\n")
    os.system('python SeqPyPlot.py -h')
    sys.exit()

else:
    if args.raw_data is not None:
        FullContainer = DataContainer.DataContainer(args)  # object containing parsed raw_data
        Analyzer = DataAnalyzer(args, FullContainer)  # object used for filtering data and writing out filtered data files

        Analyzer.seqpyfilter()  # object containing analyzed data

        DataPrinter = DataContainer.DataPrinter(args, FullContainer, Analyzer)
        DataPrinter.write_plot_data()
        DataPrinter.write_de_results()
        DataPrinter.write_filtered_data()
        Analyzer.print_analyzer_results()

        if args.ercc:
            DataPrinter.write_ercc_data()

        print(
            "Data analyzed, select genes for plotting from {}_filtered.txt and rerun the program with plot data.".format(
                args.out))
        sys.exit()

    elif args.plot_data is not None:
        if args.gene_list is not None:  # if plot data and gene list is provided

            # Create Output files
            DataContainer.PrepareOutputDirectory.make_folder(args.out)

            FullContainer = DataContainer.DataContainer(args)  # object containing parsed premade plot data
            if args.de_results is not None:  # if no results are provided, calculate them and then load them
                FullContainer.load_results()
                Analyzer = DataAnalyzer(args, FullContainer)
            else:
                Analyzer = DataAnalyzer(args, FullContainer)

            Analyzer.seqpyfilter()

            assert FullContainer.analyzed

            # Ready, get set....
            FigureList = DataContainer.MakeFigureList(args)  # object containing gene list for plotting and other attributes
            Plot_Builder = MainDataPlotter(args, Analyzer, FigureList)  # object used for generating plots
            Plot_Builder.plot_figures()

            if report:
                DataPrinter = DataContainer.DataPrinter(args, FullContainer, Analyzer)
                DataPrinter.write_plot_data()
                DataPrinter.write_de_results()
                DataPrinter.write_filtered_data()
                Analyzer.print_analyzer_results()

            if args.scatter or args.plots:
                try:
                    print "Building Scatter Plots for un-flagged data...\n"
                    Plot_Builder.make_scatter_plots()
                except IndexError:
                    print "No unflagged genes found..."
                try:
                    print "Building Scatter Plots for flagged data...\n"
                    Plot_Builder.make_scatter_plots(flagged=True)
                except IndexError:
                    print "No flagged genes found..."
            if args.bar or args.plots:
                print "Building Bar Plots...\n"
                Plot_Builder.de_bar('black')
            if args.histo or args.plots:
                print "Plotting Histograms...\n"
                Plot_Builder.plot_histograms()
            if args.log2histo or args.plots:
                print "Building Scatter Plots for collecticue log2fold data...\n"
                Plot_Builder.collective_log_plot()
                Plot_Builder.single_log_plots()
            if args.bland_alt or args.plots and args.num == 2:
                print "Building Bland-Altman Plots...\n"
                Plot_Builder.bland_altman_plot()
            if args.bland_grad or args.plots:
                print "Building Bland-Gradie Plots...\n"
                Plot_Builder.bland_gradie_plot()
                Plot_Builder.bland_gradie_plot(flagged=True)
            if args.tally:
                print "Performing Tallys...\n"
                Plot_Builder.plot_tally()
            if args.ercc:
                DataPrinter.write_ercc_data()

        else:  # IF plot data is provided without a genelist
            # perform the DE analysis using processed plot data
            FullContainer = DataContainer.DataContainer(args)
            Analyzer = DataAnalyzer(args, FullContainer)
            Analyzer.seqpyfilter()
            Plot_Builder = MainDataPlotter(args, Analyzer, None)

            if report:
                DataPrinter = DataContainer.DataPrinter(args, FullContainer, Analyzer)
                DataPrinter.write_plot_data()
                DataPrinter.write_de_results()
                DataPrinter.write_filtered_data()
                Analyzer.print_analyzer_results()

            if args.scatter or args.plots:
                try:
                    print "Building Scatter Plots for un-flagged data...\n"
                    Plot_Builder.make_scatter_plots()
                except IndexError:
                    print "No unflagged genes found..."
                try:
                    print "Building Scatter Plots for flagged data...\n"
                    Plot_Builder.make_scatter_plots(flagged=True)
                except IndexError:
                    print "No flagged genes found..."
            if args.bar or args.plots:
                print "Building Bar Plots...\n"
                Plot_Builder.de_bar('black')
            if args.histo or args.plots:
                print "Plotting Histograms...\n"
                Plot_Builder.plot_histograms()
            if args.log2histo or args.plots:
                print "Building Scatter Plots for collecticue log2fold data...\n"
                Plot_Builder.collective_log_plot()
                Plot_Builder.single_log_plots()
            if args.bland_alt or args.plots and args.num == 2:
                print "Building Bland-Altman Plots...\n"
                Plot_Builder.bland_altman_plot()
            if args.bland_grad or args.plots:
                print "Building Bland-Gradie Plots...\n"
                Plot_Builder.bland_gradie_plot()
                Plot_Builder.bland_gradie_plot(flagged=True)
            if args.tally:
                print "Performing Tallys...\n"
                Plot_Builder.plot_tally()
            if args.ercc:
                DataPrinter.write_ercc_data()

            print "Finished."
            sys.exit()
    else:
        os.system('python SeqPyPlot.py -h')
        sys.exit()

print "Cleaning up..."
print "Program Complete. View results in {}".format(args.out)