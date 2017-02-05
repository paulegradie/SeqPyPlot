import time
import os
import sys

from DataPlotter import MainDataPlotter
from DataAnalyzer import DataAnalyzer
from ArgumentCollector import Args
import DataContainer

if __name__ == '__main__':

    start = float(time.time())

    # TODO implement args.makelogs()
    # Collect args
    args = Args()
    args.make_logs()

    # if files aren't provided properly
    if args.args.raw_data == args.args.plot_data is None:
        print("Provide either a (raw_data file) or (plot_data file)")
        os.system('python SeqPyPlot.py -h')
        sys.exit()

    elif args.args.raw_data is not None and args.args.plot_data is not None:
        print ("ERROR...")
        print "\nLoad either raw_data or pre-formatted plot_data. Don't load both!\n"
        sys.exit()
    elif args.args.raw_data is not None and args.args.gene_list is not None:
        print("ERROR...")
        print "\nRun the program on raw_data, then rerun supplying the _plotter_data.txt file.\n"
        sys.exit()
    else:
        pass

    # if only raw data is supplied

    if args.args.raw_data is not None:
        FullContainer = DataContainer.DataContainer(args)  # object containing parsed raw_data
        Analyzer = DataAnalyzer(args, FullContainer)  # object used for filtering data and writing out filtered data files

        Analyzer.seqpyfilter()  # object containing analyzed data

        DataPrinter = DataContainer.DataPrinter(args, FullContainer, Analyzer)
        DataPrinter.write_plot_data()
        DataPrinter.write_de_results()
        DataPrinter.write_filtered_data()

        if args.args.ercc:
            DataPrinter.write_ercc_data()

        if args.args.tally:
            Plot_Builder = MainDataPlotter(args, Analyzer, None)  # object used for generating plots
            Plot_Builder.de_bar('black')
            Plot_Builder.plot_tally()
            Analyzer.print_de_tallies()

        print(
            "Data analyzed, select genes for plotting from {}_filtered.txt and rerun the program with plot data.".format(
                args.args.out))
        # args.make_logs()
        sys.exit()

    elif args.args.plot_data is not None:
        if args.args.gene_list is not None:  # if plot data and gene list is provided

            # Create Output files
            DataContainer.PrepareOutputDirectory.make_folder(args.args.out)

            FullContainer = DataContainer.DataContainer(args)  # object containing parsed premade plot data
            if args.args.de_results is not None:  # if no results are provided, calculate them and then load them
                FullContainer.load_results()
                Analyzer = DataAnalyzer(args, FullContainer)
            else:
                Analyzer = DataAnalyzer(args, FullContainer)

            Analyzer.seqpyfilter()

            assert FullContainer.analyzed

            # Ready, get set....
            FigureList = DataContainer.MakeFigureList(args)  # object containing gene list for plotting and other attributes
            Plot_Builder = MainDataPlotter(args, Analyzer, FigureList)  # object used for generating plots

            if args.args.report:
                DataPrinter = DataContainer.DataPrinter(args, FullContainer, Analyzer)
                DataPrinter.write_plot_data()
                DataPrinter.write_de_results()
                DataPrinter.write_filtered_data()

            if args.args.tally:
                Plot_Builder.de_bar('black')
                Plot_Builder.plot_tally()
                Analyzer.print_de_tallies()

            Plot_Builder.plot_figures()
        else:
            print "Please supply a gene list for plotting with your plot data."
            os.system('python SeqPyPlot.py -h')
            sys.exit()
    else:
        os.system('python SeqPyPlot.py -h')
        sys.exit()

    end = float(time.time())
    print("Total Run Time: {}".format((end - start) / 60.0))

