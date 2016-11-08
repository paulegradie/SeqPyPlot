from SeqPlotLib import *
import os
import time

__author__ = "Paul Gradie"


if __name__ == '__main__':

    print "SeqPlot - v0.1"

    try:

        # instantiate classes
        Argparser = ArgumentParser()
        Data_Reader = ReadNormalizedData()
        Gene_List_Parser = ParseGeneList()
        Plot_Builder = MainDataPlotter()
    
        # Set Arguments
        args = Argparser.argument_parser()
    
        time_label = Argparser.time_parser(args.time)
        condition_label = args.condition
    
        # Create Output files
        PrepareOutputDirectory(args.out)
    
        # Read Data
        plot_data = Data_Reader.read_data(args.plotter_data)  # create plot data object
    
        figure_list = Gene_List_Parser.input_list_parser(args.gene_list)  # create formatted nested figure list object
    
        # Construct Plots
        start = time.time()
        Plot_Builder.plot_figures(plot_data, time_label, condition_label, figure_list, args.out, args)
        end = time.time()
    
        print "Total Run Time: {}".format((end - start) / 60)

    except TypeError:
        print 'Parameters aren\'t set...'
        os.system('python SeqPlot.py -h')
        
