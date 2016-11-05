# SeqPyPlot

SeqPyPlot is a program for plotting time series data. Data is formatted in to a tab delimeted text file and then fed in to the main script. The script accesses code stored in SeqPyPlotLib and generates plots for the input data, saving them to .png and .svg files. 

The program was originally written to plot output from the Cufflinks suite. A pipeline was created where RNA-seq data would be mapped, quantified, and then normalized together to allow for comparisons between samples. This ouptput was then formatted and fed in to SeqPyPlot to plot expression values for lists of genes. 

This project is in early stages and currently consists of the growing SeqPyPlotLib module with useful code for plotting Time Series RNA-seq data and a run script. Several classes are defined in the module that provide any developers with tools to organize data processing and plotting in the run script.

# The run Script

In the run script, SeqPyPlot.py, you should aim to first read in any command line argumentscreate an output directory, then read in the data, and finally produce the plots.

## Usage

First make sure that all of the programs files are in the same directory. Then run the script:

>python SeqPyPlot.py -h

SeqPlot - v0.1
usage: SeqPlot.py [-h] [-t] [-o] [-c] [-n] [gene_list] [plotter_data]

Plot time series expression data without replicates.

positional arguments:
 gene_list              Single column gene list in txt file.
 plotter_data           Formatted plot-able data.

optional arguments:
 -h, --help             show this help message and exit             
 -t, --time             A comma separated list of time points.       Default: 1,2,3,4,5
 -o, --out              Output Folder Name                           Default: default_out
 -c, --condition        A comma separated list of conditions (max 2) Default: series1,series2
 -n, --num              Number of samples per plot.                  Default: 2
 
 
# Typical usage:
 
 >$ python SeqPyPlot.py -t e12,e13,e14,e15 -o Dev_genes -c Control,Mutant genelist.txt plotter_data.txt
 
 
 # Positional Arguments
These arguments need to be added in the same order as above after any options.

gene_list               A list of genes saved as a text file. One gene per row.
plotter_data            A text file containing tab delimited rows of gene expression data. See Data Formatting below.

 # Optional Arguments
 The arguments are read in as strings.
 
 -t or --time           To set this argument, write as many times as you need to match your data, and only put a comma between them.
 -o or --out            The name of the output directory AND the name of the main figure output (where all of the plots are).
 -c or --condition      To set this argument, write comma separated names for your samples. Series 1 first, Series 2 second.
 -n or --num            This is number of series you intend to plot. It defaults to 2.
 
 
 
