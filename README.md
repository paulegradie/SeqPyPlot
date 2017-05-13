
#             ***SeqPyPlot v0.2***
A tool for helping you analyze pilot data (data without replicates).

## Installation
There are no special installation requirements for this project! It is written in pure python directly from the command line (or IDE of your  choice).

## Dependencies

#### --Matplotlib
#### --Seaborn
Most machines come preinstalled with matplotlib. Its the library used for creating all of the plots so grab it if you don't have it.
Seaborn is used to styalize all of the plots. At the moment, the code breaks without it - however it would be smart of me to make the  implementation of seaborn option by adding a quick on import.

## Usage
    
#### usage: python SeqPyPlot.py [options] [-raw_data (or) -plot_data] -num 2 -time d1,d2,d3

#### Required: -raw_data or -plot_data
<pre>
-h, --help            show this help message and exit

   --------------------Required Options--------------------  
   
      -time d1,d2,d3               A comma separated list of time points.  
      -num 2                       Default: 2. Set number of plots.  
      -out Default_out             Output Folder Name.  
      -prefix SeqPyPlot_           Leading name of output file.  
      -c S1,S2                     A comma separated list of conditions (max 2).  
      -data_type htseq             Either cuffnorm, cuffdiff, deseq2, or edgeR. 
      
   --------------------Required Input Options--------------------   
   
      -raw_data None               Input file or folder.  
      -plot_data None              Formatted input data to plot.  
      -unform_plot_data            Default: False. Reorder plot data (1,1,2,2 -> 1,2,1,2).  
      -gene_list None              Single Column Gene list in txt file.  
      -input_results None          Optional. Your own flagged gene list.  
   
   --------------------Filter Options--------------------  
   
      -low 25                      Default: 2. Set the min expression value to accept.   
      -hi 1000000                  Default: 5mil. Set the max expression value to accept.   
      -dif_range 25,1000000        Default: 25-1000000. Set min difference in expression.  
      -log2 0.7                    Default: 0.7. Minimum log2fold change to accept.  
   
   --------------------Analysis Options--------------------  
   
      -find_housekeeping           Default: False. Search for housekeeping genes.  
      -num_housekeeping 3          Default: 3. Min number of housekeeping to detect.  
      -report                      Default: False. Write plot data and filter results.  
      -ercc                        Default: False. Write ERCC data to an output file.  
      -r                           Default: False. Use to remove genes not always on.  
   
   --------------------Plot Options--------------------  
      
      -plots                       Default: False. Make all plots.  
      -tally                       Default: False. Tally DE genes.  
      -bar                         Default: False. Construct bar plots.  
      -scatter                     Default: False. Construct scatter plots.  
      -scat_range 0,1200           Default: 0,5000000. Set scatter plot value range.  
      -bland_alt                   Default: False. Construct bland-altman plots.  
      -barange 0,10000             Default: 1,25000. Range for bland-altman plots.  
      -bland_log                   Default: False. Construct bland-alt-log plots.  
      -blrange 0,20                Default: 0,20. Range for Bland-alt-log plots.  
      -histo                       Default: False. Construct histogram plots.  
      -hist_range 1,1000           Default: 1.0. Lower x axis limit for histogram.  
      -log2histo                   Default: False. Construct log2fold histogram plots.  
      -svg                         Default: False. Use to svg plots.  

</pre>

## Contributing
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -am 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature`
5. Submit a pull request :D

## History and background
SeqPyPlot (rna-SEQ PYthon PLOT library) was written to provide a simple to use quick analysis tool for pilot RNA-seq data. Pilot study data is defined here as data that doesn't include replicates. A common need among developmental biologists right now is the ability to interogate pilot study data and access gene expression quickly. Usually these data are not used for publishing, but rather to guide future experiments and provide a quick peak at approximate expression levels within the target tissue. 

Quoted from the edgeR manual (a popular tool written in R used for calculating differential expression of genes:

`RNA-Seq and ChIP-Seq are still expensive technologies, so it sometimes happens that only one library can be created for each treatment condition. In these cases there are no replicate libraries from which to estimate biological variability. In this situation, the data analyst is faced with the following choices, none of which are ideal. We do not recommend any of these choices as a satisfactory alternative for biological replication. Rather, they are the best that can be done at the analysis stage [using no replicates].`

The manual goes to suggest the following option (which is admittedly not their preferred suggestion...)

`1. Be satisfied with a descriptive analysis, that might include an MDS plot and an analysis
of fold changes. Do not attempt a significance analysis. This may be the best advice.`

The motivation behind SeqPyPlot is to provide a means to analyze time series data under the assumption of no biological variability (even though we know this isn't true). With unreplicated time series data, there is no way to reasonably estimate variance within the data, since any gene could potentially change over time. SeqPyPlot therefore ignores dispersion estimates and biological variability to give a quick and dirty assessment of fold change. 

The advantage of using SeqPyPlot is twofold. First, it provides a set of optimizable parameters that can be used to filter data included in the final output. Second, it provides a very useful plotting function for rapidly visualizing gene expression. Using this program, you can create a list of flagged genes, and plot them within minutes.

## Plot Descriptions
SeqPyPlot provides a range of plots for general descriptive analysis of your data. These plots are available using the 'plot options' shown above in the usage. In this section I'll provide plot examples and their suggested interpretation.

#### Tally Plots
This plot provides feedback on your filter parameter selection. These paramters determine the number of genes that will be flagged. Possibly the most important paramter is the logfold change threshold, so by using the '-tally' option, we automatically count the number of flagged genes across a range of threshold values, while keeping the other filter parameters constant. A good value to choose for logfold change will be around the point in the plot that the line begins to steepen.

![Tally Plot](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Tally Plot")

#### Flagged Gene Bar Plots
Bar plots are provided to give information on the number of genes flagged at each stage given the current filter parameters.

![Bar Plot](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Bar Plot")

#### Filter Parameters
The following plots are intended to aid the user in finding appropriate filter parameters. The filter takes four paramters illustrated in the following image.

![Bar Plot](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Bar Plot")


#### Scatter Plots
Filter paramters

![Bar Plot](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Bar Plot")

#### Bland_altman Plots
![Bar Plot](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Bar Plot")

#### Gene Count histogram Plots
![Bar Plot](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Bar Plot")

#### Logfold change histogram Plots
![Bar Plot](https://github.com/adam-p/markdown-here/raw/master/src/common/images/icon48.png "Bar Plot")


## Credits
The current and previous versions of SeqPyPlot were built by me (Paul G). I would LOVE to add more credits! Fork this repository and help make SeqPyPlot a more useful and user friendly tool!

## License
This software is currently open source, and I'm looking to write a short paper on the utility of this software to submit to an opensource journal. IF/when that happens, citations may be in order. Maybe. Probably not worth worrying about though!




