
#             ***SeqPyPlot v0.2***
A tool for helping you analyze pilot data (data without replicates).

## SeqPyPlot is currently under a major refactorization! This will bring us to v0.4
I believe there is some potential for a bit of machine learning, or at least more intelligent analysis for this sort of dataset.
Once these ideas are incorporated (which should be much easier to accomplish after the refactor), we'll move to v1.0.


## Installation
There are no special installation requirements for this project! It is written in pure python directly from the command line (or IDE of your  choice). There are however a few dependencies. For Windows users, I strongly suggest installing the Anaconda python distribution. You can find these distributions here:

https://www.continuum.io/downloads

## Dependencies

#### --Matplotlib
#### --Seaborn
#### --Numpy

Matplotlib & Seaborn are the libraries used for creating all of the plots so grab it if you don't have it. Seaborn is used to stylize all of the plots. Numpy is used to perform some vector math during plot generation. These libraries are available via pip, however numpy could cause you some trouble. The easiest way to install all of these dependencies, as far as I'm aware, is to use one of the python distrutions such as Continuum's Anaconda distribution.

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
1. Clone it!
2. Create your feature branch: `git checkout -b my-new-feature`
3. Commit your changes: `git commit -m 'Add some feature'`
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

#### Expression Plots
The main plot time intended for SeqPyPlot users are the expression plots. These plots are new and somewhat controversial.
When we normalize RNA-seq expression data, we can use one of two methods -

1.Between-Sample normalization (used for testing for DE)
- TMM (edgeR)
- Geometic (DESeq2)

2. Within-Sample normalization (for comparing genes within a single sample).
- RPKM (single reads)
- FPKM (Paired end reads)
- CPM, TPM, Quartile, etc (Generic non genelength correcting methods)

Mixing these two forms is not well understood, so when we want to look at global quantification of RNA-seq across many samples... what do we do? If you look through the literature, you will likely not find expression plots using TMM normalized values, however we know that due to common differences in library size after library prep, RPKM/FPKM/CPM, (etc) values are not directly comparable.
The best approach we have right now is to use a relative normalization method such as TMM and explicitly say that the values shown are not absolutely quantified.

The expression plots output from SeqPyPlot use TMM normalized values and show relative difference between values across time/stages. SeqPyPlot provides a secondary measure of relative expression, akin to a scale bar, that is intended to help users interpret the relative differences between expression values in log space. SeqPyPlot calculates the mean between two measurements at a given stage and then calculates the log fold change provided by the user and plots error bars centered on the mean with the range of the logfold difference (default 0.7). The output also highlights (in yellow) any plot for a gene that is flagged at any time point.

![MDSPlot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Expression_plot.png "Expression Plot")

#### MDS Plots
After normalization, it is essential to verify that samples cluster appropriately. If you have two control samples that do not cluster together, then this indicates an inherent problem with sample composition and may invalidate your experiment.

![MDSPlot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/MDS_plot.png "MDS Plot")

#### Tally Plots
This plot provides feedback on your filter parameter selection. These paramters determine the number of genes that will be flagged. Possibly the most important paramter is the logfold change threshold, so by using the '-tally' option, we automatically count the number of flagged genes across a range of threshold values, while keeping the other filter parameters constant. A good value to choose for logfold change will be around the point in the plot that the line begins to steepen.

![Tally Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Tally_plot.png "Tally Plot")

#### Flagged Gene Bar Plots
Bar plots are provided to give information on the number of genes flagged at each stage given the current filter parameters.

![Bar Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Bar_plot.png "Bar Plot")

### Plots for Adjusting Filter Parameters
The following plots are intended to aid the user in finding appropriate filter parameters. The filter takes four parameters illustrated in the following image.

#### --Filter illustration--
![Filter Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Filter.png "Filter")

#### Scatter Plots
Scatter plots are useful for double checking that normalization was successful and that your fold change parameter is reasonable across samples.
There are two plots - flagged and un-flagged. The white line along the $f(x) = x$ vector indicates the center of the value distribution across all means, and all expressin values ar plotted against their relative mean. In other words, the mean is calculated for value pairs (control and experimental), and then each vlue is plotted against its mean. The dotted lines represents the log2fold change (from the user provided log2fold argument, default = 0.7) around the mean.

![Scatter Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Scatter_flagged_plot.png "Unflagged Scatter Plot")

![Scatter Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Scatter_unflagged_plot.png "Flaggged Scatter Plot")

#### Bland_altman Plots
Bland altman plots are useful for determining size effects on data variance. If there is a size effect, one would expect greater differences between measurements as their mean increases. This can be helpful for determining if an upper limit cutoff is necessary.
These plots are provided in logfold and in untransformed versions.

#### Untransformed
![Bland Altman Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Bland_altman_unflagged.png "Bland Altman Plot")

#### Transformed
![Bland Altman Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Bland_alt_log_unflagged.png "Bland Altman Tranformed Plot")

![Bland Altman Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Bland_altman_log_flagged_plot.png "Bland Altman Tranformed Plot")


#### Gene Count histogram Plots
RNA-seq commonly produces measurements for many loci with very low counts. It is generally important to filter out low count loci in order to not bias differential expression analysis by incorporating low count loci during normalization. With unreplicated data, it is imporrtant to filter out loci with little evidence. Gene count histograms will provide a measure of how many genes are found at a given count number, thus providing a metric for determining a lower cutoff threshold for the filter.

![Gene Count Histogram Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Gene_count_histogram_plot.png "Gene Count Histo Plot")

#### Logfold change histogram Plots
After parameter optimization, its useful to examine other characteristics of your RNA-seq data. Here we provide logfold change histograms showing flagged genes by their fold change. This can give you insight in to how many flagged genes were up-regulated vs down-regulated. Two plots are provided, per-sample plots and a combined fold change histogram plot.

#### Per Sample
![LogFold change histogram Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/logfold_histogram_plot.png "Logfold change histogram Plot")

#### Combined
![LogFold change histogram Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Logfold_change_histogram_plot2.png "Logfold change histogram Plot")

## TO DO !!!
There are a few things left to do for anyone who is interested in getting involved in this small project.

1. Fully integrate GO-term optimization
I have a clumsy script right now that will take in some gene ontology terms and and a grid of parameters and perform an exhaustive search across the parameter grid looking for those parameters that return one or more of your selected go-terms. Its VERY slow. But it works! Much of the mechanics are built in to the data container and data analyzer classes, however there are still bugs in this current version (due to to random changes made by in the moment necessity) and it needs to be integrated in to the full program architecture (its currently an add on). Or if another design might work better - for example right now its a separate script that runs the program over and over again - perhaps this is the best way to implement this tool. Either way, the script needs to be able to take command line arguments and then function properly.

2. Automatic dispersion estimation
This is the final piece to the SeqPyPlot puzzle. The limitation with unreplicated RNAseq data is there is no way to perform stats tests to find differentially expressed genes. However, we can combine the filtering provided by SeqPyPlot with dispersion estimates across a subset of genes to get a rough dispersion estimation would could then be used to test for DE genes. This would require to parts. First, we have to filter genes across control and experimental samples. Next, since we are dealing with time series data, we'd need to filter genes that are flagged across time within the control sample. Any genes that aren't flagged in either scenario are likely stable genes that we could then use for dispersion estimates (using control and mutant, and time series samples as replicates). From here, we could implement a statistical test (perhaps from edgeR) to determine differentially expressed genes.

3. GUI
I have started to build a GUI for SeqPyPlot using the Python Tkinter module. Ideally SeqPyPlot should be accessible for people who don't use command line tools. The GUI would make this program very easy for anyone to run with or without command line experience. The GUI working draft is located in SeqPyPlotLib/dev

## Credits
The current and previous versions of SeqPyPlot were built by me (Paul G). I would LOVE to add more credits! Fork this repository and help make SeqPyPlot a more useful and user friendly tool!

## License
This software is currently open source, and I'm looking to write a short paper on the utility of this software to submit to an open source journal. IF/when that happens, citations may be in order. Maybe. Probably not worth worrying about though!




