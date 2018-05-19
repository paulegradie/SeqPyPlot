
#             ***SPPLOT v0.4***
A tool for helping you analyze pilot data (data without replicates). Previously called SeqPyPlot.

## Installation

#### Dependencies
There are several dependencies required for SPPLOT. If you are a macosx user, you may already have some of these dependencies. For Windows users, I strongly suggest installing the Anaconda python distribution. You can find these distributions here:

#### Install
Use of SPLLOT requires installation. To install seqpyplot, run:

`python setup.py install`

To test the install, run the execution script using the example config and example data:

`python bin/SPPLOT.py -c bin/example_config.ini --plot --overwrite`

The examples/test_output directory should be overwritten with the same results.

## Usage
    
<pre>
$ python bin/SPPLOT.py -h
usage: SPPLOT v0.4 [-h] -c CONFIG [-o] [-i] [-p] [-t] [-u]

SPPLOT only requires a path to a configuration file. The config file is
expected to be in .ini format. See the examples/ directory for an example
config file. NYI indicates not-yet-implemented.

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIG, --config CONFIG
  -o, --overwrite       Overwrite existing output directory.
  -i, --impute          Impute missing data using neighbors: NYI
  -p, --plot            Create all plots.
  -t, --correct         Attempt to correct for heteroskedasticity: NYI
  -u, --unnorm          Do not normalize data.
</pre>

## Example config file contents
<pre>
[data_directory]
dir=C:\\Users\\paule\\code\\SeqPyPlot\\examples\\htseqCounts
output=C:\\Users\\paule\\code\\SeqPyPlot\\examples\\test_output

[data]
controls=['D1_Control.counts', 'D2_Control.counts', 'D3_Control.counts']
treated=['D1_Treated.counts', 'D2_Treated.counts', 'D3_Treated.counts']
data_type=htseq

[names]
controls=['D1_Cont', 'D2_Cont', 'D3_Cont']
treated=['D1_Treat', 'D2_Treat', 'D3_Treat']
times=['D1', 'D2', 'D3']
conditions=['Control', 'Treated']
experiment_name=Example

[params]
log2fold=2.0
low=2
hi=10000
diff=[200, 1000000]

[plot_options]
scatrange=[10, 10000]

[file_names]
genelist=C:\\Users\\paule\\code\\SeqPyPlot\\examples\\testgenelist.txt
prefix='SeqPyPlot'
</pre>

## Contributing
1. Fork it!
2. Create your feature branch: `git checkout -b my-new-feature-branch`
3. Commit your changes: `git commit -m 'Add some feature'`
4. Push to the branch: `git push origin my-new-feature-branch`
5. Submit a pull request :D

## History and background
SPPLOT (rna-Seq Python PLOT library) was written to provide a simple to use, quick analysis tool for interrogating pilot RNA-seq data. Pilot study data is defined here as data that doesn't include replicates. A common need among developmental biologists right now is the ability to interogate pilot study data and access gene expression quickly. Usually these data are not used for publishing, but rather to guide future experiments and provide a quick peak at approximate expression levels within the target tissue. 

Quoted from the edgeR manual (a popular tool written in R used for calculating differential expression of genes:

"__RNA-Seq and ChIP-Seq are still expensive technologies, so it sometimes happens that only one library can be created for each treatment condition. In these cases there are no replicate libraries from which to estimate biological variability. In this situation, the data analyst is faced with the following choices, none of which are ideal. We do not recommend any of these choices as a satisfactory alternative for biological replication. Rather, they are the best that can be done at the analysis stage [using no replicates].__```

The manual goes to suggest the following option (which is admittedly not their preferred suggestion...)

`1. Be satisfied with a descriptive analysis, that might include an MDS plot and an analysis
of fold changes. Do not attempt a significance analysis. This may be the best advice.`

The motivation behind SPPLOT is to provide a means to analyze time series data under the assumption of no biological variability (even though we know this isn't true). With unreplicated time series data, there is no way to reasonably estimate variance within the data, since any gene could potentially change over time. SPPLOT therefore ignores dispersion estimates and biological variability to give a quick and dirty assessment of fold change. 

The advantage of using SPPLOT is twofold. First, it provides a set of optimizable parameters that can be used to filter data included in the final output. Second, it provides a very useful plotting function for rapidly visualizing gene expression. Using this program, you can create a list of flagged genes, and plot them within minutes.

## Plot Descriptions
SPPLOT produces a range of plots for general descriptive analysis of your data. These plots are available using the 'plot options' shown above in the usage. This section provides plot examples along with suggested interpretation.

#### Expression Plots
The main plot time intended for SPPLOT users are the expression plots. These plots are not entirely interpretable and should be used only to generalize about the data. The rason is as follows. When we normalize RNA-seq expression data, we can use one of two methods -

1.Between-Sample normalization (used for testing for DE)
e.g.
- TMM (edgeR)
- Geometic (DESeq2)

2. Within-Sample normalization (for comparing genes within a single sample).
e.g.
- RPKM (single reads)
- FPKM (Paired end reads)
- CPM, TPM, Quartile, etc (Generic non genelength correcting methods)

Mixing these two forms is not well understood, so when we want to look at global quantification of RNA-seq across many samples... what do we do? If you look through the literature, you will likely not find expression plots using TMM normalized values, however we know that due to common differences in library size after library prep, RPKM/FPKM/CPM, (etc) values are not directly comparable.
The best approach we have right now is to use a relative normalization method such as TMM and explicitly say that the values shown are not absolutely quantified.

The expression plots output from SPPLOT use TMM normalized values and show relative difference between values across time/stages. SPPLOT provides a secondary measure of relative expression, akin to a scale bar, that is intended to help users interpret the relative differences between expression values in log space. SPPLOT calculates the mean between two measurements at a given stage and then calculates the log fold change provided by the user and plots error bars centered on the mean with the range of the logfold difference (default 0.7). The output also highlights (in yellow) any plot for a gene that is flagged at any time point.

![MDSPlot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Expression_plot.png "Expression Plot")

#### PCA Plots
After normalization, it is essential to verify that samples cluster appropriately. Typically, RNA-seq datasets from developmental time series data will contain a great deal of variation between time points, but not between samples within a time point. Of course this is not always the case. Nonetheless, observing relatively close clustering of samples from the same stage in a PCA plot is genarally desirable. If such clustering is not observed, then this may indicate a significant problem with sample composition and invalidate your experiment.

![MDSPlot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/PCA_plot.png "PCA Plot")

#### Tally Plots
This plot provides feedback on your filter parameter selection. These paramters determine the number of genes that will be flagged. Perhaps the most important parameter is the logfold change threshold. This plot shows the effect of tuning this threshold by counting the number of flagged genes across a range of threshold values, while keeping the other filter parameters constant. A good value to choose for logfold change will be around the point in the plot that the line begins to steepen.

![Tally Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Tally_plot.png "Tally Plot")

#### Flagged Gene Bar Plots
Bar plots are provided to give information on the number of genes flagged at each stage given a set of filter parameters.

![Bar Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Bar_plot.png "Bar Plot")

### Plots for Adjusting Filter Parameters
The following plots are intended to aid the user in finding appropriate filter parameters. The filter takes four parameters illustrated in the following image.

#### --Filter illustration--

The paramters used to filter genes are shown in the figure below:

![Filter Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Filter.png "Filter")

#### Scatter Plots
Scatter plots are useful for double checking that normalization was successful and that your fold change parameter is reasonable across samples.
There are two plots - flagged and un-flagged. The white line along the $f(x) = x$ vector indicates the center of the value distribution across all means, and all expressin values ar plotted against their relative mean. In other words, the mean is calculated for value pairs (control and experimental), and then each vlue is plotted against its mean. The dotted lines represents the log2fold change (from the user provided log2fold argument, default = 0.7) around the mean.

![Scatter Plot](https://github.com/paulgradie/SeqPyPlot/blob/master/examples/images/Scatter_plot.png "Scatter Plot")


## Credits
The current and previous versions of SPPLOT were built by me (Paul G). I would LOVE to add more credits! Fork this repository and help make SPPLOT a more useful and user friendly tool!

## License
Copyright 2018 Paul E. Gradie

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
