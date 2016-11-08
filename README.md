## Synopsis

SeqPyPlot is a program for plotting time series data. Data is formatted in to a tab delimeted text file and then fed in to the main script. The script accesses code stored in SeqPyPlotLib and generates plots for the input data, saving them to .png and .svg files. 

## Motivation

The program was originally written to plot output from the Cufflinks suite. A pipeline was created where RNA-seq data would be mapped, quantified, and then normalized together to allow for comparisons between samples. This ouptput was then formatted and fed in to SeqPyPlot to plot expression values for lists of genes. 

## Installation

You can either clone the repository, download the scripts, or copy the code in to new files.

## Dependencies

This program imports Numpy and MatPlotLib.

## Contributors

I am currently the only contributor on this project, but if you'd like to contribute, by all means please do. Let me know so I can include your name here if you wish.

## Code Example - Typical Usage


 >$ python SeqPyPlot.py -t e12,e13,e14,e15 -o Dev_genes -c Control,Mutant genelist.txt plotter_data.txt


Output will be a folder named "Dev_genes" with files containing plot with a maximum of six plots per file.

## Full usage

First make sure that all of the programs files are in the same directory. Then run the script:

![Alt text](images/Full_run.png?raw=true "Full Usage")
 
 
### Data Prep 

Format your data as shown below in excel, or similar program, and save as a Tab-deliminted .txt file.
![Alt text](images/DataFormat1.PNG?raw=true "Full Usage")


The gene list should be formatted in a simlar way. One gene name per line.

![Alt text](images/GeneList.PNG?raw=true "Full Usage")


## Tests and Sample Output

In the Tests folder, there is located a genelisttest.txt and a testplotdata.txt file. Using these files with the options at left to default, the following output is produced:

![Alt text](images/1_[u'Mdk',u'Mdm1',u'Mdm2',u'Mdm4',u'Mdn1',u'Mdp1'].png?raw=true "Full Usage")
![Alt text](images/2_[u'Me1',u'Me2',u'Me3',u'Mea1',u'Meaf6',u'Mecom'].png?raw=true "Full Usage")

## The run Script

If you'd like to edit the run script, SeqPyPlot.py, you should aim to first read in any command line arguments using the ArgumentParser class, then create an output directoryu using the PrepareOutputDirectory class, followed by reading in the data usnigthe  ReadNormalized Data class, and finally producing the plots using the MainDataPlotter class.
 
## License

This is free to use OpenSource software. Public Domain.
