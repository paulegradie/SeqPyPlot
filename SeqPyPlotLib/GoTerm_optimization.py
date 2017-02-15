# from __future__ import print_function
# This script is useful for optimizing the parameters for SeqPyPlot based on keywords returned from a GOTerm analysis.
# It will also produce a range of plots that are helpful for determining the most useful parameters to use.

from DataPlotter import MainDataPlotter
from DataAnalyzer import DataAnalyzer
import DataContainer
import csv
from goatools.base import download_go_basic_obo
from goatools.base import download_ncbi_associations
from goatools.obo_parser import GODag
from goatools.associations import read_ncbi_gene2go
from goatools.go_enrichment import GOEnrichmentStudy
from time import time
import sys
optimization_mode = True

#TODO replace file with script argv argument
# data = DataContainer.parse_plot_data('OVE442_Analysis_plotter_data.txt')
# analyzer = DataAnalyzer.seqpyfilter(use_iterator=None, gene_map=data)

# Create a temporary Arg class to sub for builtin argslist
class TempArgs:
    def __init__(self, plot_data):

        self.condition = None
        self.datatype = None
        self.err = None
        self.filter_results = None
        self.gene_list = None
        self.num = 2
        self.out = None
        self.raw_data = None
        self.remove = False
        self.report = False
        self.time = [1,2,3,4,5,6]
        self.ercc = False

        self.plot_data = plot_data

        self.tally = True
        # self.de_results = []

        self.prefix = None
        self.log = float()
        self.dif = float()
        self.hi = int()
        self.low = int()


class Parameters:

    def __init__(self):

        self.count = 0
        self.parameter_list = []
        self.name_list = []

        # loww = [float(x / 10.0) for x in range(10, 1000, 110)]
        # upp = [float(x / 10.0) for x in range(5000, 40000, 5000)]
        logg = [float(x) / 10.0 for x in range(6, 20, 2)]
        # diff = [float(x) / 100.0 for x in range(1600, 10000, 700)]

        # for up in upp:
        #     for low in loww:
        #         for dif in diff:
        #             for log in logg:
        #                 self.parameter_list.append((up, low, dif, log))
        #                 self.count += 1
        # for dif in diff:
        for log in logg:
            self.parameter_list.append((None, None, None, log))
            self.count += 1

        # set list of file names based on paramter list (thousands)
        self.name_list = [str(x) + '_' for x in range(len(self.parameter_list))]

def main():
    # Prep Gene Onotology useing GOaTools
    print("\nDownloading Gene Ontology database...\n")
    obo_fname = download_go_basic_obo()

    print("\nDownloading NCBI identifier associations...\n")
    gene2go = download_ncbi_associations()
    geneid2gos_mouse = read_ncbi_gene2go('gene2go', taxids=[10090])

    print("\nOrganizing Ontology...")
    obodag = GODag("go-basic.obo")


    #asign data file
    data = 'OVE442_Analysis_plotter_data.txt'

    # all of the optimization prarameters and file names tou use with temp_args
    parameters = Parameters()
    args = TempArgs(data)

    ##Establiish GeneID conversion dictionary
    TEMPContainer = DataContainer.DataContainer(args, Optimize=True)

    genes = TEMPContainer.gene_map.keys()
    conversion_table = dict()

    print "\nPreparing gene conversion list..."
    with open("gene_conversion_table.txt", 'rb') as gene_table:
        tablereader = csv.reader(gene_table, delimiter='\t')
        for row in tablereader:
            #conversiontable is symbol: entreziD
            if row[1] == 'NA':
                pass
            else:
                if row[2] not in conversion_table.keys():
                    conversion_table[str(row[2]).title()] = int(row[1])
                else:
                    pass

    print "\nDone. Performing optimization tests.\n"

    pos = 0

    param_file = open("paramaters_for_palate.txt", 'wb+')
    paramwriter = csv.writer(param_file, delimiter='\t')
    paramwriter.writerow(['Parameters: None, None, dif, log2'])
    all_params = len(parameters.parameter_list)
    cur_param = 1
    for param_set in parameters.parameter_list:
        print "Param {} of {}...".format(cur_param, all_params)


        #set Constant paramters
        args.prefix = parameters.name_list[pos] + 'optimize '

        # set variable parameters
        # args.hi = str(param_set[0])
        # args.low = str(param_set[1])
        args.log = str(param_set[3])

        args.hi = 2000000
        args.low = 25
        # args.log = 0.6
        # args.dif = str(param_set[2])
        args.dif = 65

        FullContainer = DataContainer.DataContainer(args, Optimize=True)
        Analyzer = DataAnalyzer(None, FullContainer, optimization_mode=True)
        Analyzer.seqpyfilter()
        print("Data analyzed...forming plots...")
        Plot_Builder = MainDataPlotter(args, Analyzer, None)

        Plot_Builder.de_bar('black')
        Analyzer.print_analyzer_results()
        print("Plots made, initiating GOTerm analysis...")

        # Do a GOTERM analysis for gene seet enrichment using 'gene list'

        # Convert gene_list to IDs
        gene_entrezid_list = []
        # print len(Analyzer.de_gene_list), len(conversion_table)
        for gene in Analyzer.de_gene_list:
            try:
                gene_entrezid_list.append(conversion_table[str(gene)])
            except KeyError:
                pass
                # print("Gene {} not found.").format(gene)
        #use gene_id_list for goterm analysis
        goeaobj = GOEnrichmentStudy(  # This does not containy your query genelist yet! Thats whenyou run_study
            conversion_table.values(),  # List of mouse protein-coding genes for input
            geneid2gos_mouse,       # geneid/GO associations
            obodag,                  # Ontologies
            propagate_counts=False,
            alpha=0.05,             # default significance cut-off
            methods=['fdr_bh'])     # defult multipletest correction method

        goea_results_all = goeaobj.run_study(gene_entrezid_list)
        goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < 0.05]

        keep = False
        for record in goea_results_sig:
            if 'palate' in record.name:
                print record.name
                keep = True
                print "FOUND A HIT\n"
            else:
                pass
        if keep is True:
            # print param_set
            paramwriter.writerow([param_set])
        else:
            print "Unusable parameters.\n"
        cur_param += 1

        Plot_Builder.plot_tally()

        check = True
        pos += 1
        endloop = time()

    param_file.close()

    print("Test Done")
if __name__ == "__main__":
    main()