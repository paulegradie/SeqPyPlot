import argparse
import sys


class Args:

    def __init__(self):
        self.args = self.arg_parser()
        self.args.time = self.time_parser()
        self.args.condition = self.condition_label_parser()

    def arg_parser(self):
        print("\nDefaults are shown in usage.\n")
        # type: () -> NameSpace()
        parser = argparse.ArgumentParser(description='Plot time series expression data without replicates. \n')

        # plotter args
        parser.add_argument('-time', '-t',
                            metavar='1,2,3,4,5',
                            default='None',
                            type=str,
                            dest='time',
                            help='\tA comma separated list of time points.')

        parser.add_argument('-folder_name', '-out',
                            metavar='Default_out',
                            default='Default_out',
                            type=str,
                            dest='out',
                            help='\tOutput Folder Name')

        parser.add_argument('-datatype', '-d',
                            metavar='htseq',
                            default='htseq',
                            type=str,
                            dest='datatype',
                            help='\tType of data. cuffnorm, cuffdiff, deseq2, and edgeR are supported.')

        parser.add_argument('-c', '--condition',
                            metavar='Series1,Series2',
                            default=('Series1,Series2'),
                            type=str,
                            dest='condition',
                            help='\tA comma separated list of conditions (max 2)')

        parser.add_argument('-e', '--err',
                            metavar='0.5',
                            default=None,
                            type=float,
                            dest='err',
                            help='\tDepricated. Range around the mean to visualize.')

        parser.add_argument('-fr', '--filter_results',
                            metavar='None',
                            default=None,
                            type=str,
                            dest='filter_results',
                            help='\tIf you implement your own filter. Optional.')
        ## Filter args
        parser.add_argument('-fpkm',
                            metavar='0',
                            default=1,
                            type=float,
                            dest='fpkm',
                            help='Default: 2. Set the minimum FPKM value to accept.')

        parser.add_argument('-log2',
                            '--log2fold',
                            metavar='0',
                            default=1.0,
                            type=float,
                            dest='log',
                            help='Default: 0.5. Set the minimum log2Fold change to accept.')

        parser.add_argument('-num',
                            metavar='2',
                            default=2,
                            type=int,
                            dest='num',
                            help='Default: 2. Set number of plots.')

        parser.add_argument('-r',
                            '--remove',
                            default=False,
                            action='store_true',
                            dest='remove',
                            help='Default: False. Use to remove genes not always on.')

        parser.add_argument('-tally',
                            action='store_true',
                            default=False,
                            dest='tally',
                            help='Default: False. Tally DE genes by time and sum DE count - write to separate file.')

        parser.add_argument('-report',
                            action='store_true',
                            default=False,
                            dest='report',
                            help='Default: False. Write out plot data and filter results files.')

        parser.add_argument('-ercc',
                            action='store_true',
                            default=False,
                            dest='ercc',
                            help='Default: False. Write ERCC data to an output file (Name defined by -o option.')

        parser.add_argument('-file_prefix',
                            nargs='?',
                            type=str,
                            default='',
                            metavar='SeqPyPlot_out',
                            dest='figure_name',
                            help='Leading name of output file.')

        parser.add_argument('-raw_data', '-raw',
                            nargs='?',
                            type=str,
                            default=None,
                            metavar='None',
                            dest='raw_data',
                            help='Input File (Cuffnorm output).')

        parser.add_argument('-plot_data',
                            nargs='?',
                            type=str,
                            default=None,
                            metavar='None',
                            dest='plot_data',
                            help='Formatted input data to plot')

        parser.add_argument('-gene_list',
                            nargs='?',
                            type=str,
                            default=None,
                            metavar='None',
                            dest='gene_list',
                            help='\tSingle Column Gene list in txt file.')

        parser.add_argument('-input_results',
                            nargs='?',
                            type=str,
                            default=None,
                            metavar='None',
                            dest='de_results',
                            help='\tSingle Column DE Gene list in txt file.')

        return parser.parse_args()

    @staticmethod
    def __label_parser(argument):
        # type: (character_string) -> list of strings
        try:
            parsed_list = [x for x in argument.split(',')]
            return parsed_list

        except AttributeError:
            print("The group labels have to be comma separated.")
            sys.exit()

    def time_parser(self):
        return self.__label_parser(str(self.args.time))

    def condition_label_parser(self):
        return self.__label_parser(str(self.args.condition))

    def make_logs(self):
        with open((str(self.args.out) + ".log.txt"), "w+") as logfile:
            logfile.write(
                "\nLogs for {0}:".format(str(self.args.out)))
            logfile.write(
                "\nFPKM threshold is: " + str(self.args.fpkm) + "\n")


            logfile.write(
                "\nLog2Fold threshold is: " + str(self.args.log) + "\n")


            logfile.write(
                "\nUser opted to remove non-constitutively active genes: " + str(self.args.remove) + "\n")


            logfile.write(
                "\nUser opted to remove non-constitutively active genes: \n" + str(namespace()) + "\n")
