import argparse
import sys
import os


class Args:
    def __init__(self):
        self.args = self.arg_parser()

        if self.args.time is not None:
            self.args.time = self.time_parser()
        else:
            print 'Set -time argument'
            sys.exit()


        if self.args.num == 1:
            self.args.condition = ["Gene"]
        elif self.args.num == 2:
            self.args.condition = self.condition_label_parser()
        else:
            print 'Cant support more than two yet.'

        if not os.path.isdir(self.args.out):
            os.makedirs(self.args.out)

        self.path = os.path.join(self.args.out, self.args.prefix)

    def arg_parser(self):
        # type: () -> NameSpace()
        usage = "python SeqPyPlot.py [options] [-raw_data (or) -plot_data] -num 2 -time d1,d2,d3"
        epilog = """\n
        If you are finding that your file names aren't printing in the correct order,\n
        you can rename your files with..\n
        \n
            1_[filename]\n
            2_[filename]\n

        ...to produce the correct order.\n

        order = Chronilogicaly staggered:\n

        option -num set to 1: d1, d2, d3, ...\n
        option -num set to 2: d1-ctrl, d1-exp, d2-ctrl, d2-exp, ...

        """


        parser = argparse.ArgumentParser(description='Required: -raw_data or -plot_data',
                                         prog='SeqPyPlot v0.2',
                                         usage=usage,
                                         epilog=epilog)
        # general args\
        parser.add_argument('--------------------General Options--------------------',
                            action='store_true',
                            default=False)
        parser.add_argument('-time',
                            metavar='d1,d2,d3',
                            default=None,
                            type=str,
                            dest='time',
                            help='A comma separated list of time points.')

        parser.add_argument('-num',
                            metavar='2',
                            default=2,
                            type=int,
                            dest='num',
                            help='Default: 2. Set number of plots.')

        parser.add_argument('-out',
                            metavar='Default_out',
                            default='Default_out',
                            type=str,
                            dest='out',
                            help='Output Folder Name')

        parser.add_argument('-prefix',
                            type=str,
                            default='SeqPyPlot_default',
                            metavar='SeqPyPlot_',
                            dest='prefix',
                            help='Leading name of output file.')

        parser.add_argument('-data_type',
                            metavar='htseq',
                            default='htseq',
                            type=str,
                            dest='datatype',
                            help='Either cuffnorm, cuffdiff, deseq2, or edgeR.')

        parser.add_argument('-c',
                            metavar='S1,S2',
                            default=None,
                            type=str,
                            dest='condition',
                            help='\tA comma separated list of conditions (max 2)')

        parser.add_argument('-unform_plot_data',
                            action='store_true',
                            default=False,
                            dest='unformatted_plot_data',
                            help='Default: False. Reorder plot data (1,1,2,2 -> 1,2,1,2.')


        ## Filter args

        parser.add_argument('--------------------Filter Options--------------------',
                            action='store_true',
                            default=False)
        parser.add_argument('-low',
                            metavar='0',
                            default=25,
                            type=float,
                            dest='low',
                            help='Default: 2. Set the min expression value to accept.')
        parser.add_argument('-hi',
                            metavar='5000000',
                            default=5000000,
                            type=float,
                            dest='hi',
                            help='Default: 5mil. Set the max expression value to accept.')
        parser.add_argument('-dif',
                            metavar='60',
                            default=60,
                            type=float,
                            dest='dif',
                            help='Default: 60. Set minimum difference in expression.')
        parser.add_argument('-dif_upper',
                            metavar='60',
                            default=100000,
                            type=float,
                            dest='dif_upper',
                            help='Default: 100000. Set minimum difference in expression.')
        parser.add_argument('-log2',
                            metavar='1.0',
                            default=1.0,
                            type=float,
                            dest='log',
                            help='Default: 1.0. Minimum log2 change to accept.')


        #analysis options
        parser.add_argument('--------------------Analysis Options--------------------',
                            action='store_true',
                            default=False)


        parser.add_argument('-r',
                            default=False,
                            action='store_true',
                            dest='remove',
                            help='Default: False. Use to remove genes not always on.')

        parser.add_argument('-svg',
                            default=False,
                            action='store_true',
                            dest='svg',
                            help='Default: False. Use to svg plots.')

        parser.add_argument('-tally',
                            action='store_true',
                            default=False,
                            dest='tally',
                            help='Default: False. Tally DE genes.')

        parser.add_argument('-hist_range',
                            metavar='1,1000',
                            default='1,1000',
                            type=str,
                            dest='hist_range',
                            help='Default: 1.0. Lower x axis limit for histogram.')

        parser.add_argument('-report',
                            action='store_true',
                            default=False,
                            dest='report',
                            help='Default: False. Write plot data and filter results.')

        parser.add_argument('-ercc',
                            action='store_true',
                            default=False,
                            dest='ercc',
                            help='Default: False. Write ERCC data to an output file.')

        parser.add_argument('--------------------Input Options--------------------',
                            action='store_true',
                            default=False)

        parser.add_argument('-raw_data',
                            type=str,
                            default=None,
                            metavar='None',
                            dest='raw_data',
                            help='Input file or folder.')

        parser.add_argument('-plot_data',
                            type=str,
                            default=None,
                            metavar='None',
                            dest='plot_data',
                            help='Formatted input data to plot')

        parser.add_argument('-gene_list',
                            type=str,
                            default=None,
                            metavar='None',
                            dest='gene_list',
                            help='\tSingle Column Gene list in txt file.')

        parser.add_argument('-input_results',
                            type=str,
                            default=None,
                            metavar='None',
                            dest='de_results',
                            help='Optional. Your own flagged gene list.')




        return parser.parse_args()

    @staticmethod
    def __label_parser(argument):
        # type: (character_string) -> list of strings
        try:
            parsed_list = argument.split(',')
            return parsed_list

        except AttributeError:
            print("The group labels have to be comma separated.")
            sys.exit()

    def time_parser(self):
        return self.__label_parser(str(self.args.time))

    def condition_label_parser(self):
        return self.__label_parser(str(self.args.condition))

    def make_logs(self):

        with open((str(self.path) + ".log.txt"), "w+") as logfile:

            logfile.write(
                "\nLogs for {0}:".format(str(self.args.prefix)))

            logfile.write(
                "\nExpression thresholds is: \n-Upper expression limit: {}\n-Lower expression limit: {}\n-Minimum Difference: {}\n".format(
                    self.args.hi, self.args.low, self.args.dif))

            logfile.write(
                "\nLog2Fold threshold is: " + str(self.args.log) + "\n")

args = Args()
