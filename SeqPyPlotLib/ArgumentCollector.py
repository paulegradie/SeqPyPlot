import argparse
import sys
import os
import time

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

        self.args.dif_range = self.dif_parser()

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
        parser.add_argument('--------------------Required Options--------------------',
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
        parser.add_argument('-c',
                            metavar='S1,S2',
                            default='S1,S2',
                            type=str,
                            dest='condition',
                            help='\tA comma separated list of conditions (max 2)')
        parser.add_argument('-data_type',
                            metavar='htseq',
                            default='htseq',
                            type=str,
                            dest='datatype',
                            help='Either cuffnorm, cuffdiff, deseq2, or edgeR.')

        parser.add_argument('--------------------Required Input Options--------------------',
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
        parser.add_argument('-unform_plot_data',
                            action='store_true',
                            default=False,
                            dest='unformatted_plot_data',
                            help='Default: False. Reorder plot data (1,1,2,2 -> 1,2,1,2.')
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
        ## Filter args

        parser.add_argument('--------------------Filter Options--------------------',
                            action='store_true',
                            default=False)
        parser.add_argument('-low',
                            metavar='25',
                            default=25,
                            type=float,
                            dest='low',
                            help='Default: 2. Set the min expression value to accept.')
        parser.add_argument('-hi',
                            metavar='1000000',
                            default=1000000,
                            type=float,
                            dest='hi',
                            help='Default: 5mil. Set the max expression value to accept.')
        parser.add_argument('-dif_range',
                            metavar='25,10000',
                            default='25,1000000',
                            type=str,
                            dest='dif_range',
                            help='Default: 25-1000000. Set min difference in expression.')
        parser.add_argument('-log2',
                            metavar='0.7',
                            default=0.7,
                            type=float,
                            dest='log',
                            help='Default: 0.7. Minimum log2fold change to accept.')
        #analysis options
        parser.add_argument('--------------------Analysis Options--------------------',
                            action='store_true',
                            default=False)
        parser.add_argument('-find_housekeeping',
                            action='store_true',
                            default=False,
                            dest='find_housekeeping',
                            help='Default: False. Search for housekeeping genes.')
        parser.add_argument('-num_housekeeping',
                            metavar='3',
                            default=3,
                            type=int,
                            dest='num_housekeeping',
                            help='Default: 3. Min number of Housekeeping to detect.')
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
        parser.add_argument('-r',
                            default=False,
                            action='store_true',
                            dest='remove',
                            help='Default: False. Use to remove genes not always on.')
        #analysis options
        parser.add_argument('--------------------Plot Options--------------------',
                            action='store_true',
                            default=False)
        parser.add_argument('-plots',
                            action='store_true',
                            default=False,
                            dest='plots',
                            help='Default: False. Make all plots.')
        parser.add_argument('-tally',
                            action='store_true',
                            default=False,
                            dest='tally',
                            help='Default: False. Tally DE genes.')
        parser.add_argument('-bar',
                            action='store_true',
                            default=False,
                            dest='bar',
                            help='Default: False. Construct bar plots.')
        parser.add_argument('-scatter',
                            action='store_true',
                            default=False,
                            dest='scatter',
                            help='Default: False. Construct scatter plots.')
        parser.add_argument('-scat_range',
                            metavar='0,1200',
                            default='0,1200',
                            type=str,
                            dest='scatt_range',
                            help='Default: 0,5000000. Set scatter plot value range.')
        parser.add_argument('-bland_alt',
                            action='store_true',
                            default=False,
                            dest='bland_alt',
                            help='Default: False. Construct bland-altman plots.')
        parser.add_argument('-barange',
                            metavar='0,10000',
                            default='20000',
                            type=str,
                            dest='ba_range',
                            help='Default: 1,25000. Range for BA plots.')
        parser.add_argument('-bland_log',
                            action='store_true',
                            default=False,
                            dest='bland_log',
                            help='Default: False. Construct bland-alt-log plots.')
        parser.add_argument('-blrange',
                            metavar='0,10000',
                            default='0,10000',
                            type=str,
                            dest='bl_range',
                            help='Default: 1,25000. Range for BG plots.')
        parser.add_argument('-histo',
                            action='store_true',
                            default=False,
                            dest='histo',
                            help='Default: False. Construct histogram plots.')
        parser.add_argument('-hist_range',
                            metavar='1,1000',
                            default='1,1000',
                            type=str,
                            dest='hist_range',
                            help='Default: 1.0. Lower x axis limit for histogram.')
        parser.add_argument('-log2histo',
                            action='store_true',
                            default=False,
                            dest='log2histo',
                            help='Default: False. Construct log2fold histogram plots.')
        parser.add_argument('-all',
                            action='store_true',
                            default=False,
                            dest='all',
                            help='Default: False. Make plots and tally flagged genes.')
        parser.add_argument('-svg',
                            default=False,
                            action='store_true',
                            dest='svg',
                            help='Default: False. Use to svg plots.')
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

    def dif_parser(self):
        dif_list = self.__label_parser((str(self.args.dif_range)))
        return [float(x) for x in dif_list]

    def make_logs(self):
        with open((str(self.path) + ".log.txt"), "w+") as logfile:

            logfile.write(
                "\nLogs for {0}:".format(str(self.args.prefix)))

            logfile.write(
                "\nExpression thresholds: \n-Upper expression limit: {}\n-Lower expression limit: {}\n-Difference Range: {}\n".format(
                    self.args.hi, self.args.low, self.args.dif_range))

            logfile.write(
                "\nLog2Fold threshold is: " + str(self.args.log) + "\n")

            logfile.write(
                "\nScript Paramters-----\n")
            for i in range(len(sys.argv)):

                logfile.write(sys.argv[i] + ' ')
            logfile.write('\n\n' + time.ctime())

