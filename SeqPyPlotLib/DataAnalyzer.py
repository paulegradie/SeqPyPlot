import numpy as np
import csv
import os

class DataAnalyzer(object):
    def __init__(self, args, data_container, optimization_mode=False):
        if optimization_mode is True:
            self.temp = args
            self.args = data_container.args

            self.data_frame_header = data_container.data_frame_header
            self.gene_map = data_container.gene_map

            self.de_gene_list = []  # A simple list of all DE genes
            self.de_genes_by_stage = [] # Map with stages as key and simple lists of DE genes as values

            self.prefix = None
            self.log = float()
            self.dif = float()
            self.hi = int()
            self.low = int()

        else:
            self.args = args
            self.gene_map = data_container.gene_map
            self.ercc_map = data_container.ercc_map
            self.data_frame_header = data_container.data_frame_header
            self.data_container = data_container

        self.sing_comp_map = dict()
        self.sing_comp_header = dict()

        self.de_gene_list_length = 0

        self.de_count_by_stage = dict()  # for bar graph
        self.filtered_data = dict()
        self.de_count_by_gene = dict()
        self.de_gene_list_by_stage = dict()
        self.de_gene_list = []
        # inherited

        if self.args.out is not None:
            self.path = os.path.join(self.args.out, self.args.prefix)
        else:
            self.path = os.path.join('.', self.args.prefix)

    def __getitem__(self, key):
        returnable = self.args[key]
        return returnable


    @staticmethod
    def make_comparisons(value_list):
        """O(n^2) alogirithm for making all comparisons"""

        first_half = []
        second_half = []
        while len(value_list) > 1:
            for value in value_list[1:]:
                first_half.append(value_list[0])
                second_half.append(value)
            value_list = value_list[1:]
        return first_half + second_half


    def seqpyfilter(self, use_iterator=None):

        """Running this function drops all of the outputs in to the SELF properties."""

        # args_log = self.args.log
        args = self.args
        original_map = dict()
        if self.args.num > 1:
            labels = self.args.time
        else:
            self.sing_comp_header = []
            header = self.make_comparisons(self.args.time)
            top = header[:len(header) / 2]  # split the data
            bottom = header[len(header) / 2:]
            for i in range(len(top)):
                self.sing_comp_header += [str(top[i]) + '/' + str(bottom[i])]
            labels = self.sing_comp_header

        if use_iterator is not None:
            args_log = use_iterator
            self.de_gene_list = []
            original_map = self.gene_map  # go to the end and reset gene_map to original map.

        else:
            args_log = args.log

        # self.de_gene_list_length = None  # int; Total number of de genes found
        # self.de_count_by_stage = None  # dictionary; keys = time points, values = number of de genes per time point
        # self.de_gene_list = []  # list; only de gene names
        # self.filtered_data = dict()  # dictoinary; all de genes with expression values
        # self.de_count_by_gene = dict()  # dictionary; keys = gene names, values number of time points de

        analysis_map = dict()

        if self.args.num == 1:
            for key, value in self.gene_map.items():
                analysis_map[key] = self.make_comparisons(value)

        elif self.args.num == 2:
            analysis_map = self.gene_map
        else:
            print "Can't support more than 2 series yet!"
            sys.exit()

        for key, value in sorted(analysis_map.items()):
            # print "Key, Value: ", key, value, "WEIRED"
            # continue
            # Set preconditions
            keep = False

            series1 = value[:len(value) / 2]
            series2 = value[(len(value) / 2):]

            # Main Filter
            for v in range(len(value)/2):
                sub_list = [series1[v], series2[v]]
                if None in sub_list:
                    pass
                else:
                    # Condition 1: FPKM of one or both is zero.
                    # print "SubList: ", sub_list
                    if float(sub_list[0]) == 0 or float(sub_list[1]) == 0:
                        if float(sub_list[0]) >= float(args.low) or float(sub_list[1]) >= float(args.low):
                            if abs(float(sub_list[0]) - float(sub_list[1])) >= args.dif:
                                if abs(float(sub_list[0]) - float(sub_list[1])) <= args.dif_upper:
                                    if max(float(sub_list[0]), float(sub_list[1])) <= args.hi:
                                        keep = True

                                        if use_iterator is None:  # avoid value inflation during tally routine
                                         # assign results to various maps
                                            if labels[v] not in self.de_gene_list_by_stage.keys():
                                                self.de_count_by_stage[labels[v]] = 1
                                            else:
                                                self.de_count_by_stage[labels[v]] += 1

                                            if labels[v] not in self.de_gene_list_by_stage.keys():
                                                self.de_gene_list_by_stage[labels[v]] = [key]
                                            else:
                                                self.de_gene_list_by_stage[labels[v]].append(key)

                                            if key not in self.de_count_by_gene.keys():
                                                self.de_count_by_gene[key] = 1
                                            else:
                                                self.de_count_by_gene[key] += 1


                    # Condition 2: Both are non-zero - DO THE LOG TEST
                    elif float(sub_list[0]) >= float(args.low) or float(sub_list[1]) >= float(args.low):
                        if abs(np.log2(float(sub_list[1]) / float(sub_list[0]))) >= float(args_log):
                            if abs(float(sub_list[0]) - float(sub_list[1])) >= float(args.dif):
                                if abs(float(sub_list[0]) - float(sub_list[1])) <= float(args.dif_upper):
                                    if max(float(sub_list[0]), float(sub_list[1])) <= float(args.hi):
                                        keep = True

                                        if use_iterator is None:  # avoid value inflation during tally routine

                                            # assign results to various maps
                                            if labels[v] not in self.de_gene_list_by_stage.keys():
                                                self.de_count_by_stage[labels[v]] = 1
                                            else:
                                                self.de_count_by_stage[labels[v]] += 1

                                            if labels[v] not in self.de_gene_list_by_stage.keys():
                                                self.de_gene_list_by_stage[labels[v]] = [key]
                                            else:
                                                self.de_gene_list_by_stage[labels[v]].append(key)

                                            if key not in self.de_count_by_gene.keys():
                                                self.de_count_by_gene[key] = 1
                                            else:
                                                self.de_count_by_gene[key] += 1
                    else:
                        pass
                        # print "ERROR - anomoly detected."

            if keep is True:
                # TODO consolidate these two variables to clean up
                self.de_gene_list.append(key)
                self.filtered_data[key] = value

        self.de_gene_list_length = int(len(self.de_gene_list))

        if use_iterator is not None:
            self.gene_map = original_map
        else:
            pass

        return self.de_gene_list, \
               self.de_gene_list_by_stage, \
               self.de_gene_list_length, \
               self.de_count_by_gene, \
               self.de_count_by_stage, \
               self.filtered_data

    # TODO write function to analyze control time series DE genes
    # def analyze_time_series(self, gene_map):
    #     series_1 = dict()
    #     for key, value in gene_map.items():
    #         control_data

    def print_analyzer_results(self):
        self.args = self.args

        # For log2fold parameter iterative testing
        with open(self.path + '_DE_counts_per_gene.txt', 'wb') as genecount:
            if len(self.de_count_by_gene) != 0:
                de_count_writer = csv.writer(genecount, delimiter='\t')
                de_count_writer.writerow(
                    [str("Log2threshold")] + [str('Number of DE genes detected at log2: '.format(self.args.log))])
                for key, value in sorted(self.de_count_by_gene.items()):
                    de_count_writer.writerow([key] + [int(value)])
            else:
                print "No DE genes to count for DE_counts_per_gene."

        # Data for BAR graph and Gene Ranking
        with open(self.path + '_DE_count_by_stage.txt', 'wb+') as genecount:
            if len(self.de_count_by_stage) != 0:
                de_count_writer = csv.writer(genecount, delimiter='\t')
                de_count_writer.writerow([str("Time point")] + [str('Number of DE genes detected')])
                for key, value in sorted(self.de_count_by_stage.items()):
                    de_count_writer.writerow([str(key)] + [str(value)])
            else:
                print "No DE Genes to write for DE counts by stage."

        # TODO Currently prints horizontal gene lists per stage. Use zip and collections (or Pandas) to print verticle
        # All genes DE per  Stage
        with open(self.path + '_DE_gene_lists_by_stage.txt', 'wb+') as genecount:
            if len(self.de_gene_list_by_stage) != 0:
                de_list_writer = csv.writer(genecount, delimiter='\t')
                de_list_writer.writerow(["Stage"] + ["DE_gene_list"])
                for key, value in self.de_gene_list_by_stage.items():
                    de_list_writer.writerow([key] + value)
            else:
                print "No DE Genes to write for gene_list by stage."
