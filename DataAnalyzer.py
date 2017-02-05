import numpy as np
import csv


class DataAnalyzer(object):
    def __init__(self, args, data_container):
        self.args = args
        self.de_gene_list_length = 0

        self.de_count_by_stage = dict()  # for bar graph
        self.filtered_data = dict()
        self.de_count_by_gene = dict()
        self.de_gene_list_by_stage = dict()
        self.de_gene_list = []
        self.data_container = data_container
        # inherited
        self.gene_map = data_container.gene_map
        self.ercc_map = data_container.ercc_map
        self.data_frame_header = data_container.data_frame_header

    def seqpyfilter(self, iterator=None):
        """Running this function drops all of the outputs in to the SELF properties."""
        # args_log = self.args.args.log
        args = self.args.args

        if self.args.time[0] == 'None':
            # print "FALSE"
            labels = [x for x in self.data_frame_header['Gene'] if self.data_frame_header['Gene'].index(x) % 2 != 0]
        else:
            # print "TRUE"
            labels = [i for i in self.args.time]

        if iterator is not None:
            args_log = iterator
            self.de_gene_list = []

        else:
            args_log = args.log

        # self.de_gene_list_length = None  # int; Total number of de genes found
        # self.de_count_by_stage = None  # dictionary; keys = time points, values = number of de genes per time point
        # self.de_gene_list = []  # list; only de gene names
        # self.filtered_data = dict()  # dictoinary; all de genes with expression values
        # self.de_count_by_gene = dict()  # dictionary; keys = gene names, values number of time points de

        for key, value in sorted(self.gene_map.items()):
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
                        if float(sub_list[0]) >= args.fpkm or float(sub_list[1]) >= args.fpkm:
                            if abs(float(sub_list[0]) - float(sub_list[1])) >= 3.5:
                                if max(float(sub_list[0]), float(sub_list[1])) <= 5000:
                                    keep = True

                                    # assign results to various maps
                                    if int(v) not in self.de_gene_list_by_stage.keys():
                                        self.de_count_by_stage[v] = 1
                                    else:
                                        self.de_count_by_stage[v] += 1

                                    if int(v) not in self.de_gene_list_by_stage.keys():
                                        self.de_gene_list_by_stage[v] = [key]
                                    else:
                                        self.de_gene_list_by_stage[v].append(key)

                                    if iterator is None:  # avoid value inflation during tally routine
                                        if key not in self.de_count_by_gene.keys():
                                            self.de_count_by_gene[key] = 1
                                        else:
                                            self.de_count_by_gene[key] += 1


                    # Condition 2: Both are non-zero - DO THE LOG TEST
                    elif float(sub_list[0]) >= args.fpkm or float(sub_list[1]) >= args.fpkm:
                        if abs(np.log2(float(sub_list[1]) / float(sub_list[0]))) >= args_log:
                            # print abs(np.log2(float(sub_list[1]) / float(sub_list[0])))

                            if abs(float(sub_list[0]) - float(sub_list[1])) >= 3.5:
                                if max(float(sub_list[0]), float(sub_list[1])) <= 5000:
                                    keep = True

                                    # assign results to various maps
                                    if int(v) not in self.de_gene_list_by_stage.keys():
                                        self.de_count_by_stage[v] = 1
                                    else:
                                        self.de_count_by_stage[v] += 1

                                    if int(v) not in self.de_gene_list_by_stage.keys():
                                        self.de_gene_list_by_stage[v] = [key]
                                    else:
                                        self.de_gene_list_by_stage[v].append(key)

                                    if iterator is None:  # avoid value inflation during tally routine
                                        if key not in self.de_count_by_gene.keys():
                                            self.de_count_by_gene[key] = 1
                                        else:
                                            self.de_count_by_gene[key] += 1

            if keep is True:
                # TODO consolidate these two variables to clean up
                self.de_gene_list.append(key)
                self.filtered_data[key] = value

        self.de_gene_list_length = int(len(self.de_gene_list))

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

    def print_de_tallies(self):
        self.args = self.args.args
        # For log2fold parameter iterative testing
        with open(self.args.out + '_DE_counts_per_gene.txt', 'wb') as genecount:
            if len(self.de_count_by_gene) != 0:
                de_count_writer = csv.writer(genecount, delimiter='\t')
                de_count_writer.writerow(
                    [str("Log2threshold")] + [str('Number of DE genes detected at log2: '.format(self.args.log))])
                for key, value in sorted(self.de_count_by_gene.items()):
                    de_count_writer.writerow([key] + [int(value)])
            else:
                print "No DE genes to count for DE_counts_per_gene."

        # Data for BAR graph and Gene Ranking
        with open(self.args.out + '_DE_count_by_stage.txt', 'wb+') as genecount:
            if len(self.de_count_by_stage) != 0:
                de_count_writer = csv.writer(genecount, delimiter='\t')
                de_count_writer.writerow([str("Time point")] + [str('Number of DE genes detected')])
                for key, value in sorted(self.de_count_by_stage.items()):
                    de_count_writer.writerow([int(key)] + [int(value)])
            else:
                print "No DE Genes to write for DE counts by stage."

        # TODO Currently prints horizontal gene lists per stage. Use zip and collections (or Pandas) to print verticle
        # All genes DE per  Stage
        with open(self.args.out + '_DE_gene_lists_by_stage.txt', 'wb+') as genecount:
            if len(self.de_gene_list_by_stage) != 0:
                de_list_writer = csv.writer(genecount, delimiter='\t')
                de_list_writer.writerow(["Stage"] + ["DE_gene_list"])
                for key, value in self.de_gene_list_by_stage.items():
                    de_list_writer.writerow([key] + value)
            else:
                print "No DE Genes to write for gene_list by stage."
