import numpy as np
import csv
import os

class DataAnalyzer(object):
    def __init__(self, args, data_container, optimization_mode=False):
        """
        :param args:
        :param data_container:
        :param optimization_mode: Set usage to
        """

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
        self.unflagged_genes = dict()
        self.de_count_by_gene = dict()
        self.de_gene_list_by_stage = dict()
        self.de_gene_list = []
        self.sing_time_series_data = dict()

        self.foldchange_map = dict()
        self.housekeeping_dict = dict()
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
        analysis_map = dict()

        #in order to analyze housekeeping genes, data must be fed in as single series to avoid ambiguiety
        if self.args.num == 1:
            self.sing_comp_header = []
            header = self.make_comparisons(self.args.time)
            top = header[:len(header) / 2]  # split the data
            bottom = header[len(header) / 2:]
            self.housekeeping_dict = dict()

            for i in range(len(top)):
                self.sing_comp_header += [str(top[i]) + '/' + str(bottom[i])]

            labels = self.sing_comp_header

            for key, value in self.gene_map.items():
                analysis_map[key] = self.make_comparisons(value)
            self.sing_time_series_data = analysis_map

            if use_iterator is None and self.args.find_housekeeping:
                print "Searching for Housekeeping Genes...\n"
                std_dev_marker = 0.0
                floor = 5500

                while len(self.housekeeping_dict) < self.args.num_housekeeping:
                    std_dev_marker += 10.0
                    for key, value in self.gene_map.items():
                        vals = np.asarray([float(x) if x != 0 else 0 for x in value])
                        if np.mean(vals) >= floor:
                            std_dev = np.std(vals)
                            if std_dev <= std_dev_marker:
                                if key not in self.housekeeping_dict.keys():
                                    self.housekeeping_dict[key] = value + [round(std_dev, ndigits=3)] + [floor]
                                    print "Gene Found: ", key, '\t--Mean:\t', np.mean(vals), '\t--stdDev:\t', round(std_dev, ndigits=3), "\tFloor:\t", floor
                                else:
                                    pass
                        else:
                            pass
                    if std_dev_marker >= 50.0:
                        std_dev_marker = 0.0
                        floor -= 250
                    if floor < 500:
                        print "Fewer than {} acceptable Housekeeping Genes Found...\n".format(self.args.num_housekeeping)
                        break

                if len(self.housekeeping_dict) > 0:
                    with open(self.path + '_Housekeeping_Genes.txt', 'wb+') as housekeeping:
                        self.housekeeping_dict["Gene"] = self.args.time + ["Std_Dev"] + ["Floor"]
                        hkwriter = csv.writer(housekeeping, delimiter='\t')
                        hkwriter.writerow(['Gene'] + self.housekeeping_dict['Gene'])
                        for key, value in sorted(self.housekeeping_dict.items(), key=lambda k: k[1][-1]):
                            if key == 'Gene':
                                pass
                            else:
                                hkwriter.writerow([key] + value)

        elif self.args.num == 2:
            labels = self.args.time
            analysis_map = self.gene_map

        else:
            print "Doesn't suppport num > 2  -- DA, line 107"
            sys.exit()



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


        self.foldchange_map["Gene"] = labels

        for key, value in sorted(analysis_map.items()):
            if key not in self.foldchange_map.keys():
                self.foldchange_map[key] = [0] * (len(value) / 2)
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
                    # Condition 1: Either is zero.
                    if float(sub_list[0]) == 0 or float(sub_list[1]) == 0:
                        if float(sub_list[0]) >= float(args.low) or float(sub_list[1]) >= float(args.low):
                            if abs(float(sub_list[0]) - float(sub_list[1])) >= float(args.dif_range[0]):
                                if abs(float(sub_list[0]) - float(sub_list[1])) <= float(args.dif_range[1]):
                                    if max(float(sub_list[0]), float(sub_list[1])) <= float(args.hi):
                                        keep = True

                                        if use_iterator is None:  # avoid value inflation during tally routine
                                         # assign results to various maps
                                            if float(sub_list[0]) > float(sub_list[1]):
                                                self.foldchange_map[key][v] = "Turned OFF"
                                            elif float(sub_list[0]) < float(sub_list[1]):
                                                self.foldchange_map[key][v] = "Turned ON"
                                            else:
                                                self.foldchange_map[key][v] = "No Change"

                                            if labels[v] not in self.de_gene_list_by_stage.keys():
                                                try:
                                                    self.de_count_by_stage[labels[v]] = 1
                                                except IndexError:
                                                    print "\nCheck -time argument.\n"
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
                            if abs(float(sub_list[0]) - float(sub_list[1])) >= float(args.dif_range[0]):
                                if abs(float(sub_list[0]) - float(sub_list[1])) <= float(args.dif_range[1]):
                                    if max(float(sub_list[0]), float(sub_list[1])) <= float(args.hi):
                                        keep = True

                                        if use_iterator is None:  # avoid value inflation during tally routine

                                            # assign results to various maps
                                            self.foldchange_map[key][v] = round(np.log2(float(sub_list[1]) / float(sub_list[0])), ndigits=3)
                                            try:
                                                if labels[v] not in self.de_gene_list_by_stage.keys():
                                                    self.de_count_by_stage[labels[v]] = 1
                                                else:
                                                    self.de_count_by_stage[labels[v]] += 1
                                            except IndexError:
                                                print "\nCheck -time argument.\n"

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

            if keep is True:
                # TODO consolidate these two variables to clean up
                self.de_gene_list.append(key)
                self.filtered_data[key] = value
            else:
                self.unflagged_genes[key] = value

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
               self.filtered_data, \
               self.unflagged_genes

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
                for key, value in sorted(self.de_gene_list_by_stage.items()):
                    de_list_writer.writerow([key] + value)
            else:
                print "No flagged genes to write for gene_list by stage."

        with open(self.path + '_logFoldChange_matrix.txt', 'wb+') as logmatrix:
            if len(self.foldchange_map) != 0:
                foldwriter = csv.writer(logmatrix, delimiter='\t')
                foldwriter.writerow(['Gene'] + self.foldchange_map['Gene'])
                for key, value in sorted(self.foldchange_map.items()):
                    if key == 'Gene':
                        pass
                    else:
                        foldwriter.writerow([key] + value)
            else:
                print "No flagged genes."

