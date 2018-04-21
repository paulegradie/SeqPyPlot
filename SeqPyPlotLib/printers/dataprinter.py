

#TODO update this to print pandas dataframe results
class DataPrinter:

    """
    This should have a data container and analyzer instance passed upon instantiation.

    Available upon passing container object:
        self.args = args
        self.datatype = self.args.datatype
        self.data_frame_header = dict()
        self.gene_map = dict()
        self.ercc_map = dict()
        self.gene_map_const = dict()
        self.de_gene_list = []
        self.de_count_by_stage = None
        self.de_gene_list = None
        self.analyzed = self.analyzed()
        self.filtered_data = dict()
    """
    def __init__(self, args, container, analyzer):

        """
        :param args: Args object
        :param container: Data Container Object
        :param analyzer:  Data analyzer Object
        """
        self.container = container
        self.analyzer = analyzer
        self.filtered_data = analyzer.filtered_data
        self.args = args

        if self.args.out is not None:
            self.path = os.path.join(self.args.out, self.args.prefix)
        else:
            self.path = os.path.join('.', self.args.prefix)

    def write_ercc_data(self):

        with open((self.path + ".all_ercc.txt"), "wb+") as ercc:
            ercc_data_writer = csv.writer(ercc, delimiter='\t')
            for key, value in self.container.data_frame_header.items():
                ercc_data_writer.writerow(["ERCC"] + value)
            for key, value in self.container.ercc_map.items():
                ercc_data_writer.writerow([key] + value)

    def write_plot_data(self):
        with open(self.path + '_plotter_data.txt', 'wb+') as plot:
            plot_data_writer = csv.writer(plot, delimiter='\t')

            for key, value in self.container.data_frame_header.items():
                plot_data_writer.writerow([key] + value)

            for key, value in sorted(self.container.gene_map.items()):
                plot_data_writer.writerow([key] + value)

    def write_filtered_data(self):
        with open((self.path + '_filtered.txt'), "wb+") as filt:

            filtered_data_writer = csv.writer(filt, delimiter='\t')
            for key, value in self.container.data_frame_header.items():
                filtered_data_writer.writerow([key] + value)

            for key, value in self.filtered_data.items():
                filtered_data_writer.writerow([key] + value)

    def write_de_results(self):
        with open(self.path + "_de_gene_list.txt", 'wb+') as results_out:
            result_writer = csv.writer(results_out, delimiter='\n')
            for de_gene in sorted(self.analyzer.filtered_data):
                result_writer.writerow([de_gene])