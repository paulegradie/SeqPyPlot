
# TODO update this function to read the data in to pd dataframe
class CuffNormParser:
    """
    A class for parsing cufflnorm outputs and reading them in to the standard
    pandas df datafram structure used by the DataContainer
    """

    def parse_cuffnorm(self, infile):
            """
            This function expects output from Cuffnorm. If the ERCC option is set, it will also return any ERCCs.
                Returns: dictionaries with Gene:[ExpressionData]
            """
            self.data_frame_header = dict()
            self.gene_map = dict()
            self.ercc_map = dict()
            self.gene_map_const = dict()
            # open the files one at a time from within the loop
            try:
                with open(infile, 'rb') as dataFile:

                    data_frame = csv.reader(dataFile, delimiter='\t')

                    # for each row, ignore the data_frame_header row, then execute the following code
                    for gene in data_frame:

                        zero = False
                        count = 9
                        if "tracking_id" in str(gene[0]):  # generate data_frame_header
                            self.data_frame_header["Gene"] = []
                            for column_header in enumerate(gene):
                                if column_header[0] >= count:
                                    self.data_frame_header['Gene'] += [column_header[1]]
                                    count += 4

                        elif "ERCC-" in gene[3]:  # collect ERCC data if ERCC option set
                            self.ercc_map[gene[3]] = []
                            for column in enumerate(gene):
                                if column[0] == count:
                                    self.ercc_map[gene[3]] += [column[1]]
                                    count += 4

                        else:  # Collect all data

                            self.gene_map[gene[3]] = []
                            for column in enumerate(gene):
                                if int(column[0]) == count:
                                    #if value is less than 0.01, change it to zero...
                                    if column[1] <= 1.0:
                                        self.gene_map[gene[3]] += [0]
                                    else:
                                        self.gene_map[gene[3]] += [column[1]]
                                    if str(column[1]) <= 0.5:
                                        zero = True
                                    count += 4

                            if zero is False:
                                self.gene_map_const[gene[3]] = self.gene_map[gene[3]]
            except IOError:
                print "/nAre you loading cuffnorm data? - Try resetting data_type argument."

            self.reorder(self.data_frame_header)
            self.reorder(self.gene_map)
            self.reorder(self.ercc_map)

            return self.gene_map, self.ercc_map, self.data_frame_header