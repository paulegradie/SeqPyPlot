import sys

# TODO convert this to a pandas upload
def parse_plot_data(self, datafile):
    """This reads a preformatted data file output from SeqPyPlot."""
    # type: (input normalized data file) -> dictionary, key=gene name, value=normalized expression data
    # convert plotter data to a dictionary for quick access
    gene_map = dict()
    data_frame_header = dict()
    ercc_map = dict()

    if os.path.exists(datafile):
        with open(datafile, 'rb') as datafile:
            data_reader = csv.reader(datafile, delimiter='\t')

            header = True
            for row in data_reader:
                if header:
                    data_frame_header[row[0]] = row[1:]
                    header = False
                else:
                    temprow = row[1:]
                    finalrow = []
                    for i in temprow:
                        if i == '':
                            finalrow.append(None)
                        elif i < 1.0:
                            finalrow.append(0.0)
                        else:
                            finalrow.append(i)
                    gene_map[row[0].capitalize()] = finalrow
            if self.args.num == 1:
                pass
            elif self.args.num == 2:
                if self.args.unformatted_plot_data:
                    self.reorder(data_frame_header)
                    # print self.data_frame_header
                    self.reorder(gene_map)
                    self.reorder(ercc_map)
                    for key, value in gene_map.items():
                        series1 = value[:len(value) / 2]  # split the data
                        series2 = value[len(value) / 2:]
                        gene_map[key] = self.__average_flanking__(series1) + self.__average_flanking__(series2)
            else:
                print "num == more than 2 - does not support. DataConta. line 347"
                sys.exit()
    else:
        print "Couldn't open *_plotter_file.txt"
        sys.exit()

    return gene_map, ercc_map, data_frame_header

